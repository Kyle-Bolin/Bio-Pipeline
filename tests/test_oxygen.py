"""Oxygen-sensitive enzyme warnings, offline on the bundled models."""

import cobra
import pytest
from cobra.flux_analysis import pfba

from cellfactory import analyze, load_model, oxygen_sensitive_reactions, oxygen_sensitivity_check, set_medium
from cellfactory.cli import main
from cellfactory.oxygen import FluxSet, oxygen_caveats
from cellfactory.simulate import _biomass_id


@pytest.fixture
def core():
    return load_model("e_coli_core")


def _state(model, objective=None, knockouts=(), kind="state", label="this flux state"):
    objective = objective or _biomass_id(model)
    with model:
        for rid in knockouts:
            model.reactions.get_by_id(rid).knock_out()
        model.objective = objective
        fluxes = pfba(model).fluxes
    return FluxSet(label, fluxes, objective, tuple(knockouts), kind=kind)


def test_matches_by_name_where_the_model_has_no_ec_annotations(core):
    table = oxygen_sensitive_reactions(core)
    assert list(table["reaction"]) == ["PFL"]
    assert table.iloc[0]["matched_by"] == "enzyme name"
    assert "glycyl radical" in table.iloc[0]["enzyme"]


def test_repaired_enzymes_are_not_flagged():
    # B12 diol/glycerol dehydratases are inactivated in use and repaired by a reactivase; aerobic 1,3-propanediol
    # production with them is a commercial process, so they are deliberately out of the table.
    listed = set(oxygen_sensitive_reactions(load_model("iYS1720"))["reaction"])
    assert "PPDD" not in listed and "GLYCDH" not in listed


def test_matches_by_ec_number_and_lists_closed_reactions_too():
    table = oxygen_sensitive_reactions(load_model("iJO1366")).set_index("reaction")
    assert table.loc["PFL", "matched_by"] == "EC 2.3.1.54"
    assert table.loc["OBTFL", "matched_by"] == "EC 2.3.1.54"
    assert table.loc["RNTR1c2", "enzyme"].startswith("class III")
    assert {"HYD1pp", "HYD2pp", "HYD3pp"} <= set(table.index)
    # Formate-hydrogen lyase is listed, but iJO1366 ships it closed, so it can never be flagged.
    assert table.loc["FHL", "open"] is False or not table.loc["FHL", "open"]


@pytest.mark.parametrize("model_id,reaction", [
    ("e_coli_core", "PPS"),  # phosphoenolpyruvate synthase: "pyruvate synthase" is a substring, not a match
    ("e_coli_core", "NADTRHD"),  # transhydrogenase, not a hydrogenase
    ("e_coli_core", "GAPD"),  # a dehydrogenase
    ("iJO1366", "PDH"),
    ("iJO1366", "THD2pp"),
])
def test_lookalike_names_are_not_flagged(model_id, reaction):
    assert reaction not in set(oxygen_sensitive_reactions(load_model(model_id))["reaction"])


def test_anaerobic_flux_state_is_not_checked_at_all(core):
    set_medium(core, "glucose", 10.0, aerobic=False)
    state = _state(core)
    assert state.fluxes["PFL"] > 1.0  # the cell really is fermenting through PFL
    result = oxygen_sensitivity_check(core, [state])
    assert not result.flagged and result.skipped == ["this flux state"] and result.checked == []
    assert oxygen_caveats(result) == []


def test_aerobic_flux_state_through_pfl_is_flagged_as_load_bearing(core):
    # Deleting pyruvate dehydrogenase is what sends an aerobic cell's acetyl-CoA through PFL in this model.
    set_medium(core, "glucose", 10.0, aerobic=True)
    state = _state(core, knockouts=["PDH"], label="the cell's flux state at maximum growth")
    assert state.fluxes["PFL"] > 1.0
    result = oxygen_sensitivity_check(core, [state])
    assert result.flagged and result.checked == ["the cell's flux state at maximum growth"]
    row = result.table.iloc[0]
    assert row["reaction"] == "PFL"
    assert row["verdict"] in ("load-bearing", "infeasible without it")
    assert row["objective_without"] < row["objective"]
    caveats = oxygen_caveats(result)
    assert len(caveats) == 1 and "PFL" in caveats[0] and "no bounds were changed" in caveats[0].lower()


def test_a_design_that_only_works_through_pfl_is_called_unbuildable(core):
    set_medium(core, "glucose", 10.0, aerobic=True)
    state = _state(core, knockouts=["PDH"], kind="design", label="design 1 (PDH)")
    result = oxygen_sensitivity_check(core, [state])
    (caveat,) = oxygen_caveats(result)
    assert "not buildable" in caveat and "design 1 (PDH)" in caveat


def test_several_affected_designs_give_one_caveat(core):
    set_medium(core, "glucose", 10.0, aerobic=True)
    states = [_state(core, knockouts=["PDH"], kind="design", label=f"design {i} (PDH)") for i in (1, 2, 3)]
    caveats = oxygen_caveats(oxygen_sensitivity_check(core, states))
    assert len(caveats) == 1 and "3 knockout design(s)" in caveats[0]


def test_flux_that_costs_nothing_is_called_avoidable():
    # iJO1366 has both the aerobic (class I, RNDR*) and the anaerobic (class III, RNTR*c2) ribonucleotide
    # reductase, and parsimonious FBA picks the class III route in air anyway. Growth does not depend on it,
    # so the verdict is "avoidable".
    model = load_model("iJO1366")
    set_medium(model, "glucose", 10.0, aerobic=True)
    result = oxygen_sensitivity_check(model, [_state(model)])
    rnr = result.table[result.table["reaction"].str.startswith("RNTR")]
    assert not rnr.empty and set(rnr["verdict"]) == {"avoidable"}
    assert (rnr["cost"] < 1e-3).all()
    caveat = next(c for c in oxygen_caveats(result) if "RNTR1c2" in c)
    assert "does not depend on it" in caveat


def test_the_check_never_changes_the_model(core):
    set_medium(core, "glucose", 10.0, aerobic=True)
    before = {r.id: r.bounds for r in core.reactions}
    medium_before = dict(core.medium)
    growth_before = core.slim_optimize()
    oxygen_sensitivity_check(core, [_state(core, knockouts=["PDH"]), _state(core)])
    assert {r.id: r.bounds for r in core.reactions} == before
    assert dict(core.medium) == medium_before
    assert core.slim_optimize() == pytest.approx(growth_before)


def test_a_model_without_oxygen_sensitive_reactions_is_reported_empty(core):
    core.remove_reactions([core.reactions.PFL])
    result = oxygen_sensitivity_check(core, [_state(core)])
    assert result.candidates.empty and not result.flagged and oxygen_caveats(result) == []


def test_analyze_skips_anaerobic_runs(core):
    anaerobic = analyze(core.copy(), "glucose", "ethanol", aerobic=False, qa_gate=False)
    assert "Oxygen-sensitive enzymes" not in anaerobic.checks
    assert not [c for c in anaerobic.caveats if "aerobic cell" in c]


def test_analyze_wires_the_check_into_an_aerobic_run(core):
    # With PDH gone, the aerobic cell's own best flux state goes through PFL, so the pipeline must flag it
    # without anyone building a FluxSet by hand.
    core.reactions.PDH.knock_out()
    report = analyze(core, "glucose", "succinate", qa_gate=False)
    table = report.checks["Oxygen-sensitive enzymes"]
    assert "PFL" in set(table["reaction"])
    assert set(table["where"]) <= {"the cell's flux state at maximum growth", "the max-yield route to EX_succ_e"}
    assert [c for c in report.caveats if "PFL" in c]


def test_analyze_flags_a_design_the_cell_could_not_run(core):
    pytest.importorskip("straindesign")
    report = analyze(core, "glucose", "succinate", qa_gate=False, design="optknock", design_knockouts_max=3)
    designs = next(df for name, df in report.checks.items() if name.startswith("Knockout designs"))
    assert not designs.empty
    table = report.checks["Oxygen-sensitive enzymes"]
    flagged = table[table["kind"] == "design"]
    # The designs delete PDH and route acetyl-CoA through PFL, which is PR #17's artifact.
    assert not flagged.empty and "PFL" in set(flagged["reaction"])
    assert len(set(flagged["where"])) <= 5  # MAX_OXYGEN_DESIGNS
    assert [c for c in report.caveats if "not buildable" in c]


def test_the_resolve_keeps_the_bounds_the_flux_state_was_solved_under(core):
    # Maximum succinate with growth floored, which is how the tools build a production flux state. Without the
    # floor the re-solve answers an easier question and PFL looks avoidable, though with the floor kept there is
    # no solution at all without it.
    from cellfactory import resolve_exchange
    from cellfactory.tools import _fluxes

    set_medium(core, "glucose", 10.0, aerobic=True)
    core.reactions.PDH.knock_out()
    target = resolve_exchange(core, "succinate")
    biomass = _biomass_id(core)
    fluxes = _fluxes(core, target, 0.3)
    assert fluxes["PFL"] > 1.0
    floored = FluxSet("a floored production state", fluxes, target, lower_bounds={biomass: 0.3})
    assert oxygen_sensitivity_check(core, [floored]).table.iloc[0]["verdict"] == "infeasible without it"
    # The same state judged without the floor is the bug the review caught: it reads as harmless.
    unfloored = FluxSet("the same state, floor forgotten", fluxes, target)
    assert oxygen_sensitivity_check(core, [unfloored]).table.iloc[0]["verdict"] == "avoidable"


def test_analyze_can_turn_the_check_off(core, monkeypatch):
    called = []
    monkeypatch.setattr("cellfactory.pipeline.oxygen_sensitivity_check",
                        lambda *a, **k: called.append(a) or pytest.fail("should not run"))
    analyze(core, "glucose", "succinate", qa_gate=False, oxygen_check=False)
    assert called == []


def test_cli_no_oxygen_flag(core, monkeypatch, capsys):
    seen = {}

    def fake_analyze(model, substrate, product, **kwargs):
        seen.update(kwargs)
        return analyze(model, substrate, product, **{**kwargs, "qa_gate": False})

    monkeypatch.setattr("cellfactory.cli.analyze", fake_analyze)
    assert main(["run", "-p", "succinate", "--no-oxygen", "--no-qa"]) == 0
    assert seen["oxygen_check"] is False
    assert main(["run", "-p", "succinate", "--no-qa"]) == 0
    assert seen["oxygen_check"] is True


def test_oxygen_uptake_decides_not_the_medium(core):
    # Oxygen is in the medium, but this flux state doesn't use any, so the enzymes are fine.
    set_medium(core, "glucose", 10.0, aerobic=True)
    with core:
        core.reactions.EX_o2_e.bounds = (0.0, 0.0)
        core.reactions.PDH.knock_out()
        fluxes = pfba(core).fluxes
    state = FluxSet("an oxygen-free flux state", fluxes, _biomass_id(core))
    assert fluxes["PFL"] > 1.0
    result = oxygen_sensitivity_check(core, [state])
    assert not result.flagged and result.skipped == ["an oxygen-free flux state"]


def test_report_text_and_json_carry_the_check(core):
    report = analyze(core, "glucose", "succinate", qa_gate=False)
    report.checks["Oxygen-sensitive enzymes"] = oxygen_sensitivity_check(
        core, [_state(core, knockouts=["PDH"])]).table
    assert "Oxygen-sensitive enzymes" in report.to_text()
    assert isinstance(cobra.io.json.json.dumps(report.to_dict(), default=float), str)


def test_flux_distribution_tool_works_without_a_product():
    # Regression: the flux set's objective used to be the label "growth", which cobra rejects as a reaction id,
    # so the whole tool call failed as soon as anything was flagged.
    from cellfactory import tools

    out = tools.flux_distribution(model="iJO1366")
    assert out["growth"] > 0.5
    assert out["oxygen_warnings"]
    assert not [c for c in out["caveats"] if "did not run" in c]


def test_flux_distribution_tool_survives_a_broken_check(monkeypatch):
    from cellfactory import tools

    def boom(*args, **kwargs):
        raise RuntimeError("solver exploded")

    monkeypatch.setattr("cellfactory.tools.oxygen_sensitivity_check", boom)
    out = tools.flux_distribution(model="e_coli_core")
    assert out["growth"] > 0.5 and out["oxygen_warnings"] == []
    assert [c for c in out["caveats"] if "did not run" in c]

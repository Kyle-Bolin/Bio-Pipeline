"""Max-min Driving Force: the linear program on hand-made chains; eQuilibrator itself runs offline."""

import json
import os
import sys

import cobra
import pytest

from cellfactory import MDFResult, ThermoUnavailable, analyze, load_model, max_min_driving_force, product_route_fluxes
from cellfactory import mdf as mdf_module
from cellfactory import thermodynamics
from cellfactory.cli import main

RT = 8.314462618e-3 * 298.15
SPAN = RT * 9.2103404  # ln(1e-2 / 1e-6): how far one metabolite's concentration can move ΔG'


def chain_model() -> cobra.Model:
    """A -> B -> C, A -> D (made by group contribution), a hydrolysis with water, and a transport step."""
    model = cobra.Model("chain")
    a, b, c, d, w = (cobra.Metabolite(f"{x}_c", compartment="c") for x in ("a", "b", "c", "d", "h2o"))
    w.formula = "H2O"
    c_e = cobra.Metabolite("c_e", compartment="e")
    for rid, mets in {"S1": {a: -1, b: 1}, "S2": {b: -1, c: 1}, "GC": {a: -1, d: 1}, "HYD": {d: -1, w: -1, c: 1},
                      "CTEX": {c: -1, c_e: 1}}.items():
        rxn = cobra.Reaction(rid, lower_bound=-1000, upper_bound=1000)
        model.add_reactions([rxn])
        rxn.add_metabolites(mets)
    return model


def estimator(dg):
    def estimate(rxn):
        if rxn.id not in dg:
            raise LookupError("not in eQuilibrator's compound set")
        value, sigma, method = dg[rxn.id]
        return value, sigma, method, {"h2o_c"}
    return estimate


UPHILL = {"S1": (20.0, 1.0, "component contribution"), "S2": (20.0, 1.0, "component contribution")}
DOWNHILL = {"S1": (-10.0, 1.0, "component contribution"), "S2": (-10.0, 1.0, "component contribution")}


def test_each_step_feasible_alone_but_not_together():
    # +20 kJ/mol is within the 22.8 kJ/mol a concentration ratio can overturn, so each step alone can run...
    model = chain_model()
    for step in ("S1", "S2"):
        assert max_min_driving_force(model, {step: 1.0}, estimator=estimator(UPHILL)).mdf > 0
    # ...but together they climb 40 kJ/mol with only one span of A-to-C ratio to pay for it.
    result = max_min_driving_force(model, {"S1": 1.0, "S2": 1.0}, estimator=estimator(UPHILL))
    assert result.mdf == pytest.approx(-(40 - SPAN) / 2, abs=1e-3)
    assert result.mdf_optimistic == pytest.approx(-(40 - 2 * 1.96 - SPAN) / 2, abs=1e-3)
    assert result.verdict == "infeasible"
    assert set(result.bottlenecks) == {"S1", "S2"}


def test_downhill_chain_and_concentrations():
    result = max_min_driving_force(chain_model(), {"S1": 1.0, "S2": 1.0}, estimator=estimator(DOWNHILL))
    assert result.mdf == pytest.approx((20 + SPAN) / 2, abs=1e-3)
    assert result.verdict == "feasible"
    conc = result.concentrations.set_index("metabolite")["concentration"]
    assert conc["a_c"] == pytest.approx(1e-2) and conc["c_c"] == pytest.approx(1e-6)
    steps = result.reactions.set_index("reaction")
    assert steps.loc["S1", "driving_force"] == pytest.approx(result.mdf, abs=1e-6)


def test_reverse_flux_uses_the_reverse_direction():
    # Running the uphill chain backwards is the downhill chain.
    result = max_min_driving_force(chain_model(), {"S1": -2.0, "S2": -2.0}, estimator=estimator(UPHILL))
    assert result.mdf == pytest.approx((40 + SPAN) / 2, abs=1e-3)


def test_low_driving_force_tier_follows_temperature():
    # Pinned concentrations make the MDF exactly 6 kJ/mol: above RT ln 10 at 25 °C (5.7), below it at 60 °C (6.4).
    dg = {"S1": (-6.0, 0.1, "component contribution")}
    pinned = {"a_c": (1e-3, 1e-3), "b_c": (1e-3, 1e-3)}
    cool = max_min_driving_force(chain_model(), {"S1": 1.0}, estimator=estimator(dg), bounds=pinned)
    warm = max_min_driving_force(chain_model(), {"S1": 1.0}, estimator=estimator(dg), bounds=pinned,
                                 conditions=thermodynamics.Conditions(temperature=333.15))
    assert cool.mdf == pytest.approx(6.0) and warm.mdf == pytest.approx(6.0)
    assert cool.verdict == "feasible" and warm.verdict == "low driving force"


def test_failed_lp_becomes_a_caveat(monkeypatch):
    _needs_equilibrator()

    def broken(*args, **kwargs):
        raise RuntimeError("MDF linear program did not solve (infeasible)")

    monkeypatch.setattr(mdf_module, "_solve", broken)
    report = analyze("e_coli_core", "glucose", "ethanol", aerobic=False, qa_gate=False, online_ids=False, mdf=True)
    assert any("Max-min driving force not computed: MDF linear program did not solve" in c for c in report.caveats)


def test_low_driving_force_tier():
    dg = {"S1": (SPAN - 2.0, 0.1, "component contribution")}  # uphill, but a full A:B ratio leaves 2 kJ/mol downhill
    result = max_min_driving_force(chain_model(), {"S1": 1.0}, estimator=estimator(dg))
    assert result.mdf == pytest.approx(2.0, abs=1e-3)
    assert result.verdict == "low driving force"


def test_bounds_override_pins_a_metabolite():
    bounds = {"a_c": (1e-3, 1e-3), "c_c": (1e-3, 1e-3)}  # same concentration at both ends: no ratio to exploit
    result = max_min_driving_force(chain_model(), {"S1": 1.0, "S2": 1.0}, estimator=estimator(DOWNHILL),
                                   bounds=bounds)
    assert result.mdf == pytest.approx(10.0, abs=1e-3)


def test_water_has_no_concentration():
    dg = {"HYD": (-5.0, 0.1, "component contribution")}
    result = max_min_driving_force(chain_model(), {"HYD": 1.0}, estimator=estimator(dg))
    assert result.mdf == pytest.approx(5.0 + SPAN, abs=1e-3)  # d and c vary, water does not
    assert "h2o_c" not in set(result.concentrations["metabolite"])


def test_group_contribution_bottleneck_is_uncertain_not_a_kill():
    dg = {"GC": (40.0, 1.0, "group contribution")}
    result = max_min_driving_force(chain_model(), {"GC": 1.0}, estimator=estimator(dg))
    assert result.mdf_optimistic < 0
    assert result.verdict == "uncertain"
    assert "group-contribution" in result.reason


def test_steps_without_estimate_are_listed_and_caveated():
    result = max_min_driving_force(chain_model(), {"S1": 1.0, "S2": 1.0, "CTEX": 1.0, "GC": 0.0},
                                   estimator=estimator(DOWNHILL))
    assert result.excluded == ["CTEX"]  # GC carries no flux, so it is not part of the route
    assert result.reactions.set_index("reaction").loc["CTEX", "note"].startswith("transport")
    notes = mdf_module.mdf_caveats(result, "the route")
    assert any("1 of 3 steps" in n and "can only be lower" in n for n in notes)


def test_nothing_estimated():
    result = max_min_driving_force(chain_model(), {"CTEX": 1.0}, estimator=estimator({}))
    assert result.verdict == "not computed" and result.mdf is None
    assert "not computed" in mdf_module.mdf_caveats(result, "the route")[0]


def test_route_fluxes_for_core_ethanol():
    model = load_model("e_coli_core")
    model.reactions.EX_o2_e.lower_bound = 0
    fluxes = product_route_fluxes(model, "EX_etoh_e")
    assert fluxes["EX_etoh_e"] == pytest.approx(20.0, abs=1e-3)  # glucose -> 2 ethanol at 10 glucose
    assert fluxes["PFK"] > 0 and "Biomass_Ecoli_core" not in fluxes.index
    assert model.objective.expression.has(model.reactions.Biomass_Ecoli_core.forward_variable)  # left unchanged


def test_route_fluxes_refuses_an_impossible_product():
    model = load_model("e_coli_core")
    model.reactions.EX_glc__D_e.lower_bound = 0
    with pytest.raises(ValueError, match="can't make"):
        product_route_fluxes(model, "EX_etoh_e")


def test_missing_equilibrator_becomes_a_caveat(monkeypatch, capsys):
    monkeypatch.setitem(sys.modules, "equilibrator_api", None)
    thermodynamics._equilibrator.cache_clear()
    try:
        report = analyze("e_coli_core", "glucose", "ethanol", qa_gate=False, online_ids=False, mdf=True)
        assert "Max-min driving force" not in report.checks
        assert any("Max-min driving force not computed: equilibrator-api is not installed" in c
                   for c in report.caveats)
        assert main(["run", "--product", "ethanol", "--no-qa", "--offline", "--mdf"]) == 0
        assert "Max-min driving force not computed" in capsys.readouterr().out
    finally:
        thermodynamics._equilibrator.cache_clear()


# ---- eQuilibrator itself, from the data it ships with ----

def _needs_equilibrator():
    try:
        thermodynamics._equilibrator(False)
    except ThermoUnavailable as err:
        (pytest.fail if os.environ.get("CI") else pytest.skip)(f"eQuilibrator not available: {err}")


GLYCOLYSIS = {"PGI": 1, "PFK": 1, "FBA": 1, "TPI": 1, "GAPD": 2, "PGK": -2, "PGM": -2, "ENO": 2, "PYK": 2}


def test_glycolysis_runs_forward_but_not_backward():
    _needs_equilibrator()
    model = load_model("e_coli_core")
    table = thermodynamics.reaction_thermodynamics(model, list(GLYCOLYSIS))
    assert (table["thermo_direction"] == "either").all()  # no single step is one-way on its own
    forward = max_min_driving_force(model, GLYCOLYSIS, table=table)
    backward = max_min_driving_force(model, {r: -v for r, v in GLYCOLYSIS.items()}, table=table)
    assert isinstance(forward, MDFResult) and forward.verdict == "feasible"
    assert backward.verdict == "infeasible" and "PYK" in backward.bottlenecks


def test_analyze_reports_mdf():
    _needs_equilibrator()
    report = analyze("e_coli_core", "glucose", "ethanol", aerobic=False, qa_gate=False, online_ids=False, mdf=True)
    steps = report.checks["Max-min driving force"].set_index("reaction")
    assert steps.loc["EX_etoh_e", "note"] == "boundary reaction"
    assert steps["bottleneck"].any()
    assert any(c.startswith("Max-min driving force of the zero-growth max-yield route to EX_etoh_e") and "feasible" in c
               for c in report.caveats)


def test_analyze_reuses_the_thermo_table(monkeypatch):
    _needs_equilibrator()
    calls = []
    real = thermodynamics.reaction_thermodynamics

    def counting(*args, **kwargs):
        calls.append(args[1] if len(args) > 1 else kwargs.get("reactions"))
        return real(*args, **kwargs)

    monkeypatch.setattr("cellfactory.pipeline.reaction_thermodynamics", counting)
    monkeypatch.setattr(mdf_module, "reaction_thermodynamics", counting)
    analyze("e_coli_core", "glucose", "ethanol", aerobic=False, qa_gate=False, online_ids=False, thermo=True, mdf=True)
    assert calls == [None]  # the whole-model table, once; MDF takes its ΔG'° from there


def test_cli_flag(capsys):
    _needs_equilibrator()
    assert main(["run", "--product", "ethanol", "--anaerobic", "--no-qa", "--offline", "--mdf", "--json", "-"]) == 0
    out = json.loads(capsys.readouterr().out)
    rows = {r["reaction"]: r for r in out["checks"]["Max-min driving force"]}
    assert rows["PFK"]["driving_force"] > 0
    assert any(c.startswith("Max-min driving force") for c in out["caveats"])

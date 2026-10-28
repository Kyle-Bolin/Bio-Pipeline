"""Reaction thermodynamics: the direction rules run on a fake estimator; eQuilibrator itself runs offline."""

import json
import os
import socket
import sys

import cobra
import pytest

from cellfactory import ThermoUnavailable, analyze, apply_directionality, load_model, reaction_thermodynamics
from cellfactory import thermodynamics
from cellfactory.cli import main


def toy_model() -> cobra.Model:
    """Single-compartment reactions with ΔG'° chosen to land on each verdict, plus ones that must be skipped."""
    model = cobra.Model("toy")
    a, b, c, w = (cobra.Metabolite(f"{x}_c", compartment="c") for x in ("a", "b", "c", "h2o"))
    a_e = cobra.Metabolite("a_e", compartment="e")
    reactions = {
        "DOWNHILL": ({a: -1, b: 1}, (-1000, 1000)),  # reversible in the model, one-way by physics
        "NEAR_EQ": ({a: -1, b: 1}, (-1000, 1000)),
        "UPHILL_FORCED": ({b: -1, c: 1}, (0, 1000)),  # model only allows the impossible direction
        "UPHILL_REVERSE": ({b: -1, c: 1}, (-1000, 0)),  # model already runs it the possible way
        "HYDROLYSIS": ({a: -1, w: -1, c: 1}, (-1000, 1000)),  # water's concentration must not count
        "GROUPS": ({a: -1, c: 1}, (-1000, 1000)),  # one-way, but only by a group-contribution estimate
        "UNKNOWN": ({c: -1, b: 1}, (-1000, 1000)),
        "FUZZY": ({a: -1, c: 1}, (-1000, 1000)),
        "ATEX": ({a_e: -1, a: 1}, (-1000, 1000)),
    }
    for rid, (mets, bounds) in reactions.items():
        rxn = cobra.Reaction(rid, lower_bound=bounds[0], upper_bound=bounds[1])
        model.add_reactions([rxn])
        rxn.add_metabolites(mets)
    model.add_boundary(a_e, type="exchange")
    return model


DG = {  # reaction -> (dg0_prime, sigma, method)
    "DOWNHILL": (-40.0, 1.0, "component contribution"),
    "NEAR_EQ": (5.0, 1.0, "component contribution"),
    "UPHILL_FORCED": (40.0, 1.0, "component contribution"),
    "UPHILL_REVERSE": (40.0, 1.0, "component contribution"),
    "HYDROLYSIS": (-30.0, 1.0, "component contribution"),
    "GROUPS": (-60.0, 3.0, "group contribution"),
    "FUZZY": (-40.0, 5000.0, "group contribution"),
}


def fake_estimator(rxn):
    if rxn.id not in DG:
        raise LookupError("not in eQuilibrator's compound set: c_c")
    return (*DG[rxn.id], {"h2o_c"})


@pytest.fixture
def table():
    return reaction_thermodynamics(toy_model(), estimator=fake_estimator).set_index("reaction")


def test_directions_and_verdicts(table):
    # 1 µM-10 mM spans ln(1e4) * RT = 22.8 kJ/mol per metabolite at 25 °C, so a 1-in/1-out
    # reaction needs |ΔG'°| > 22.8 kJ/mol plus 1.96 sigma to be one-way.
    assert table.loc["DOWNHILL", "thermo_direction"] == "forward"
    assert table.loc["DOWNHILL", "verdict"] == "can tighten"
    assert table.loc["NEAR_EQ", "thermo_direction"] == "either"
    assert table.loc["NEAR_EQ", "verdict"] == "agree"
    assert table.loc["UPHILL_FORCED", "thermo_direction"] == "reverse"
    assert table.loc["UPHILL_FORCED", "verdict"] == "CONFLICT"
    assert table.loc["UPHILL_REVERSE", "verdict"] == "agree"
    assert table.loc["DOWNHILL", "dg_prime_max"] == pytest.approx(-40.0 + 8.314462618e-3 * 298.15 * 9.2103, abs=0.01)


def test_water_has_no_concentration_range(table):
    # Only a and c vary, so the range is the same width as for a 1-in/1-out reaction.
    row = table.loc["HYDROLYSIS"]
    assert row["dg_prime_max"] - row["dg_prime_min"] == pytest.approx(2 * 22.84, abs=0.05)
    assert row["thermo_direction"] == "forward"


def test_reactions_without_estimate_say_why(table):
    assert table.loc["ATEX", "note"].startswith("transport")
    assert table.loc["EX_a_e", "note"] == "boundary reaction"
    assert "not in eQuilibrator" in table.loc["UNKNOWN", "note"]
    assert table.loc["FUZZY", "note"] == "estimate too uncertain to decide a direction"
    for rid in ("ATEX", "EX_a_e", "UNKNOWN", "FUZZY"):
        assert table.loc[rid, "verdict"] == ""


def test_apply_only_narrows_and_skips_conflicts_and_group_estimates(table):
    model = toy_model()
    before = {r.id: r.bounds for r in model.reactions}
    changed = apply_directionality(model, table.reset_index())
    assert set(changed["reaction"]) == {"DOWNHILL", "HYDROLYSIS"}
    assert model.reactions.DOWNHILL.bounds == (0, 1000)
    assert model.reactions.UPHILL_FORCED.bounds == before["UPHILL_FORCED"]  # conflicts are for a person to judge
    assert model.reactions.GROUPS.bounds == before["GROUPS"]
    for rxn in model.reactions:
        lo, hi = before[rxn.id]
        assert lo <= rxn.lower_bound and rxn.upper_bound <= hi

    apply_directionality(model, table.reset_index(), group_contribution=True)
    assert model.reactions.GROUPS.bounds == (0, 1000)


def test_reverse_tightening_closes_upper_bound():
    model = toy_model()
    model.reactions.UPHILL_FORCED.bounds = (-1000, 1000)
    apply_directionality(model, reaction_thermodynamics(model, estimator=fake_estimator))
    assert model.reactions.UPHILL_FORCED.bounds == (-1000, 0)


def test_caveats_name_conflicts(table):
    notes = thermodynamics.thermo_caveats(table.reset_index(), tightened=2)
    assert "for 6 of 9 reactions" in notes[0]
    assert any("UPHILL_FORCED" in n for n in notes)
    assert any("2 reversible reaction(s) were made one-way" in n for n in notes)


def test_missing_equilibrator_becomes_a_caveat(monkeypatch):
    def unavailable(download):
        raise ThermoUnavailable("equilibrator-api is not installed (equilibrator_api)")

    monkeypatch.setattr(thermodynamics, "_equilibrator", unavailable)
    report = analyze("e_coli_core", "glucose", "ethanol", qa_gate=False, online_ids=False, thermo_bounds=True)
    assert "Reaction thermodynamics" not in report.checks
    assert any("not checked: equilibrator-api is not installed" in c for c in report.caveats)
    assert report.production.growth > 0.8


def test_missing_package_degrades_without_traceback(monkeypatch, capsys):
    # The real failure, not a mocked one: equilibrator_api can't be imported.
    monkeypatch.setitem(sys.modules, "equilibrator_api", None)
    thermodynamics._equilibrator.cache_clear()
    try:
        with pytest.raises(ThermoUnavailable, match="not installed"):
            reaction_thermodynamics(load_model("e_coli_core"), reactions=["PGI"])
        report = analyze("e_coli_core", "glucose", "ethanol", qa_gate=False, online_ids=False, thermo=True)
        assert any("not checked: equilibrator-api is not installed" in c for c in report.caveats)
        assert main(["run", "--product", "ethanol", "--no-qa", "--offline", "--thermo"]) == 0
        assert "Reaction thermodynamics not checked" in capsys.readouterr().out
    finally:
        thermodynamics._equilibrator.cache_clear()


# ---- eQuilibrator itself, from the data it ships with ----

def _needs_equilibrator():
    try:
        thermodynamics._equilibrator(False)
    except ThermoUnavailable as err:
        # CI installs eQuilibrator, so there a missing one is a broken step, not a reason to go green.
        (pytest.fail if os.environ.get("CI") else pytest.skip)(f"eQuilibrator not available: {err}")


@pytest.fixture(scope="module")
def core_table():
    _needs_equilibrator()
    return reaction_thermodynamics(load_model("e_coli_core")).set_index("reaction")


def test_runs_without_network(monkeypatch):
    _needs_equilibrator()
    thermodynamics._equilibrator.cache_clear()

    def refuse(*args, **kwargs):
        raise AssertionError("eQuilibrator must not touch the network unless download=True")

    monkeypatch.setattr(socket.socket, "connect", refuse)
    monkeypatch.setattr(socket, "create_connection", refuse)
    try:
        table = reaction_thermodynamics(load_model("e_coli_core"), reactions=["PGI", "PFK"])
    finally:
        thermodynamics._equilibrator.cache_clear()
    assert table["dg0_prime"].notna().all()


def test_core_model_estimates(core_table):
    t = core_table
    assert t.loc["PFK", "method"] == "component contribution"
    assert t.loc["FBA", "dg0_prime"] > 0 and t.loc["FBA", "thermo_direction"] == "either"
    assert t.loc["ATPM", "thermo_direction"] == "forward"
    assert t.loc["GLCpts", "note"].startswith("transport")
    assert t.loc["Biomass_Ecoli_core", "note"] == "biomass pseudo-reaction"
    # The core model's fumarate reductase uses ubiquinol; eQuilibrator's quinone estimate makes it look impossible.
    assert t.loc["FRD7", "verdict"] == "CONFLICT" and t.loc["FRD7", "method"] == "group contribution"
    # The exact count follows eQuilibrator's bundled compound set; 49 of 74 with 0.8.1.
    assert t["thermo_direction"].notna().sum() >= 40


def test_conditions_change_the_estimate():
    _needs_equilibrator()
    model = load_model("e_coli_core")
    # Glyceraldehyde-3-phosphate dehydrogenase releases a proton, so lower pH makes it less favourable.
    at = {ph: reaction_thermodynamics(model, ["GAPD"], conditions=thermodynamics.Conditions(p_h=ph))
          .loc[0, "dg0_prime"] for ph in (6.5, 7.5)}
    assert at[6.5] > at[7.5]


def test_analyze_tightens_malate_oxidase_in_ijo1366():
    _needs_equilibrator()
    model = load_model("iJO1366")
    assert model.reactions.MOX.lower_bound < 0
    report = analyze(model, "glucose", "ethanol", qa_gate=False, online_ids=False, thermo_bounds=True)
    assert model.reactions.MOX.bounds == (0, 1000)
    table = report.checks["Reaction thermodynamics"].set_index("reaction")
    assert table.loc["MOX", "verdict"] == "can tighten"
    assert any("were made one-way" in c for c in report.caveats)
    assert report.production.growth == pytest.approx(0.982, abs=0.01)


def test_cli_flag(capsys):
    _needs_equilibrator()
    assert main(["run", "--product", "ethanol", "--no-qa", "--offline", "--thermo", "--json", "-"]) == 0
    out = json.loads(capsys.readouterr().out)
    rows = {r["reaction"]: r for r in out["checks"]["Reaction thermodynamics"]}
    assert rows["PFK"]["thermo_direction"] == "either"
    assert any("eQuilibrator" in c for c in out["caveats"])

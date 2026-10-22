"""Any substrate, any organism, any product."""

import json

import pytest

from cellfactory import analyze, load_model, resolve_exchange, set_medium
from cellfactory.cli import main


@pytest.fixture
def core():
    return load_model("e_coli_core")


@pytest.mark.parametrize("query", ["ethanol", "Ethanol", "etoh", "etoh_e", "EX_etoh_e"])
def test_resolve_exchange_accepts_names_and_ids(core, query):
    assert resolve_exchange(core, query) == "EX_etoh_e"


def test_unknown_compound_suggests_matches(core):
    with pytest.raises(KeyError, match="Close matches"):
        resolve_exchange(core, "gluc")  # not exact, but "glucose" names exist


def test_internal_product_gets_demand_reaction(core):
    rxn_id = resolve_exchange(core, "atp_c", create_demand=True)
    assert rxn_id == "DM_atp_c"


def test_set_medium_swaps_the_carbon_source(core):
    set_medium(core, "fructose", uptake=8)
    assert core.medium.get("EX_fru_e") == 8
    assert "EX_glc__D_e" not in core.medium
    assert core.slim_optimize() > 0.1


def test_xylose_to_ethanol_theoretical_yield():
    # 3 xylose (C5) -> 5 ethanol (C2) + 5 CO2, so 5/3 mol/mol
    report = analyze(load_model("iJO1366"), "xylose", "ethanol", aerobic=False)
    assert report.production.max_theoretical_yield == pytest.approx(5 / 3, abs=1e-6)


def test_other_organism_salmonella_glycerol():
    report = analyze("iYS1720", "glycerol", "ethanol", aerobic=False)
    assert report.model_id == "iYS1720" and report.substrate == "EX_glyc_e"
    assert report.production.growth > 0


def test_analyze_with_knockouts_and_checks(core):
    report = analyze(core, "glucose", "succinate", aerobic=False, max_knockouts=1, run_checks=True)
    assert report.knockouts.steps[0].knocked_out == []
    assert report.checks["Free-energy checks"]["ok"].all()
    assert "Predicted essential genes" in report.checks


def test_cli_writes_json(tmp_path, capsys):
    out = tmp_path / "report.json"
    assert main(["run", "-p", "ethanol", "--anaerobic", "--json", str(out)]) == 0
    assert "guaranteed yield" in capsys.readouterr().out
    data = json.loads(out.read_text())
    assert data["product"] == "EX_etoh_e"
    assert data["production"]["max_theoretical_yield"] == pytest.approx(2.0, abs=1e-6)


def test_cli_prints_json_to_stdout(capsys):
    assert main(["run", "-p", "ethanol", "--anaerobic", "--json", "-"]) == 0
    assert json.loads(capsys.readouterr().out)["product"] == "EX_etoh_e"


def test_cli_reports_missing_compound(capsys):
    assert main(["run", "-p", "lycopene"]) == 1
    assert "not a metabolite" in capsys.readouterr().err

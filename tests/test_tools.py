"""The agent-facing tools (tools.py) and the MCP server that exposes them."""

import asyncio
import inspect
import json

import pytest

import cellfactory
from cellfactory import knockouts, pipeline, tools

from .test_kinetics import BRENDA as BRENDA_SAMPLE  # synthetic BRENDA fixture, so no network and no licence


def test_simulate_production_returns_json_safe_yields():
    out = tools.simulate_production("ethanol", aerobic=False)
    json.dumps(out, allow_nan=False)  # must survive strict JSON
    assert out["product"] == "EX_etoh_e"
    assert out["production"]["max_theoretical_yield"] == pytest.approx(2.0, abs=1e-3)
    assert 0 < out["production"]["guaranteed_yield"] < 2


def test_calls_do_not_leak_state_between_each_other():
    tools.simulate_production("atp_c")  # adds a demand reaction inside the call
    model, _ = tools._model("e_coli_core")
    assert "DM_atp_c" not in model.reactions
    assert model.medium.get("EX_o2_e") == 1000.0
    anaerobic = tools.simulate_production("ethanol", aerobic=False)["production"]["growth"]
    aerobic = tools.simulate_production("ethanol")["production"]["growth"]
    assert aerobic > anaerobic


def test_search_knockouts_summarises_the_design():
    out = tools.search_knockouts("succinate", aerobic=False, max_knockouts=2)
    summary = out["summary"]
    assert summary["improved"]
    assert summary["best_guaranteed_yield"] > summary["wild_type_guaranteed_yield"]
    assert summary["genes"]


def test_search_metabolites_lists_exchanges():
    out = tools.search_metabolites("succin")
    exchanges = {m["exchange"] for m in out["metabolites"]}
    assert "EX_succ_e" in exchanges


def test_envelope_and_checks():
    points = tools.production_envelope("ethanol", aerobic=False, points=5)["points"]
    assert len(points) == 5 and points[0]["growth"] == 0
    checks = tools.check_model()
    assert checks["ok"]
    assert checks["predicted_essential_genes"]["total"] > 0


def test_pathway_adds_a_non_native_product():
    names = {p["name"] for p in tools.list_pathways()["pathways"]}
    assert names == {"vanillin", "lycopene"}
    out = tools.simulate_production("lycopene", model="iJO1366", pathway="lycopene")
    json.dumps(out, allow_nan=False)
    assert out["pathway"]["steps"] == {"HET_CRTE": "added", "HET_CRTB": "added", "HET_CRTI": "added"}
    # Regression value from this code (iJO1366, aerobic glucose), not a literature number; it matches
    # tests/test_heterologous.py, which explains it from the MEP pathway's ~0.8 IPP per glucose.
    assert out["production"]["max_theoretical_yield"] == pytest.approx(0.102, abs=0.005)
    model, _ = tools._model("iJO1366")
    assert "HET_CRTE" not in model.reactions and "lycop_c" not in model.metabolites
    env = tools.production_envelope("vanillin", model="iJO1366", pathway="vanillin", points=3)
    assert env["points"][0]["flux_maximum"] > 0
    assert env["pathway"]["steps"]["HET_OMT"] == "added"


def test_knockout_searches_run_on_the_added_pathway(monkeypatch):
    # A real lycopene search on iJO1366 takes minutes (MOMA ~2 min per round), so stub the search
    # algorithms and check the wiring: the pathway is in the model and its genes are protected.
    seen = {}

    def fake_moma(model, target, substrate, exclude_genes, **kwargs):
        seen["moma"] = (target, "HET_CRTI" in model.reactions, set(exclude_genes))
        return [knockouts.MomaStep([], [], 0.9, 0.01, 0.001)]

    def fake_greedy(model, target, *args, **kwargs):
        seen["fba"] = (target, "HET_CRTI" in model.reactions)
        return knockouts.SearchResult(target=target, steps=[knockouts.KnockoutStep([], [], 0.9, 0.0)])

    monkeypatch.setattr(tools, "moma_knockout_search", fake_moma)
    monkeypatch.setattr(pipeline, "greedy_knockout_search", fake_greedy)
    out = tools.search_knockouts_moma("lycopene", model="iJO1366", pathway="lycopene", max_knockouts=1)
    assert seen["moma"] == ("DM_lycop_c", True, {"crtB", "crtE", "crtI"})
    assert out["pathway"]["product_metabolite"] == "lycop_c"
    out = tools.search_knockouts("lycopene", model="iJO1366", pathway="lycopene", max_knockouts=1)
    assert seen["fba"] == ("DM_lycop_c", True)
    assert out["pathway"]["genes"] == ["crtB", "crtE", "crtI"]


def test_moma_knockouts_raise_the_mutant_yield():
    out = tools.search_knockouts_moma("succinate", max_knockouts=2)
    json.dumps(out, allow_nan=False)
    summary = out["summary"]
    assert summary["improved"] and summary["best_yield"] > summary["parent_yield"]
    assert out["steps"][0]["knocked_out"] == []
    assert all(s["growth"] >= 0.3 * out["steps"][0]["growth"] - 1e-9 for s in out["steps"])


def test_resolve_compound_reports_how_it_matched():
    out = tools.resolve_compound("kegg:C00031", online=False)
    json.dumps(out, allow_nan=False)
    assert out["ok"] and out["model_metabolite"] == "glc__D"
    assert out["matched_via"] == "annotation kegg.compound"
    assert out["metabolites"][0]["exchange"] == "EX_glc__D_e"
    miss = tools.resolve_compound("glucoze", online=False)
    assert not miss["ok"] and "glc__D" in miss["suggestions"][0] and "not a metabolite" in miss["error"]


def test_flux_distribution_shows_the_fermentation_route():
    out = tools.flux_distribution("ethanol", aerobic=False, limit=10)
    json.dumps(out, allow_nan=False)
    assert out["product_flux"] == pytest.approx(20.0, abs=1e-3)  # 2 ethanol per glucose at zero growth
    assert {r["reaction"] for r in out["uptake"]} == {"EX_glc__D_e"}
    assert {"EX_etoh_e", "EX_co2_e"} <= {r["reaction"] for r in out["secretion"]}
    assert "ALCD2x" in {r["reaction"] for r in out["reactions"]}
    growth = tools.flux_distribution()
    assert growth["growth"] == pytest.approx(0.874, abs=1e-3) and growth["product_flux"] is None


def test_compare_conditions_runs_a_grid_and_keeps_failures_as_rows():
    out = tools.compare_conditions("ethanol", "glucose, glycerol", oxygen="both")
    json.dumps(out, allow_nan=False)
    assert len(out["rows"]) == 4
    glycerol = [r for r in out["rows"] if r["substrate"] == "glycerol"]
    assert all(r["error"] for r in glycerol)  # e_coli_core has no glycerol
    assert out["best"]["substrate"] == "glucose" and out["best"]["aerobic"] is False
    assert out["best"]["max_carbon_yield"] == pytest.approx(2 / 3, abs=1e-3)


def test_qa_model_returns_a_verdict():
    out = tools.qa_model()
    json.dumps(out, allow_nan=False)
    assert out["ok"] and out["engine"]
    assert out["hard_failures"]["total"] == 0


def test_check_thermodynamics_flags_the_route_not_the_whole_model():
    pytest.importorskip("equilibrator_api")
    growth = tools.check_thermodynamics()
    json.dumps(growth, allow_nan=False)
    assert growth["summary"]["route_feasible"]
    assert growth["summary"]["reactions_checked"] < 95  # only reactions carrying flux
    # Anaerobic succinate runs fumarate reductase (FRD7), whose ubiquinone ΔG'° is a known
    # group-contribution error (see tests/benchmarks/cases.py): flagged, and labelled as such.
    succ = tools.check_thermodynamics("succinate", aerobic=False, min_growth=0.05)
    rows = succ["against_thermodynamics"]["rows"]
    assert [r["reaction"] for r in rows] == ["FRD7"]
    assert rows[0]["method"] == "group contribution" and rows[0]["flux"] > 0
    assert succ["summary"]["against_thermodynamics_group_contribution"] == 1


def test_simulate_production_thermo_bounds_keeps_output_small():
    pytest.importorskip("equilibrator_api")
    out = tools.simulate_production("ethanol", aerobic=False, thermo_bounds=True)
    assert isinstance(out["checks"]["Reaction thermodynamics"], str)
    assert any("eQuilibrator" in c for c in out["caveats"])


@pytest.mark.parametrize("call, message", [
    (lambda: tools.simulate_production("unobtainium"), "not a metabolite"),
    (lambda: tools.simulate_production("ethanol", substrate="ethanol"), "same compound"),
    (lambda: tools.search_knockouts("ethanol", max_knockouts=0), "between 1 and 10"),
    (lambda: tools.describe_model("no/such/model.xml"), "Could not load"),
    (lambda: tools.simulate_production("taxol", pathway="taxol"), "Unknown pathway"),
    (lambda: tools.simulate_production("vanillin", pathway="vanillin"), "no 3-dehydroshikimate"),
    (lambda: tools.search_knockouts_moma("ethanol", max_knockouts=9), "between 1 and 5"),
    (lambda: tools.compare_conditions("ethanol", oxygen="sometimes"), "oxygen must be"),
    (lambda: tools.compare_conditions("ethanol", substrates=" , "), "between 1 and 8"),
    (lambda: tools.flux_distribution(limit=0), "between 1 and 200"),
    (lambda: tools.flux_distribution("ethanol", min_growth=5.0), "No feasible flux"),
    (lambda: tools.check_enzyme_kinetics(limit=0), "between 1 and 500"),
    (lambda: tools.design_knockouts("ethanol", method="optgene"), "Unknown method"),
    (lambda: tools.design_knockouts("ethanol", max_knockouts=0), "between 1 and 10"),
    (lambda: tools.design_knockouts("ethanol", max_designs=0), "between 1 and 20"),
    (lambda: tools.design_knockouts("ethanol", time_limit=0), "positive number of seconds"),
    (lambda: tools.design_knockouts("ethanol", substrate="ethanol"), "same compound"),
    (lambda: tools.check_driving_force("ethanol", substrate="ethanol"), "same compound"),
])
def test_bad_input_raises_a_helpful_error(call, message):
    with pytest.raises(tools.ToolInputError, match=message):
        call()


def test_every_public_function_is_exposed_or_excluded():
    """Fails when a new public cellfactory function has no MCP tool and no stated reason.

    Fix it by adding a tool to tools.TOOLS and an entry in tools.COVERAGE,
    or an entry in tools.NOT_EXPOSED explaining why agents don't need it.
    """
    names = {fn.__name__ for fn in tools.TOOLS}
    missing = set(cellfactory.__all__) - set(tools.COVERAGE) - set(tools.NOT_EXPOSED)
    assert not missing, f"Not exposed as MCP tools and not in NOT_EXPOSED: {sorted(missing)}"
    assert set(tools.COVERAGE.values()) <= names


def test_tools_take_only_json_friendly_arguments():
    for fn in tools.TOOLS:
        assert fn.__doc__, f"{fn.__name__} needs a docstring: it is the tool description agents see"
        for param in inspect.signature(fn).parameters.values():
            assert param.annotation in ("str", "int", "float", "bool"), (fn.__name__, param.name)


def test_mcp_server_lists_and_calls_tools():
    pytest.importorskip("mcp")
    from cellfactory.mcp_server import build_server

    async def go():
        server = build_server()
        listed = {t.name: t for t in await server.list_tools()}
        assert set(listed) == {fn.__name__ for fn in tools.TOOLS}
        assert "product" in listed["simulate_production"].input_schema["properties"]
        result = await server.call_tool("simulate_production", {"product": "ethanol", "aerobic": False})
        assert result.structured_content["product"] == "EX_etoh_e"

    asyncio.run(go())


def test_mcp_server_reports_bad_input_as_tool_error():
    pytest.importorskip("mcp")
    from cellfactory.mcp_server import ToolError, build_server

    with pytest.raises(ToolError, match="not a metabolite"):
        asyncio.run(build_server().call_tool("simulate_production", {"product": "unobtainium"}))


def test_check_driving_force_screens_the_whole_route():
    pytest.importorskip("equilibrator_api")
    out = tools.check_driving_force("ethanol", aerobic=False)
    json.dumps(out, allow_nan=False)
    assert out["route_for"] == "EX_etoh_e"
    # Anaerobic ethanol is how the cell makes a living, so the whole route must be able to run at once.
    assert out["verdict"] in ("feasible", "low driving force") and out["summary"]["route_feasible"]
    assert out["mdf"] is not None and out["mdf_optimistic"] >= out["mdf"]
    assert out["summary"]["estimated"] > 0 and out["bottlenecks"]
    assert out["steps"]["rows"][0]["reaction"] in out["bottlenecks"]  # least downhill step first
    assert out["concentrations"]["total"] > 0
    assert "EX_etoh_e" not in {r["reaction"] for r in out["steps"]["rows"]}  # exchanges aren't route steps


def test_check_enzyme_kinetics_reads_brenda_and_reports_the_source(monkeypatch):
    monkeypatch.setenv("CELLFACTORY_BRENDA", str(BRENDA_SAMPLE))
    out = tools.check_enzyme_kinetics("ethanol", aerobic=False, online=False)
    json.dumps(out, allow_nan=False)
    assert out["sources"] == ["BRENDA"] and out["organism"] == "Escherichia coli"
    rows = {r["reaction"]: r for r in out["reactions"]["rows"]}
    pfk = rows["PFK"]  # EC 2.7.1.11 is in the fixture, and e_coli_core borrows its ECs from iJO1366
    assert pfk["value_type"] == "measured" and pfk["organism_match"] == "host" and pfk["source"] == "BRENDA"
    assert pfk["kcat_per_s"] > 0 and pfk["flux"] > 0
    assert out["summary"]["measured"] >= 1 and out["summary"]["measured_in_host"] >= 1
    assert out["summary"]["enzyme_steps"] <= out["summary"]["reactions"]
    assert any("test-tube values" in c for c in out["caveats"])
    # Without a source there is nothing to look up, and that is said rather than guessed.
    monkeypatch.delenv("CELLFACTORY_BRENDA")
    blind = tools.check_enzyme_kinetics("ethanol", aerobic=False, online=False)
    assert blind["sources"] == [] and blind["summary"]["measured"] == 0
    assert any("no BRENDA file" in c for c in blind["caveats"])


def test_design_knockouts_finds_a_combination_and_checks_it_with_fba():
    pytest.importorskip("straindesign")
    out = tools.design_knockouts("succinate", aerobic=False, method="optknock", max_knockouts=3, time_limit=120)
    json.dumps(out, allow_nan=False)
    assert out["status"] == "optimal" and out["method"] == "OptKnock" and out["solver"]
    s = out["summary"]
    assert s["designs"] and s["improved"]
    assert s["best_guaranteed_yield"] > out["wild_type_guaranteed_yield"]
    best = out["designs"]["rows"][0]
    assert best["rank"] == 1 and best["n_knockouts"] >= 1 and best["genes"]
    assert best["growth"] >= out["min_growth"] and best["coupling"] in ("strong", "at max growth")
    model, _ = tools._model("e_coli_core")
    assert model.reactions.get_by_id(best["knockouts"].split(", ")[0]).bounds != (0.0, 0.0)  # model left alone

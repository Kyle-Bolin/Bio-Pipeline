"""The agent-facing tools (tools.py) and the MCP server that exposes them."""

import asyncio
import inspect
import json

import pytest

import cellfactory
from cellfactory import knockouts, pipeline, tools


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


@pytest.mark.parametrize("call, message", [
    (lambda: tools.simulate_production("unobtainium"), "not a metabolite"),
    (lambda: tools.simulate_production("ethanol", substrate="ethanol"), "same compound"),
    (lambda: tools.search_knockouts("ethanol", max_knockouts=0), "between 1 and 10"),
    (lambda: tools.describe_model("no/such/model.xml"), "Could not load"),
    (lambda: tools.simulate_production("taxol", pathway="taxol"), "Unknown pathway"),
    (lambda: tools.simulate_production("vanillin", pathway="vanillin"), "no 3-dehydroshikimate"),
    (lambda: tools.search_knockouts_moma("ethanol", max_knockouts=9), "between 1 and 5"),
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

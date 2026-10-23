"""Identifier resolution: database ids -> model metabolites. Offline; MetaNetX is always mocked."""

import io
import json
import urllib.error

import pytest

from cellfactory import analyze, load_model, parse_identifier, resolve_exchange, resolve_metabolite
from cellfactory import identifiers
from cellfactory.cli import main


@pytest.fixture(scope="module")
def _core():
    return load_model("e_coli_core")


@pytest.fixture
def core(_core):
    return _core.copy()


@pytest.fixture(autouse=True)
def isolated_cache(tmp_path, monkeypatch):
    monkeypatch.setenv("CELLFACTORY_CACHE", str(tmp_path))
    monkeypatch.delenv("CELLFACTORY_OFFLINE", raising=False)


def no_network(*args, **kwargs):
    raise AssertionError("MetaNetX must not be contacted here")


@pytest.mark.parametrize("query, expected", [
    ("kegg:C00031", ("kegg.compound", "C00031")),
    ("C00031", ("kegg.compound", "C00031")),
    ("CHEBI:4167", ("chebi", "4167")),
    ("chebi:4167", ("chebi", "4167")),
    ("MNXM41", ("metanetx.chemical", "MNXM41")),
    ("mnx:mnxm41", ("metanetx.chemical", "MNXM41")),
    ("HMDB0000122", ("hmdb", "122")),
    ("hmdb:HMDB00122", ("hmdb", "122")),  # old 5-digit and new 7-digit HMDB ids are the same compound
    ("metacyc:GLC", ("biocyc", "GLC")),
    ("seed:cpd00027", ("seed.compound", "cpd00027")),
    ("bigg:glc__D", ("bigg.metabolite", "glc__D")),
    ("WQZGKKKJIJFFOK-GASJEMHNSA-N", ("inchi_key", "WQZGKKKJIJFFOK-GASJEMHNSA-N")),
    ("glucose", None),
    ("glc__D", None),
    ("foo:bar", None),
])
def test_parse_identifier(query, expected):
    assert parse_identifier(query) == expected


@pytest.mark.parametrize("query", ["kegg:C00469", "C00469", "CHEBI:16236", "chebi:16236", "bigg:etoh",
                                   "seed:cpd00363", "hmdb:HMDB0000108"])
def test_database_ids_resolve_offline_through_annotations(core, query):
    res = resolve_metabolite(core, query, lookup=no_network)
    assert res.ok and res.model_id == "etoh" and res.via.startswith("annotation")
    assert res.caveat is None
    assert resolve_exchange(core, query) == "EX_etoh_e"


@pytest.mark.parametrize("query, via", [("glucose", "alias"), ("glc__D_e", "model id"), ("D-Glucose", "name")])
def test_plain_queries_still_resolve_first(core, query, via):
    res = resolve_metabolite(core, query, lookup=no_network)
    assert (res.model_id, res.via) == ("glc__D", via)


def test_shared_annotation_is_refused_not_guessed():
    # iJO1366 annotates the D-galactonate and D-gluconate forms of 2-dehydro-3-deoxy-... with the same ChEBI ids.
    model = load_model("iJO1366")
    res = resolve_metabolite(model, "CHEBI:57990", lookup=no_network)
    assert not res.ok
    assert res.candidates == ["2ddglcn", "2dh3dgal"]
    with pytest.raises(KeyError, match="matches several metabolites.*2ddglcn, 2dh3dgal"):
        resolve_exchange(model, "CHEBI:57990", create_demand=True)


def test_metanetx_fills_a_gap_in_the_annotations(core):
    calls = []

    def lookup(ns, local):
        calls.append((ns, local))
        return {"mnx": ["MNXM303"], "bigg": ["etoh"]}

    res = resolve_metabolite(core, "metacyc:ETOH-NOT-ANNOTATED", lookup=lookup)
    assert calls == [("biocyc", "ETOH-NOT-ANNOTATED")]
    assert res.ok and res.model_id == "etoh" and res.via == "MetaNetX MNXM303"
    assert "MetaNetX" in res.caveat


def test_metanetx_translation_absent_from_model(core):
    res = resolve_metabolite(core, "kegg:C99999", lookup=lambda ns, local: {"mnx": ["MNXM1"], "bigg": ["nope"]})
    assert not res.ok and "no metabolite in the model carries" in res.note


def test_offline_never_calls_metanetx(core):
    res = resolve_metabolite(core, "kegg:C99999", online=False, lookup=no_network)
    assert not res.ok and "lookup was off" in res.note


def test_unreachable_metanetx_degrades_with_a_message(core, monkeypatch):
    seen = {}

    def refuse(req, timeout):
        seen["timeout"] = timeout
        raise urllib.error.URLError("connection refused")

    monkeypatch.setattr(identifiers.urllib.request, "urlopen", refuse)
    with pytest.raises(KeyError, match="MetaNetX could not be reached"):
        resolve_exchange(core, "kegg:C99999")
    assert seen["timeout"] == identifiers.TIMEOUT  # a bounded wait, never a hang


def test_offline_env_var_skips_the_network(core, monkeypatch):
    monkeypatch.setenv("CELLFACTORY_OFFLINE", "1")
    monkeypatch.setattr(identifiers.urllib.request, "urlopen", no_network)
    res = resolve_metabolite(core, "kegg:C99999")
    assert not res.ok and "CELLFACTORY_OFFLINE" in res.note


class _Response(io.BytesIO):
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def test_metanetx_answer_is_parsed_and_cached(monkeypatch):
    answer = {"results": {"bindings": [
        {"chem": {"value": "https://rdf.metanetx.org/chem/MNXM303"},
         "bigg": {"value": "https://identifiers.org/bigg.metabolite:etoh"}},
    ]}}
    calls = []

    def fake(req, timeout):
        calls.append(req.full_url)
        return _Response(json.dumps(answer).encode())

    monkeypatch.setattr(identifiers.urllib.request, "urlopen", fake)
    assert identifiers.metanetx_lookup("kegg.compound", "C00469") == {"mnx": ["MNXM303"], "bigg": ["etoh"]}
    assert "C00469" in calls[0]
    monkeypatch.setattr(identifiers.urllib.request, "urlopen", no_network)
    assert identifiers.metanetx_lookup("kegg.compound", "C00469") == {"mnx": ["MNXM303"], "bigg": ["etoh"]}


def test_garbled_metanetx_answer_is_unavailable(monkeypatch):
    monkeypatch.setattr(identifiers.urllib.request, "urlopen", lambda req, timeout: _Response(b"<html>oops"))
    with pytest.raises(identifiers.MetaNetXUnavailable):
        identifiers.metanetx_lookup("kegg.compound", "C00469")


def test_near_miss_suggestions_catch_typos(core):
    res = resolve_metabolite(core, "acetaldehyd", lookup=no_network)
    assert not res.ok and any(s.startswith("acald") for s in res.suggestions)
    res = resolve_metabolite(core, "Ethanoll", lookup=no_network)
    assert any(s.startswith("etoh") for s in res.suggestions)


def test_analyze_records_how_ids_were_matched(core):
    report = analyze(core, "kegg:C00031", "CHEBI:16236", aerobic=False, qa_gate=False, online_ids=False)
    assert (report.substrate, report.product) == ("EX_glc__D_e", "EX_etoh_e")
    table = report.checks["Identifier resolution"].set_index("role")
    assert table.loc["substrate", "model_metabolite"] == "glc__D"
    assert table.loc["product", "matched_via"] == "annotation chebi"
    assert report.production.max_theoretical_yield == pytest.approx(2.0, abs=0.01)


def test_analyze_caveats_a_metanetx_match(core, monkeypatch):
    monkeypatch.setattr(identifiers, "metanetx_lookup", lambda ns, local: {"mnx": ["MNXM303"], "bigg": ["etoh"]})
    report = analyze(core, "glucose", "metacyc:SOMETHING-ELSE", aerobic=False, qa_gate=False)
    assert report.product == "EX_etoh_e"
    assert any("MetaNetX" in c for c in report.caveats)


def test_cli_offline_flag(capsys, monkeypatch):
    monkeypatch.setattr(identifiers.urllib.request, "urlopen", no_network)
    assert main(["run", "-p", "kegg:C99999", "--no-qa", "--offline"]) == 1
    assert "lookup was off" in capsys.readouterr().err
    assert main(["run", "-p", "kegg:C00469", "--anaerobic", "--no-qa", "--offline"]) == 0
    assert "annotation kegg.compound" in capsys.readouterr().out

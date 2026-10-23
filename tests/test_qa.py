"""The MEMOTE model QA gate (qa.py). Offline: MEMOTE runs locally, and is faked where it isn't installed."""

import json
import os
import subprocess
import sys
from pathlib import Path

import cobra
import pytest

from cellfactory import ModelQAError, analyze, load_model, memote_qa, qa
from cellfactory.cli import main as cli

FIXTURES = Path(__file__).parent / "data" / "memote"
NO_PYTHON = "/nonexistent/python"  # makes the gate fall back as if MEMOTE were missing


def fixture(name: str) -> dict:
    return json.loads((FIXTURES / f"{name}.json").read_text())


def free_atp(model: cobra.Model) -> cobra.Model:
    """A balanced but impossible reaction: ADP + Pi -> ATP with no energy source."""
    rxn = cobra.Reaction("FREE_ATP", lower_bound=0, upper_bound=1000)
    model.add_reactions([rxn])
    rxn.build_reaction_from_string("adp_c + pi_c + h_c --> atp_c + h2o_c")
    return model


@pytest.fixture(autouse=True)
def _private_cache(tmp_path, monkeypatch):
    monkeypatch.setenv("CELLFACTORY_CACHE", str(tmp_path / "cache"))
    monkeypatch.delenv("CELLFACTORY_MEMOTE_PYTHON", raising=False)


def _fake_memote(monkeypatch, raw: dict) -> list:
    calls = []

    def fake(model, tests, python, timeout):
        calls.append(tests)
        return raw
    monkeypatch.setattr(qa, "_run_memote", fake)
    return calls


MEMOTE_PYTHON = os.environ.get("CELLFACTORY_MEMOTE_PYTHON") or sys.executable  # read before the fixture clears it


def memote_installed() -> bool:
    try:
        return subprocess.run([MEMOTE_PYTHON, "-c", "import memote"], capture_output=True).returncode == 0
    except OSError:
        return False


# --- classifying MEMOTE's results (recorded runs of MEMOTE 0.17.0) ---

def test_core_model_passes():
    table = qa.classify(fixture("e_coli_core"))
    assert not (table["gate"] == "hard fail").any()
    assert "skipped" not in set(table["memote_result"])


def test_free_atp_reaction_is_a_hard_fail():
    table = qa.classify(fixture("e_coli_core_free_atp"))
    hard = table[table["gate"] == "hard fail"]
    assert hard["parameter"].str.contains("atp_c").any()
    assert hard["note"].str.contains("FREE_ATP").any()  # the note names the reactions in the cycle


def test_cycles_that_need_curator_closed_reactions_only_warn():
    """MEMOTE flags iJO1366 for ATP etc. after reopening e.g. superoxide dismutase; the model as shipped is fine."""
    raw = fixture("iJO1366")
    table = qa.classify(raw)
    egc = table[table["test"] == "test_detect_energy_generating_cycles"]
    flagged = egc[egc["memote_result"] == "failed"]
    assert len(flagged) == 6 and set(flagged["gate"]) == {"warning"}
    assert not (table["gate"] == "hard fail").any()


def test_unrechecked_cycle_and_inconsistency_are_hard_fails():
    raw = fixture("iJO1366")
    raw["egc_shipped_bounds"].pop("MNXM3")
    raw["tests"]["test_stoichiometric_consistency"]["result"] = "failed"
    hard = qa.classify(raw).query("gate == 'hard fail'")
    assert set(hard["test"]) == {"test_detect_energy_generating_cycles", "test_stoichiometric_consistency"}
    assert hard["note"].str.contains("not re-checked").any()


# --- memote_qa: caching and degrading ---

def test_results_are_cached_per_model_checksum(monkeypatch):
    calls = _fake_memote(monkeypatch, fixture("e_coli_core"))
    model = load_model("e_coli_core")
    first, second = memote_qa(model), memote_qa(model)
    assert len(calls) == 1 and calls[0] == qa.FAST_TESTS
    assert not first.cached and second.cached and second.ok
    model.reactions.PFK.upper_bound = 0
    assert memote_qa(model).checksum != first.checksum and len(calls) == 2
    memote_qa(model, use_cache=False)
    assert len(calls) == 3


def test_missing_memote_falls_back_to_builtin_checks():
    result = memote_qa(load_model("e_coli_core"), python=NO_PYTHON)
    assert result.ok and "MEMOTE unavailable" in result.engine
    assert "MEMOTE could not run" in result.caveats[0] and "stoichiometric consistency" in result.caveats[0]


def test_fallback_still_catches_free_atp():
    result = memote_qa(free_atp(load_model("e_coli_core")), python=NO_PYTHON)
    assert not result.ok
    with pytest.raises(ModelQAError, match="ATP made from nothing"):
        result.raise_if_failed()


def test_timeout_degrades_instead_of_hanging(monkeypatch):
    def slow(*args, **kwargs):
        raise subprocess.TimeoutExpired(args[0], kwargs.get("timeout"))
    monkeypatch.setattr(qa.subprocess, "run", slow)
    result = memote_qa(load_model("e_coli_core"), timeout=1)
    assert result.ok and "did not finish within 1 s" in result.caveats[0]


@pytest.mark.parametrize("damage", [
    lambda raw: raw["tests"].clear(),
    lambda raw: raw["tests"].pop("test_stoichiometric_consistency"),
    lambda raw: raw["tests"].pop("test_detect_energy_generating_cycles"),
    lambda raw: raw["tests"]["test_stoichiometric_consistency"].update(result="skipped"),
], ids=["no results", "no consistency test", "no cycle test", "consistency skipped"])
def test_incomplete_memote_results_fall_back_instead_of_passing(monkeypatch, damage):
    raw = fixture("e_coli_core")
    damage(raw)
    _fake_memote(monkeypatch, raw)
    result = memote_qa(load_model("e_coli_core"))
    assert "MEMOTE unavailable" in result.engine
    assert "MEMOTE could not run" in result.caveats[0]
    assert set(result.table["test"]) == {"free_energy_checks"}


def test_incomplete_cache_entry_is_rerun(monkeypatch):
    model = load_model("e_coli_core")
    cache = qa._default_cache()
    cache.mkdir(parents=True)
    (cache / f"{qa.model_checksum(model)}-fast-v{qa.RUNNER_VERSION}.json").write_text('{"tests": {}}')
    calls = _fake_memote(monkeypatch, fixture("e_coli_core"))
    result = memote_qa(model)
    assert len(calls) == 1 and not result.cached and result.engine.startswith("memote")


def test_fallback_warns_immediately(caplog):
    memote_qa(load_model("e_coli_core"), python=NO_PYTHON)
    assert "MEMOTE could not run" in caplog.text


# --- the gate in the pipeline and CLI ---

def test_analyze_stops_on_a_hard_failure(monkeypatch):
    _fake_memote(monkeypatch, fixture("e_coli_core_free_atp"))
    with pytest.raises(ModelQAError, match="atp_c made from nothing"):
        analyze(free_atp(load_model("e_coli_core")), "glucose", "ethanol")


def test_analyze_reports_qa_table_and_caveats(monkeypatch):
    raw = fixture("e_coli_core")
    raw["tests"]["test_gene_product_annotation_presence"]["result"] = "failed"
    _fake_memote(monkeypatch, raw)
    report = analyze("e_coli_core", "glucose", "ethanol", aerobic=False)
    assert any(name.startswith("Model QA gate (memote") for name in report.checks)
    assert "Presence of Gene Annotation" in report.caveats[0]
    assert "Caveats:" in report.to_text()
    assert report.to_dict()["caveats"] == report.caveats


def test_cli_qa_flag_exits_with_the_reason(tmp_path, capsys, monkeypatch):
    path = tmp_path / "broken.json"
    cobra.io.save_json_model(free_atp(load_model("e_coli_core")), str(path))
    monkeypatch.setenv("CELLFACTORY_MEMOTE_PYTHON", NO_PYTHON)
    assert cli(["run", "--model", str(path), "--product", "ethanol"]) == 1  # the gate is on by default
    assert "failed the model QA gate" in capsys.readouterr().err
    assert cli(["run", "--model", str(path), "--product", "ethanol", "--no-qa"]) == 0


# --- against real MEMOTE, where it is installed ---

def require_memote():
    """Skip locally without MEMOTE; in CI, where a step installs it, a missing MEMOTE is a failure."""
    if not memote_installed():
        message = f"MEMOTE not installed for {MEMOTE_PYTHON} (see CELLFACTORY_MEMOTE_PYTHON)"
        if os.environ.get("CI"):
            pytest.fail(message)
        pytest.skip(message)


def test_real_memote_passes_core_and_catches_free_atp():
    require_memote()
    good = memote_qa(load_model("e_coli_core"), python=MEMOTE_PYTHON)
    assert good.ok and good.engine.startswith("memote ")
    bad = memote_qa(free_atp(load_model("e_coli_core")), python=MEMOTE_PYTHON)
    assert not bad.ok
    assert bad.hard_failures["note"].str.contains("FREE_ATP").any()

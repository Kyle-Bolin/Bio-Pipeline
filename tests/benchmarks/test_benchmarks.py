"""Run every case in cases.py. ``pytest tests/benchmarks -v`` prints one line per known result."""

from __future__ import annotations

import os
import warnings
from functools import lru_cache
from pathlib import Path

import pytest
from cobra.flux_analysis import single_gene_deletion

import cobra

from cellfactory import load_model, memote_qa, production_at_max_growth, resolve_exchange, resolve_metabolite, set_medium

from .cases import CASES, ESSENTIALITY, IDENTIFIERS, KEIO, QA, Case, EssentialityCase, IdCase, QACase


# Folder of pre-downloaded SBML files (<id>.xml), checked before BiGG. CI fills
# and caches it so the BiGG cases run without hitting BiGG on every PR.
MODEL_DIR = os.environ.get("BENCHMARK_MODEL_DIR")


@lru_cache(maxsize=None)
def _cached(model_id: str):
    if MODEL_DIR and (Path(MODEL_DIR) / f"{model_id}.xml").exists():
        return load_model(Path(MODEL_DIR) / f"{model_id}.xml")
    try:
        return load_model(model_id)
    except RuntimeError as err:  # BiGG unreachable (offline, firewalled CI): skip rather than fail
        return err


def model_or_skip(model_id: str):
    model = _cached(model_id)
    if isinstance(model, Exception):
        pytest.skip(f"{model_id} not available offline: {model}")
    return model.copy()


def _param(case, case_id):
    marks = [pytest.mark.xfail(strict=True, reason=case.known_failure)] if case.known_failure else []
    return pytest.param(case, id=case_id, marks=marks)


def _measure(case: Case) -> float:
    model = model_or_skip(case.model)
    substrate = set_medium(model, case.substrate, 10.0, case.aerobic)
    for rxn_id in case.knockouts:
        model.reactions.get_by_id(rxn_id).knock_out()
    if case.metric == "growth":
        return model.slim_optimize(error_value=0.0)
    result = production_at_max_growth(model, resolve_exchange(model, case.product), substrate)
    return {"max_yield": result.max_theoretical_yield, "guaranteed_yield": result.guaranteed_yield}[case.metric]


@pytest.mark.parametrize("case", [_param(c, c.id) for c in CASES])
def test_known_result(case: Case):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Solver status is 'infeasible'")
        value = _measure(case)
    detail = f"{case.id}: got {value:.4f}, expected {case.expected} (source: {case.source})"
    if isinstance(case.expected, tuple):
        low, high = case.expected
        assert low <= value <= high, detail
    else:
        assert value == pytest.approx(case.expected, abs=case.tol), detail


@pytest.mark.parametrize("case", [_param(c, f"{c.model} {c.name}") for c in ESSENTIALITY])
def test_gene_essentiality(case: EssentialityCase):
    model = model_or_skip(case.model)
    set_medium(model, "glucose", 10.0, aerobic=True)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Solver status is 'infeasible'")
        growth = single_gene_deletion(model, gene_list=[case.gene])["growth"].fillna(0.0).iloc[0]
    predicted_essential = growth < 0.01
    assert predicted_essential == case.essential, (
        f"{case.name} ({case.gene}): model growth {growth:.3f}/h after deletion, "
        f"but it is {'essential' if case.essential else 'non-essential'} in the lab ({KEIO})"
    )


@pytest.mark.parametrize("case", [_param(c, c.id) for c in QA])
def test_model_qa(case: QACase, tmp_path):
    model = model_or_skip(case.model)
    for i, equation in enumerate(case.add_reactions):
        rxn = cobra.Reaction(f"PLANTED_DEFECT_{i}", lower_bound=0, upper_bound=1000)
        model.add_reactions([rxn])
        rxn.build_reaction_from_string(equation)
    result = memote_qa(model, cache_dir=tmp_path)
    if not result.engine.startswith("memote"):
        # CI installs MEMOTE, so there a missing one is a broken step, not a reason to go green.
        (pytest.fail if os.environ.get("CI") else pytest.skip)(f"MEMOTE not available: {result.caveats[0]}")
    if case.expect_hard_fail:
        assert case.expect_hard_fail in set(result.hard_failures["test"]), (
            f"{case.id}: expected {case.expect_hard_fail} to fail\n{result.table.to_string(index=False)}")
    if case.metric == "gate_ok":
        value = result.ok
    else:
        egc = result.table[result.table["test"] == "test_detect_energy_generating_cycles"]
        value = bool((egc["memote_result"] == "passed").all())
    detail = result.table[result.table["memote_result"] != "passed"].to_string(index=False)
    assert value == case.expected, f"{case.id}: got {value} (source: {case.source})\n{detail}"


@pytest.mark.parametrize("case", [_param(c, c.id) for c in IDENTIFIERS])
def test_identifier_resolution(case: IdCase):
    model = model_or_skip(case.model)
    res = resolve_metabolite(model, case.query, online=False)
    got = tuple(res.candidates) if res.candidates else res.model_id
    assert got == case.expected, f"{case.id}: {case.query} -> {got or res.note}, expected {case.expected} ({case.source})"

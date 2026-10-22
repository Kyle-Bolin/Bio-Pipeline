"""Model QA gate: vet a metabolic model with MEMOTE before simulating it.

MEMOTE (Lieven et al. 2020, Nat Biotechnol 38:272, doi:10.1038/s41587-020-0446-y)
is the community test suite for genome-scale models. ``memote_qa`` runs a fast
subset of it (or the whole suite) and sorts each result into:

* **hard fail**: the model can't be trusted for any yield. Stoichiometric
  inconsistency (some reaction creates or destroys mass) or an energy-generating
  cycle (ATP, NADH, ... made from nothing) that the model can actually run.
* **warning**: worth knowing, not a reason to stop: unbalanced reactions,
  metabolites without formulas, missing database annotations.

MEMOTE checks energy-generating cycles after resetting every reaction to its
reversibility, which reopens reactions the curators closed on purpose. iJO1366
"fails" that way for ATP, GTP and four more, and has no such cycle with its own
bounds. So every cycle MEMOTE finds is re-run with the model's shipped
directions; only one that survives is a hard fail.

MEMOTE pins older pandas and pytest than the rest of the pipeline, so it runs
in a subprocess, by default with this Python. Point ``CELLFACTORY_MEMOTE_PYTHON``
at a separate virtualenv's python to keep it out of this one. Results are
cached per model checksum (``CELLFACTORY_CACHE``). Without MEMOTE the gate falls back to the built-in
``free_energy_checks`` and says, in its caveats, what that leaves unchecked.
"""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

import cobra
import pandas as pd

from .checks import free_energy_checks

HARD = {"test_stoichiometric_consistency", "test_detect_energy_generating_cycles"}
WARN = {
    "test_reaction_mass_balance", "test_reaction_charge_balance",
    "test_metabolites_formula_presence", "test_metabolites_charge_presence",
    "test_metabolite_annotation_presence", "test_reaction_annotation_presence",
    "test_gene_product_annotation_presence",
}
FAST_TESTS = sorted(HARD | WARN)  # a few seconds on e_coli_core, ~20 s on iJO1366

RUNNER = Path(__file__).with_name("_memote_runner.py")
RUNNER_VERSION = 2  # bump when the runner's output changes, so old cache entries are ignored
COLUMNS = ["test", "parameter", "title", "memote_result", "metric", "gate", "note"]


class ModelQAError(ValueError):
    """The model failed a hard QA check; simulating it would give untrustworthy numbers."""


class MemoteUnavailable(RuntimeError):
    """MEMOTE could not be run: not installed, timed out, or crashed."""


@dataclass
class QAResult:
    model_id: str
    checksum: str
    engine: str  # "memote 0.17.0" or "built-in free_energy_checks (MEMOTE unavailable)"
    table: pd.DataFrame
    caveats: list[str] = field(default_factory=list)
    cached: bool = False

    @property
    def hard_failures(self) -> pd.DataFrame:
        return self.table[self.table["gate"] == "hard fail"]

    @property
    def ok(self) -> bool:
        return self.hard_failures.empty

    def raise_if_failed(self) -> None:
        if self.ok:
            return
        lines = [f"{r.test}{f' [{r.parameter}]' if r.parameter else ''}: {r.note}"
                 for r in self.hard_failures.itertuples()]
        raise ModelQAError(f"{self.model_id} failed the model QA gate ({self.engine}):\n  " + "\n  ".join(lines)
                           + "\nFix the model, or run without the QA gate if you know why this is acceptable.")


def _default_cache() -> Path:
    root = os.environ.get("CELLFACTORY_CACHE") or Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache"),
                                                       "cellfactory")
    return Path(root) / "memote"


def model_checksum(model: cobra.Model) -> str:
    """SHA-256 of the model's full JSON form: same reactions, bounds and annotations -> same checksum."""
    payload = json.dumps(cobra.io.model_to_dict(model), sort_keys=True, default=str)
    return hashlib.sha256(payload.encode()).hexdigest()


def _run_memote(model: cobra.Model, tests: list[str] | None, python: str, timeout: float) -> dict:
    with tempfile.TemporaryDirectory() as tmp:
        model_file, out_file = Path(tmp, "model.json"), Path(tmp, "result.json")
        cobra.io.save_json_model(model, str(model_file))
        cmd = [python, "-I", str(RUNNER), str(model_file), str(out_file)]
        for test in tests or []:
            cmd += ["--test", test]
        env = {**os.environ, "DISABLE_PANDERA_IMPORT_WARNING": "True"}
        try:
            proc = subprocess.run(cmd, cwd=tmp, env=env, capture_output=True, text=True, timeout=timeout)
        except FileNotFoundError as err:
            raise MemoteUnavailable(f"no Python at {python!r} ({err})") from err
        except subprocess.TimeoutExpired as err:
            raise MemoteUnavailable(f"MEMOTE did not finish within {timeout:.0f} s") from err
        if proc.returncode == 3:
            raise MemoteUnavailable(proc.stderr.strip() or "MEMOTE is not installed")
        if proc.returncode != 0 or not out_file.exists():
            tail = (proc.stderr or proc.stdout).strip().splitlines()[-3:]
            raise MemoteUnavailable(f"MEMOTE crashed (exit {proc.returncode}): {' / '.join(tail)}")
        return json.loads(out_file.read_text())


def _short(items, n: int = 5) -> str:
    items = list(items)
    return ", ".join(map(str, items[:n])) + (f" and {len(items) - n} more" if len(items) > n else "")


def classify(raw: dict) -> pd.DataFrame:
    """Turn the runner's JSON into one row per MEMOTE result, with this pipeline's verdict in ``gate``.

    Skipped results are dropped (MEMOTE skips e.g. a cycle test for a cofactor
    the model doesn't have). Failures of tests outside HARD and WARN, which only
    run with ``memote_qa(full=True)``, are reported as "info".
    """
    recheck, met_ids = raw.get("egc_shipped_bounds", {}), raw.get("egc_metabolites", {})
    rows = []
    for test, res in sorted(raw["tests"].items()):
        outcomes = res["result"] if isinstance(res["result"], dict) else {"": res["result"]}
        metrics = res.get("metric") if isinstance(res.get("metric"), dict) else {"": res.get("metric")}
        for param, outcome in outcomes.items():
            if outcome == "skipped":
                continue
            metric = metrics.get(param)
            gate, note = "pass", ""
            if outcome in ("failed", "error"):
                gate = "hard fail" if test in HARD else "warning" if test in WARN else "info"
                note = (res.get("message", {}).get(param) if param else res.get("message")) or ""
                note = " ".join(str(note).split("\n")[0].split())
            if test == "test_stoichiometric_consistency" and gate == "hard fail":
                note = "some reaction creates or destroys mass, so no set of positive molecular masses balances every reaction"
            if test == "test_detect_energy_generating_cycles" and outcome == "failed":
                met = met_ids.get(param, param)
                cycle = recheck.get(param)
                if cycle == []:
                    gate = "warning"
                    note = (f"{met} from nothing only after MEMOTE reopens reactions the model closes on purpose; "
                            "no such cycle with the model's own bounds")
                elif cycle:
                    note = f"{met} made from nothing with the model's own bounds, via {_short(cycle)}"
                else:
                    note = f"{met} made from nothing (not re-checked with the model's own bounds)"
            if test == "test_detect_energy_generating_cycles":
                param = f"{param} ({met_ids[param]})" if param in met_ids else param
            rows.append((test, param, res.get("title", test), outcome, metric, gate, note))
    return pd.DataFrame(rows, columns=COLUMNS)


def _fallback(model: cobra.Model, reason: str) -> tuple[pd.DataFrame, list[str]]:
    energy = free_energy_checks(model)
    rows = [("free_energy_checks", r.check, r.check, "passed" if r.ok else "failed", r.value,
             "pass" if r.ok else "hard fail", "" if r.ok else f"{r.value:.3g} with nothing taken up")
            for r in energy.itertuples()]
    caveats = [f"MEMOTE could not run ({reason}). Only the built-in free-energy checks were applied: "
               "stoichiometric consistency, cycles for cofactors other than ATP, and annotations are unchecked. "
               "Install it with `pip install cellfactory[checks]` or set CELLFACTORY_MEMOTE_PYTHON."]
    return pd.DataFrame(rows, columns=COLUMNS), caveats


def memote_qa(
    model: cobra.Model,
    full: bool = False,
    python: str | None = None,
    cache_dir: str | Path | None = None,
    use_cache: bool = True,
    timeout: float = 1800.0,
) -> QAResult:
    """Vet ``model`` with MEMOTE. Returns a QAResult; call ``raise_if_failed()`` to stop on hard failures.

    ``full=True`` runs the whole MEMOTE suite (minutes on a genome-scale model)
    instead of the fast subset that the gate decision rests on. ``python`` is
    the interpreter with MEMOTE installed (default: ``CELLFACTORY_MEMOTE_PYTHON``,
    else this one). Results are cached in ``cache_dir`` (default
    ``$CELLFACTORY_CACHE/memote`` or ``~/.cache/cellfactory/memote``) unless ``use_cache=False``.
    The model is not modified.
    """
    python = python or os.environ.get("CELLFACTORY_MEMOTE_PYTHON") or sys.executable
    checksum = model_checksum(model)
    tests = None if full else FAST_TESTS
    cache_file = None
    if use_cache:
        name = f"{checksum}-{'full' if full else 'fast'}-v{RUNNER_VERSION}.json"
        cache_file = Path(cache_dir or _default_cache()) / name
    raw, cached = None, False
    if cache_file is not None and cache_file.exists():
        try:
            raw, cached = json.loads(cache_file.read_text()), True
        except (OSError, ValueError):
            raw = None  # unreadable cache entry: run again
    if raw is None:
        try:
            raw = _run_memote(model, tests, python, timeout)
        except MemoteUnavailable as err:
            table, caveats = _fallback(model, str(err))
            return QAResult(model.id, checksum, "built-in free_energy_checks (MEMOTE unavailable)", table, caveats)
        if cache_file is not None:
            try:
                cache_file.parent.mkdir(parents=True, exist_ok=True)
                cache_file.write_text(json.dumps(raw))
            except OSError:
                pass  # a read-only cache shouldn't stop the analysis
    table = classify(raw)
    caveats = []
    warned = table[table["gate"] == "warning"]
    if not warned.empty:
        caveats.append(f"MEMOTE warnings on {model.id}: " + "; ".join(sorted(set(warned["title"]))) + ".")
    return QAResult(model.id, checksum, f"memote {raw.get('memote_version', '?')}", table, caveats, cached)

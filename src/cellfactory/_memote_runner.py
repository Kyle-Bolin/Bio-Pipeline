"""Run MEMOTE on one model file and write the results as JSON.

Standalone on purpose: ``qa.py`` runs this file with whichever Python has
MEMOTE installed (``CELLFACTORY_MEMOTE_PYTHON``), which may be a separate
virtualenv, since MEMOTE pins older pandas and pytest than the pipeline uses.
So it imports only MEMOTE, cobra and the standard library, never cellfactory.

Usage: python _memote_runner.py MODEL.json OUT.json [--test NAME ...] [--solver-timeout S]
Exit code 3 means MEMOTE is not installed in this interpreter.
"""

from __future__ import annotations

import argparse
import io
import json
import sys
from contextlib import redirect_stdout


def _keep_shipped_directions(model):
    """Like MEMOTE's ``close_boundaries_sensibly``, but a reaction the model ships closed stays closed.

    MEMOTE resets every reaction to (0, 1) or (-1, 1) from its reversibility,
    which reopens reactions curators switched off on purpose (iJO1366 closes
    catalase and superoxide dismutase, for example) and turns reverse-only
    reactions forward. That finds cycles the model as shipped cannot run.
    """
    for rxn in model.reactions:
        lower, upper = rxn.bounds
        rxn.bounds = (-1 if lower < 0 else 0, 1 if upper > 0 else 0)
    for rxn in model.boundary:
        rxn.bounds = (0, 0)


def _egc_with_shipped_bounds(model, met_id, consistency, helpers):
    original = helpers.close_boundaries_sensibly
    helpers.close_boundaries_sensibly = _keep_shipped_directions
    try:
        with model:
            return consistency.detect_energy_generating_cycles(model, met_id)
    finally:
        helpers.close_boundaries_sensibly = original


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("model")
    parser.add_argument("out")
    parser.add_argument("--test", action="append", default=[], help="run only these MEMOTE tests (repeatable)")
    parser.add_argument("--solver-timeout", type=int, default=60)
    args = parser.parse_args(argv)
    try:
        import memote
        from memote.suite import api
        from memote.support import consistency, helpers
        from memote.utils import jsonify
    except ImportError as err:
        print(f"MEMOTE is not installed in {sys.executable} ({err})", file=sys.stderr)
        return 3
    import cobra

    model = cobra.io.load_json_model(args.model)
    with redirect_stdout(io.StringIO()):  # MEMOTE's pytest run prints a full test log
        _, result = api.test_model(model, results=True, exclusive=args.test or None,
                                   pytest_args=["-q", "-p", "no:cacheprovider"],
                                   solver_timeout=args.solver_timeout)
    tests = json.loads(jsonify(result))["tests"]

    shipped = cobra.io.load_json_model(args.model)
    recheck, met_ids = {}, {}
    egc = tests.get("test_detect_energy_generating_cycles")
    if egc and isinstance(egc.get("result"), dict):
        main_comp = helpers.find_compartment_id_in_model(shipped, "c")
        for mnx_id, outcome in egc["result"].items():
            if outcome == "skipped":  # the model has no such cofactor
                continue
            met_ids[mnx_id] = helpers.find_met_in_model(shipped, mnx_id, main_comp)[0].id
            if outcome == "failed":
                recheck[mnx_id] = _egc_with_shipped_bounds(shipped, mnx_id, consistency, helpers)

    with open(args.out, "w", encoding="utf-8") as fh:
        json.dump({"memote_version": memote.__version__, "tests": tests,
                   "egc_shipped_bounds": recheck, "egc_metabolites": met_ids}, fh)
    return 0


if __name__ == "__main__":
    sys.exit(main())

"""Command line: cellfactory run --model iML1515 --substrate glucose --product ethanol"""

from __future__ import annotations

import argparse
import json
import sys

from .models import BUNDLED, find_metabolites, load_model
from .pipeline import analyze


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="cellfactory", description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    run = sub.add_parser("run", help="simulate substrate -> product in an organism")
    run.add_argument("--model", "-m", default="e_coli_core",
                     help=f"BiGG id or path to an SBML/JSON/MAT file. Offline: {', '.join(BUNDLED)}")
    run.add_argument("--substrate", "-s", default="glucose", help="carbon source: name, metabolite id or exchange id")
    run.add_argument("--product", "-p", required=True, help="target compound: name, metabolite id or exchange id")
    run.add_argument("--uptake", type=float, default=10.0, help="max substrate uptake, mmol/gDW/h (default 10)")
    run.add_argument("--anaerobic", action="store_true", help="remove oxygen from the medium")
    run.add_argument("--knockouts", "-k", type=int, default=0, help="run a greedy knockout search up to N deletions")
    run.add_argument("--min-growth", type=float, default=0.05, help="reject designs growing slower (1/h)")
    run.add_argument("--checks", action="store_true", help="also run model sanity checks")
    run.add_argument("--json", metavar="FILE", help="write the full report as JSON to FILE, or to stdout (instead of text) with -")

    find = sub.add_parser("find", help="search a model's metabolites by name or id")
    find.add_argument("query")
    find.add_argument("--model", "-m", default="e_coli_core")

    args = parser.parse_args(argv)
    try:
        model = load_model(args.model)
        if args.command == "find":
            for met in find_metabolites(model, args.query):
                exchange = next((r.id for r in met.reactions if r in model.exchanges), "")
                print(f"{met.id:20} {met.formula or '':14} {exchange:16} {met.name}")
            return 0
        report = analyze(model, args.substrate, args.product, uptake=args.uptake, aerobic=not args.anaerobic,
                         max_knockouts=args.knockouts, min_growth=args.min_growth, run_checks=args.checks)
    except (KeyError, ValueError, RuntimeError) as err:
        print(f"error: {err.args[0] if err.args else err}", file=sys.stderr)
        return 1
    if args.json == "-":
        json.dump(report.to_dict(), sys.stdout, indent=2, default=float)
        print()
        return 0
    print(report.to_text())
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(report.to_dict(), fh, indent=2, default=float)
    return 0


if __name__ == "__main__":
    sys.exit(main())

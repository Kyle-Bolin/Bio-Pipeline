"""Showcase: vanillin and lycopene from glucose, E. coli vs baker's yeast.

Reproduces every number and chart on the results page:

    python scripts/showcase.py                  # ~15 min, writes docs/showcase/
    python scripts/showcase.py --knockouts 0    # skip the knockout search (~1 min)

Models: E. coli K-12 iJO1366 (bundled with COBRApy) and the consensus yeast
model yeast-GEM, pinned to release v9.1.1 and downloaded from GitHub on first
run. Neither organism makes vanillin or lycopene, so the heterologous pathways
in ``cellfactory.heterologous`` are added first.
"""

from __future__ import annotations

import argparse
import datetime
import json
import urllib.request
from pathlib import Path

import cobra
from cobra.flux_analysis import production_envelope

from cellfactory import analyze, free_energy_checks, load_model, max_theoretical_yield, set_medium
from cellfactory.checks import ECOLI_ESSENTIAL
from cellfactory.heterologous import MET, PATHWAYS, add_pathway, find_or_create
from cellfactory.knockouts import moma_knockout_search

GLUCOSE_MW = 180.16
UPTAKE = 10.0  # mmol glucose / gDW / h, the same for both organisms so yields compare fairly
YEAST_GEM_TAG = "v9.1.1"
YEAST_GEM_URL = f"https://raw.githubusercontent.com/SysBioChalmers/yeast-GEM/{YEAST_GEM_TAG}/model/yeast-GEM.xml"
CACHE = Path.home() / ".cache" / "cellfactory"

ORGANISMS = {
    "ecoli": {"label": "E. coli (iJO1366)", "maintenance": "ATPM"},
    # yeast-GEM's non-growth-associated maintenance reaction
    "yeast": {"label": "S. cerevisiae (yeast-GEM v9.1.1)", "maintenance": "r_4046"},
}

# Precursor each pathway starts from, to explain differences between organisms.
PRECURSOR = {"vanillin": "3dhsk", "lycopene": "ipdp"}

# Plausibility band for biomass yield on glucose (g dry cells / g glucose) under fully aerobic,
# glucose-limited growth. Yeast: 0.50 g/g measured by Verduyn et al. 1991 (Antonie van Leeuwenhoek
# 59:49, doi:10.1007/BF00582119). E. coli: typically ~0.4-0.55 g/g. FBA assumes perfect efficiency,
# so it should land at or a little above the measured value.
BIOMASS_YIELD_BAND = (0.4, 0.65)


def load_yeast(path: str | None = None) -> cobra.Model:
    if path is None:
        CACHE.mkdir(parents=True, exist_ok=True)
        cached = CACHE / f"yeast-GEM-{YEAST_GEM_TAG}.xml"
        if not cached.exists():
            print(f"downloading yeast-GEM {YEAST_GEM_TAG} ...")
            urllib.request.urlretrieve(YEAST_GEM_URL, cached)
        path = str(cached)
    model = load_model(path)
    # yeast-GEM uses its own ids; give glucose and oxygen exchanges the BiGG ids cellfactory expects.
    model.reactions.r_1714.id = "EX_glc__D_e"
    model.reactions.r_1992.id = "EX_o2_e"
    model.repair()
    return model


# The bundled iJO1366 and yeast-GEM carry locus tags with no gene names, so name the genes
# the search actually reports. Names from EcoCyc (E. coli) and SGD (yeast).
GENE_NAMES = {
    "b1761": "gdhA", "b0452": "tesB", "b2241": "glpA", "b3236": "mdh", "b0871": "poxB", "b0529": "folD", "b3956": "ppc", "b2297": "pta", "b1380": "ldhA",
    "b0114": "aceE", "b0115": "aceF", "b2029": "gnd", "b1852": "zwf", "b4025": "pgi",
    "b2914": "rpiA", "b0720": "gltA", "b1136": "icd", "b2779": "eno", "b1702": "pykF",
    "YOR375C": "GDH1", "YAL062W": "GDH3", "YDL066W": "IDP1", "YNL037C": "IDH1",
    "YKL029C": "MAE1", "YKR080W": "MTD1", "YIL125W": "KGD1", "YLR044C": "PDC1",
    "YOR374W": "ALD4", "YPL061W": "ALD6", "YBR196C": "PGI1", "YHR183W": "ZWF1",
    "YOR374W-A": "ALD4", "YDR035W": "ARO3", "YBR249C": "ARO4",
}


def gene_label(model: cobra.Model, gene_id: str) -> str:
    name = model.genes.get_by_id(gene_id).name or GENE_NAMES.get(gene_id, "")
    return f"{name} ({gene_id})" if name and name != gene_id else gene_id


def add_sah_hydrolase(model: cobra.Model) -> None:
    """Yeast's SAH1: S-adenosyl-homocysteine + H2O -> adenosine + homocysteine.

    E. coli lacks it. It recycles SAH (left over after every SAM methylation, such as the
    last vanillin step) with Pfs + LuxS instead, which throws the ribose away as the
    signalling molecule AI-2: five carbons lost per vanillin.
    """
    rxn = cobra.Reaction("HET_SAHH", name="S-adenosyl-L-homocysteine hydrolase (yeast SAH1)", upper_bound=1000.0)
    m = model.metabolites
    rxn.add_metabolites({m.ahcys_c: -1, m.h2o_c: -1, m.adn_c: 1, m.hcys__L_c: 1})
    model.add_reactions([rxn])
    rxn.gene_reaction_rule = "SAH1"


# Engineering ideas the model can test: (organism, product) -> (label, model edit)
WHAT_IFS = {("ecoli", "vanillin"): ("add yeast SAH1 (SAH hydrolase)", add_sah_hydrolase)}


def run_case(base: cobra.Model, organism: str, product: str, knockouts: int) -> dict:
    model = base.copy()
    pathway = PATHWAYS[product]
    steps = add_pathway(model, pathway)
    product_met = steps.pop("product")
    report = analyze(model, "EX_glc__D_e", product_met, uptake=UPTAKE, aerobic=True)
    p = report.production
    target = report.product
    n_carbon = model.metabolites.get_by_id(product_met).elements["C"]

    env = production_envelope(model, reactions=[model.reactions.get_by_id(_biomass(model))],
                              objective=target, points=15)
    envelope = [{"growth": float(r[_biomass(model)]), "max_flux": float(r["flux_maximum"])}
                for _, r in env.iterrows()]

    # Precursor supply: the most precursor the cell could make per glucose.
    with model:
        pre = find_or_create(model, PRECURSOR[product])[0]
        dm = model.add_boundary(pre, type="demand")
        precursor_yield = max_theoretical_yield(model, dm.id, "EX_glc__D_e")

    what_if = None
    if (organism, product) in WHAT_IFS:
        label, edit = WHAT_IFS[organism, product]
        with model:
            edit(model)
            what_if = {"label": label, "max_yield_mol": max_theoretical_yield(model, target, "EX_glc__D_e")}

    ko = []
    if knockouts:
        heterologous = {g.id for s in pathway.steps for g in model.genes if g.id in s.genes.split()}
        for s in moma_knockout_search(model, target, max_knockouts=knockouts, exclude_genes=heterologous):
            ko.append({"genes": [gene_label(model, g) for g in s.knocked_out],
                       "reactions": [f"{r} ({model.reactions.get_by_id(r).name})" for r in s.reactions],
                       "growth": s.growth, "product_yield": s.product_yield})

    added = [model.reactions.get_by_id(r) for r, how in steps.items() if how == "added"]
    return {
        "organism": organism, "product": product,
        "pathway_steps": steps,
        "pathway_mass_balanced": all(not r.check_mass_balance() for r in added),
        "growth": p.growth,
        "max_yield_mol": p.max_theoretical_yield,
        "max_yield_g": p.max_theoretical_yield * pathway.molar_mass / GLUCOSE_MW,
        "max_yield_carbon": p.max_theoretical_yield * n_carbon / 6,
        "guaranteed_yield_at_max_growth": p.guaranteed_yield,
        "precursor": MET[PRECURSOR[product]].name,
        "precursor_max_yield_mol": precursor_yield,
        "envelope": envelope, "uptake": UPTAKE,
        "what_if": what_if,
        "knockouts": ko,
    }


def flag_known_failures(results: dict) -> dict:
    """Mark predicted deletions that the lab data says are wrong.

    iJO1366 lets the cell live without PEP carboxylase (ppc) through a bypass real
    E. coli can't use, so any design resting on that deletion is a model artefact.
    """
    for case in results["cases"]:
        for step in case["knockouts"]:
            step["known_failures"] = [] if case["organism"] != "ecoli" else [
                ECOLI_ESSENTIAL[b] for g in step["genes"]
                if (b := g.split("(")[-1].strip(")")) in ECOLI_ESSENTIAL
            ]
    return results


def sanity(model: cobra.Model, organism: str) -> dict:
    set_medium(model, "EX_glc__D_e", UPTAKE, aerobic=True)
    growth = model.slim_optimize()
    biomass_yield = growth / (UPTAKE * GLUCOSE_MW / 1000)
    fe = free_energy_checks(model, maintenance=ORGANISMS[organism]["maintenance"])
    return {
        "growth": growth,
        "biomass_yield_g_per_g": biomass_yield,
        "biomass_yield_plausible": BIOMASS_YIELD_BAND[0] <= biomass_yield <= BIOMASS_YIELD_BAND[1],
        "free_energy": fe.to_dict(orient="records"),
    }


def _biomass(model: cobra.Model) -> str:
    return next(r.id for r in model.reactions if r.objective_coefficient)


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default="docs/showcase", help="output folder")
    ap.add_argument("--knockouts", type=int, default=2, help="max genes to delete in the MOMA search (0 = skip)")
    ap.add_argument("--yeast-model", help="path to yeast-GEM.xml (default: download the pinned release)")
    args = ap.parse_args(argv)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    models = {"ecoli": load_model("iJO1366"), "yeast": load_yeast(args.yeast_model)}
    results = {"organisms": {k: v["label"] for k, v in ORGANISMS.items()},
               "settings": {"glucose_uptake": UPTAKE, "aerobic": True, "max_knockouts": args.knockouts,
                            "yeast_gem": YEAST_GEM_TAG, "cobra": cobra.__version__},
               "generated": datetime.date.today().isoformat(), "sanity": {}, "cases": []}
    for org, model in models.items():
        results["sanity"][org] = sanity(model.copy(), org)
        for product in PATHWAYS:
            print(f"{org} -> {product}")
            case = run_case(model, org, product, args.knockouts)
            results["cases"].append(case)
            print(f"  max yield {case['max_yield_mol']:.3f} mol/mol, knockouts: {case['knockouts'][-1:] }")
    flag_known_failures(results)
    (out / "results.json").write_text(json.dumps(results, indent=2, default=float))
    from showcase_charts import make_charts  # noqa: E402 (sibling module; keeps plotting out of the science)
    make_charts(results, out)
    print(f"wrote {out}/ (open index.html for the results page)")


if __name__ == "__main__":
    main()

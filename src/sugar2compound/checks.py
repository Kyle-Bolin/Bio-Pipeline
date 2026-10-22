"""Sanity checks that catch common model quirks before you trust a prediction.

Models are hand-curated and imperfect. Before believing "knock out X to get 2x
yield", check that the model behaves like a real cell in cases where we already
know the answer.
"""

from __future__ import annotations

import cobra
import pandas as pd
from cobra.flux_analysis import find_essential_genes, single_gene_deletion

from .models import _carbon_uptakes

# E. coli genes with well-established behaviour on glucose minimal medium, from
# the Keio single-gene knockout collection (Baba et al. 2006, Mol Syst Biol) and
# growth screens on minimal media (Joyce et al. 2006, J Bacteriol).
# Keys are E. coli "b-numbers", the locus tags models use as gene ids.
# For another organism, pass your own lists (e.g. from a published knockout screen).
ECOLI_ESSENTIAL = {
    "b1779": "gapA: glyceraldehyde-3-phosphate dehydrogenase (glycolysis)",
    "b2926": "pgk: phosphoglycerate kinase (glycolysis)",
    "b2779": "eno: enolase (glycolysis)",
    "b0720": "gltA: citrate synthase (no glutamate without it)",
    "b3956": "ppc: PEP carboxylase (refills the TCA cycle on glucose)",
    "b0004": "thrC: threonine synthase (makes the amino acid threonine)",
    "b3960": "argH: argininosuccinate lyase (makes arginine)",
    "b4177": "purA: adenylosuccinate synthase (makes AMP for DNA/RNA)",
    "b1281": "pyrF: OMP decarboxylase (makes pyrimidines for DNA/RNA)",
}
ECOLI_NONESSENTIAL = {
    "b4025": "pgi: glucose-6-phosphate isomerase (pentose phosphate route bypasses it)",
    "b3916": "pfkA: phosphofructokinase (isozyme pfkB covers)",
    "b1852": "zwf: glucose-6-phosphate dehydrogenase",
    "b1380": "ldhA: lactate dehydrogenase (fermentation side product)",
    "b2297": "pta: phosphotransacetylase (acetate production)",
    "b2296": "ackA: acetate kinase (acetate production)",
    "b1241": "adhE: alcohol dehydrogenase (ethanol production)",
    "b1101": "ptsG: glucose transporter (other transporters take over)",
    "b0723": "sdhA: succinate dehydrogenase subunit",
}


def gene_essentiality_check(
    model: cobra.Model,
    known_essential: dict[str, str] | None = None,
    known_nonessential: dict[str, str] | None = None,
    threshold: float = 0.01,
) -> pd.DataFrame:
    """Compare predicted vs known essentiality for genes with lab data.

    The known lists map gene id -> description and default to the E. coli
    lists above. A gene is predicted essential if deleting it drops growth below
    ``threshold`` (1/h). Genes absent from the model (e.g. the 137-gene core
    model lacks amino acid synthesis) are reported as "not in model".
    """
    if known_essential is None and known_nonessential is None:
        known_essential, known_nonessential = ECOLI_ESSENTIAL, ECOLI_NONESSENTIAL
    known = {**{g: (True, d) for g, d in (known_essential or {}).items()},
             **{g: (False, d) for g, d in (known_nonessential or {}).items()}}
    in_model = [g for g in known if g in model.genes]
    growth = {}
    if in_model:
        deletions = single_gene_deletion(model, gene_list=in_model)
        for ids, value in zip(deletions["ids"], deletions["growth"]):
            (gene,) = ids
            growth[gene] = 0.0 if pd.isna(value) else value
    rows = []
    for gene, (essential, description) in known.items():
        if gene not in growth:
            rows.append((gene, description, essential, None, None, "not in model"))
            continue
        predicted = growth[gene] < threshold
        rows.append((gene, description, essential, growth[gene], predicted,
                     "agree" if predicted == essential else "MISMATCH"))
    return pd.DataFrame(rows, columns=["gene", "description", "known_essential",
                                       "predicted_growth", "predicted_essential", "verdict"])


def predicted_essential_genes(model: cobra.Model, threshold: float = 0.01) -> pd.DataFrame:
    """Every gene the model says the cell can't live without, in the current medium.

    Works for any organism. Without lab data to compare against, use it as a
    reading list: a surprising entry (or a surprising absence) is a lead.
    """
    genes = find_essential_genes(model, threshold=threshold)
    return pd.DataFrame(sorted(((g.id, g.name, ", ".join(r.id for r in g.reactions)) for g in genes)),
                        columns=["gene", "name", "reactions"])


def free_energy_checks(model: cobra.Model, maintenance: str = "ATPM") -> pd.DataFrame:
    """Physics checks: no growth or ATP should come from nothing.

    1. With every uptake closed, the cell must not grow.
    2. With every uptake closed, the ATP maintenance reaction (``maintenance``,
       usually ATPM in BiGG models) must be 0.
       A model that makes ATP from nothing has a thermodynamically impossible
       loop, and it will over-predict every yield.
    3. With every carbon source closed but salts open, there is no growth.
    """
    rows = []
    with model:
        for ex in model.exchanges:
            ex.lower_bound = 0.0
        rows.append(("growth with no uptake at all", model.slim_optimize(error_value=0.0)))
        if maintenance in model.reactions:
            atpm = model.reactions.get_by_id(maintenance)
            atpm.lower_bound = 0.0
            model.objective = atpm
            rows.append(("ATP made from nothing", model.slim_optimize(error_value=0.0)))
    with model:
        medium = model.medium
        for rxn_id in _carbon_uptakes(model, min_uptake=0.0):
            medium.pop(rxn_id)
        model.medium = medium
        rows.append(("growth with no carbon source", model.slim_optimize(error_value=0.0)))
    df = pd.DataFrame(rows, columns=["check", "value"])
    df["ok"] = df["value"].abs() < 1e-6
    return df


def mass_balance_problems(model: cobra.Model) -> pd.DataFrame:
    """Non-boundary reactions whose atoms or charge don't balance.

    Some imbalances are deliberate (biomass, or pseudo-reactions); others are
    curation errors that let a model create atoms out of thin air.
    """
    rows = []
    for rxn in model.reactions:
        if rxn.boundary or rxn.objective_coefficient:
            continue
        imbalance = rxn.check_mass_balance()
        if imbalance:
            rows.append((rxn.id, rxn.name, dict(imbalance)))
    return pd.DataFrame(rows, columns=["reaction", "name", "imbalance"])

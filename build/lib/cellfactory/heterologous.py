"""Add non-native ("heterologous") pathways to a model.

E. coli and yeast make neither vanillin nor lycopene. To produce them, metabolic
engineers insert genes from other organisms (plants, fungi, other bacteria) that
encode the missing enzymes. In a model that means adding the missing reactions,
attached to the organism's own precursor metabolites.

Models name metabolites differently: BiGG models (iJO1366, iML1515) use ids like
``frdp_c``, while yeast-GEM uses ``s_0190`` with the name "farnesyl diphosphate".
Each metabolite here therefore carries both a BiGG id and a name, and we look up
whichever the model uses. Metabolites the model lacks (the new intermediates and
the product) are created, with formula and charge so mass balance can be checked.
"""

from __future__ import annotations

from dataclasses import dataclass

import cobra


@dataclass(frozen=True)
class Met:
    bigg: str
    name: str
    formula: str
    charge: int


@dataclass(frozen=True)
class Step:
    id: str
    name: str
    genes: str  # gene_reaction_rule; source organism in the pathway description
    stoich: dict[str, float]  # metabolite key -> coefficient (negative = consumed)


@dataclass(frozen=True)
class Pathway:
    product: str  # metabolite key of the final product
    description: str
    steps: tuple[Step, ...]
    molar_mass: float  # g/mol, to report g product per g substrate


MET = {m.bigg: m for m in [
    Met("h", "H+", "H", 1),
    Met("h2o", "H2O", "H2O", 0),
    Met("o2", "oxygen", "O2", 0),
    Met("atp", "ATP", "C10H12N5O13P3", -4),
    Met("amp", "AMP", "C10H12N5O7P", -2),
    Met("ppi", "diphosphate", "HO7P2", -3),
    Met("nadph", "NADPH", "C21H26N7O17P3", -4),
    Met("nadp", "NADP(+)", "C21H25N7O17P3", -3),
    Met("amet", "S-adenosyl-L-methionine", "C15H23N6O5S", 1),
    Met("ahcys", "S-adenosyl-L-homocysteine", "C14H20N6O5S", 0),
    # Vanillin route (from the shikimate pathway, which makes aromatic amino acids)
    Met("3dhsk", "3-dehydroshikimate", "C7H7O5", -1),
    Met("34dhbz", "protocatechuate", "C7H5O4", -1),
    Met("34dhbald", "protocatechualdehyde", "C7H6O3", 0),
    Met("vanln", "vanillin", "C8H8O3", 0),
    # Lycopene route (from the isoprenoid precursors IPP and FPP)
    Met("ipdp", "isopentenyl diphosphate", "C5H9O7P2", -3),
    Met("frdp", "farnesyl diphosphate", "C15H25O7P2", -3),
    Met("ggdp", "geranylgeranyl diphosphate", "C20H33O7P2", -3),
    Met("phyto", "phytoene", "C40H64", 0),
    Met("lycop", "lycopene", "C40H56", 0),
]}

VANILLIN = Pathway(
    product="vanln",
    description=(
        "De novo vanillin from glucose (Hansen et al. 2009, Appl Environ Microbiol 75:2765): "
        "3-dehydroshikimate dehydratase (aroZ, Podospora pauciseta), aromatic carboxylic acid reductase "
        "(car, Nocardia iowensis, activated by a phosphopantetheinyl transferase) and catechol "
        "O-methyltransferase (hsOMT, human)."
    ),
    steps=(
        Step("HET_3DSD", "3-dehydroshikimate dehydratase", "aroZ", {"3dhsk": -1, "34dhbz": 1, "h2o": 1}),
        Step("HET_ACAR", "aromatic carboxylic acid reductase", "car",
             {"34dhbz": -1, "atp": -1, "nadph": -1, "h": -1, "34dhbald": 1, "amp": 1, "ppi": 1, "nadp": 1}),
        Step("HET_OMT", "catechol O-methyltransferase", "hsOMT",
             {"34dhbald": -1, "amet": -1, "vanln": 1, "ahcys": 1, "h": 1}),
    ),
    molar_mass=152.15,
)

LYCOPENE = Pathway(
    product="lycop",
    description=(
        "Lycopene from the cell's own isoprenoid precursors via the Pantoea crtEBI genes (Misawa et al. "
        "1990, J Bacteriol 172:6704): GGPP synthase (crtE), phytoene synthase (crtB) and phytoene "
        "desaturase (crtI). crtI is modelled with O2 as the final electron acceptor (4 desaturations)."
    ),
    steps=(
        Step("HET_CRTE", "geranylgeranyl diphosphate synthase", "crtE", {"frdp": -1, "ipdp": -1, "ggdp": 1, "ppi": 1}),
        Step("HET_CRTB", "phytoene synthase", "crtB", {"ggdp": -2, "phyto": 1, "ppi": 2}),
        Step("HET_CRTI", "phytoene desaturase", "crtI", {"phyto": -1, "o2": -2, "lycop": 1, "h2o": 4}),
    ),
    molar_mass=536.87,
)

PATHWAYS = {"vanillin": VANILLIN, "lycopene": LYCOPENE}


def find_or_create(model: cobra.Model, key: str, compartment: str = "c") -> tuple[cobra.Metabolite, bool]:
    """The model's metabolite for ``key`` (by BiGG id, then by name), creating it if absent.

    Returns (metabolite, created).
    """
    spec = MET[key]
    bigg_id = f"{spec.bigg}_{compartment}"
    if bigg_id in model.metabolites:
        return model.metabolites.get_by_id(bigg_id), False
    hits = [m for m in model.metabolites
            if m.compartment == compartment and (m.name or "").lower() == spec.name.lower()]
    if hits:
        return hits[0], False
    met = cobra.Metabolite(bigg_id, formula=spec.formula, name=spec.name, compartment=compartment, charge=spec.charge)
    model.add_metabolites([met])
    return met, True


def _native(model: cobra.Model, step: Step, mets: dict[str, cobra.Metabolite]) -> cobra.Reaction | None:
    """A model reaction doing the same chemistry as ``step`` (same substrates and products), if any."""
    consumed = {mets[k] for k, c in step.stoich.items() if c < 0 and k not in ("h", "h2o")}
    produced = {mets[k] for k, c in step.stoich.items() if c > 0 and k not in ("h", "h2o")}
    for rxn in next(iter(produced)).reactions:
        if consumed <= set(rxn.reactants) and produced <= set(rxn.products):
            return rxn
    return None


def add_pathway(model: cobra.Model, pathway: Pathway | str, compartment: str = "c") -> dict[str, str]:
    """Insert ``pathway`` into ``model`` and return what happened to each step.

    A step the organism can already do natively (e.g. yeast's own GGPP synthase,
    BTS1) is skipped. Returns {step id: "added" | "native (<reaction id>)"} and
    the product's metabolite id under the key "product".
    """
    if isinstance(pathway, str):
        pathway = PATHWAYS[pathway]
    outcome: dict[str, str] = {}
    for step in pathway.steps:
        mets = {k: find_or_create(model, k, compartment)[0] for k in step.stoich}
        native = _native(model, step, mets)
        if native is not None:
            outcome[step.id] = f"native ({native.id})"
            continue
        rxn = cobra.Reaction(step.id, name=step.name, lower_bound=0.0, upper_bound=1000.0)
        rxn.add_metabolites({mets[k]: c for k, c in step.stoich.items()})
        model.add_reactions([rxn])
        rxn.gene_reaction_rule = step.genes
        outcome[step.id] = "added"
    outcome["product"] = find_or_create(model, pathway.product, compartment)[0].id
    return outcome

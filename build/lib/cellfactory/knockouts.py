"""A simple greedy knockout search.

Idea: deleting a gene removes the reaction(s) its enzyme catalyses. If we delete
the routes that compete with our product (e.g. lactate and acetate), the cell
may be forced to send more carbon to the product just to keep growing.

Real tools (OptKnock, OptGene, cameo, StrainDesign) solve this more cleverly.
Greedy search is the "brute force with a heuristic" version: each round, try
every candidate knockout, keep the one that helps most, repeat.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field

import cobra
from cobra.flux_analysis import pfba

from .simulate import DEFAULT_SUBSTRATE, production_at_max_growth

# Reactions we never knock out: exchanges/demands/sinks are model plumbing, not enzymes.
_BOUNDARY_PREFIXES = ("EX_", "DM_", "SK_")


@dataclass
class KnockoutStep:
    knocked_out: list[str]
    genes: list[str]
    growth: float
    guaranteed_yield: float


@dataclass
class SearchResult:
    target: str
    steps: list[KnockoutStep] = field(default_factory=list)

    def to_frame(self):
        import pandas as pd

        return pd.DataFrame([s.__dict__ for s in self.steps])


def candidate_reactions(model: cobra.Model, active_only: bool = True, exclude_transport: bool = True) -> list[str]:
    """Reactions worth trying to knock out.

    Only reactions with genes (otherwise there is nothing to delete in the lab),
    excluding boundary reactions and the biomass objective. With ``active_only``,
    keep just reactions carrying flux in a parsimonious FBA solution. That is a
    heuristic to make genome-scale searches fast; it can miss some designs.

    ``exclude_transport`` skips reactions that move a molecule between
    compartments (e.g. across a membrane). Models often give transporters
    vague or catch-all gene assignments, so "knock out the hydrogen porin"
    looks great on paper and is meaningless in the lab.
    """
    active = None
    if active_only:
        fluxes = pfba(model).fluxes
        active = set(fluxes[fluxes.abs() > 1e-6].index)
    out = []
    for rxn in model.reactions:
        if not rxn.genes or rxn.id.startswith(_BOUNDARY_PREFIXES) or rxn.objective_coefficient:
            continue
        if active is not None and rxn.id not in active:
            continue
        if exclude_transport and len(rxn.compartments) > 1:
            continue
        out.append(rxn.id)
    return out


def greedy_knockout_search(
    model: cobra.Model,
    target: str,
    substrate: str = DEFAULT_SUBSTRATE,
    max_knockouts: int = 3,
    min_growth: float = 0.05,
    candidates: list[str] | None = None,
) -> SearchResult:
    """Greedily pick up to ``max_knockouts`` reaction deletions that raise guaranteed yield.

    ``min_growth`` (1/h) rejects designs that barely grow; a strain that grows at
    0.001/h is useless in a fermenter even if its yield is perfect.
    """
    result = SearchResult(target)
    knocked: list[str] = []
    with model, warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Solver status is 'infeasible'")  # lethal knockouts
        base = production_at_max_growth(model, target, substrate)
        result.steps.append(KnockoutStep([], [], base.growth, base.guaranteed_yield))
        pool = list(candidates) if candidates is not None else candidate_reactions(model)
        best_yield = base.guaranteed_yield

        for _ in range(max_knockouts):
            best = None
            for rxn_id in pool:
                if rxn_id in knocked:
                    continue
                with model:
                    model.reactions.get_by_id(rxn_id).knock_out()
                    res = production_at_max_growth(model, target, substrate)
                if res.growth >= min_growth and res.guaranteed_yield > best_yield + 1e-6:
                    best_yield = res.guaranteed_yield
                    best = (rxn_id, res)
            if best is None:
                break  # no single extra knockout helps any more
            rxn_id, res = best
            model.reactions.get_by_id(rxn_id).knock_out()
            knocked.append(rxn_id)
            genes = sorted({g.id for r in knocked for g in model.reactions.get_by_id(r).genes})
            result.steps.append(KnockoutStep(list(knocked), genes, res.growth, res.guaranteed_yield))
    return result


@dataclass
class MomaStep:
    knocked_out: list[str]  # gene ids
    reactions: list[str]  # reactions those genes disable
    growth: float
    product_flux: float  # mmol/gDW/h
    product_yield: float  # mol product / mol substrate


def _disabled_reactions(model: cobra.Model, gene_id: str) -> list[str]:
    with model:
        model.genes.get_by_id(gene_id).knock_out()
        return [r.id for r in model.genes.get_by_id(gene_id).reactions if r.bounds == (0.0, 0.0)]


def moma_knockout_search(
    model: cobra.Model,
    target: str,
    substrate: str = DEFAULT_SUBSTRATE,
    reference_fraction: float = 0.01,
    max_knockouts: int = 2,
    min_growth_fraction: float = 0.3,
    exclude_genes: set[str] | None = None,
) -> list[MomaStep]:
    """Greedy gene-knockout search scored with linear MOMA instead of FBA.

    For products the cell has no reason to make (most non-native compounds), no
    single deletion forces production, so the FBA search above finds nothing,
    and in FBA a deletion can never raise the *maximum* yield either. MOMA
    (minimisation of metabolic adjustment, Segrè, Vitkup & Church 2002) instead
    assumes a freshly made mutant keeps its fluxes as close as possible to its
    parent's. Starting from a parent that sends ``reference_fraction`` of the
    maximum possible product flux to the product, it predicts how a deletion
    reroutes carbon, which is how Alper et al. 2005 picked lycopene targets.

    Candidates are genes of reactions active in the parent (no transporters).
    A mutant must keep ``min_growth_fraction`` of the parent's growth rate.
    """
    from cobra.flux_analysis.moma import add_moma

    exclude_genes = exclude_genes or set()
    product = model.reactions.get_by_id(target)
    biomass = next(r.id for r in model.reactions if r.objective_coefficient)
    with model, warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Solver status is 'infeasible'")
        with model:
            model.objective = target
            max_flux = model.slim_optimize()
        with model:
            product.lower_bound = reference_fraction * max_flux
            parent = pfba(model)

        def step(genes: list[str], reactions: list[str], sol) -> MomaStep:
            uptake = -sol.fluxes[substrate]
            flux = sol.fluxes[target]
            return MomaStep(list(genes), list(reactions), float(sol.fluxes[biomass]), float(flux),
                            float(flux / uptake) if uptake > 1e-9 else 0.0)

        steps = [step([], [], parent)]
        active = {r for r, v in parent.fluxes.items() if abs(v) > 1e-6}
        pool = {}
        for rxn_id in active:
            rxn = model.reactions.get_by_id(rxn_id)
            if rxn.boundary or rxn.objective_coefficient or len(rxn.compartments) > 1:
                continue
            for gene in rxn.genes:
                if gene.id not in exclude_genes and gene.id not in pool:
                    pool[gene.id] = None
        pool = {g: r for g in pool if (r := _disabled_reactions(model, g))}

        add_moma(model, parent, linear=True)
        knocked: list[str] = []
        for _ in range(max_knockouts):
            best = None
            for gene_id in pool:
                if gene_id in knocked:
                    continue
                with model:
                    for g in [*knocked, gene_id]:
                        model.genes.get_by_id(g).knock_out()
                    sol = model.optimize()
                    if sol.status != "optimal":
                        continue
                    cand = step([*knocked, gene_id], [r for g in [*knocked, gene_id] for r in pool[g]], sol)
                if cand.growth >= min_growth_fraction * steps[0].growth and cand.product_yield > (
                        best or steps[-1]).product_yield * (1 + 1e-4):
                    best = cand
            if best is None:
                break
            knocked = best.knocked_out
            steps.append(best)
    return steps

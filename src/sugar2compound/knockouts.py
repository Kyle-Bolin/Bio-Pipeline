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

from .simulate import production_at_max_growth

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
        base = production_at_max_growth(model, target)
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
                    res = production_at_max_growth(model, target)
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

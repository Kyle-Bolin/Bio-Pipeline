"""Flux balance analysis (FBA) helpers.

FBA in one paragraph: let S be the stoichiometric matrix (rows = metabolites,
columns = reactions) and v the vector of reaction rates ("fluxes"). At steady
state nothing accumulates, so S v = 0. Each flux has bounds lb <= v <= ub
(e.g. glucose uptake <= 10). Among all v that satisfy those constraints, FBA
picks the one maximising an objective, usually the biomass reaction (growth).
That is a linear program, solved in milliseconds.
"""

from __future__ import annotations

from dataclasses import dataclass

import cobra
import pandas as pd
from cobra.flux_analysis import production_envelope

from .models import GLUCOSE_EXCHANGE


@dataclass
class ProductionResult:
    target: str
    growth: float  # 1/h: biomass flux at the optimum
    target_flux_min: float  # mmol/gDW/h: least product the cell is *forced* to make at max growth
    target_flux_max: float  # mmol/gDW/h: most product it *could* make at max growth
    glucose_uptake: float  # mmol/gDW/h
    max_theoretical_yield: float  # mol product / mol glucose if growth were zero

    @property
    def guaranteed_yield(self) -> float:
        """mol product per mol glucose that is guaranteed when the cell grows as fast as it can."""
        return self.target_flux_min / self.glucose_uptake if self.glucose_uptake else 0.0


def glucose_uptake(model: cobra.Model, solution: cobra.Solution | None = None) -> float:
    solution = solution or model.optimize()
    return -solution.fluxes[GLUCOSE_EXCHANGE]  # uptake is negative flux on an exchange


def max_theoretical_yield(model: cobra.Model, target: str) -> float:
    """Best possible mol product / mol glucose, ignoring the cell's need to grow."""
    with model:  # changes inside a `with model:` block are undone on exit
        model.objective = target
        sol = model.optimize()
        if sol.status != "optimal" or sol.objective_value <= 0:
            return 0.0
        return sol.objective_value / glucose_uptake(model, sol)


def production_at_max_growth(model: cobra.Model, target: str, fraction_of_optimum: float = 0.999) -> ProductionResult:
    """Grow as fast as possible, then ask how much ``target`` the cell must / could make.

    A cell evolves to grow, not to help us, so the *minimum* product flux at max
    growth is the pessimistic, realistic number. Strain designers call a strain
    "growth-coupled" when that minimum is above zero: the cell can only grow
    fast by making our product.
    """
    sol = model.optimize()
    if sol.status != "optimal" or sol.objective_value < 1e-6:
        return ProductionResult(target, 0.0, 0.0, 0.0, 0.0, 0.0)
    growth = sol.objective_value
    uptake = glucose_uptake(model, sol)
    fva = cobra.flux_analysis.flux_variability_analysis(
        model, reaction_list=[target], fraction_of_optimum=fraction_of_optimum
    )
    return ProductionResult(
        target=target,
        growth=growth,
        target_flux_min=max(fva.at[target, "minimum"], 0.0),
        target_flux_max=fva.at[target, "maximum"],
        glucose_uptake=uptake,
        max_theoretical_yield=max_theoretical_yield(model, target),
    )


def envelope(model: cobra.Model, target: str, points: int = 20) -> pd.DataFrame:
    """Trade-off curve between growth and product flux (the "production envelope").

    Each row is a growth rate; the columns give the min and max product flux
    that are possible at that growth rate.
    """
    return production_envelope(model, reactions=[model.reactions.get_by_id(_biomass_id(model))],
                               objective=target, points=points)


def _biomass_id(model: cobra.Model) -> str:
    (rxn,) = [r for r in model.reactions if r.objective_coefficient != 0]
    return rxn.id

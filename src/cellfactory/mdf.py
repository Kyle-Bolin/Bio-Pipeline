"""Max-min Driving Force (MDF): can a whole route run forward at once?

``thermodynamics.reaction_thermodynamics`` asks of each reaction on its own
whether *some* concentrations let it run. A route can pass that test step by
step and still be impossible: two neighbouring steps may each need the shared
intermediate at opposite ends of its range. MDF (Noor et al. 2014, PLoS Comput
Biol 10:e1003483, doi:10.1371/journal.pcbi.1003483) asks the joint question:

    choose one concentration per metabolite, within 1 µM-10 mM, so that the
    least downhill step of the route is as downhill as possible.

That is a linear program in ln(concentration). Its optimum B is the MDF, in
kJ/mol. B <= 0 means no concentration profile in range lets every step run in
the direction the route uses: a hard kill. A small positive B means the route
can run, but its bottleneck step sits close to equilibrium, where the enzyme
wastes much of its effort on the backward reaction (at a driving force of
RT ln 10, about 5.7 kJ/mol, one backward turnover for every ten forward).

The route is a flux distribution: reaction id -> flux, the sign giving the
direction. ``product_route_fluxes`` picks the obvious one for a product: the
parsimonious flux distribution at maximum product yield. Steps without a ΔG'°
estimate (transport, boundary, biomass, compounds eQuilibrator lacks) are left
out of the program and listed. Leaving out a step removes a constraint, so a
positive MDF over the remaining steps is an upper bound on the route's real
MDF, while a non-positive one is a kill that stays valid.

This is a simpler version of equilibrator-pathway's MDF, which also takes the
correlations between ΔG'° errors into account (two reactions sharing a
compound share its error). Here each reaction's error is treated alone: the
"optimistic" MDF lets every ΔG'° sit 1.96 standard errors in the route's
favour at once, which is more generous than the joint uncertainty allows. That
generosity is deliberate; it only makes the kill harder to reach.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import cobra
import pandas as pd

from .thermodynamics import CONCENTRATION_BOUNDS, R, Z95, Conditions, reaction_thermodynamics

LOW_DRIVING_FORCE = R * 298.15 * math.log(10)  # kJ/mol; below this the bottleneck runs backwards >1 time in 10
FLUX_TOL = 1e-6  # mmol/gDW/h; smaller fluxes are solver noise, not part of the route
ROUTE_COLUMNS = ["reaction", "name", "flux", "dg0_prime", "dg0_sigma", "method", "driving_force", "bottleneck", "note"]


@dataclass
class MDFResult:
    """Outcome of ``max_min_driving_force``. Driving forces are -ΔG' in the direction the route uses, kJ/mol."""

    mdf: float | None  # with eQuilibrator's point estimates
    mdf_optimistic: float | None  # with every ΔG'° shifted 1.96 sigma in the route's favour
    verdict: str  # "feasible", "low driving force", "uncertain", "infeasible" or "not computed"
    reactions: pd.DataFrame  # one row per route step, ROUTE_COLUMNS
    concentrations: pd.DataFrame = field(default_factory=pd.DataFrame)  # metabolite, concentration (M) at the MDF
    reason: str = ""

    @property
    def bottlenecks(self) -> list[str]:
        """Steps that hold the MDF down: making any of them more favourable would raise it."""
        return list(self.reactions.loc[self.reactions["bottleneck"].astype(bool), "reaction"])

    @property
    def excluded(self) -> list[str]:
        return list(self.reactions.loc[self.reactions["dg0_prime"].isna(), "reaction"])


def _is_fixed(met: cobra.Metabolite) -> bool:
    """Water's activity is 1 and H+ is set by pH: neither has a concentration to choose."""
    base = met.id[: -len(met.compartment) - 1] if met.compartment and met.id.endswith("_" + met.compartment) else met.id
    return base in ("h2o", "h") or (met.formula, met.charge) in (("H2O", 0), ("H", 1))


def product_route_fluxes(model: cobra.Model, product: str) -> pd.Series:
    """The reactions a cell would use to make ``product`` (an exchange or demand id) at maximum yield.

    Maximises the product flux under the model's current medium and bounds,
    then picks the distribution with the least total flux (parsimonious FBA),
    which drops pointless loops. Returns the non-zero fluxes. The model is not
    changed.
    """
    message = f"{model.id} can't make {product} under these conditions, so there is no route to screen."
    with model:
        model.objective = product
        try:
            solution = cobra.flux_analysis.pfba(model)
        except cobra.exceptions.OptimizationError as err:
            raise ValueError(message) from err
    if solution.status != "optimal" or solution.fluxes[product] < FLUX_TOL:
        raise ValueError(message)
    fluxes = solution.fluxes
    return fluxes[fluxes.abs() > FLUX_TOL]


def max_min_driving_force(
    model: cobra.Model,
    fluxes: pd.Series | dict[str, float],
    table: pd.DataFrame | None = None,
    conditions: Conditions = Conditions(),
    concentrations: tuple[float, float] = CONCENTRATION_BOUNDS,
    bounds: dict[str, tuple[float, float]] | None = None,
    estimator=None,
) -> MDFResult:
    """Max-min Driving Force of the route given by ``fluxes`` (reaction id -> flux; the sign is the direction).

    ``table`` is a ``reaction_thermodynamics`` table to take ΔG'° from; without
    one, the route's reactions are estimated now (``estimator`` as there).
    ``concentrations`` bounds every metabolite (molar); ``bounds`` overrides it
    per metabolite id, e.g. to pin a measured cofactor level.
    Raises ``ThermoUnavailable`` if it has to estimate and eQuilibrator can't run.
    """
    fluxes = pd.Series(fluxes, dtype=float)
    fluxes = fluxes[fluxes.abs() > FLUX_TOL]
    if table is None:
        table = reaction_thermodynamics(model, list(fluxes.index), conditions=conditions,
                                        concentrations=concentrations, estimator=estimator)
    thermo = table.set_index("reaction")
    rows, steps = [], []
    for rid, flux in fluxes.items():
        rxn = model.reactions.get_by_id(rid)
        row = thermo.loc[rid] if rid in thermo.index else None
        if row is None or pd.isna(row["dg0_prime"]) or pd.isna(row["dg_prime_min"]):
            note = "not in the thermodynamics table" if row is None else (row["note"] or "no usable ΔG'° estimate")
            rows.append([rid, rxn.name, flux, None, None, None, None, False, note])
            continue
        rows.append([rid, rxn.name, flux, row["dg0_prime"], row["dg0_sigma"], row["method"], None, False, ""])
        steps.append((len(rows) - 1, rxn, 1.0 if flux > 0 else -1.0, row["dg0_prime"], row["dg0_sigma"]))
    route = pd.DataFrame(rows, columns=ROUTE_COLUMNS)
    if not steps:
        return MDFResult(None, None, "not computed", route, reason="no step of the route has a ΔG'° estimate")

    rt = R * conditions.temperature
    best = _solve(model, steps, rt, concentrations, bounds or {}, margin=0.0)
    optimistic = _solve(model, steps, rt, concentrations, bounds or {}, margin=Z95)
    mdf, ln_c, forces, duals = best
    for (i, *_), force, dual in zip(steps, forces, duals):
        route.at[i, "driving_force"] = force
        route.at[i, "bottleneck"] = dual > 1e-9
    conc = pd.DataFrame([(m, math.exp(x)) for m, x in ln_c.items()], columns=["metabolite", "concentration"])
    result = MDFResult(mdf, optimistic[0], "", route, conc)
    result.verdict, result.reason = _verdict(result)
    return result


def _solve(model, steps, rt, concentrations, bounds, margin):
    """maximise B subject to  -d_i (ΔG'°_i - d_i margin sigma_i + RT S_i . ln c) >= B  for every step i."""
    interface = model.solver.interface
    lp = interface.Model()
    b = interface.Variable("B", lb=-1e4, ub=1e4)  # bounded so a route of only fixed metabolites can't be unbounded
    ln_c, constraints = {}, []
    for _, rxn, *_ in steps:
        for met in rxn.metabolites:
            if met.id not in ln_c and not _is_fixed(met):
                low, high = bounds.get(met.id, concentrations)
                ln_c[met.id] = interface.Variable(f"lnc_{len(ln_c)}", lb=math.log(low), ub=math.log(high))
    for k, (_, rxn, direction, dg0, sigma) in enumerate(steps):
        dg = dg0 - direction * margin * sigma
        expr = sum(-direction * rt * coeff * ln_c[met.id] for met, coeff in rxn.metabolites.items() if met.id in ln_c)
        # driving force = -d (dg + RT S.lnc) >= B   <=>   expr - B >= d dg
        constraints.append(interface.Constraint(expr - b, lb=direction * dg, name=f"step_{k}"))
    lp.add([b, *ln_c.values(), *constraints])
    lp.objective = interface.Objective(b, direction="max")
    if lp.optimize() != "optimal":
        raise RuntimeError(f"MDF linear program did not solve ({lp.status})")
    values = {m: v.primal for m, v in ln_c.items()}
    # Reported at the point estimates, whatever margin the program used.
    forces = [-direction * (dg0 + rt * sum(coeff * values[met.id] for met, coeff in rxn.metabolites.items()
                                           if met.id in values))
              for _, rxn, direction, dg0, _ in steps]
    duals = [abs(c.dual or 0.0) for c in constraints]
    return b.primal, values, forces, duals


def _verdict(result: MDFResult) -> tuple[str, str]:
    groups = result.reactions.loc[result.reactions["bottleneck"].astype(bool), "method"].eq("group contribution").any()
    if result.mdf_optimistic <= 0 and not groups:
        return "infeasible", ("no concentrations between the bounds let every step run in the route's direction, "
                              "even giving each ΔG'° estimate the benefit of its error")
    if result.mdf <= 0:
        return "uncertain", ("the route is infeasible with eQuilibrator's best estimates but not beyond their error"
                             + (", and a bottleneck step rests on a group-contribution estimate" if groups else ""))
    if result.mdf < LOW_DRIVING_FORCE:
        return "low driving force", (f"feasible, but the bottleneck step is within {LOW_DRIVING_FORCE:.1f} kJ/mol of "
                                     "equilibrium, so its enzyme spends much of its effort running backwards")
    return "feasible", "every step can be kept comfortably downhill at once"


def mdf_caveats(result: MDFResult, label: str) -> list[str]:
    """Plain-words caveats for the report."""
    r = result.reactions
    if result.verdict == "not computed":
        return [f"Max-min driving force of {label} not computed: {result.reason}."]
    estimated = int(r["dg0_prime"].notna().sum())
    head = (f"Max-min driving force of {label}: {result.mdf:.1f} kJ/mol "
            f"({result.mdf_optimistic:.1f} allowing for estimate error), {result.verdict}: {result.reason}.")
    if result.bottlenecks:
        head += f" Bottleneck: {', '.join(result.bottlenecks[:8])}."
    out = [head]
    if result.verdict in ("infeasible", "uncertain"):
        out.append(f"Yields computed for {label} may not be reachable: FBA used a set of directions that "
                   "thermodynamics says cannot all hold at once at 1 µM-10 mM.")
    if estimated < len(r):
        out.append(f"{len(r) - estimated} of {len(r)} steps in {label} (transport, exchange, biomass or compounds "
                   "without an estimate) were left out of the MDF, so the true value can only be lower.")
    return out

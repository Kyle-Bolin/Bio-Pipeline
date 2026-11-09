"""Knockout design with StrainDesign: OptKnock, RobustKnock and minimal cut sets.

``knockouts.greedy_knockout_search`` tries one deletion at a time and keeps the
best. That misses designs where no single deletion helps but a combination
does: deleting one competing route just sends the carbon down another one, and
only closing both forces it to the product. StrainDesign (Schneider, Bekiaris &
Klamt 2022, Bioinformatics 38:4981, doi:10.1093/bioinformatics/btac632) instead
writes the whole question as one mixed-integer linear program (MILP): pick at
most N reactions to delete, binary on/off variables, such that the cell's best
growth forces product out. A MILP solver searches every combination at once.

Three formulations, from optimistic to strict:

* ``optknock`` (Burgard, Pharkya & Maranas 2003): maximise product flux at the
  mutant's maximum growth, *assuming the cell picks the flux state we like best*
  among the equally fast ones. Its designs can be "coupled" only on paper.
* ``robustknock`` (Tepper & Shlomi 2010): maximise the *minimum* product flux at
  maximum growth, i.e. the guaranteed yield this package reports everywhere.
* ``mcs`` (minimal cut sets; von Kamp & Klamt 2017): the fewest deletions after
  which every flux state with growth >= ``min_growth`` makes at least
  ``min_yield`` mol product per mol substrate, at any growth rate, not just the
  fastest. This is "strong" growth coupling.

Whatever the solver claims, every design is re-checked here with the same FBA
the rest of the report uses (``simulate.production_at_max_growth``), and the
table reports those numbers, not the MILP's. Each design is also translated to
genes: the smallest set of genes whose deletion removes each reaction, and any
*other* reactions those genes take down with them (shared enzyme subunits),
with the yield re-checked at the gene level.

Solvers: StrainDesign supports Gurobi, CPLEX, SCIP and GLPK. SCIP (open source)
comes with StrainDesign itself. GLPK is never picked automatically: it has no
indicator constraints and stalls on genome-scale problems; pass
``solver="glpk"`` to use it anyway on a small model.
"""

from __future__ import annotations

import contextlib
import itertools
import logging
import math
import time
import warnings
from dataclasses import dataclass, field

import cobra
import pandas as pd

from .knockouts import SearchResult, candidate_reactions
from .simulate import DEFAULT_SUBSTRATE, _biomass_id, max_theoretical_yield, production_at_max_growth

METHODS = {"optknock": "OptKnock", "robustknock": "RobustKnock", "mcs": "Minimal cut sets"}
SOLVER_PREFERENCE = ("gurobi", "cplex", "scip")  # GLPK only on request
FLUX_TOL = 1e-6
DEFAULT_MIN_YIELD_FRACTION = 0.3  # mcs: demand at least 30% of the chemical ceiling unless told otherwise
BIG_M = 1000.0  # retry value when indicator constraints run into numerical trouble
COLUMNS = ["rank", "knockouts", "n_knockouts", "growth", "guaranteed_yield", "best_case_yield", "coupling",
           "genes", "also_disabled", "gene_level_growth", "gene_level_guaranteed_yield"]


class DesignUnavailable(RuntimeError):
    """StrainDesign or a suitable MILP solver is not installed."""


@dataclass
class DesignResult:
    """Outcome of ``design_knockouts``. Yields are mol product / mol substrate, re-checked with FBA."""

    method: str
    target: str
    solver: str | None
    status: str  # "optimal", "time limit", "infeasible" (no design within the limits) or "error"
    designs: pd.DataFrame = field(default_factory=lambda: pd.DataFrame(columns=COLUMNS))
    wild_type_yield: float = 0.0  # guaranteed yield before any deletion
    min_growth: float = 0.0
    min_yield: float | None = None  # mcs only
    notes: list[str] = field(default_factory=list)

    @property
    def best(self) -> pd.Series | None:
        return None if self.designs.empty else self.designs.iloc[0]


def milp_solver(preferred: str | None = None) -> str:
    """The MILP solver to hand StrainDesign: ``preferred`` if given and installed, else the best installed one."""
    sd = _straindesign()
    available = set(getattr(sd, "avail_solvers", set()))
    if preferred:
        if preferred not in available:
            raise DesignUnavailable(f"solver {preferred!r} is not installed (installed: "
                                    f"{', '.join(sorted(available)) or 'none'})")
        return preferred
    for name in SOLVER_PREFERENCE:
        if name in available:
            return name
    raise DesignUnavailable("no MILP solver suited to strain design is installed (need Gurobi, CPLEX or SCIP; "
                            "`pip install pyscipopt` gives SCIP). GLPK is installed but stalls on these problems; "
                            "pass solver='glpk' to use it on a small model anyway")


def _straindesign():
    try:
        import straindesign
    except ImportError as err:
        raise DesignUnavailable("StrainDesign is not installed (pip install 'cellfactory[design]')") from err
    return straindesign


def _not_straindesign(record: logging.LogRecord) -> bool:
    return "straindesign" not in record.pathname.replace("\\", "/").split("/")


@contextlib.contextmanager
def _quiet():
    """StrainDesign reports progress through the root logger at WARNING; keep its records out of the output.

    Only records logged from StrainDesign's own files are dropped, so other threads (the API server) keep theirs.
    """
    logging.root.addFilter(_not_straindesign)
    try:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="Solver status is 'infeasible'")
            yield
    finally:
        logging.root.removeFilter(_not_straindesign)


def _term(coef: float, rxn_id: str) -> str:
    return f"{coef:.6g} {rxn_id}"


def _modules(sd, model: cobra.Model, method: str, target: str, substrate: str, biomass: str, min_growth: float,
             min_yield: float | None):
    growth = f"{biomass} >= {min_growth:.6g}"
    if method == "mcs":
        # Forbid every state that grows yet yields less than min_yield: product + min_yield * EX_substrate <= 0,
        # since substrate uptake is negative flux on its exchange. Then make sure growth stays possible.
        weak = f"{target} + {_term(min_yield, substrate)} <= 0"
        return [sd.SDModule(model, sd.names.SUPPRESS, constraints=[weak, growth]),
                sd.SDModule(model, sd.names.PROTECT, constraints=[growth])]
    kind = sd.names.OPTKNOCK if method == "optknock" else sd.names.ROBUSTKNOCK
    return [sd.SDModule(model, kind, inner_objective=biomass, outer_objective=target, constraints=[growth])]


def _solve(sd, model, modules, ko_cost, max_knockouts, max_designs, solver, time_limit, big_m=None):
    kwargs = dict(sd_modules=modules, max_cost=max_knockouts, max_solutions=max_designs,
                  solution_approach=sd.names.BEST, solver=solver, ko_cost=ko_cost, time_limit=time_limit)
    if big_m:
        kwargs["M"] = big_m
    with _quiet():
        return sd.compute_strain_designs(model, **kwargs)


def _gene_set(rxn: cobra.Reaction, max_size: int = 6, max_tries: int = 5000) -> list[str]:
    """Fewest genes whose deletion switches ``rxn`` off, preferring genes that serve fewest other reactions."""
    genes = sorted(rxn.genes, key=lambda g: (len(g.reactions), g.id))
    tries = 0
    for size in range(1, min(len(genes), max_size) + 1):
        for combo in itertools.combinations(genes, size):  # sorted, so the first hit has the fewest side effects
            tries += 1
            if not rxn.gpr.eval({g.id for g in combo}):
                return sorted(g.id for g in combo)
            if tries > max_tries:
                return sorted(g.id for g in genes)
    return sorted(g.id for g in genes)


def _forced_flux(model: cobra.Model, target: str, biomass: str, min_growth: float) -> float:
    """Least product flux over every state growing at >= min_growth (not just the fastest)."""
    with model:
        model.reactions.get_by_id(biomass).lower_bound = min_growth
        model.objective = model.reactions.get_by_id(target)
        model.objective_direction = "min"
        value = model.slim_optimize(error_value=math.nan)
    return value


def _check(model: cobra.Model, knockouts: list[str], target: str, substrate: str, biomass: str,
           min_growth: float) -> dict:
    row = {"knockouts": ", ".join(knockouts), "n_knockouts": len(knockouts)}
    with model:
        for rid in knockouts:
            model.reactions.get_by_id(rid).knock_out()
        res = production_at_max_growth(model, target, substrate)
        forced = _forced_flux(model, target, biomass, min_growth) if res.growth >= min_growth else math.nan
    row.update(growth=res.growth, guaranteed_yield=res.guaranteed_yield, best_case_yield=res.best_case_yield)
    if res.growth < min_growth:
        row["coupling"] = "too slow"
    elif not math.isnan(forced) and forced > FLUX_TOL:
        row["coupling"] = "strong"  # product forced at every growth rate >= min_growth
    elif res.guaranteed_yield > FLUX_TOL:
        row["coupling"] = "at max growth"
    else:
        row["coupling"] = "none"
    geneless = [rid for rid in knockouts if not model.reactions.get_by_id(rid).genes]
    genes = sorted({g for rid in knockouts for g in _gene_set(model.reactions.get_by_id(rid))})
    row["genes"] = ", ".join(genes)
    if geneless:
        # No gene to delete for these reactions, so the design can't be built by deleting genes.
        row.update(also_disabled=f"(no genes: {', '.join(geneless)})", gene_level_growth=math.nan,
                   gene_level_guaranteed_yield=math.nan)
        return row
    wild_bounds = {r.id: r.bounds for g in genes for r in model.genes.get_by_id(g).reactions}
    with model:
        for gid in genes:
            model.genes.get_by_id(gid).knock_out()
        # Only reactions the gene deletion actually shut: ones already closed in the model are not side effects.
        lost = sorted({rid for rid, bounds in wild_bounds.items()
                       if rid not in knockouts and bounds != (0.0, 0.0)
                       and model.reactions.get_by_id(rid).bounds == (0.0, 0.0)})
        gres = production_at_max_growth(model, target, substrate) if lost else res
    row["also_disabled"] = ", ".join(lost)
    row.update(gene_level_growth=gres.growth, gene_level_guaranteed_yield=gres.guaranteed_yield)
    return row


def design_knockouts(
    model: cobra.Model,
    target: str,
    substrate: str = DEFAULT_SUBSTRATE,
    method: str = "optknock",
    max_knockouts: int = 3,
    min_growth: float = 0.05,
    min_yield: float | None = None,
    max_designs: int = 5,
    candidates: list[str] | None = None,
    solver: str | None = None,
    time_limit: float = 120.0,
) -> DesignResult:
    """Find up to ``max_knockouts`` reaction deletions that make the cell make ``target`` to grow.

    ``method`` is ``optknock``, ``robustknock`` or ``mcs`` (see the module
    docstring). ``min_growth`` (1/h) is the slowest growth a design may leave.
    ``min_yield`` (mol/mol, ``mcs`` only) defaults to 30% of the maximum
    theoretical yield. ``candidates`` are the reactions allowed to be deleted;
    by default every reaction with genes except exchanges, biomass and
    transporters (``knockouts.candidate_reactions`` with ``active_only=False``,
    since the MILP is not limited to what the wild type uses).

    The MILP is given ``max_designs`` solutions to find and ``time_limit``
    seconds in all (a big-M retry after a solver failure gets what is left); StrainDesign expands each solution into its equivalent variants
    (deleting either of two reactions in series), so more rows may come back.
    Rows are sorted by FBA-checked guaranteed yield at gene level (rows whose
    reactions have no genes, so no gene-level value, last). Raises
    ``DesignUnavailable`` when StrainDesign or a MILP solver is missing; a
    solver failure is returned as ``status="error"`` with the reason in notes.
    """
    if method not in METHODS:
        raise ValueError(f"unknown design method {method!r}; choose one of {', '.join(METHODS)}")
    if max_knockouts < 1:
        raise ValueError("max_knockouts must be at least 1")
    sd = _straindesign()
    solver = milp_solver(solver)
    biomass = _biomass_id(model)
    for rid in (target, substrate):
        if rid not in model.reactions:
            raise KeyError(f"{rid} is not a reaction of {model.id}")
    with _quiet():
        wild = production_at_max_growth(model, target, substrate)
        if method == "mcs" and min_yield is None:
            min_yield = DEFAULT_MIN_YIELD_FRACTION * max_theoretical_yield(model, target, substrate)
    result = DesignResult(METHODS[method], target, solver, "optimal", wild_type_yield=wild.guaranteed_yield,
                          min_growth=min_growth, min_yield=min_yield if method == "mcs" else None)
    if wild.growth < min_growth:
        result.status = "infeasible"
        result.notes.append(f"the unmodified model grows at {wild.growth:.3f}/h, below the {min_growth}/h minimum")
        return result
    if method == "mcs" and not min_yield:
        result.status = "infeasible"
        result.notes.append(f"{target} cannot be made from {substrate} at all, so there is no yield to couple")
        return result
    pool = candidates if candidates is not None else candidate_reactions(model, active_only=False)
    ko_cost = {rid: 1 for rid in pool if rid in model.reactions and rid not in (target, substrate, biomass)}
    work = model.copy()  # StrainDesign compresses and edits the network it is given
    modules = _modules(sd, work, method, target, substrate, biomass, min_growth, min_yield)
    deadline = time.monotonic() + time_limit
    try:
        sols = _solve(sd, work, modules, ko_cost, max_knockouts, max_designs, solver, time_limit)
    except (ValueError, KeyError, TypeError):
        raise  # a bad problem (say an id StrainDesign can't parse), not solver trouble: retrying won't help
    except Exception as first:  # noqa: BLE001 - solver failures surface as arbitrary errors inside StrainDesign
        # Indicator constraints sometimes hit numerical trouble (seen with SCIP + RobustKnock); big-M usually works.
        remaining = deadline - time.monotonic()
        if remaining < 5:
            result.status = "error"
            result.notes.append(f"the {solver} MILP failed ({type(first).__name__}: {first}) with no time left in "
                                f"its {time_limit:g} s limit to retry")
            return result
        try:
            work = model.copy()
            modules = _modules(sd, work, method, target, substrate, biomass, min_growth, min_yield)
            sols = _solve(sd, work, modules, ko_cost, max_knockouts, max_designs, solver, remaining, BIG_M)
        except Exception as second:  # noqa: BLE001
            result.status = "error"
            result.notes.append(f"the {solver} MILP failed ({type(first).__name__}: {first}; retried with big-M: "
                                f"{type(second).__name__}: {second}). A commercial solver or a smaller model may help")
            return result
        result.notes.append(f"{solver} failed with indicator constraints, so the MILP was re-run with big-M "
                            f"constraints (M = {BIG_M:g}), which can cut off the best design: treat these as good "
                            "designs, not necessarily optimal ones")
    status = str(getattr(sols, "status", "optimal"))
    if status.startswith("time_limit"):
        result.status = "time limit"
        result.notes.append(f"the MILP hit its {time_limit:g} s time limit, so better designs may exist")
    elif status == "infeasible":
        result.status = "infeasible"
    elif status != "optimal":
        result.status = "error"
        result.notes.append(f"StrainDesign returned status {status!r}")
    seen, rows = set(), []
    with _quiet():
        for design in getattr(sols, "reaction_sd", None) or []:
            kos = sorted(rid for rid, v in design.items() if v < 0 and rid in model.reactions)
            if tuple(kos) in seen:
                continue
            seen.add(tuple(kos))
            rows.append(_check(model, kos, target, substrate, biomass, min_growth))
    if rows:
        # Rank by what deleting the genes would do, since genes, not reactions, are what a lab deletes.
        df = pd.DataFrame(rows).sort_values(["gene_level_guaranteed_yield", "guaranteed_yield", "n_knockouts", "growth"],
                                            ascending=[False, False, True, False], kind="stable", na_position="last")
        df.insert(0, "rank", range(1, len(df) + 1))
        result.designs = df[COLUMNS].reset_index(drop=True)
    elif result.status == "optimal":
        result.status = "infeasible"
    return result


def design_caveats(result: DesignResult, greedy: SearchResult | None = None) -> list[str]:
    """Plain-words caveats for the report."""
    label = f"Knockout design ({result.method})"
    limits = f"with up to the allowed deletions and growth >= {result.min_growth}/h"
    if result.status == "error":
        return [f"{label} not computed: {'; '.join(result.notes)}."]
    out = []
    best = result.best
    if best is None:
        why = (f"no design makes every growing state yield >= {result.min_yield:.2f} mol/mol" if result.min_yield
               else f"no design forces more product than the wild type")
        out.append(f"{label}: {why} {limits}" + (f" ({'; '.join(result.notes)})." if result.notes else "."))
        return out
    out.append(f"{label}, solved with {result.solver}: best design deletes {best['knockouts'] or 'nothing'} "
               f"(genes {best['genes'] or 'none'}): guaranteed yield {best['guaranteed_yield']:.3f} mol/mol at "
               f"{best['growth']:.3f}/h, coupling {best['coupling']}; wild type {result.wild_type_yield:.3f}.")
    out += [f"{label}: {n}." for n in result.notes]
    optimistic = result.designs[result.designs["coupling"].isin(["none", "too slow"])]
    if len(optimistic):
        out.append(f"{len(optimistic)} of {len(result.designs)} {result.method} designs force no product once "
                   "re-checked: the MILP assumed the cell would choose the product among equally fast flux states.")
    geneless = result.designs[result.designs["gene_level_guaranteed_yield"].isna()]
    if len(geneless):
        out.append(f"{len(geneless)} design(s) delete reactions that have no genes in the model, so they can't be built "
                   "by deleting genes; their gene-level columns are empty.")
    shifted = result.designs[(result.designs["also_disabled"] != "") & result.designs["gene_level_growth"].notna()]
    if len(shifted):
        dead = shifted[shifted["gene_level_growth"] < result.min_growth]
        out.append(f"In {len(shifted)} design(s) the genes to delete also disable other reactions (isozymes shared "
                   "with another reaction, or shared enzyme subunits); see the gene-level columns"
                   + (f", and {len(dead)} of them stop growth once done with genes." if len(dead) else "."))
    if greedy is not None and greedy.steps:
        g = greedy.steps[-1].guaranteed_yield
        if best["guaranteed_yield"] > g + 1e-3:
            out.append(f"The {result.method} design beats the greedy search ({g:.3f} mol/mol): it found a "
                       "combination of deletions where no single deletion helps on its own.")
        elif g > best["guaranteed_yield"] + 1e-3:
            out.append(f"The greedy search found a better guaranteed yield ({g:.3f} mol/mol) than the "
                       f"{result.method} design; the MILP may have hit its limits or solved a different objective.")
    out.append("Knockout designs assume the mutant grows as fast as its network allows; freshly built strains "
               "usually need adaptive evolution to get there, and regulation, enzyme levels and toxicity are not "
               "modelled.")
    return out

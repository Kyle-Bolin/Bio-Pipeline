"""cellfactory as a set of JSON-in, JSON-out tools, for AI agents and scripts.

Each function takes plain arguments (strings, numbers, booleans), returns a
JSON-safe dict, and raises ``ToolInputError`` with a message that says what to
try instead. ``mcp_server.py`` exposes every function in ``TOOLS`` over MCP;
nothing here depends on MCP, so the same calls work from Python or a notebook.

Loaded models are cached per process and every call runs inside ``with model:``,
so one call's medium or demand reaction never leaks into the next.
"""

from __future__ import annotations

import math
import threading
import warnings
from typing import Any

import cobra
import pandas as pd

from .checks import free_energy_checks, mass_balance_problems, predicted_essential_genes
from .heterologous import MET, PATHWAYS, add_pathway, find_or_create
from .knockouts import candidate_reactions, greedy_knockout_search, moma_knockout_search
from .models import BUNDLED, find_metabolites, load_model, resolve_exchange, set_medium
from .pipeline import analyze
from .precursors import fate_to_text, nearest_precursors, product_fate, routes_to_text
from .simulate import _biomass_id, envelope

MAX_ROWS = 50  # cap on list-shaped results, so a genome-scale answer fits in an agent's context

MODEL_NOTES = {
    "e_coli_core": "E. coli central metabolism only (95 reactions). Fast; use it to try things out.",
    "iJO1366": "E. coli K-12 MG1655, genome scale (2,583 reactions, 2011).",
    "iYS1720": "Salmonella Typhimurium LT2, genome scale (3,357 reactions).",
}


class ToolInputError(ValueError):
    """Bad input the caller can fix: unknown model, unknown compound, no growth, ..."""


_cache: dict[str, cobra.Model] = {}
_locks: dict[str, threading.Lock] = {}
_cache_lock = threading.Lock()


def _model(model_id: str) -> tuple[cobra.Model, threading.Lock]:
    with _cache_lock:
        if model_id not in _cache:
            try:
                _cache[model_id] = load_model(model_id)
            except (RuntimeError, OSError) as err:
                reason = str(err).rstrip(".")
                raise ToolInputError(f"Could not load model {model_id!r}: {reason}. "
                                     "Call list_models for ids that work offline.") from err
            _locks[model_id] = threading.Lock()
        return _cache[model_id], _locks[model_id]


def _clean(value: Any) -> Any:
    """Make a value JSON-safe: NaN/inf become None, numpy and pandas types become Python ones."""
    if isinstance(value, dict):
        return {str(k): _clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_clean(v) for v in value]
    if isinstance(value, pd.DataFrame):
        return _clean(value.to_dict(orient="records"))
    if hasattr(value, "item") and not isinstance(value, (str, bytes)):
        value = value.item()  # numpy scalar
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _rows(df: pd.DataFrame, limit: int = MAX_ROWS) -> dict:
    return {"total": len(df), "truncated": len(df) > limit, "rows": _clean(df.head(limit))}


def _run(model_id: str, fn):
    """Run ``fn(model)`` on the cached model with every change rolled back afterwards."""
    model, lock = _model(model_id)
    with lock, model, warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Solver status is 'infeasible'")
        try:
            return fn(model)
        except KeyError as err:
            raise ToolInputError(err.args[0] if err.args else str(err)) from err
        except ValueError as err:
            raise ToolInputError(str(err)) from err


def _add_pathway(m: cobra.Model, pathway: str) -> dict | None:
    """Insert the named heterologous pathway into ``m`` (inside ``_run``, so it is rolled back).

    Raises ToolInputError for an unknown name, or when the host lacks the
    precursor the pathway starts from (the product could never be made).
    """
    if not pathway:
        return None
    spec = PATHWAYS.get(pathway.lower())
    if spec is None:
        raise ToolInputError(f"Unknown pathway {pathway!r}. Choose one of: {', '.join(PATHWAYS)} "
                             "(list_pathways describes them).")
    made = {k for step in spec.steps for k, c in step.stoich.items() if c > 0}
    for k in sorted({k for step in spec.steps for k, c in step.stoich.items() if c < 0} - made):
        if find_or_create(m, k)[1]:
            raise ToolInputError(f"{m.id} has no {MET[k].name} for the {pathway} pathway to start from. "
                                 "Use a genome-scale model such as iJO1366.")
    steps = add_pathway(m, spec)
    product = steps.pop("product")
    return {"name": pathway.lower(), "product_metabolite": product, "steps": steps,
            "genes": sorted({g for s in spec.steps for g in s.genes.split()})}


def list_models() -> dict[str, Any]:
    """Models that load offline, and how to use any other one.

    Any BiGG id (http://bigg.ucsd.edu/models, e.g. iML1515 for E. coli, iMM904
    for yeast) also works if the machine can reach BiGG; so does a path to an
    SBML (.xml), JSON or MATLAB (.mat) model file.
    """
    return {
        "bundled": [{"id": k, "description": MODEL_NOTES.get(k, "")} for k in BUNDLED],
        "other_sources": "Any BiGG model id (downloaded and cached on first use), or a path to an "
                         ".xml/.sbml/.json/.mat model file.",
    }


def describe_model(model: str = "e_coli_core") -> dict[str, Any]:
    """Size, growth objective, default medium and default growth rate of a model."""
    def go(m: cobra.Model) -> dict:
        objective = [r.id for r in m.reactions if r.objective_coefficient]
        return {
            "model": m.id,
            "name": m.name,
            "reactions": len(m.reactions),
            "metabolites": len(m.metabolites),
            "genes": len(m.genes),
            "compartments": dict(m.compartments),
            "objective": objective,
            "medium": dict(sorted(m.medium.items())),
            "growth_rate_default_medium": m.slim_optimize(error_value=0.0),
        }
    return _clean(_run(model, go))


def list_pathways() -> dict[str, Any]:
    """Non-native pathways that can be added to a host with the ``pathway`` argument of other tools.

    E. coli and yeast make neither vanillin nor lycopene. Passing
    ``pathway="vanillin"`` to simulate_production, search_knockouts,
    search_knockouts_moma or production_envelope inserts the missing enzymes
    (from plants, fungi or other bacteria) before simulating; then pass the
    pathway's product name (e.g. "vanillin") as ``product``. Steps the host
    already has natively are reused, not duplicated. These pathways need a
    genome-scale model (iJO1366): e_coli_core lacks their precursors.
    """
    rows = []
    for name, p in PATHWAYS.items():
        rows.append({
            "name": name,
            "product": MET[p.product].name,
            "molar_mass_g_per_mol": p.molar_mass,
            "description": p.description,
            "steps": [{"reaction": s.id, "enzyme": s.name, "genes": s.genes,
                       "equation": " + ".join(f"{-c:g} {MET[k].name}" for k, c in s.stoich.items() if c < 0)
                                   + " -> " + " + ".join(f"{c:g} {MET[k].name}" for k, c in s.stoich.items() if c > 0)}
                      for s in p.steps],
        })
    return {"pathways": rows}


def search_metabolites(query: str, model: str = "e_coli_core", limit: int = 25) -> dict[str, Any]:
    """Find metabolite ids by name or id fragment (e.g. "succin"), before calling the other tools.

    ``exchange`` is the reaction that moves the compound in or out of the cell;
    a compound with one can be used as substrate or product directly. One
    without can still be a product (a demand reaction is added for it).
    """
    def go(m: cobra.Model) -> dict:
        hits = find_metabolites(m, query)
        rows = [{"id": met.id, "name": met.name, "formula": met.formula, "compartment": met.compartment,
                 "exchange": next((r.id for r in met.reactions if r in m.exchanges), None)} for met in hits]
        return {"model": m.id, "query": query, "total": len(rows), "truncated": len(rows) > limit,
                "metabolites": rows[:limit]}
    return _clean(_run(model, go))


def simulate_production(product: str, substrate: str = "glucose", model: str = "e_coli_core",
                        uptake: float = 10.0, aerobic: bool = True, pathway: str = "") -> dict[str, Any]:
    """Feed ``substrate`` to ``model`` and report growth and ``product`` yield.

    Substrate and product accept everyday names ("glucose", "ethanol"), metabolite
    ids ("etoh") or exchange ids ("EX_etoh_e"). Yields are mol product per mol
    substrate. ``guaranteed_yield`` (least product made at max growth) is the
    realistic number; ``max_theoretical_yield`` ignores growth and is the ceiling.

    For a compound the organism cannot make, set ``pathway`` to a name from
    list_pathways (e.g. "lycopene", with model "iJO1366") to add the missing
    enzymes first; ``pathway`` in the result says which steps were added and
    which the host already had.
    """
    def go(m: cobra.Model) -> dict:
        added = _add_pathway(m, pathway)
        out = analyze(m, substrate, product, uptake=uptake, aerobic=aerobic).to_dict()
        if added:
            out["pathway"] = added
        return out
    return _clean(_run(model, go))


def search_knockouts(product: str, substrate: str = "glucose", model: str = "e_coli_core",
                     max_knockouts: int = 3, min_growth: float = 0.05,
                     uptake: float = 10.0, aerobic: bool = True, pathway: str = "") -> dict[str, Any]:
    """Greedy search for up to ``max_knockouts`` reaction deletions that raise guaranteed yield.

    Returns one step per deletion with the reactions removed, the genes to
    delete in the lab, growth (1/h) and guaranteed yield. Takes seconds on
    e_coli_core and minutes on genome-scale models. Designs growing slower than
    ``min_growth`` are rejected. Verify a design with check_model before trusting it.
    ``pathway`` adds a non-native pathway first, as in simulate_production.

    For a product the cell has no reason to make (most non-native compounds) this
    FBA search usually finds nothing; try search_knockouts_moma instead.
    """
    if not 1 <= max_knockouts <= 10:
        raise ToolInputError("max_knockouts must be between 1 and 10.")

    def go(m: cobra.Model) -> dict:
        added = _add_pathway(m, pathway)
        report = analyze(m, substrate, product, uptake=uptake, aerobic=aerobic,
                         max_knockouts=max_knockouts, min_growth=min_growth)
        out = report.to_dict()
        if added:
            out["pathway"] = added
        steps = out["knockouts"]
        out["summary"] = {
            "wild_type_guaranteed_yield": steps[0]["guaranteed_yield"],
            "best_guaranteed_yield": steps[-1]["guaranteed_yield"],
            "knockouts": steps[-1]["knocked_out"],
            "genes": steps[-1]["genes"],
            "improved": len(steps) > 1,
        }
        return out
    return _clean(_run(model, go))


def search_knockouts_moma(product: str, substrate: str = "glucose", model: str = "e_coli_core",
                          pathway: str = "", max_knockouts: int = 2, min_growth_fraction: float = 0.3,
                          reference_fraction: float = 0.01, uptake: float = 10.0,
                          aerobic: bool = True) -> dict[str, Any]:
    """Greedy gene-deletion search scored with linear MOMA, for products the cell has no reason to make.

    search_knockouts (FBA) assumes the mutant re-optimises for growth, so for a
    non-native product no deletion forces production and it finds nothing.
    MOMA (Segre, Vitkup & Church 2002) instead assumes a fresh mutant keeps its
    fluxes as close as possible to its parent's, which predicts how a deletion
    reroutes carbon toward the product; Alper et al. 2005 picked lycopene
    targets this way. The parent sends ``reference_fraction`` of the maximum
    possible product flux to the product. Mutants must keep
    ``min_growth_fraction`` of the parent's growth.

    Returns one step per deleted gene: the genes, the reactions they disable,
    growth (1/h), product flux (mmol/gDW/h) and product yield (mol/mol
    substrate). Step 0 is the parent. Genes of an added ``pathway`` are never
    deleted. Slow: seconds on e_coli_core, around a minute per knockout round
    on iJO1366, so keep ``max_knockouts`` at 2 or 3. MOMA yields are
    predictions of a fresh mutant, not a growth-coupled guarantee.
    """
    if not 1 <= max_knockouts <= 5:
        raise ToolInputError("max_knockouts must be between 1 and 5.")
    if not 0 < reference_fraction < 1:
        raise ToolInputError("reference_fraction must be between 0 and 1 (e.g. 0.01).")
    if not 0 <= min_growth_fraction < 1:
        raise ToolInputError("min_growth_fraction must be between 0 and 1 (e.g. 0.3).")

    def go(m: cobra.Model) -> dict:
        added = _add_pathway(m, pathway)
        substrate_ex = set_medium(m, substrate, uptake, aerobic)
        target = resolve_exchange(m, product, create_demand=True)
        if target == substrate_ex:
            raise ToolInputError("Substrate and product are the same compound.")
        if m.slim_optimize(error_value=0.0) < 1e-6:
            raise ToolInputError(f"{m.id} can't grow on {substrate_ex} under these conditions. "
                                 "Try the other oxygen setting or another substrate.")
        with m:
            m.objective = target
            if m.slim_optimize(error_value=0.0) < 1e-9:
                raise ToolInputError(f"{m.id} cannot make {target} at all, so no deletion can help. "
                                     "Add a pathway (list_pathways) or check the compound with simulate_production.")
        exclude = {g for g in (added or {}).get("genes", []) if g in m.genes}
        steps = moma_knockout_search(m, target, substrate=substrate_ex, reference_fraction=reference_fraction,
                                     max_knockouts=max_knockouts, min_growth_fraction=min_growth_fraction,
                                     exclude_genes=exclude)
        rows = [{"knocked_out": s.knocked_out,
                 "gene_names": [m.genes.get_by_id(g).name or g for g in s.knocked_out],
                 "reactions": s.reactions, "growth": s.growth, "product_flux": s.product_flux,
                 "product_yield": s.product_yield} for s in steps]
        out = {
            "model": m.id, "substrate": substrate_ex, "product": target,
            "units": "growth 1/h, product_flux mmol/gDW/h, product_yield mol/mol substrate",
            "steps": rows,
            "summary": {"parent_yield": rows[0]["product_yield"], "best_yield": rows[-1]["product_yield"],
                        "genes": rows[-1]["knocked_out"], "improved": len(rows) > 1},
            "note": "MOMA predicts a freshly made mutant's fluxes, not an evolved one; it is not growth coupling.",
        }
        if added:
            out["pathway"] = added
        return out
    return _clean(_run(model, go))


def production_envelope(product: str, substrate: str = "glucose", model: str = "e_coli_core",
                        uptake: float = 10.0, aerobic: bool = True, points: int = 20,
                        pathway: str = "") -> dict[str, Any]:
    """Trade-off between growth and product flux: for each growth rate, the min and max product flux.

    ``pathway`` adds a non-native pathway first, as in simulate_production.
    """
    def go(m: cobra.Model) -> dict:
        _add_pathway(m, pathway)
        set_medium(m, substrate, uptake, aerobic)
        target = resolve_exchange(m, product, create_demand=True)
        df = envelope(m, target, points=points)
        df = df.rename(columns={_biomass_id(m): "growth"})[["growth", "flux_minimum", "flux_maximum"]]
        return {"model": m.id, "product": target, "units": "growth 1/h, flux mmol/gDW/h", "points": _clean(df)}
    return _clean(_run(model, go))


def check_model(model: str = "e_coli_core", substrate: str = "glucose",
                aerobic: bool = True) -> dict[str, Any]:
    """Sanity checks before trusting a prediction: free-energy loops, mass balance, essential genes.

    ``ok`` is False if the model grows or makes ATP from nothing, which inflates
    every yield. Mass-imbalanced reactions and predicted essential genes are
    listed for review (the first 50 of each).
    """
    def go(m: cobra.Model) -> dict:
        set_medium(m, substrate, aerobic=aerobic)
        energy = free_energy_checks(m)
        return {
            "model": m.id,
            "ok": bool(energy["ok"].all()),
            "free_energy_checks": _clean(energy),
            "mass_imbalanced_reactions": _rows(mass_balance_problems(m)),
            "predicted_essential_genes": _rows(predicted_essential_genes(m)),
        }
    return _clean(_run(model, go))


def knockout_candidates(model: str = "e_coli_core", substrate: str = "glucose",
                        aerobic: bool = True) -> dict[str, Any]:
    """Reactions worth trying to delete: gene-associated, not transport or boundary, active under these conditions."""
    def go(m: cobra.Model) -> dict:
        set_medium(m, substrate, aerobic=aerobic)
        ids = candidate_reactions(m)
        rows = [{"reaction": r, "name": m.reactions.get_by_id(r).name,
                 "genes": sorted(g.id for g in m.reactions.get_by_id(r).genes)} for r in ids]
        return {"model": m.id, "total": len(rows), "truncated": len(rows) > MAX_ROWS * 4,
                "candidates": rows[:MAX_ROWS * 4]}
    return _clean(_run(model, go))


def find_precursors(target: str, model: str = "e_coli_core", donor: str = "iJO1366",
                    substrate: str = "glucose", uptake: float = 10.0, aerobic: bool = True,
                    max_new_enzymes: int = 4, max_total_steps: int = 8, top: int = 5) -> dict[str, Any]:
    """Which metabolites the host already makes are closest to ``target``, and what to clone in.

    Retrobiosynthesis: search the chemistry backwards from the product until it
    reaches something the organism already produces, minimising the number of
    enzymes you would have to add. Use it when you know the compound you want
    but not the pathway.

    ``donor`` is where the non-native reactions come from: another model id
    (``iJO1366`` is a good donor for ``e_coli_core``), or an empty string to
    search only the host's own chemistry. Host and donor must use the same
    identifier namespace; both BiGG is fine.

    Each route lists its steps in order, which of them the host already has,
    and the donor's genes for the ones it does not, so the new reactions are a
    cloning shortlist. ``precursor_max_flux`` is how much of the starting
    metabolite the host can make under this medium.

    Raising ``max_new_enzymes`` finds longer pathways but costs time; a real
    pathway worth building is usually under five steps. A step count is not a
    difficulty ranking: one step can be an enzyme that will not fold in E. coli.
    """
    if not 0 <= max_new_enzymes <= 15:
        raise ToolInputError("max_new_enzymes must be between 0 and 15.")
    if not 1 <= max_total_steps <= 20:
        raise ToolInputError("max_total_steps must be between 1 and 20.")
    donor_model = _model(donor)[0] if donor else None
    null_donor = "(host only)"

    def go(m: cobra.Model) -> dict:
        routes = nearest_precursors(
            m, target, donor_model, substrate=substrate, uptake=uptake, aerobic=aerobic,
            max_heterologous_steps=max_new_enzymes, max_total_steps=max_total_steps, top=top,
        )
        return {
            "model": m.id,
            "donor": donor or null_donor,
            "target": target,
            "conditions": {"substrate": substrate, "uptake": uptake, "aerobic": aerobic},
            "total": len(routes),
            "routes": [r.to_dict() for r in routes],
            "summary": routes_to_text(routes, limit=top),
            "note": "Model prediction. Step count is not difficulty, and the search assumes every "
                    "step's other inputs come from the host.",
        }
    return _clean(_run(model, go))


def check_product_fate(target: str, model: str = "e_coli_core", donor: str = "",
                       substrate: str = "glucose", uptake: float = 10.0,
                       aerobic: bool = True) -> dict[str, Any]:
    """What the host would do to ``target`` once it makes it: eat it, export it, or react with it.

    Two halves. ``consumed_by`` is solid model output: every host reaction that
    consumes the product, whether it carries flux while the cell grows, the
    genes behind it, and whether deleting it leaves the cell viable.
    ``protective_knockouts`` is that list filtered to the safe, gene-associated,
    non-transport ones.

    The rest is chemistry the model cannot see, and is labelled as heuristic.
    ``charge_note`` reads the product's charge at cytosolic pH, which decides
    whether it can leave the cell on its own or needs a transporter.
    ``reactive_pools_the_host_makes`` is measured by FBA: thiols, free amines
    and reactive oxygen species the cell is actually producing, any of which
    can attack a product without an enzyme. ``functional_group_hints`` are
    guessed from the metabolite's name, since these models carry no structures.

    Set ``donor`` to a model id when asking about a product the host cannot
    make yet, so the compound can be looked up. A target the host model does
    not contain is the good case: no native enzyme touches it.
    """
    donor_model = _model(donor)[0] if donor else None

    def go(m: cobra.Model) -> dict:
        fate = product_fate(m, target, donor_model, substrate=substrate, uptake=uptake,
                            aerobic=aerobic)
        return {"model": m.id, **fate, "summary": fate_to_text(fate)}
    return _clean(_run(model, go))


# Every function exposed over MCP. Add new tools here; tests/test_tools.py fails
# if a public cellfactory function is neither covered by a tool nor listed in
# NOT_EXPOSED with a reason.
TOOLS = [
    list_models,
    describe_model,
    list_pathways,
    search_metabolites,
    simulate_production,
    search_knockouts,
    search_knockouts_moma,
    production_envelope,
    check_model,
    knockout_candidates,
    find_precursors,
    check_product_fate,
]

# Public functions in cellfactory.__all__ and which tool covers each, or why none does.
COVERAGE = {
    "Report": "simulate_production",
    "analyze": "simulate_production",
    "candidate_reactions": "knockout_candidates",
    "Route": "find_precursors",
    "Step": "find_precursors",
    "envelope": "production_envelope",
    "fate_to_text": "check_product_fate",
    "find_metabolites": "search_metabolites",
    "free_energy_checks": "check_model",
    "greedy_knockout_search": "search_knockouts",
    "load_model": "describe_model",
    "mass_balance_problems": "check_model",
    "nearest_precursors": "find_precursors",
    "product_fate": "check_product_fate",
    "max_theoretical_yield": "simulate_production",
    "predicted_essential_genes": "check_model",
    "production_at_max_growth": "simulate_production",
    "resolve_exchange": "search_metabolites",
    "routes_to_text": "find_precursors",
    "set_medium": "simulate_production",
}
NOT_EXPOSED = {
    "gene_essentiality_check": "Only meaningful for E. coli with its built-in lab data; needs gene lists "
                               "as input. Expose once it takes another organism's screen data.",
    "set_glucose_medium": "Shortcut for set_medium; every tool already takes a substrate.",
}

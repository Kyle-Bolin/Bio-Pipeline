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
import os
import threading
import time
import warnings
from typing import Any

import cobra
import pandas as pd
from cobra.flux_analysis import pfba

from .checks import free_energy_checks, mass_balance_problems, predicted_essential_genes
from .heterologous import MET, PATHWAYS, add_pathway, find_or_create
from .identifiers import resolve_metabolite
from .kinetics import KineticsUnavailable
from .kinetics import enzyme_kinetics as _enzyme_kinetics
from .kinetics import kinetics_caveats, sabiork_lookup
from .knockouts import candidate_reactions, greedy_knockout_search, moma_knockout_search
from .mdf import max_min_driving_force as _max_min_driving_force
from .mdf import mdf_caveats
from .models import BUNDLED, find_metabolites, load_model, resolve_exchange, set_medium
from .pipeline import analyze
from .precursors import fate_to_text, nearest_precursors, product_fate, routes_to_text
from .qa import memote_qa
from .simulate import _biomass_id, envelope
from .strain_design import METHODS as DESIGN_METHODS
from .strain_design import DesignUnavailable, design_caveats
from .strain_design import design_knockouts as _design_knockouts
from .thermodynamics import Conditions, ThermoUnavailable, reaction_thermodynamics, thermo_caveats

MAX_ROWS = 50  # cap on list-shaped results, so a genome-scale answer fits in an agent's context
MAX_TIME_LIMIT = 3600  # seconds a single solve may hold a model's lock
SABIO_BUDGET = 60.0  # seconds of SABIO-RK lookups per call; the rest become the "partial" note

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
                        uptake: float = 10.0, aerobic: bool = True, pathway: str = "",
                        thermo_bounds: bool = False) -> dict[str, Any]:
    """Feed ``substrate`` to ``model`` and report growth and ``product`` yield.

    Substrate and product accept everyday names ("glucose", "ethanol"), metabolite
    ids ("etoh") or exchange ids ("EX_etoh_e"). Yields are mol product per mol
    substrate. ``guaranteed_yield`` (least product made at max growth) is the
    realistic number; ``max_theoretical_yield`` ignores growth and is the ceiling.

    For a compound the organism cannot make, set ``pathway`` to a name from
    list_pathways (e.g. "lycopene", with model "iJO1366") to add the missing
    enzymes first; ``pathway`` in the result says which steps were added and
    which the host already had.

    ``thermo_bounds=True`` first makes reversible reactions one-way where
    eQuilibrator's ΔG'° estimate rules a direction out (needs the [thermo]
    extra; see check_thermodynamics), so the yield cannot lean on a reaction
    running backwards. Compare with and without to see whether it matters.
    """
    def go(m: cobra.Model) -> dict:
        added = _add_pathway(m, pathway)
        out = analyze(m, substrate, product, uptake=uptake, aerobic=aerobic,
                      thermo_bounds=thermo_bounds).to_dict()
        thermo = out["checks"].pop("Reaction thermodynamics", None)
        if thermo is not None:  # one row per reaction: thousands on a genome-scale model
            out["checks"]["Reaction thermodynamics"] = (
                f"{len(thermo)} reactions estimated; call check_thermodynamics for the ones that matter.")
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
    deleted. Slow: under a second on e_coli_core, about two minutes per
    knockout round on iJO1366 (lycopene, measured), so keep ``max_knockouts`` at 2 or 3. MOMA yields are
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
        added = _add_pathway(m, pathway)
        set_medium(m, substrate, uptake, aerobic)
        target = resolve_exchange(m, product, create_demand=True)
        df = envelope(m, target, points=points)
        df = df.rename(columns={_biomass_id(m): "growth"})[["growth", "flux_minimum", "flux_maximum"]]
        out = {"model": m.id, "product": target, "units": "growth 1/h, flux mmol/gDW/h", "points": _clean(df)}
        if added:
            out["pathway"] = added
        return out
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

def resolve_compound(query: str, model: str = "e_coli_core", online: bool = True) -> dict[str, Any]:
    """Which model metabolite a name or database id refers to, and how the match was made.

    Accepts model ids (``etoh``, ``etoh_c``), everyday names ("glucose"), exact
    metabolite names, and database ids: ``kegg:C00031`` or bare ``C00031``,
    ``CHEBI:4167``, ``MNXM41``, ``metacyc:GLC``, ``seed:cpd00027``,
    ``HMDB0000122``, an InChIKey. Matching tries the model's own ids and names,
    then its database annotations, then (if ``online``) a MetaNetX
    cross-reference.

    ``ok`` says whether exactly one compound matched. On a miss, ``suggestions``
    lists near misses; when an id fits several model compounds (often
    stereoisomers or charge forms), ``candidates`` lists them. A ``caveat`` means
    the match came through MetaNetX and is worth double-checking. The other
    tools already accept these ids as substrate or product; call this one to see
    the match before trusting it.
    """
    def go(m: cobra.Model) -> dict:
        res = resolve_metabolite(m, query, online=online)
        mets = [{"id": met.id, "name": met.name, "formula": met.formula, "charge": met.charge,
                 "compartment": met.compartment,
                 "exchange": next((r.id for r in met.reactions if r in m.exchanges), None)}
                for met in res.metabolites]
        return {"model": m.id, "query": query, "ok": res.ok, "model_metabolite": res.model_id or None,
                "matched_via": res.via or None, "metabolites": mets, "candidates": res.candidates,
                "suggestions": res.suggestions, "caveat": res.caveat, "note": res.note or None,
                "error": None if res.ok else res.error(m.id)}
    return _clean(_run(model, go))


def _fluxes(m: cobra.Model, objective: str, min_growth: float) -> pd.Series:
    """Parsimonious FBA fluxes for growth, or for ``objective`` (a reaction id) with growth held at ``min_growth``."""
    with m:
        if objective != "growth":
            biomass = m.reactions.get_by_id(_biomass_id(m))
            biomass.lower_bound = min_growth
            m.objective = objective
        try:
            flux = pfba(m).fluxes
        except cobra.exceptions.Infeasible as err:
            raise ToolInputError(f"No feasible flux distribution for {objective} at growth >= {min_growth}/h. "
                                 "Lower min_growth or check the medium.") from err
    return flux.where(flux.abs() > 1e-9, 0.0)  # solver noise like -3e-15 reads as a backwards flux


def check_thermodynamics(product: str = "", substrate: str = "glucose", model: str = "e_coli_core",
                         uptake: float = 10.0, aerobic: bool = True, pathway: str = "",
                         min_growth: float = 0.0, all_reactions: bool = False) -> dict[str, Any]:
    """Does a flux route ask any reaction to run in a direction physics rules out?

    Estimates ΔG'° (kJ/mol) with eQuilibrator (needs the [thermo] extra, offline)
    and, across metabolite concentrations of 1 µM to 10 mM plus a 95% margin for
    the estimate's error, decides whether each reaction can only go forward,
    only in reverse, or either way.

    Without ``product`` it checks the route the cell uses to grow; with one, the
    route that makes the most ``product`` while growing at least ``min_growth``
    (1/h). The default 0 is a cell that has stopped growing (the stoichiometric
    best case); try about 0.05 for the route a growing production strain uses. ``pathway`` adds a non-native pathway first, as in simulate_production.
    Flux comes from parsimonious FBA, the shortest route.

    ``against_thermodynamics`` is the headline: active reactions carrying flux
    in a direction ruled out, i.e. the route needs something impossible. Rows
    whose ``method`` is "group contribution" rest on a weaker estimate and can
    be the estimate's fault rather than the model's. ``conflicts`` lists
    reactions the model only allows in the impossible direction, and
    ``can_tighten`` reversible ones that could be made one-way (pass
    thermo_bounds=true to simulate_production to do that). Transport, boundary
    and biomass reactions are skipped. By default only reactions on the route
    are estimated (seconds); ``all_reactions`` checks the whole model (about a
    minute on genome scale).
    """
    def go(m: cobra.Model) -> dict:
        added = _add_pathway(m, pathway)
        substrate_ex = set_medium(m, substrate, uptake, aerobic)
        target = resolve_exchange(m, product, create_demand=True) if product else "growth"
        if target == substrate_ex:
            raise ToolInputError("Substrate and product are the same compound.")
        if m.slim_optimize(error_value=0.0) < 1e-6:
            raise ToolInputError(f"{m.id} can't grow on {substrate_ex} under these conditions. "
                                 "Try the other oxygen setting or another substrate.")
        flux = _fluxes(m, target, min_growth)
        active = flux[flux.abs() > 1e-6]
        try:
            table = reaction_thermodynamics(m, reactions=None if all_reactions else list(active.index))
        except ThermoUnavailable as err:
            raise ToolInputError(f"Thermodynamics unavailable: {err}. Install the [thermo] extra "
                                 "(pip install 'cellfactory[thermo]').") from err
        table["flux"] = table["reaction"].map(flux).fillna(0.0)
        wrong = table[((table["thermo_direction"] == "forward") & (table["flux"] < -1e-6))
                      | ((table["thermo_direction"] == "reverse") & (table["flux"] > 1e-6))]
        cols = ["reaction", "name", "flux", "dg0_prime", "dg0_sigma", "dg_prime_min", "dg_prime_max",
                "thermo_direction", "model_direction", "method"]
        estimated = table["thermo_direction"].notna()
        out = {
            "model": m.id, "substrate": substrate_ex, "route_for": target,
            "units": "kJ/mol; flux mmol/gDW/h",
            "conditions": str(Conditions()),
            "summary": {"reactions_checked": len(table), "estimated": int(estimated.sum()),
                        "active_reactions": len(active),
                        "against_thermodynamics": len(wrong),
                        "against_thermodynamics_group_contribution":
                            int((wrong["method"] == "group contribution").sum()),
                        "conflicts": int((table["verdict"] == "CONFLICT").sum()),
                        "can_tighten": int((table["verdict"] == "can tighten").sum()),
                        "route_feasible": wrong.empty},
            "against_thermodynamics": _rows(wrong[cols]),
            "conflicts": _rows(table[table["verdict"] == "CONFLICT"][cols]),
            "can_tighten": _rows(table[table["verdict"] == "can tighten"][cols]),
            "unestimated": _rows(table[~estimated & ~table["note"].isin(["boundary reaction", "biomass pseudo-reaction"])]
                                 [["reaction", "name", "note"]], limit=20),
            "caveats": thermo_caveats(table),
            "note": "A direction the estimate allows is not a direction the enzyme runs fast; "
                    "this rules things out, it does not rule them in.",
        }
        if added:
            out["pathway"] = added
        return out
    return _clean(_run(model, go))


def qa_model(model: str = "e_coli_core", full: bool = False) -> dict[str, Any]:
    """Vet a model with MEMOTE before trusting anything it predicts.

    Hard failures (mass created from nothing, ATP or other energy carriers
    regenerated with no input) mean every yield from the model is suspect;
    simulate_production refuses such a model. Warnings and info rows are worth
    reading but do not block. ``full`` runs the whole MEMOTE suite (minutes on a
    genome-scale model) instead of the fast subset the gate uses. Results are
    cached per model checksum. If MEMOTE isn't installed, the built-in
    free-energy checks run instead and ``engine`` says so.
    """
    def go(m: cobra.Model) -> dict:
        res = memote_qa(m, full=full)
        table = res.table
        return {
            "model": m.id, "ok": res.ok, "engine": res.engine, "cached": res.cached, "checksum": res.checksum,
            "counts": {k: int(v) for k, v in table["gate"].value_counts().items()},
            "hard_failures": _rows(res.hard_failures),
            "warnings": _rows(table[table["gate"] == "warning"]),
            "info": _rows(table[table["gate"] == "info"], limit=20),
            "caveats": res.caveats,
        }
    return _clean(_run(model, go))


def flux_distribution(product: str = "", substrate: str = "glucose", model: str = "e_coli_core",
                      uptake: float = 10.0, aerobic: bool = True, pathway: str = "",
                      min_growth: float = 0.0, limit: int = 30) -> dict[str, Any]:
    """Where the carbon goes: one flux distribution, its exchanges and its busiest reactions.

    Without ``product``: the cell growing as fast as it can. With one: making
    as much ``product`` as possible while growing at least ``min_growth`` (1/h);
    the default 0 is a cell that has stopped growing, so try about 0.05 for a
    growing production strain.
    Parsimonious FBA picks the distribution with the least total flux, the usual
    stand-in for "the cell doesn't make enzymes it doesn't need"; the optimum is
    often not unique, so treat other equally good routes as possible.

    ``uptake`` and ``secretion`` are the exchange fluxes (what goes in and comes
    out, mmol/gDW/h). ``reactions`` are the ``limit`` internal reactions with the
    largest flux, with their equations and genes; a negative flux runs the
    reaction right to left. ``pathway`` adds a non-native pathway first.
    """
    if not 1 <= limit <= 200:
        raise ToolInputError("limit must be between 1 and 200.")

    def go(m: cobra.Model) -> dict:
        added = _add_pathway(m, pathway)
        substrate_ex = set_medium(m, substrate, uptake, aerobic)
        target = resolve_exchange(m, product, create_demand=True) if product else "growth"
        if target == substrate_ex:
            raise ToolInputError("Substrate and product are the same compound.")
        if m.slim_optimize(error_value=0.0) < 1e-6:
            raise ToolInputError(f"{m.id} can't grow on {substrate_ex} under these conditions. "
                                 "Try the other oxygen setting or another substrate.")
        flux = _fluxes(m, target, min_growth)
        boundary = {r.id for r in m.boundary}
        ex = flux[flux.index.isin(boundary) & (flux.abs() > 1e-6)]
        internal = flux[~flux.index.isin(boundary) & (flux.abs() > 1e-6)]
        top = internal.reindex(internal.abs().sort_values(ascending=False).index)

        def ex_rows(series: pd.Series) -> list[dict]:
            return [{"reaction": r, "metabolite": next(iter(m.reactions.get_by_id(r).metabolites)).id,
                     "name": m.reactions.get_by_id(r).name, "flux": v} for r, v in series.items()]

        biomass = _biomass_id(m)
        out = {
            "model": m.id, "substrate": substrate_ex, "optimised_for": target,
            "units": "mmol/gDW/h, uptake negative; growth 1/h",
            "growth": flux[biomass],
            "product_flux": flux[target] if product else None,
            "active_reactions": len(internal),
            "uptake": ex_rows(ex[ex < 0].sort_values()),
            "secretion": ex_rows(ex[ex > 0].sort_values(ascending=False)),
            "reactions": [{"reaction": r, "name": m.reactions.get_by_id(r).name, "flux": v,
                           "equation": m.reactions.get_by_id(r).build_reaction_string(use_metabolite_names=True),
                           "genes": sorted(g.id for g in m.reactions.get_by_id(r).genes)}
                          for r, v in top.head(limit).items() if r != biomass],
            "truncated": len(top) > limit,
        }
        if added:
            out["pathway"] = added
        return out
    return _clean(_run(model, go))


def _split(text: str) -> list[str]:
    return [t.strip() for t in text.split(",") if t.strip()]


def compare_conditions(product: str, substrates: str = "glucose", model: str = "e_coli_core",
                       oxygen: str = "both", uptake: float = 10.0, pathway: str = "") -> dict[str, Any]:
    """Yield of ``product`` across several substrates and oxygen settings, one row each.

    ``substrates`` is a comma-separated list ("glucose, glycerol, xylose");
    ``oxygen`` is "aerobic", "anaerobic" or "both". Each row is one
    simulate_production run: growth, guaranteed and theoretical yield (mol/mol
    substrate) and carbon yield (fraction of substrate carbon ending in product),
    which is the fair comparison between substrates of different size. A
    condition the cell cannot grow under comes back with ``error`` set instead of
    numbers. ``best`` is the row with the highest guaranteed carbon yield.
    At most 8 substrates.
    """
    subs = _split(substrates)
    if not 1 <= len(subs) <= 8:
        raise ToolInputError("Give between 1 and 8 comma-separated substrates.")
    settings = {"aerobic": [True], "anaerobic": [False], "both": [True, False]}.get(oxygen.lower())
    if settings is None:
        raise ToolInputError('oxygen must be "aerobic", "anaerobic" or "both".')

    def go(m: cobra.Model) -> dict:
        added = _add_pathway(m, pathway)
        rows = []
        for sub in subs:
            for aerobic in settings:
                row = {"substrate": sub, "aerobic": aerobic}
                try:
                    with m:
                        r = analyze(m, sub, product, uptake=uptake, aerobic=aerobic)
                    p = r.production
                    row.update(substrate_exchange=r.substrate, product=r.product, growth=p.growth,
                               guaranteed_yield=p.guaranteed_yield, best_case_yield=p.best_case_yield,
                               max_theoretical_yield=p.max_theoretical_yield,
                               guaranteed_carbon_yield=p.guaranteed_yield * r.carbon_ratio,
                               max_carbon_yield=p.max_theoretical_yield * r.carbon_ratio, error=None)
                except (KeyError, ValueError) as err:
                    row["error"] = err.args[0] if err.args else str(err)
                rows.append(row)
        ok = [r for r in rows if r["error"] is None]
        best = max(ok, key=lambda r: r["guaranteed_carbon_yield"]) if ok else None
        out = {"model": m.id, "product": product,
               "units": "growth 1/h; yields mol/mol substrate; carbon yields fraction of substrate carbon",
               "rows": rows, "best": best}
        if added:
            out["pathway"] = added
        return out
    return _clean(_run(model, go))


def _route(m: cobra.Model, product: str, substrate: str, uptake: float, aerobic: bool,
           min_growth: float) -> tuple[str, str, pd.Series]:
    """Set the medium, then (substrate exchange, route target, pFBA fluxes of the route's internal reactions)."""
    substrate_ex = set_medium(m, substrate, uptake, aerobic)
    target = resolve_exchange(m, product, create_demand=True) if product else "growth"
    if target == substrate_ex:
        raise ToolInputError("Substrate and product are the same compound.")
    if m.slim_optimize(error_value=0.0) < 1e-6:
        raise ToolInputError(f"{m.id} can't grow on {substrate_ex} under these conditions. "
                             "Try the other oxygen setting or another substrate.")
    flux = _fluxes(m, target, min_growth)
    skip = {r.id for r in m.boundary} | {_biomass_id(m)}
    return substrate_ex, target, flux[(flux.abs() > 1e-6) & ~flux.index.isin(skip)]


def check_driving_force(product: str = "", substrate: str = "glucose", model: str = "e_coli_core",
                        uptake: float = 10.0, aerobic: bool = True, pathway: str = "",
                        min_growth: float = 0.0) -> dict[str, Any]:
    """Can every step of a route run forward at once? The Max-min Driving Force (MDF) of the whole route.

    check_thermodynamics asks of each reaction alone whether some
    concentrations let it run. A route can pass that test step by step and
    still be impossible: two neighbouring steps may need the shared
    intermediate at opposite ends of its range. MDF (Noor et al. 2014) picks
    one concentration per metabolite, within 1 µM-10 mM, so that the least
    downhill step is as downhill as possible; that least driving force, in
    kJ/mol, is the MDF.

    ``verdict`` is "infeasible" (no concentrations in range let every step run
    the route's way, even giving each ΔG'° estimate the benefit of its error),
    "uncertain" (infeasible on the best estimates but not beyond their error),
    "low driving force" (feasible, but the bottleneck sits within about
    5.7 kJ/mol of equilibrium, where more than a tenth of its forward
    turnovers are undone by the reverse reaction, so the step needs more
    enzyme for the same net flux) or "feasible". ``bottlenecks`` are the steps holding the MDF
    down. Steps with no ΔG'° estimate are left out and listed, which can only
    make the true MDF lower, and cofactor ratios such as ATP:ADP are left free,
    so a positive MDF is optimistic.

    The route is picked as in check_thermodynamics: without ``product``, the one
    the cell uses to grow; with one, the route making the most product while
    growing at least ``min_growth`` (1/h, 0 being the stoichiometric best case).
    Needs the [thermo] extra (eQuilibrator, which works offline).
    """
    def go(m: cobra.Model) -> dict:
        added = _add_pathway(m, pathway)
        substrate_ex, target, route = _route(m, product, substrate, uptake, aerobic, min_growth)
        if route.empty:
            raise ToolInputError(f"The route for {target} has no internal reactions to screen.")
        try:
            res = _max_min_driving_force(m, route)
        except ThermoUnavailable as err:
            raise ToolInputError(f"Thermodynamics unavailable: {err}. Install the [thermo] extra "
                                 "(pip install 'cellfactory[thermo]').") from err
        except RuntimeError as err:  # the linear program itself failed
            raise ToolInputError(f"The MDF linear program did not solve: {err}.") from err
        steps = res.reactions.sort_values("driving_force", na_position="last")
        conc = res.concentrations
        out = {
            "model": m.id, "substrate": substrate_ex, "route_for": target,
            "units": "kJ/mol (driving force is -ΔG' the way the route runs); flux mmol/gDW/h; "
                     "concentration mol/L",
            "conditions": str(Conditions()),
            "mdf": res.mdf, "mdf_optimistic": res.mdf_optimistic,
            "verdict": res.verdict, "reason": res.reason, "bottlenecks": res.bottlenecks,
            "summary": {"steps": len(res.reactions), "estimated": int(res.reactions["dg0_prime"].notna().sum()),
                        "excluded": len(res.excluded), "bottlenecks": len(res.bottlenecks),
                        "route_feasible": res.verdict in ("feasible", "low driving force")},
            "steps": _rows(steps, limit=MAX_ROWS * 4),
            "concentrations": _rows(conc.sort_values("concentration", ascending=False) if len(conc) else conc,
                                    limit=MAX_ROWS * 4),
            "caveats": mdf_caveats(res, f"the route for {target}"),
        }
        if added:
            out["pathway"] = added
        return out
    return _clean(_run(model, go))


def _budgeted_lookup(budget: float = SABIO_BUDGET):
    """sabiork_lookup, but the whole call gets ``budget`` seconds.

    Each request already has its own timeout, yet a server that answers slowly
    rather than failing could keep a route's few dozen EC numbers going for
    minutes while this call holds the model's lock. Past the budget the lookup
    raises the same error an unreachable server does, which kinetics turns into
    its "SABIO-RK (partial)" note: the EC numbers already fetched are kept.
    """
    deadline = time.monotonic() + budget

    def lookup(ec: str):
        if time.monotonic() > deadline:
            raise KineticsUnavailable(f"the {budget:g} s budget for SABIO-RK lookups ran out")
        return sabiork_lookup(ec)
    return lookup


def check_enzyme_kinetics(product: str = "", substrate: str = "glucose", model: str = "e_coli_core",
                          uptake: float = 10.0, aerobic: bool = True, pathway: str = "",
                          min_growth: float = 0.0, organism: str = "", online: bool = True,
                          limit: int = 50) -> dict[str, Any]:
    """Turnover number (kcat, 1/s) and Km (mM) of each enzyme on a route, with where every number came from.

    FBA says which reactions run, not how much enzyme running them takes. Flux
    divided by kcat is the least enzyme a step needs, so a slow enzyme on a
    busy step is where a real strain is likely to stall. Reactions are matched
    to enzyme data by EC number; a model with no EC annotations (e_coli_core)
    borrows them from iJO1366 by reaction id.

    Sources, most trusted first: BRENDA's flat file when the server has one
    (its licence forbids bundling it, so set ``CELLFACTORY_BRENDA`` to the
    path), then SABIO-RK's web service (``online``; best effort, and offline or
    blocked it becomes a note rather than an error), and finally a predictor's
    kcat table (DLKcat, TurNuP) if the server sets ``CELLFACTORY_KCAT_PREDICTIONS``
    to one, used only where nothing was measured. ``value_type`` per row is
    "measured", "predicted", "missing", or "no enzyme" for spontaneous,
    diffusion and pseudo-reactions. ``organism_match`` says whether the kcat
    was measured in the host, the same genus or another organism; ``organism``
    names the host species ("Escherichia coli") and is guessed from well-known
    model ids when empty.

    The route is picked as in check_thermodynamics: without ``product``, the one
    the cell uses to grow; with one, the route making the most product while
    growing at least ``min_growth`` (1/h). Rows come back busiest first.
    Measured kcats are test-tube values; inside a cell enzymes usually run
    below them. One SABIO-RK query per EC number the first time (cached on
    disk afterwards), and they are given a minute in all: EC numbers past that
    are left unlooked-up and said so in the notes, rather than holding the
    model for a slow server.
    """
    if not 1 <= limit <= 500:
        raise ToolInputError("limit must be between 1 and 500.")

    def go(m: cobra.Model) -> dict:
        added = _add_pathway(m, pathway)
        substrate_ex, target, route = _route(m, product, substrate, uptake, aerobic, min_growth)
        if route.empty:
            raise ToolInputError(f"The route for {target} has no internal reactions to look up.")
        predictions = os.environ.get("CELLFACTORY_KCAT_PREDICTIONS") or None
        try:
            res = _enzyme_kinetics(m, reactions=route, organism=organism or None, online=online,
                                   predictions=predictions, lookup=_budgeted_lookup() if online else None)
        except ValueError as err:
            if predictions is None:
                raise
            raise ToolInputError(f"CELLFACTORY_KCAT_PREDICTIONS ({predictions}) could not be read: {err}") from err
        table = res.table.copy()
        if not table.empty:
            table.insert(1, "flux", table["reaction"].map(route))
            table = table.sort_values("flux", key=lambda f: f.abs(), ascending=False)
        enz = table[table["value_type"] != "no enzyme"] if not table.empty else table
        counts = enz["value_type"].value_counts() if len(enz) else {}
        out = {
            "model": m.id, "substrate": substrate_ex, "route_for": target, "organism": res.organism,
            "units": "kcat 1/s; Km mM; flux mmol/gDW/h",
            "sources": res.sources or [],
            "summary": {"reactions": len(table), "enzyme_steps": len(enz),
                        "measured": int(counts.get("measured", 0)),
                        "predicted": int(counts.get("predicted", 0)),
                        "missing": int(counts.get("missing", 0)),
                        "measured_in_host": int((enz["organism_match"] == "host").sum()) if len(enz) else 0},
            "reactions": _rows(table, limit=limit),
            "caveats": kinetics_caveats(res, f"the route for {target}"),
        }
        if added:
            out["pathway"] = added
        return out
    return _clean(_run(model, go))


def design_knockouts(product: str, substrate: str = "glucose", model: str = "e_coli_core",
                     method: str = "optknock", max_knockouts: int = 3, min_growth: float = 0.05,
                     min_yield: float = 0.0, max_designs: int = 5, time_limit: float = 120.0,
                     uptake: float = 10.0, aerobic: bool = True, pathway: str = "") -> dict[str, Any]:
    """Solve for a *combination* of deletions that forces the product out, with StrainDesign (MILP).

    search_knockouts deletes one reaction at a time and keeps the best, so it
    misses designs where no single deletion helps but two together do: closing
    one competing route just sends the carbon down another. StrainDesign
    (Schneider et al. 2022) writes the whole question as one mixed-integer
    linear program and searches every combination at once.

    ``method`` picks how strict the coupling has to be:

    * ``optknock``: most product at the mutant's fastest growth, assuming the
      cell picks the flux state we like among the equally fast ones, so a
      design can be coupled only on paper.
    * ``robustknock``: most *guaranteed* product at fastest growth, the number
      the rest of this package reports.
    * ``mcs`` (minimal cut sets): the fewest deletions after which *every*
      state growing at least ``min_growth`` yields at least ``min_yield``
      mol/mol (default: 30% of the theoretical maximum). This is strong growth
      coupling.

    Every design is re-checked with the same FBA as the other tools and
    translated to genes: ``genes`` is the smallest gene set whose deletion
    removes the reactions, ``also_disabled`` the other reactions those genes
    take down with them, and the gene-level columns the yield once the genes,
    not the reactions, are gone. ``coupling`` is "strong" (forced at every
    growth rate), "at max growth", "none" or "too slow". Rows are ranked by
    gene-level guaranteed yield.

    Needs the [design] extra (StrainDesign with a MILP solver). ``time_limit``
    (seconds, at most an hour) caps the solver; a design that hits the limit is
    a good one, not provably the best. The solve holds the model while it runs,
    so a genome-scale model is better given a smaller candidate problem than an
    open-ended limit.
    """
    if method not in DESIGN_METHODS:
        raise ToolInputError(f"Unknown method {method!r}. Choose one of: {', '.join(DESIGN_METHODS)}.")
    if not 1 <= max_knockouts <= 10:
        raise ToolInputError("max_knockouts must be between 1 and 10.")
    if not 1 <= max_designs <= 20:
        raise ToolInputError("max_designs must be between 1 and 20.")
    if not 0 < time_limit <= MAX_TIME_LIMIT:
        raise ToolInputError(f"time_limit must be between 0 and {MAX_TIME_LIMIT} seconds. "
                             "The solve holds the model while it runs, so it cannot be left open-ended.")

    def go(m: cobra.Model) -> dict:
        added = _add_pathway(m, pathway)
        substrate_ex = set_medium(m, substrate, uptake, aerobic)
        target = resolve_exchange(m, product, create_demand=True)
        if target == substrate_ex:
            raise ToolInputError("Substrate and product are the same compound.")
        if m.slim_optimize(error_value=0.0) < 1e-6:
            raise ToolInputError(f"{m.id} can't grow on {substrate_ex} under these conditions. "
                                 "Try the other oxygen setting or another substrate.")
        try:
            res = _design_knockouts(m, target, substrate_ex, method=method, max_knockouts=max_knockouts,
                                    min_growth=min_growth, min_yield=min_yield or None,
                                    max_designs=max_designs, time_limit=time_limit)
        except DesignUnavailable as err:
            raise ToolInputError(f"Knockout design unavailable: {err}. Install the [design] extra "
                                 "(pip install 'cellfactory[design]').") from err
        best = res.best
        out = {
            "model": m.id, "substrate": substrate_ex, "product": target,
            "method": res.method, "solver": res.solver, "status": res.status,
            "units": "growth 1/h; yields mol/mol substrate",
            "wild_type_guaranteed_yield": res.wild_type_yield,
            "min_growth": res.min_growth, "min_yield": res.min_yield,
            "summary": {"designs": len(res.designs),
                        "best_guaranteed_yield": None if best is None else best["guaranteed_yield"],
                        "best_knockouts": None if best is None else best["knockouts"],
                        "best_genes": None if best is None else best["genes"],
                        "best_coupling": None if best is None else best["coupling"],
                        "improved": best is not None and best["guaranteed_yield"] > res.wild_type_yield + 1e-6},
            "designs": _rows(res.designs, limit=MAX_ROWS),
            "notes": res.notes,
            "caveats": design_caveats(res),
        }
        if added:
            out["pathway"] = added
        return out
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
    resolve_compound,
    check_thermodynamics,
    check_driving_force,
    check_enzyme_kinetics,
    design_knockouts,
    qa_model,
    flux_distribution,
    compare_conditions,
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
    "resolve_metabolite": "resolve_compound",
    "Resolution": "resolve_compound",
    "reaction_thermodynamics": "check_thermodynamics",
    "ThermoUnavailable": "check_thermodynamics",
    "apply_directionality": "simulate_production",
    "memote_qa": "qa_model",
    "QAResult": "qa_model",
    "enzyme_kinetics": "check_enzyme_kinetics",
    "KineticsResult": "check_enzyme_kinetics",
    "max_min_driving_force": "check_driving_force",
    "MDFResult": "check_driving_force",
    "design_knockouts": "design_knockouts",
    "DesignResult": "design_knockouts",
    "DesignUnavailable": "design_knockouts",
}
NOT_EXPOSED = {
    "load_kcat_predictions": "Reads a local DLKcat/TurNuP output file for enzyme_kinetics; the file lives on the "
                             "server, so check_enzyme_kinetics reads its path from CELLFACTORY_KCAT_PREDICTIONS "
                             "rather than taking a path from the caller.",
    "gene_essentiality_check": "Only meaningful for E. coli with its built-in lab data; needs gene lists "
                               "as input. Expose once it takes another organism's screen data.",
    "set_glucose_medium": "Shortcut for set_medium; every tool already takes a substrate.",
    "parse_identifier": "String helper behind resolve_metabolite; no use to an agent on its own.",
    "product_route_fluxes": "Picks the zero-growth max-yield route for max_min_driving_force; "
                            "check_driving_force and check_enzyme_kinetics pick the route themselves, from the same "
                            "parsimonious FBA as check_thermodynamics, so min_growth means the same in all three.",
    "ModelQAError": "Raised by analyze's QA gate; simulate_production and search_knockouts turn it into a "
                    "ToolInputError naming the failed check.",
}

"""Nearest native precursor search: the fewest new enzymes that reach a target compound.

The question this answers is "I want the cell to make Z, what does it already
make that is closest to Z?". Metabolic engineers call it retrobiosynthesis.
Rather than designing a pathway by hand, you search the chemistry backwards
from the product until you land on something the host organism already has.

How it works:

* **The graph.** Nodes are metabolites, edges are reactions. Edges come from
  two places: the *host* model, which are reactions the organism already runs
  and so cost nothing, and a *donor* reaction set, which are reactions it does
  not have and so cost one enzyme each to clone in. Any other COBRApy model
  works as the donor (``iJO1366`` is a good donor for ``e_coli_core``), which
  keeps the whole search offline. A MetaNetX, Rhea or ATLAS adapter can later
  supply donor reactions through the same interface.
* **Hub metabolites are cut out.** Water, protons, ATP, CoA and the other
  cofactors take part in hundreds of reactions, so leaving them in makes every
  compound two steps from every other compound through ATP. They are dropped
  as connections (the reactions still appear in the route, they just cannot be
  *entered* through a cofactor).
* **The search.** Dijkstra backwards from the target, with cost
  ``(heterologous steps, total steps)``, stopping at the first metabolite the
  host can actually make.
* **"Already makes it" means FBA says so.** A metabolite that is listed in the
  host model but cannot be produced on your carbon source is not a usable
  starting point, so each candidate is checked by maximising a demand reaction
  on it under the chosen medium.

Caveats worth saying out loud. The graph connects a reaction's substrates to
its products without tracking which atoms actually move, so a route can in
principle hop through a pair that shares no carbon; the cofactor filter removes
most of these but not all. Step count is not difficulty: one step can be a
cytochrome P450 that will not fold in *E. coli*. And matching metabolites
between host and donor is done on identifiers, so both must use the same
namespace (BiGG with BiGG). Cross-namespace mapping is a separate job, tracked
as MetaNetX identifier resolution in the pathway-tools backlog.
"""

from __future__ import annotations

import heapq
import logging
from dataclasses import asdict, dataclass, field

import cobra

from .models import ALIASES, load_model, set_medium

log = logging.getLogger(__name__)

# Cofactors and inorganics, by BiGG id without the compartment suffix. These are
# carriers, not precursors: reaching ATP tells you nothing about reaching a
# terpene. Acetyl-CoA, pyruvate and the other real building blocks are kept.
CURRENCY = frozenset("""
    h h2o o2 co2 hco3 pi ppi pppi nh4 so4 so3 h2s h2o2 o2s photon
    atp adp amp gtp gdp gmp ctp cdp cmp utp udp ump itp idp imp datp dadp damp
    nad nadh nadp nadph fad fadh2 fmn fmnh2 coa acp
    thf mlthf methf 10fthf 5mthf dhf
    q8 q8h2 mqn8 mql8 2dmmq8 2dmmql8 trdox trdrd grxox grxrd
    na1 k cl ca2 mg2 mn2 fe2 fe3 zn2 cu cu2 ni2 cobalt2 mobd cd2 hg2 ag tungs sel
""".split())

# A reaction with more participants than this is a lumped pseudo-reaction
# (biomass, an ATP maintenance sink): real chemistry, no. It would connect
# everything to everything.
MAX_PARTICIPANTS = 20

# BiGG's stand-in gene for spontaneous or unassigned reactions. There is nothing
# to delete in the lab, so it must never appear in a knockout shortlist.
PLACEHOLDER_GENES = frozenset({"s0001"})


def _is_transport(rxn: cobra.Reaction) -> bool:
    """True if the reaction only moves a compound between compartments."""
    by_compartment: dict[str, set[str]] = {}
    for met in rxn.metabolites:
        by_compartment.setdefault(met.compartment, set()).add(_base_id(met.id))
    if len(by_compartment) < 2:
        return False
    shared = set.intersection(*by_compartment.values())
    return bool(shared)


@dataclass(frozen=True)
class Step:
    """One reaction in a route, running in the direction that reaches the target."""

    reaction: str
    name: str
    equation: str
    native: bool  # the host already has this reaction
    genes: tuple[str, ...]  # for a non-native step, the donor's genes for it
    substrate: str  # the metabolite entered from
    product: str  # the metabolite reached
    also_needs: tuple[str, ...] = ()  # the step's other inputs, all of them host-made


@dataclass
class Route:
    """A way to reach ``target`` from a metabolite the host already makes."""

    target: str
    target_name: str
    precursor: str
    precursor_name: str
    heterologous_steps: int  # enzymes you would have to add
    total_steps: int
    new_reactions: list[str]
    new_genes: list[str]  # donor genes for the new reactions, as a cloning shortlist
    steps: list[Step]
    precursor_max_flux: float | None = None  # mmol/gDW/h the host can make of it, or None if unchecked
    carbon_delta: int | None = None  # target carbons minus precursor carbons

    def to_dict(self) -> dict:
        out = asdict(self)
        out["steps"] = [{**asdict(s), "genes": list(s.genes)} for s in self.steps]
        return out


@dataclass
class _Edge:
    reaction: cobra.Reaction
    native: bool
    substrate: str
    product: str
    co_substrates: tuple[str, ...] = ()  # the reaction's other inputs in this direction


@dataclass
class _Graph:
    """Reverse adjacency: ``producers[m]`` are the edges that make ``m``."""

    producers: dict[str, list[_Edge]] = field(default_factory=dict)
    metabolites: dict[str, cobra.Metabolite] = field(default_factory=dict)
    host_metabolites: set[str] = field(default_factory=set)


def _base_id(met_id: str) -> str:
    """``akg_c`` -> ``akg``. BiGG ids end in a one- or two-letter compartment tag."""
    head, _, tail = met_id.rpartition("_")
    return head if head and len(tail) <= 2 else met_id


def _usable_reactions(model: cobra.Model) -> list[cobra.Reaction]:
    return [
        r for r in model.reactions
        if not r.boundary
        and len(r.metabolites) <= MAX_PARTICIPANTS
        and not r.objective_coefficient
        and "BIOMASS" not in r.id.upper()
    ]


def _hub_metabolites(reactions: list[cobra.Reaction], max_degree: int | None) -> set[str]:
    """Metabolite ids too well connected to be informative, plus the cofactor list."""
    degree: dict[str, int] = {}
    for rxn in reactions:
        for met in rxn.metabolites:
            degree[met.id] = degree.get(met.id, 0) + 1
    if max_degree is None:
        max_degree = max(25, round(0.05 * len(reactions)))
    return {m for m, d in degree.items() if d > max_degree or _base_id(m) in CURRENCY}


def build_graph(host: cobra.Model, donor: cobra.Model | None = None,
                max_degree: int | None = None) -> _Graph:
    """Reverse reaction graph over the host plus whatever the donor adds.

    Edges run product -> substrate, because the search runs backwards from the
    target. Hub metabolites are kept out of the edges but stay in the printed
    reaction equations.
    """
    host_reactions = _usable_reactions(host)
    host_reaction_ids = {r.id for r in host.reactions}
    donor_reactions = [r for r in _usable_reactions(donor) if r.id not in host_reaction_ids] if donor else []

    graph = _Graph(host_metabolites={m.id for m in host.metabolites})
    hubs = _hub_metabolites(host_reactions + donor_reactions, max_degree)

    for reactions, native in ((host_reactions, True), (donor_reactions, False)):
        for rxn in reactions:
            consumed = [m for m, coeff in rxn.metabolites.items() if coeff < 0]
            produced = [m for m, coeff in rxn.metabolites.items() if coeff > 0]
            directions = []
            if rxn.upper_bound > 0:
                directions.append((consumed, produced))
            if rxn.lower_bound < 0:
                directions.append((produced, consumed))
            for substrates, products in directions:
                real_inputs = [m.id for m in substrates if m.id not in hubs]
                for product in products:
                    if product.id in hubs:
                        continue
                    for substrate in substrates:
                        if substrate.id in hubs or substrate.id == product.id:
                            continue
                        graph.metabolites.setdefault(substrate.id, substrate)
                        graph.metabolites.setdefault(product.id, product)
                        graph.producers.setdefault(product.id, []).append(_Edge(
                            rxn, native, substrate.id, product.id,
                            co_substrates=tuple(m for m in real_inputs if m != substrate.id),
                        ))
    return graph


def resolve_metabolite(graph: _Graph, query: str) -> str:
    """Turn a user's word for a compound into a metabolite id present in the graph.

    Accepts an exact id (``lyco_c``), an id without its compartment (``lyco``),
    an everyday name (``lycopene``) or an exact metabolite name. Cytosolic
    metabolites win ties, since that is where a new pathway is expressed.
    """
    candidates = [query, ALIASES.get(query.lower(), query)]
    for candidate in candidates:
        if candidate in graph.metabolites:
            return candidate
        matches = [m for m in graph.metabolites if _base_id(m) == candidate]
        if matches:
            return min(matches, key=lambda m: (not m.endswith("_c"), m))
    lowered = query.lower()
    by_name = [m for m, met in graph.metabolites.items() if (met.name or "").lower() == lowered]
    if by_name:
        return min(by_name, key=lambda m: (not m.endswith("_c"), m))
    near = sorted({met.name for m, met in graph.metabolites.items()
                   if lowered in m.lower() or lowered in (met.name or "").lower()})[:8]
    hint = ("Did you mean: " + ", ".join(near) + "?") if near else (
        "It is in neither the host nor the donor reaction set. Try a larger donor model.")
    raise KeyError(f"{query!r} is not a metabolite in the reaction graph. {hint}")


class _Producibility:
    """Can the host actually make this metabolite, on this medium? Cached, FBA-backed.

    The obvious test — hang a demand reaction on the metabolite and maximise it
    — is wrong for anything carrying a conserved moiety. Draining succinyl-CoA
    also drains the cell's CoA, which nothing in a central-metabolism model
    replaces, so the demand test reports zero for a metabolite the cell plainly
    makes. Instead we ask the steady-state question: can any reaction that
    produces this metabolite carry flux? Steady state keeps the CoA pool
    balanced, so the moiety never leaks.
    """

    TOL = 1e-6
    MAX_REACTIONS_PER_METABOLITE = 25  # bound the cost for hub-ish metabolites

    def __init__(self, model: cobra.Model, enabled: bool = True, max_checks: int = 4000):
        self.model = model
        self.enabled = enabled
        self.max_checks = max_checks
        self.checks = 0
        self._can_make: dict[str, float | None] = {}
        self._max_rate: dict[str, float | None] = {}

    def flux(self, met_id: str) -> float | None:
        """Most mmol/gDW/h of ``met_id`` the host can make, or None if not checked.

        Unlike :meth:`is_native` this scans every producing reaction, so the
        number does not depend on which one happened to be tried first. Only
        worth paying for the routes actually reported.
        """
        if met_id not in self._max_rate:
            self._max_rate[met_id] = self._measure(met_id, stop_early=False)
        return self._max_rate[met_id]

    def _measure(self, met_id: str, stop_early: bool = True) -> float | None:
        if not self.enabled or self.checks >= self.max_checks:
            return None
        if met_id not in self.model.metabolites:
            return 0.0
        met = self.model.metabolites.get_by_id(met_id)
        makers = [r for r in met.reactions if not r.boundary and not r.objective_coefficient]
        best = 0.0
        for rxn in makers[:self.MAX_REACTIONS_PER_METABOLITE]:
            coeff = rxn.metabolites[met]
            if coeff > 0 and rxn.upper_bound <= 0:
                continue
            if coeff < 0 and rxn.lower_bound >= 0:
                continue  # the reaction can only consume it
            self.checks += 1
            with self.model:
                self.model.objective = rxn
                self.model.objective.direction = "max" if coeff > 0 else "min"
                value = self.model.slim_optimize()
            if value is None or value != value:  # infeasible, or NaN
                continue
            rate = abs(coeff) * (value if coeff > 0 else -value)
            best = max(best, rate)
            if stop_early and best > self.TOL:
                break
        return best

    def is_native(self, met_id: str) -> bool:
        """True if the host can make it at all. Stops at the first reaction that works.

        A metabolite left unchecked because the budget ran out counts as native
        if the host models it, which is the lenient direction: better to offer
        a precursor the user can rule out than to hide a real one.
        """
        if met_id not in self.model.metabolites:
            return False
        if met_id not in self._can_make:
            self._can_make[met_id] = self._measure(met_id, stop_early=True)
        rate = self._can_make[met_id]
        return True if rate is None else rate > self.TOL


def _carbons(met: cobra.Metabolite) -> int | None:
    try:
        return int(met.elements.get("C", 0)) or None
    except (TypeError, ValueError):  # a formula the model never filled in
        return None


def _route(graph: _Graph, target: str, precursor: str,
           came_from: dict[str, tuple[str, _Edge]]) -> Route:
    """Walk the predecessor chain from precursor back out to the target."""
    steps: list[Step] = []
    node = precursor
    while node != target:
        product, edge = came_from[node]
        rxn = edge.reaction
        steps.append(Step(
            reaction=rxn.id,
            name=rxn.name or rxn.id,
            equation=rxn.build_reaction_string(),
            native=edge.native,
            genes=tuple(sorted(g.id for g in rxn.genes)) if not edge.native else (),
            substrate=node,
            product=product,
            also_needs=edge.co_substrates,
        ))
        node = product

    new = [s for s in steps if not s.native]
    target_c = _carbons(graph.metabolites[target])
    precursor_c = _carbons(graph.metabolites[precursor])
    return Route(
        target=target,
        target_name=graph.metabolites[target].name or target,
        precursor=precursor,
        precursor_name=graph.metabolites[precursor].name or precursor,
        heterologous_steps=len(new),
        total_steps=len(steps),
        new_reactions=[s.reaction for s in new],
        new_genes=sorted({g for s in new for g in s.genes}),
        steps=steps,
        carbon_delta=(target_c - precursor_c) if target_c and precursor_c else None,
    )


def nearest_precursors(
    host: cobra.Model,
    target: str,
    donor: cobra.Model | str | None = "iJO1366",
    *,
    substrate: str = "glucose",
    uptake: float = 10.0,
    aerobic: bool = True,
    max_heterologous_steps: int = 4,
    max_total_steps: int = 8,
    top: int = 10,
    max_degree: int | None = None,
    verify_with_fba: bool = True,
) -> list[Route]:
    """Metabolites the host already makes, ranked by how few new enzymes reach ``target``.

    ``donor`` supplies the reactions the host does not have: another COBRApy
    model, a model id or path that :func:`cellfactory.load_model` understands,
    or None to search the host's own chemistry only.

    The search is run under ``substrate``/``uptake``/``aerobic``, because
    whether a precursor is available depends on what the cell is being fed.
    Both models must use the same identifier namespace.

    Returns at most ``top`` routes, best first: fewest heterologous steps, then
    most of the precursor available, then closest in carbon count.
    """
    if max_heterologous_steps < 0 or max_total_steps < 1:
        raise ValueError("max_heterologous_steps must be >= 0 and max_total_steps >= 1.")
    if isinstance(donor, str):
        donor = load_model(donor)

    with host:
        set_medium(host, substrate, uptake, aerobic)
        graph = build_graph(host, donor, max_degree=max_degree)
        target_id = resolve_metabolite(graph, target)
        native = _Producibility(host, enabled=verify_with_fba)

        if native.is_native(target_id):
            route = Route(
                target=target_id,
                target_name=graph.metabolites[target_id].name or target_id,
                precursor=target_id,
                precursor_name=graph.metabolites[target_id].name or target_id,
                heterologous_steps=0,
                total_steps=0,
                new_reactions=[],
                new_genes=[],
                steps=[],
                precursor_max_flux=native.flux(target_id),
                carbon_delta=0,
            )
            return [route]

        routes = _search(graph, target_id, native, max_heterologous_steps, max_total_steps, top)
        # Measured inside the block: how much of a precursor the host can make
        # depends on the medium, which is rolled back when the block exits.
        for route in routes:
            route.precursor_max_flux = native.flux(route.precursor)

    routes.sort(key=lambda r: (r.heterologous_steps, -(r.precursor_max_flux or 0.0),
                               abs(r.carbon_delta or 99), r.precursor))
    return routes[:top]


def _search(graph: _Graph, target: str, native: _Producibility,
            max_het: int, max_total: int, top: int) -> list[Route]:
    """Dijkstra backwards from ``target``, stopping each branch at a native metabolite."""
    came_from: dict[str, tuple[str, _Edge]] = {}
    best: dict[str, tuple[int, int]] = {target: (0, 0)}
    queue: list[tuple[int, int, str]] = [(0, 0, target)]
    found: list[Route] = []

    while queue:
        het, total, node = heapq.heappop(queue)
        if (het, total) > best.get(node, (het, total)):
            continue  # a cheaper path to this node was already expanded
        if node != target and native.is_native(node):
            found.append(_route(graph, target, node, came_from))
            if len(found) >= top * 3:  # enough to rank well without exploring the whole network
                break
            continue  # a route ends here; going further back would not improve it
        if total >= max_total:
            continue
        arrived_by = came_from[node][1].reaction.id if node in came_from else None
        for edge in graph.producers.get(node, ()):
            cost = (het + (0 if edge.native else 1), total + 1)
            if cost[0] > max_het or cost[1] > max_total:
                continue
            # A reversible reaction can otherwise be used in both directions in
            # a row, which "makes" the target out of itself.
            if edge.reaction.id == arrived_by:
                continue
            # Every other input of the step has to be something the host
            # already makes, or the step cannot run once it is cloned in.
            if not all(native.is_native(m) for m in edge.co_substrates):
                continue
            if cost < best.get(edge.substrate, (max_het + 1, max_total + 1)):
                best[edge.substrate] = cost
                came_from[edge.substrate] = (node, edge)
                heapq.heappush(queue, (cost[0], cost[1], edge.substrate))
    return found


# Pools a cell holds that attack a product chemically, with no enzyme involved.
# Each is checked against what the host actually produces before it is reported.
REACTIVE_POOLS = {
    "gthrd": ("thiol", "Reduced glutathione. E. coli holds it at 5-10 mM precisely to mop up "
                       "electrophiles, so an electrophilic product walks into a wall of it and "
                       "leaves as a glutathione conjugate."),
    "coa": ("thiol", "Free coenzyme A. Another cytosolic thiol that adds to electrophilic centres."),
    "cys__L": ("thiol", "Free cysteine, and the cysteine side chains all over the proteome."),
    "nh4": ("amine", "Free ammonium. Condenses with aldehydes and ketones to imines."),
    "lys__L": ("amine", "Free lysine, and the lysine side chains on every protein. Aldehydes form "
                        "Schiff bases with these, which is a large part of why aldehyde products "
                        "are hard to accumulate."),
    "h2o2": ("oxidant", "Hydrogen peroxide from aerobic metabolism. Oxidises conjugated systems "
                        "and thiols."),
    "o2s": ("oxidant", "Superoxide. Same story, faster."),
}

# Functional groups that make a product chemically fragile. BiGG models carry no
# structures, so this reads the metabolite's name: a hint to check, not a finding.
GROUP_HINTS = (
    (("aldehyde", "anal", "-al "), "aldehyde",
     "Aldehydes condense with amines to Schiff bases and oxidise readily. The most reactive "
     "product class there is; expect losses unless it is consumed or exported fast."),
    (("lactone", "olide"), "lactone",
     "Lactones hydrolyse on their own at neutral to alkaline pH, and cells are full of "
     "promiscuous esterases that speed it up."),
    (("ester",), "ester", "Same hydrolysis risk as a lactone, enzymatic and spontaneous."),
    (("epoxide", "oxirane"), "epoxide",
     "A strong electrophile. Glutathione will find it."),
    (("quinone",), "quinone", "Redox active and electrophilic; conjugates with thiols."),
    (("carotene", "lycopene", "phytoene", "retinal"), "polyene",
     "Long conjugated chains are oxidised by reactive oxygen species. This is a known reason "
     "carotenoid titers fall off under strong aeration."),
)


def _charge_note(charge: int | None) -> str:
    if charge is None:
        return "The model does not give this metabolite a charge, so nothing can be said about it."
    if charge <= -1:
        return (
            f"Charge {charge} at cytosolic pH, so the product exists as an anion, a salt rather than "
            "the free acid. That is mostly good: an ion does not dissolve into the membrane, so it "
            "avoids the solvent toxicity that limits ethanol and butanol, and neutralised fermentation "
            "is how organic acids are produced industrially. The cost is that an ion cannot diffuse "
            "out either, so it accumulates inside until it stalls the pathway or poisons the cell. "
            "Expect to need a transporter, and expect to pay downstream to acidify the salt back."
        )
    if charge >= 1:
        return (
            f"Charge {charge} at cytosolic pH, so the product is a cation and cannot cross the "
            "membrane unaided. Plan on a transporter."
        )
    return (
        "Neutral at cytosolic pH, so the product can cross the membrane on its own and may not need "
        "a transporter. The flip side is that a neutral lipophilic compound partitions into the "
        "membrane and disorders it, which is the toxicity that caps alcohol and solvent titers."
    )


def product_fate(
    host: cobra.Model,
    target: str,
    donor: cobra.Model | str | None = None,
    *,
    substrate: str = "glucose",
    uptake: float = 10.0,
    aerobic: bool = True,
    min_growth_fraction: float = 0.1,
) -> dict:
    """What the host would do to ``target`` once it is made: eat it, or react with it.

    Two halves, and they are not equally solid.

    The enzymatic half comes out of the model and is as good as the model is:
    every host reaction that consumes the target, the genes behind it, whether
    it carries flux while the cell grows, and whether deleting it still leaves
    the cell viable. Those deletions are the ones that protect your product.

    The chemical half is heuristic and labelled as such. It reports the
    product's charge at cytosolic pH, which decides whether it can leave the
    cell at all, and which reactive pools the host is actually producing under
    these conditions, checked by FBA rather than assumed. Functional-group
    risks are guessed from the metabolite's name, because genome-scale models
    carry formulas but not structures.

    ``donor`` is only needed for a product the host cannot make yet: without it
    the target's name cannot be looked up. A target absent from the host model
    is the good case, and says so — no native enzyme touches it.
    """
    if isinstance(donor, str):
        donor = load_model(donor)
    with host:
        set_medium(host, substrate, uptake, aerobic)
        graph = build_graph(host, donor)
        graph.metabolites.update({m.id: m for m in host.metabolites})
        target_id = resolve_metabolite(graph, target)
        met = graph.metabolites[target_id]
        native = _Producibility(host)

        consumers = _consumers(host, target_id, min_growth_fraction) if target_id in host.metabolites else []
        pools = []
        for base, (kind, note) in REACTIVE_POOLS.items():
            for met_id in (f"{base}_c", f"{base}_e"):
                flux = native.flux(met_id)
                if flux and flux > _Producibility.TOL:
                    pools.append({"metabolite": met_id, "kind": kind,
                                  "max_flux": round(flux, 4), "why_it_matters": note})
                    break

    name = (met.name or target_id).lower()
    groups = [{"group": group, "why_it_matters": note}
              for keywords, group, note in GROUP_HINTS if any(k in name for k in keywords)]

    return {
        "target": target_id,
        "target_name": met.name or target_id,
        "formula": met.formula,
        "charge": met.charge,
        "charge_note": _charge_note(met.charge),
        "in_host_model": target_id in host.metabolites,
        "consumed_by": consumers,
        "protective_knockouts": [c["reaction"] for c in consumers if c["safe_to_delete"]],
        "reactive_pools_the_host_makes": pools,
        "functional_group_hints": groups,
        "caveats": [
            "Reactions and knockout effects come from the model. Reactive pools are measured by FBA, "
            "but whether they actually attack your product is chemistry the model knows nothing about.",
            "Functional groups are guessed from the metabolite's name; genome-scale models carry no "
            "structures. Treat them as things to check, not findings.",
            "An acid meeting a base is not destruction. It is neutralisation, and the product survives "
            "as its salt. Real losses come from enzymes eating it, or from the reactive pools above.",
        ],
    }


def _consumers(host: cobra.Model, target_id: str, min_growth_fraction: float) -> list[dict]:
    """Host reactions that consume the target, with flux, genes and knockout safety."""
    met = host.metabolites.get_by_id(target_id)
    solution = host.optimize()
    wild_type_growth = solution.objective_value if solution.status == "optimal" else 0.0
    rows = []
    for rxn in met.reactions:
        if rxn.objective_coefficient:
            continue
        coeff = rxn.metabolites[met]
        consumes = (coeff < 0 and rxn.upper_bound > 0) or (coeff > 0 and rxn.lower_bound < 0)
        if not consumes:
            continue
        flux = float(solution.fluxes.get(rxn.id, 0.0)) if solution.status == "optimal" else 0.0
        eating = (flux > 1e-6 and coeff < 0) or (flux < -1e-6 and coeff > 0)
        with host:
            rxn.knock_out()
            growth_without = max(host.slim_optimize(error_value=0.0), 0.0)
        keeps_growing = (
            wild_type_growth <= 0 or growth_without >= min_growth_fraction * wild_type_growth
        )
        partners = sorted(
            m.id for m, c in rxn.metabolites.items()
            if m.id != target_id and ((c < 0) == (coeff < 0)) and _base_id(m.id) not in CURRENCY
        )
        genes = sorted(g.id for g in rxn.genes if g.id not in PLACEHOLDER_GENES)
        rows.append({
            "reaction": rxn.id,
            "name": rxn.name or rxn.id,
            "equation": rxn.build_reaction_string(),
            "boundary": rxn.boundary,
            "transport": _is_transport(rxn),
            "genes": genes,
            "reacts_with": partners,
            "flux_at_max_growth": round(flux, 6),
            "active": bool(eating),
            "growth_if_deleted": round(growth_without, 6),
            # Deleting a transporter is as likely to trap the product inside the
            # cell as to protect it, so it is never offered as protective.
            "safe_to_delete": bool(
                keeps_growing and genes and not rxn.boundary and not _is_transport(rxn)
            ),
        })
    rows.sort(key=lambda r: (not r["active"], not r["safe_to_delete"], r["reaction"]))
    return rows


def fate_to_text(fate: dict, limit: int = 8) -> str:
    """Human-readable summary of what :func:`product_fate` found."""
    lines = [f"Fate of {fate['target_name']} ({fate['target']}) in the host:", "", fate["charge_note"], ""]
    if not fate["in_host_model"]:
        lines.append("The host model has no enzyme that touches this compound at all, so nothing "
                     "native is going to eat it. That is the good case.")
    elif not fate["consumed_by"]:
        lines.append("No host reaction consumes it.")
    else:
        lines.append(f"{len(fate['consumed_by'])} host reaction(s) consume it:")
        for row in fate["consumed_by"][:limit]:
            tag = "ACTIVE" if row["active"] else "idle  "
            safe = "deletable" if row["safe_to_delete"] else "do not delete"
            partners = f" with {', '.join(row['reacts_with'])}" if row["reacts_with"] else ""
            genes = f" genes: {', '.join(row['genes'])}" if row["genes"] else ""
            lines.append(f"  [{tag}] {row['reaction']}{partners} — {safe}.{genes}")
        if fate["protective_knockouts"]:
            lines += ["", "Knockouts that would protect the product: "
                          + ", ".join(fate["protective_knockouts"])]
    if fate["functional_group_hints"]:
        lines += ["", "Chemically fragile groups suggested by the name (check these):"]
        lines += [f"  {g['group']}: {g['why_it_matters']}" for g in fate["functional_group_hints"]]
    if fate["reactive_pools_the_host_makes"]:
        lines += ["", "Reactive pools the host is producing under these conditions:"]
        lines += [f"  {p['metabolite']} ({p['kind']}): {p['why_it_matters']}"
                  for p in fate["reactive_pools_the_host_makes"]]
    return "\n".join(lines)


def routes_to_text(routes: list[Route], limit: int = 5) -> str:
    """Human-readable summary of what :func:`nearest_precursors` found."""
    if not routes:
        return "No native precursor found within the step limits. Try a larger donor model or more steps."
    lines = [f"Nearest native precursors for {routes[0].target_name} ({routes[0].target}):"]
    for i, route in enumerate(routes[:limit], 1):
        flux = "unchecked" if route.precursor_max_flux is None else f"{route.precursor_max_flux:.2f} mmol/gDW/h"
        lines.append(
            f"{i}. {route.precursor_name} ({route.precursor}) — "
            f"{route.heterologous_steps} new enzyme(s), host can make {flux}"
        )
        for step in route.steps:
            tag = "have" if step.native else "ADD "
            genes = f"  genes: {', '.join(step.genes)}" if step.genes else ""
            lines.append(f"     [{tag}] {step.reaction}: {step.equation}{genes}")
    return "\n".join(lines)

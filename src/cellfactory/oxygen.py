"""Oxygen-sensitive enzymes: warn when an aerobic prediction leans on one.

Some enzymes only work when there is no oxygen. The best known is *pyruvate
formate-lyase* (PFL), which carries a glycyl radical on its own protein
backbone: oxygen cleaves that radical and the enzyme is finished, which is why
*E. coli* ferments pyruvate to formate only when air runs out. Other examples
are the iron-sulfur pyruvate:flavodoxin/ferredoxin oxidoreductases, the Ni-Fe
and Fe-Fe hydrogenases, nitrogenase, the B12 diol dehydratases and the class
III (anaerobic) ribonucleotide reductase.

Genome-scale models mostly do not know any of this. A BiGG model lists PFL as a
reaction with bounds (0, 1000) whatever the medium, so FBA may happily run it in
fully aerobic conditions. That is not a harmless detail: PR #17's OptKnock run
produced a "best" aerobic succinate design that deleted pyruvate dehydrogenase
and pushed acetyl-CoA through PFL instead — a strain nobody could ever build.

This module **warns and changes nothing**. Reaction bounds and the medium are
left exactly as the model ships them (Ferrous's decision, 2026-11-11), because
closing reactions silently moves every number in the report. Instead each
suspect reaction is matched by EC number and enzyme name (not by hard-coded
BiGG ids, which differ between models), checked against the flux distributions
the pipeline already computed, and reported with one extra piece of evidence:
closing it and re-solving, to say whether the prediction actually *needs* the
anaerobic enzyme or merely happened to pick it out of several equally good
routes.
"""

from __future__ import annotations

import re
from collections.abc import Mapping
from dataclasses import dataclass, field

import cobra
import pandas as pd

from .models import OXYGEN_EXCHANGE

FLUX_TOL = 1e-6
OXYGEN_TOL = 1e-6  # below this oxygen uptake the condition is anaerobic, whatever the medium says
COST_TOL = 1e-3  # relative drop in the objective below which closing a reaction counts as free
REPORT_COST = 0.005  # only quantify the drop in words above this; below it, solver noise dominates
COLUMNS = ["where", "kind", "reaction", "name", "enzyme", "flux", "matched_by", "objective", "objective_without",
           "cost", "verdict", "why"]


@dataclass(frozen=True)
class SensitiveEnzyme:
    """A class of enzymes that cannot work in the presence of oxygen."""

    key: str
    label: str
    why: str  # what oxygen does to it, in one clause
    ec_prefixes: tuple[str, ...] = ()
    name_patterns: tuple[str, ...] = ()
    exclude_patterns: tuple[str, ...] = ()


# Sources for each class are in docs/learning/07-oxygen-sensitive-enzymes.md.
# EC prefixes match a reaction's `ec-code` annotation; name patterns are regexes against the reaction
# name, for models (such as e_coli_core) that carry no EC annotations at all.
SENSITIVE_ENZYMES: tuple[SensitiveEnzyme, ...] = (
    SensitiveEnzyme(
        "glycyl_radical_lyase",
        "glycyl radical formate-lyase (PFL family)",
        "oxygen cleaves the enzyme's glycyl radical and destroys the protein, so the cell only makes the "
        "active enzyme when oxygen is absent",
        ec_prefixes=("2.3.1.54",),
        name_patterns=(r"formate[ -]lyase",),
        exclude_patterns=(r"hydrogen[ -]lyase",),  # formate-hydrogen lyase is the hydrogenase class below
    ),
    SensitiveEnzyme(
        "oxo_acid_oxidoreductase",
        "pyruvate/2-oxoacid:flavodoxin or ferredoxin oxidoreductase",
        "its iron-sulfur clusters and radical intermediate are oxidised by oxygen",
        ec_prefixes=("1.2.7.",),
        # "(?<![a-z])" keeps phosphoenol*pyruvate* synthase (PPS), an ordinary aerobic enzyme, out.
        name_patterns=(r"(?<![a-z])pyruvate synthase",
                       r"(pyruvate|oxoglutarate|ketoglutarate).*(flavodoxin|ferredoxin)",
                       r"(flavodoxin|ferredoxin).*oxidoreductase"),
    ),
    SensitiveEnzyme(
        "hydrogenase",
        "hydrogenase / formate-hydrogen lyase",
        "the Ni-Fe or Fe-Fe active site is inhibited and in several cases damaged by oxygen, and the cell "
        "builds these enzymes anaerobically",
        ec_prefixes=("1.12.",),
        # The word on its own, so de*hydrogenase* and trans*hydrogenase* (both aerobically fine) don't match.
        name_patterns=(r"(?:^|[\s(\-])hydrogenase", r"hydrogen[ -]lyase"),
        # A lumped phosphonate reaction that happens to release H2: not a hydrogenase, and we make no claim about it.
        exclude_patterns=(r"phosphonate",),
    ),
    SensitiveEnzyme(
        "class_iii_rnr",
        "class III (anaerobic) ribonucleotide reductase",
        "it also uses a glycyl radical, which oxygen destroys; aerobic cells use the class I enzyme instead",
        ec_prefixes=("1.17.4.2",),
        name_patterns=(r"ribonucleoside-triphosphate reductase",),
    ),
    SensitiveEnzyme(
        "nitrogenase",
        "nitrogenase",
        "its iron-molybdenum cofactor is irreversibly damaged by oxygen",
        ec_prefixes=("1.18.6.", "1.19.6."),
        name_patterns=(r"nitrogenase",),
    ),
    SensitiveEnzyme(
        "b12_diol_dehydratase",
        "B12-dependent diol/glycerol dehydratase",
        "the radical its B12 cofactor generates is quenched by oxygen",
        ec_prefixes=("4.2.1.28", "4.2.1.30"),
        name_patterns=(r"(propanediol|glycerol) dehydratase",),
    ),
)


@dataclass
class FluxSet:
    """One flux distribution to check, with enough context to re-solve it.

    ``objective`` is the reaction that was maximised (biomass, or a product
    exchange); ``knockouts`` are reactions that were closed first, so a design's
    flux state can be re-solved the same way it was produced.
    """

    label: str  # how to name this flux state to a reader
    fluxes: Mapping[str, float] | pd.Series
    objective: str
    knockouts: tuple[str, ...] = ()
    kind: str = "state"  # "design" for a knockout design's flux state, which is judged more harshly


@dataclass
class OxygenResult:
    """What ``oxygen_sensitivity_check`` found. ``table`` has one row per (flux state, reaction)."""

    table: pd.DataFrame = field(default_factory=lambda: pd.DataFrame(columns=COLUMNS))
    candidates: pd.DataFrame = field(default_factory=pd.DataFrame)  # the model's oxygen-sensitive reactions
    checked: list[str] = field(default_factory=list)  # labels of the flux states actually checked (aerobic ones)
    skipped: list[str] = field(default_factory=list)  # labels skipped because no oxygen was taken up

    @property
    def flagged(self) -> bool:
        return not self.table.empty


def _ec_codes(reaction: cobra.Reaction) -> list[str]:
    value = reaction.annotation.get("ec-code")
    if value is None:
        return []
    return [str(v) for v in (value if isinstance(value, list) else [value])]


def _match(reaction: cobra.Reaction) -> tuple[SensitiveEnzyme, str] | None:
    """The enzyme class this reaction belongs to, and how it was matched."""
    name = reaction.name or ""
    codes = _ec_codes(reaction)
    for enzyme in SENSITIVE_ENZYMES:
        if any(re.search(p, name, re.I) for p in enzyme.exclude_patterns):
            continue
        hit = next((c for c in codes if any(c.startswith(p) for p in enzyme.ec_prefixes)), None)
        if hit is not None:
            return enzyme, f"EC {hit}"
        if any(re.search(p, name, re.I) for p in enzyme.name_patterns):
            return enzyme, "enzyme name"
    return None


def oxygen_sensitive_reactions(model: cobra.Model) -> pd.DataFrame:
    """Every reaction in ``model`` catalysed by an enzyme oxygen would disable.

    Matched by the reaction's ``ec-code`` annotation first, then by its name, so
    the same table works for an annotated genome-scale model and for
    e_coli_core, which has no EC annotations. Columns: reaction, name, enzyme,
    matched_by, why, open (whether the model lets it carry flux at all).
    """
    rows = []
    for rxn in model.reactions:
        found = _match(rxn)
        if found is None:
            continue
        enzyme, matched_by = found
        rows.append({"reaction": rxn.id, "name": rxn.name, "enzyme": enzyme.label, "matched_by": matched_by,
                     "why": enzyme.why, "open": rxn.bounds != (0.0, 0.0)})
    return pd.DataFrame(rows, columns=["reaction", "name", "enzyme", "matched_by", "why", "open"])


def _oxygen_uptake(fluxes: Mapping[str, float] | pd.Series, oxygen_exchange: str) -> float:
    try:
        value = fluxes[oxygen_exchange]
    except (KeyError, IndexError):
        return 0.0
    return max(-float(value), 0.0)  # uptake is negative flux on an exchange


def _objective_without(model: cobra.Model, flux_set: FluxSet, reaction: str) -> float:
    """The best achievable objective of ``flux_set`` with ``reaction`` closed."""
    with model:
        for rid in flux_set.knockouts:
            if rid in model.reactions:
                model.reactions.get_by_id(rid).knock_out()
        model.reactions.get_by_id(reaction).knock_out()
        model.objective = flux_set.objective
        return model.slim_optimize(error_value=0.0)


def oxygen_sensitivity_check(
    model: cobra.Model,
    flux_sets: list[FluxSet],
    oxygen_exchange: str = OXYGEN_EXCHANGE,
    flux_tol: float = FLUX_TOL,
) -> OxygenResult:
    """Flag oxygen-sensitive reactions carrying flux in an aerobic flux state.

    Each :class:`FluxSet` is checked only if the cell actually takes up oxygen in
    it: a flux state with no oxygen uptake is anaerobic, whatever the medium
    allows, and these enzymes are then perfectly normal. For each flagged
    reaction the same objective is re-optimised with that one reaction closed, and
    the ``verdict`` says which of two very different problems you have:

    * ``avoidable`` — the objective is unchanged without it, so the solver
      merely picked an anaerobic enzyme out of several equally good routes. The
      number stands; the flux picture is misleading.
    * ``load-bearing`` — the objective drops, so this prediction genuinely
      depends on an enzyme that cannot run in air. A strain design in this
      state is not buildable; a growth or yield number is overstated.
    * ``infeasible without it`` — nothing at all is possible without the
      reaction. For the cell's own flux state that usually means the model draws
      that step only as its anaerobic enzyme (a gap in the model); for a
      knockout design it means the design is not buildable aerobically.

    The model's bounds and medium are never changed.
    """
    candidates = oxygen_sensitive_reactions(model)
    result = OxygenResult(candidates=candidates)
    open_candidates = candidates[candidates["open"]] if not candidates.empty else candidates
    if open_candidates.empty:
        return result
    rows = []
    for flux_set in flux_sets:
        if _oxygen_uptake(flux_set.fluxes, oxygen_exchange) <= OXYGEN_TOL:
            result.skipped.append(flux_set.label)
            continue
        result.checked.append(flux_set.label)
        objective_value = float(flux_set.fluxes[flux_set.objective]) if flux_set.objective in flux_set.fluxes else 0.0
        for row in open_candidates.itertuples():
            flux = float(flux_set.fluxes.get(row.reaction, 0.0))
            if abs(flux) <= flux_tol:
                continue
            without = _objective_without(model, flux_set, row.reaction)
            if without <= flux_tol and objective_value > flux_tol:
                verdict, cost = "infeasible without it", 1.0
            elif objective_value <= flux_tol:
                verdict, cost = "load-bearing", float("nan")
            else:
                cost = max(objective_value - without, 0.0) / objective_value
                verdict = "avoidable" if cost <= COST_TOL else "load-bearing"
            rows.append({"where": flux_set.label, "kind": flux_set.kind, "reaction": row.reaction, "name": row.name,
                         "enzyme": row.enzyme, "flux": flux, "matched_by": row.matched_by,
                         "objective": objective_value, "objective_without": without, "cost": cost,
                         "verdict": verdict, "why": row.why})
    if rows:
        result.table = pd.DataFrame(rows, columns=COLUMNS)
    return result


def _join(values) -> str:
    return ", ".join(dict.fromkeys(str(v) for v in values))


def _reactions_phrase(group: pd.DataFrame) -> str:
    reactions = _join(f"{r.reaction} ({r.name})" for r in group.itertuples())
    return f"{reactions}: {_join(group['enzyme'])}, which oxygen disables"


def _design_caveats(designs: pd.DataFrame) -> list[str]:
    """One caveat per verdict across all affected designs, so five sibling designs don't give five caveats."""
    out = []
    dead = designs[designs["verdict"].isin(["load-bearing", "infeasible without it"])]
    if not dead.empty:
        labels = _join(dead["where"])
        out.append(f"Aerobically, {len(dead['where'].unique())} knockout design(s) only work through an enzyme "
                   f"oxygen destroys, so they are not buildable as they stand: {labels}, each needing "
                   f"{_reactions_phrase(dead)}. Deleting those reactions in the design's model leaves it unable to "
                   "reach the growth it was designed for. Either design without oxygen, or close those reactions "
                   "before designing and re-run. This check only warns: no bounds were changed.")
    weak = designs[designs["verdict"] == "avoidable"]
    weak = weak[~weak["where"].isin(dead["where"])] if not dead.empty else weak
    if not weak.empty:
        out.append(f"{len(weak['where'].unique())} knockout design(s) route some flux through "
                   f"{_reactions_phrase(weak)}, though they do not depend on it: {_join(weak['where'])}. The "
                   "yields stand; the flux picture does not.")
    return out


def oxygen_caveats(result: OxygenResult) -> list[str]:
    """Plain-words caveats for the report: what an aerobic cell could not actually do."""
    if not result.flagged:
        return []
    out = []
    table = result.table
    states = table[table["kind"] != "design"]
    for verdict in ("load-bearing", "infeasible without it", "avoidable"):
        part = states[states["verdict"] == verdict]
        for where, group in part.groupby("where", sort=False):
            if verdict == "load-bearing":
                worst = group["cost"].max()
                drop = f", costing up to {worst:.0%} of it" if pd.notna(worst) and worst >= REPORT_COST else ""
                out.append(f"With oxygen present, {where} needs {_reactions_phrase(group)}. Closing those "
                           f"reactions lowers what this flux state can reach{drop}, so a real aerobic cell would "
                           "not reach these numbers. Re-run without oxygen, or close those reactions yourself, to "
                           "see the aerobic number. This check only warns: no bounds were changed.")
            elif verdict == "infeasible without it":
                out.append(f"With oxygen present, {where} cannot work at all without {_reactions_phrase(group)}, "
                           "and the model offers no other route for that step. That is usually a gap in the model "
                           "— a step drawn only as its anaerobic enzyme — rather than a wrong prediction, but "
                           "every number from this flux state rests on an enzyme the cell could not run in air. "
                           "This check only warns: no bounds were changed.")
            else:
                out.append(f"In {where} the solver routes flux through {_reactions_phrase(group)}. The result does "
                           "not depend on it (closing those reactions costs nothing), so the yield and growth "
                           "numbers stand, but the flux picture shows a route an aerobic cell has no working "
                           "enzyme for. This check only warns: no bounds were changed.")
    out += _design_caveats(table[table["kind"] == "design"])
    return out

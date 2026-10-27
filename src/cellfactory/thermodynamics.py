"""Reaction thermodynamics: ΔG'° from eQuilibrator, and which directions it rules out.

Not to be confused with ``checks.free_energy_checks``, which computes no free
energies at all: it checks that the model can't make ATP or grow from nothing.

eQuilibrator (Beber et al. 2022, Nucleic Acids Res 50:D603, doi:10.1093/nar/gkab1106)
estimates each reaction's standard transformed Gibbs energy ΔG'° with the
component-contribution method (Noor et al. 2013, PLoS Comput Biol 9:e1003098,
doi:10.1371/journal.pcbi.1003098), adjusted for pH, ionic strength, pMg and
temperature. ``reaction_thermodynamics`` estimates it for every reaction in a
model and asks one question per reaction: across every metabolite
concentration a cell could plausibly hold (1 µM to 10 mM by default, the
bounds of Noor et al. 2014's MDF), plus a 95% margin for the estimate's own
error, can ΔG' still take either sign? If not, one direction is impossible.

The rule is deliberately conservative: a 2-in/2-out reaction needs |ΔG'°| above
about 45 kJ/mol before concentrations can't flip it, so most reactions come out
"either direction", including ones that are irreversible in a living cell for
kinetic or regulatory reasons (PFK, pyruvate kinase). What it does catch is a
model letting a reaction run a direction physics forbids.

``apply_directionality`` closes the forbidden direction of reactions the model
left reversible. It never opens a bound, and never touches a reaction whose
model bounds disagree with the estimate: those are reported as conflicts for a
person to look at.

Transport reactions (more than one compartment) are skipped: their ΔG' depends
on membrane potential and the pH difference across the membrane, which a
single-compartment estimate doesn't have. So are boundary and biomass reactions.

Runs offline: equilibrator-api 0.8.1 ships a compound set covering BiGG and KEGG
ids, and the predictor's parameters. Metabolites outside that set are reported
as unmapped. ``download=True`` lets eQuilibrator fetch its full ~1.3 GB compound
database from Zenodo instead, once, the first time a lookup needs it.
"""

from __future__ import annotations

import logging
import math
import os
from dataclasses import dataclass
from functools import lru_cache

import cobra
import pandas as pd

log = logging.getLogger(__name__)

R = 8.314462618e-3  # kJ/(mol K)
CONCENTRATION_BOUNDS = (1e-6, 1e-2)  # M; Noor et al. 2014, PLoS Comput Biol 10:e1003483 (MDF)
Z95 = 1.96  # the estimate's standard error, scaled to a 95% margin
UNRELIABLE_SIGMA = 100.0  # kJ/mol; above this the estimate can't decide any direction (eQuilibrator puts 1e5 on unknowns)
# Water's activity is 1 and H+ is set by pH, so neither has a concentration to vary.
FIXED = {"XLYOFNOQVPJJNP-UHFFFAOYSA-N", "GPRLSGONYQIRFK-UHFFFAOYSA-N"}  # InChIKeys of H2O, H+
ANNOTATIONS = [("bigg.metabolite", "bigg.metabolite"), ("kegg.compound", "kegg"), ("chebi", "chebi"),
               ("metanetx.chemical", "metanetx.chemical")]  # model annotation key -> eQuilibrator namespace
COLUMNS = ["reaction", "name", "dg0_prime", "dg0_sigma", "dg_prime_min", "dg_prime_max", "method",
           "thermo_direction", "model_direction", "verdict", "note"]


class ThermoUnavailable(RuntimeError):
    """eQuilibrator can't be used: not installed, or its data isn't available offline."""


@dataclass(frozen=True)
class Conditions:
    p_h: float = 7.5  # E. coli keeps its cytoplasm near 7.5
    ionic_strength: float = 0.25  # M
    p_mg: float = 3.0
    temperature: float = 298.15  # K; eQuilibrator's estimates are fitted at 25 °C

    def __str__(self) -> str:
        return f"pH {self.p_h}, ionic strength {self.ionic_strength} M, pMg {self.p_mg}, {self.temperature} K"


def _offline() -> bool:
    return os.environ.get("CELLFACTORY_OFFLINE", "").lower() not in ("", "0", "false", "no")


@lru_cache(maxsize=2)
def _equilibrator(download: bool):
    """A ComponentContribution object; a few seconds to build, so built once per process."""
    try:
        import equilibrator_cache_data
        from component_contribution.parameters import CCModelParameters
        from equilibrator_api import ComponentContribution
        from equilibrator_cache.api import create_compound_cache_from_sqlite_file, create_compound_cache_from_zenodo
    except ImportError as err:
        raise ThermoUnavailable(f"equilibrator-api is not installed ({err.name}); "
                                "pip install 'cellfactory[thermo]'") from err
    if download:
        ccache = create_compound_cache_from_zenodo(tier="auto")
    else:
        # The plain mini cache never escalates, so a lookup outside it is a miss, not a 1.3 GB download.
        if not equilibrator_cache_data.is_available() or not CCModelParameters.is_bundled_available():
            raise ThermoUnavailable("this equilibrator-api install carries no bundled compound data or parameters, "
                                    "and downloading them is off")
        ccache = create_compound_cache_from_sqlite_file(str(equilibrator_cache_data.MINI_COMPOUND_CACHE_PATH))
    return ComponentContribution(ccache=ccache)


def _compound(cc, met: cobra.Metabolite):
    base = met.id[: -len(met.compartment) - 1] if met.compartment and met.id.endswith("_" + met.compartment) else met.id
    queries = [f"bigg.metabolite:{base}"]
    for key, namespace in ANNOTATIONS:
        values = met.annotation.get(key, [])
        for value in [values] if isinstance(values, str) else values:
            queries.append(f"{namespace}:{value.removeprefix('CHEBI:')}" if namespace == "chebi" else f"{namespace}:{value}")
    for query in dict.fromkeys(queries):
        try:
            compound = cc.get_compound(query)
        except Exception:  # an odd annotation value is a miss, not a crash
            continue
        if compound is not None:
            return compound
    return None


def _model_direction(rxn: cobra.Reaction) -> str:
    fwd, rev = rxn.upper_bound > 0, rxn.lower_bound < 0
    return {(True, True): "reversible", (True, False): "forward", (False, True): "reverse"}.get((fwd, rev), "blocked")


def _verdict(thermo: str, model: str) -> str:
    if thermo == "either" or model == "blocked":
        return "agree"
    if model == "reversible":
        return "can tighten"
    return "agree" if model == thermo else "CONFLICT"


def reaction_thermodynamics(
    model: cobra.Model,
    reactions: list[str] | None = None,
    conditions: Conditions = Conditions(),
    concentrations: tuple[float, float] = CONCENTRATION_BOUNDS,
    download: bool = False,
    estimator=None,
) -> pd.DataFrame:
    """ΔG'° and the thermodynamically possible direction of each reaction, in kJ/mol.

    Columns: ``dg0_prime`` and its standard error ``dg0_sigma``; ``dg_prime_min``
    / ``dg_prime_max``, the range ΔG' can take over ``concentrations`` (molar);
    ``thermo_direction`` ("forward", "reverse" or "either", after a 95% margin
    for the estimate's error); ``model_direction`` from the current bounds; and
    ``verdict``: "agree", "can tighten" (model allows a direction physics rules
    out, but also the possible one) or "CONFLICT" (model allows only the
    impossible direction). Reactions without an estimate have a ``note`` saying
    why and an empty verdict.

    ``estimator`` is for tests: a callable taking a cobra.Reaction and returning
    ``(dg0_prime, sigma, method, fixed)``, where ``fixed`` holds the ids of
    metabolites with no concentration to vary (water, H+), or raising
    ``LookupError`` with a note.
    Raises ``ThermoUnavailable`` if eQuilibrator can't run.
    """
    if estimator is None:
        estimator = _EquilibratorEstimator(conditions, download=download and not _offline())
    low, high = (math.log(c) for c in concentrations)
    rt = R * conditions.temperature
    rows = []
    for rxn in (model.reactions if reactions is None else [model.reactions.get_by_id(r) for r in reactions]):
        model_dir = _model_direction(rxn)
        skip = ("boundary reaction" if rxn.boundary
                else "biomass pseudo-reaction" if rxn.objective_coefficient or "biomass" in rxn.id.lower()
                else "transport: depends on membrane potential and pH gradient" if len(rxn.compartments) > 1
                else None)
        if skip:
            rows.append((rxn.id, rxn.name, None, None, None, None, None, None, model_dir, "", skip))
            continue
        try:
            dg0, sigma, method, fixed = estimator(rxn)
        except LookupError as err:
            rows.append((rxn.id, rxn.name, None, None, None, None, None, None, model_dir, "", str(err)))
            continue
        if sigma > UNRELIABLE_SIGMA:
            rows.append((rxn.id, rxn.name, dg0, sigma, None, None, method, None, model_dir, "",
                         "estimate too uncertain to decide a direction"))
            continue
        # Products at their lowest and substrates at their highest give the most negative ΔG', and vice versa.
        varying = {m: c for m, c in rxn.metabolites.items() if m.id not in fixed}
        dg_min = dg0 + rt * sum(c * (low if c > 0 else high) for c in varying.values())
        dg_max = dg0 + rt * sum(c * (high if c > 0 else low) for c in varying.values())
        thermo = ("forward" if dg_max + Z95 * sigma < 0 else "reverse" if dg_min - Z95 * sigma > 0 else "either")
        rows.append((rxn.id, rxn.name, dg0, sigma, dg_min, dg_max, method, thermo, model_dir,
                     _verdict(thermo, model_dir), ""))
    return pd.DataFrame(rows, columns=COLUMNS)


class _EquilibratorEstimator:
    def __init__(self, conditions: Conditions, download: bool):
        from equilibrator_api import Q_, Reaction

        self.cc, self.Q_, self.Reaction = _equilibrator(download), Q_, Reaction
        self.conditions = conditions

    def __call__(self, rxn: cobra.Reaction):
        cc, c = self.cc, self.conditions
        # The ComponentContribution object is shared, so set the conditions on every call.
        cc.p_h, cc.ionic_strength = self.Q_(c.p_h), self.Q_(c.ionic_strength, "M")
        cc.p_mg, cc.temperature = self.Q_(c.p_mg), self.Q_(c.temperature, "K")
        sparse, missing, unknown, fixed = {}, [], [], set()
        for met, coeff in rxn.metabolites.items():
            compound = _compound(cc, met)
            if compound is None:
                missing.append(met.id)
                continue
            if not compound.atom_bag:
                unknown.append(met.id)
            if compound.inchi_key in FIXED:
                fixed.add(met.id)
            sparse[compound] = sparse.get(compound, 0.0) + coeff
        if missing:
            raise LookupError(f"not in eQuilibrator's compound set: {', '.join(missing)}")
        # Lipids, quinone pools and other lumped compounds are often in the set by name but with no structure.
        if unknown:
            raise LookupError(f"no chemical structure in eQuilibrator for: {', '.join(unknown)}")
        reaction = self.Reaction({k: v for k, v in sparse.items() if v})
        if not reaction.is_balanced():
            raise LookupError("unbalanced as eQuilibrator reads it; not estimated")
        dg = cc.standard_dg_prime(reaction)
        method = "group contribution" if cc.is_using_group_contribution(reaction) else "component contribution"
        return dg.value.m_as("kJ/mol"), dg.error.m_as("kJ/mol"), method, fixed


def apply_directionality(model: cobra.Model, table: pd.DataFrame, group_contribution: bool = False) -> pd.DataFrame:
    """Close the impossible direction of every "can tighten" reaction in ``table``. Returns the changed rows.

    Only ever narrows a bound to 0. CONFLICT rows are left alone: the model's
    curators may know something (a coupled transport, a different compartment
    pH) the estimate doesn't, so a person should look before anything changes.

    Estimates that fall back on group contribution are skipped unless
    ``group_contribution=True``: they can be wrong by more than their stated
    error. eQuilibrator puts ubiquinone reduction by succinate (SUCDi) at about
    -74 kJ/mol where the quinone's redox potential says a few kJ/mol (see
    tests/benchmarks/cases.py).
    """
    changed = table[table["verdict"] == "can tighten"]
    if not group_contribution:
        changed = changed[changed["method"] != "group contribution"]
    for row in changed.itertuples():
        rxn = model.reactions.get_by_id(row.reaction)
        if row.thermo_direction == "forward":
            rxn.lower_bound = 0.0
        else:
            rxn.upper_bound = 0.0
    return changed


def thermo_caveats(table: pd.DataFrame, conditions: Conditions = Conditions(), tightened: int | None = None) -> list[str]:
    """Plain-words caveats for the report."""
    estimated = table["thermo_direction"].notna()
    eligible = table["note"].ne("boundary reaction") & table["note"].ne("biomass pseudo-reaction")
    unmapped = int((eligible & ~estimated & ~table["note"].str.startswith("transport")).sum())
    transport = int(table["note"].str.startswith("transport").sum())
    out = [f"Reaction directions checked against eQuilibrator ΔG'° estimates ({conditions}) for "
           f"{int(estimated.sum())} of {int(eligible.sum())} reactions. {transport} transport and {unmapped} other "
           "reactions had no usable estimate and keep the model's directions unchecked."]
    conflicts = table[table["verdict"] == "CONFLICT"]
    if len(conflicts):
        shaky = conflicts["method"] == "group contribution"
        out.append(f"{len(conflicts)} reaction(s) can run only in a direction thermodynamics rules out at 1 µM-10 mM: "
                   f"{', '.join(conflicts['reaction'][:10])}. Check them before trusting yields that use them"
                   + (f"; {int(shaky.sum())} of them use group-contribution estimates, so the error may be in the "
                      "estimate rather than the model." if shaky.any() else "."))
    gc = int((table["method"] == "group contribution").sum())
    if gc:
        out.append(f"{gc} ΔG'° estimate(s) rest on group contribution (a compound with no measured data), "
                   "which is less accurate; see dg0_sigma.")
    if tightened:
        out.append(f"{tightened} reversible reaction(s) were made one-way from these estimates before simulating.")
    return out

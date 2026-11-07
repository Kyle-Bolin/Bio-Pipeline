"""Loading metabolic models and setting up growth conditions.

A *genome-scale metabolic model* (GEM) is a list of every chemical reaction an
organism's enzymes can catalyse, plus which genes encode those enzymes. COBRApy
reads them from SBML/JSON/MAT files; BiGG (http://bigg.ucsd.edu) hosts curated
ones for E. coli, yeast, Salmonella, human cells and ~100 more organisms.
"""

from __future__ import annotations

import logging
from pathlib import Path

import cobra
from cobra.io import load_json_model, load_matlab_model, read_sbml_model
from cobra.io import load_model as _cobra_load_model

from .identifiers import ALIASES, resolve_metabolite  # noqa: F401  (ALIASES re-exported)

log = logging.getLogger(__name__)

OXYGEN_EXCHANGE = "EX_o2_e"

# Models COBRApy ships with, so they work offline. Keys are BiGG ids.
BUNDLED = {
    "e_coli_core": "textbook",  # E. coli, 95 reactions: central metabolism only. Fast, great for learning.
    "iJO1366": "iJO1366",  # E. coli K-12 MG1655, 2583 reactions (2011).
    "iYS1720": "salmonella",  # Salmonella Typhimurium LT2, 3357 reactions.
}

_READERS = {".xml": read_sbml_model, ".sbml": read_sbml_model, ".json": load_json_model, ".mat": load_matlab_model}


def load_model(model: str | Path = "e_coli_core", fallback: str | None = None) -> cobra.Model:
    """Load a model from a file path, a bundled model, or BiGG.

    ``model`` can be:

    * a path to an SBML (.xml/.sbml), JSON or MATLAB (.mat) file, for models
      from anywhere (BiGG, BioModels, a paper's supplement, your own);
    * a bundled BiGG id (see ``BUNDLED``), loaded from disk;
    * any other BiGG id, e.g. ``iML1515`` (E. coli) or ``iMM904`` (yeast),
      downloaded from BiGG and cached by COBRApy.

    If a download fails and ``fallback`` is set, that model is loaded instead,
    so notebooks still run offline. Leave it unset unless the fallback is the
    same organism: silently swapping yeast for E. coli would be worse than an error.
    """
    path = Path(model)
    if path.suffix.lower() in _READERS:
        return _READERS[path.suffix.lower()](str(path))
    if model in BUNDLED:
        return _cobra_load_model(BUNDLED[model])
    try:
        return _cobra_load_model(str(model))
    except Exception as err:  # network errors surface as several exception types
        if fallback is None:
            raise RuntimeError(
                f"Could not load {model!r} from BiGG ({err}). Check the id at http://bigg.ucsd.edu/models, "
                "or download the SBML file and pass its path instead."
            ) from err
        log.warning("Could not fetch %s from BiGG (%s); using %s instead.", model, err, fallback)
        return load_model(fallback)


def find_metabolites(model: cobra.Model, query: str) -> list[cobra.Metabolite]:
    """Metabolites whose id or name matches ``query`` (case-insensitive substring)."""
    q = query.lower()
    return [m for m in model.metabolites if q in m.id.lower() or q in (m.name or "").lower()]


def resolve_exchange(model: cobra.Model, compound: str, create_demand: bool = False, online: bool = True) -> str:
    """Turn a compound into the id of the reaction that moves it in or out of the cell.

    ``compound`` may be an exchange reaction id (``EX_etoh_e``), a metabolite id
    with or without compartment (``etoh_e``, ``etoh``), an exact metabolite
    name (``Ethanol``), or a database id such as ``kegg:C00031``, ``CHEBI:4167``
    or ``MNXM41`` (see ``identifiers.resolve_metabolite``; ``online=False``
    keeps that lookup to the model's own annotations). Models represent "the
    cell takes up / secretes X" as an exchange reaction on the extracellular
    metabolite, so that is what we look for.

    Some compounds the cell makes have no exchange (they stay inside the cell,
    like the pigment lycopene). With ``create_demand=True`` we add a *demand*
    reaction that drains the cytosolic metabolite, which is how you model
    "accumulate product inside the cell".
    """
    if compound in model.reactions:
        return compound
    res = resolve_metabolite(model, compound, online=online)
    if not res.ok:
        raise KeyError(res.error(model.id))
    mets = res.metabolites
    for met in sorted(mets, key=lambda m: m.compartment != "e"):  # prefer extracellular
        for rxn in met.reactions:
            if rxn in model.exchanges:
                return rxn.id
    if create_demand:
        internal = next((m for m in mets if m.compartment == "c"), mets[0])
        return model.add_boundary(internal, type="demand").id
    raise KeyError(f"{compound!r} has no exchange reaction in {model.id}; pass create_demand=True to add one.")


def _carbon_uptakes(model: cobra.Model, min_uptake: float) -> list[str]:
    out = []
    for rxn_id, uptake in model.medium.items():
        (met,) = model.reactions.get_by_id(rxn_id).metabolites
        if met.elements.get("C", 0) and met.id.rsplit("_", 1)[0] not in ("co2", "hco3") and uptake >= min_uptake:
            out.append(rxn_id)
    return out


def set_medium(
    model: cobra.Model,
    substrate: str = "EX_glc__D_e",
    uptake: float = 10.0,
    aerobic: bool = True,
    oxygen: str = OXYGEN_EXCHANGE,
) -> str:
    """Feed the cell ``substrate`` as its main carbon source. Returns the exchange id.

    Uptake rates are in mmol per gram of dry cell weight per hour (mmol/gDW/h);
    10 is typical for E. coli on glucose. BiGG models ship with a minimal medium
    (salts + one sugar). We close every major carbon uptake in it (>= 1
    mmol/gDW/h, ignoring CO2) and open the new substrate. Trace carbon
    nutrients such as vitamins stay on, since some organisms need them.

    ``aerobic=False`` removes oxygen. Fermentation products such as ethanol are
    mostly made without it.
    """
    substrate = resolve_exchange(model, substrate)
    medium = model.medium
    for rxn_id in _carbon_uptakes(model, min_uptake=1.0):
        medium.pop(rxn_id)
    medium[substrate] = uptake
    if aerobic:
        if oxygen not in model.reactions:
            raise KeyError(f"No oxygen exchange {oxygen!r} in {model.id}; pass oxygen=<exchange id>.")
        medium[oxygen] = 1000.0
    else:
        medium.pop(oxygen, None)
    model.medium = medium
    return substrate


def set_glucose_medium(model: cobra.Model, glucose_uptake: float = 10.0, aerobic: bool = True) -> None:
    """Shortcut for ``set_medium`` with glucose."""
    set_medium(model, "EX_glc__D_e", glucose_uptake, aerobic)

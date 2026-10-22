"""Loading metabolic models and setting up growth conditions.

A *genome-scale metabolic model* (GEM) is a list of every chemical reaction an
organism's enzymes can catalyse, plus which genes encode those enzymes. COBRApy
reads them from SBML/JSON files; BiGG (http://bigg.ucsd.edu) hosts curated ones.
"""

from __future__ import annotations

import logging

import cobra
from cobra.io import load_model as _cobra_load_model

log = logging.getLogger(__name__)

GLUCOSE_EXCHANGE = "EX_glc__D_e"
OXYGEN_EXCHANGE = "EX_o2_e"

# Models COBRApy ships with, so they work offline.
BUNDLED = {
    "e_coli_core": "textbook",  # 95 reactions: central metabolism only. Fast, great for learning.
    "iJO1366": "iJO1366",  # 2583 reactions: full E. coli K-12 MG1655 (2011).
}


def load_model(model_id: str = "e_coli_core", fallback: str | None = "iJO1366") -> cobra.Model:
    """Load a model by BiGG id.

    Bundled models load from disk. Anything else (e.g. ``iML1515``, the 2017
    successor to iJO1366) is downloaded from BiGG and cached by COBRApy. If the
    download fails and ``fallback`` is set, the bundled fallback is used instead
    so the notebooks still run without internet.
    """
    if model_id in BUNDLED:
        return _cobra_load_model(BUNDLED[model_id])
    try:
        return _cobra_load_model(model_id)
    except Exception as err:  # network errors surface as several exception types
        if fallback is None:
            raise
        log.warning("Could not fetch %s from BiGG (%s); using bundled %s instead.", model_id, err, fallback)
        return load_model(fallback, fallback=None)


def set_glucose_medium(model: cobra.Model, glucose_uptake: float = 10.0, aerobic: bool = True) -> None:
    """Minimal medium: glucose as the only carbon source, plus salts.

    Uptake rates are in mmol per gram of dry cell weight per hour (mmol/gDW/h).
    10 mmol/gDW/h of glucose is a typical measured value for E. coli.

    The model's default medium is already glucose minimal (M9-like); we just
    set the glucose rate and switch oxygen on or off. Fermentation products such
    as ethanol are mostly made *without* oxygen, which is why ``aerobic`` matters.
    """
    medium = model.medium
    medium[GLUCOSE_EXCHANGE] = glucose_uptake
    if aerobic:
        medium[OXYGEN_EXCHANGE] = 1000.0
    else:
        medium.pop(OXYGEN_EXCHANGE, None)
    model.medium = medium

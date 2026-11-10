"""Browser front end for cellfactory: static files served by the REST API.

The source is a React + TypeScript app in ``frontend/`` at the repository root.
``npm run build`` there writes the bundle into ``static/``, which is committed,
so installing the package and running ``cellfactory-api`` needs no Node. The
bundle makes no external requests, so it works offline. It talks to the REST
API over ``GET /api/tools``, ``POST /api/tools/<name>`` (one endpoint per
function in ``cellfactory.tools.TOOLS``) and ``/api/jobs``, so a server only
needs to mount ``STATIC_DIR`` at ``/`` and the API at ``/api``.
"""

from __future__ import annotations

from pathlib import Path

STATIC_DIR = Path(__file__).parent / "static"

# Tools the page calls by name; mirrors TOOL_NAMES in frontend/src/api.ts.
# tests/test_web.py checks the two lists agree and that each is still in TOOLS.
TOOLS_USED = (
    "list_models", "describe_model", "list_pathways", "search_metabolites",
    "simulate_production", "production_envelope", "search_knockouts",
    "search_knockouts_moma", "find_precursors", "check_product_fate", "check_model",
    "compare_conditions", "flux_distribution", "check_thermodynamics", "check_driving_force",
    "check_enzyme_kinetics", "design_knockouts", "qa_model",
)

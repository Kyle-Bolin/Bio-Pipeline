"""Browser front end for cellfactory: static files served by the REST API.

The page in ``static/`` is plain HTML, CSS and JavaScript with no build step and
no external requests, so it works offline. It talks to the REST API over
``GET /api/tools`` and ``POST /api/tools/<name>`` (one endpoint per function in
``cellfactory.tools.TOOLS``), so a server only needs to mount ``STATIC_DIR`` at
``/`` and the API at ``/api``.
"""

from __future__ import annotations

from pathlib import Path

STATIC_DIR = Path(__file__).parent / "static"

# Tools the page calls by name. tests/test_web.py checks each is still in TOOLS.
TOOLS_USED = (
    "list_models", "describe_model", "list_pathways", "search_metabolites",
    "simulate_production", "production_envelope", "search_knockouts",
    "search_knockouts_moma", "find_precursors", "check_product_fate", "check_model",
)

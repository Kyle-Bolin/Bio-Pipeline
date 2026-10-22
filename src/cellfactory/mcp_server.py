"""MCP server: lets Claude and other AI agents call cellfactory as tools.

Run ``cellfactory-mcp`` (stdio). To add it to Claude Code:

    claude mcp add cellfactory -- cellfactory-mcp

The tools themselves live in ``tools.py``; this file only registers them.
"""

from __future__ import annotations

import functools
import logging

try:  # mcp >= 2
    from mcp.server.mcpserver import MCPServer
    from mcp.server.mcpserver.exceptions import ToolError
except ImportError:  # mcp 1.x
    from mcp.server.fastmcp import FastMCP as MCPServer
    from mcp.server.fastmcp.exceptions import ToolError

from .tools import TOOLS, ToolInputError

INSTRUCTIONS = """\
cellfactory simulates a microbe turning a substrate (e.g. glucose) into a product (e.g. ethanol) with
flux balance analysis on genome-scale metabolic models.

Typical workflow:
1. list_models, then describe_model to pick an organism. e_coli_core is small and fast; iJO1366 is full E. coli.
2. search_metabolites if unsure what a compound is called in the model.
3. simulate_production for growth and yield. Try aerobic=false for fermentation products (ethanol, lactate, succinate).
4. search_knockouts to find gene deletions that raise the guaranteed yield (slow on genome-scale models).
5. check_model before trusting a striking result.

Yields are mol product per mol substrate. guaranteed_yield is the realistic number; max_theoretical_yield is a
ceiling that ignores growth. These are model predictions, not lab results: they ignore regulation, kinetics
and toxicity.
"""


def _tool(fn):
    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
        try:
            return fn(*args, **kwargs)
        except ToolInputError as err:
            raise ToolError(str(err)) from err
    return wrapper


def build_server() -> MCPServer:
    server = MCPServer("cellfactory", instructions=INSTRUCTIONS, log_level="WARNING")
    for fn in TOOLS:
        server.tool()(_tool(fn))
    return server


def main() -> None:
    logging.getLogger("cobra").setLevel(logging.WARNING)  # its INFO chatter floods the client's server log
    build_server().run()


if __name__ == "__main__":
    main()

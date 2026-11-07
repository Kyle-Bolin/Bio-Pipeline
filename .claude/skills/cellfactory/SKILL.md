---
name: cellfactory
description: Simulate a microbe turning a substrate into a product, find the pathway to a compound, and design gene knockouts with the cellfactory MCP tools. Use for questions like "how much ethanol can E. coli make from xylose", "what is the shortest route from E. coli metabolism to this compound", or "which deletions raise succinate yield".
---

# Using cellfactory

The `cellfactory` MCP server (`cellfactory-mcp`, see `.mcp.json`) runs flux balance analysis on
genome-scale metabolic models. These steps keep answers fast and honest.

## Workflow

1. **Pick the model.** Start on `e_coli_core` (seconds per call) to explore, then confirm on a
   genome-scale model (`iJO1366`, or `iML1515` if BiGG is reachable). `list_models` shows what
   loads offline.
2. **Resolve names first.** If a compound isn't an everyday name (glucose, ethanol, succinate...),
   call `search_metabolites` and pass the id it returns. For a database id (`kegg:C00469`, `CHEBI:16236`,
   `MNXM303`), `resolve_compound` shows which model metabolite it maps to and how; read its `caveat` before
   trusting a MetaNetX-only match. A product with no `exchange` still works:
   a demand reaction is added, which models the compound accumulating inside the cell.
3. **Set oxygen deliberately.** Fermentation products (ethanol, lactate, succinate, formate) need
   `aerobic: false`; biomass-linked products and most aerobic pathways need `true`. If unsure, run both.
4. **If the pathway is unknown, search for it first.** `find_precursors` takes a target compound and
   returns the metabolites the host already makes that are fewest enzymes away, with the donor's genes
   for each step. Raise `max_new_enzymes` only as far as needed; a route worth building is usually under
   five steps, and the step count is not a difficulty ranking. For vanillin and lycopene there are ready-made
   pathways: `list_pathways` describes them, and passing `pathway: "vanillin"` (on `iJO1366`; `e_coli_core`
   lacks the precursors) to `simulate_production`, `search_knockouts`, `search_knockouts_moma` or
   `production_envelope` adds the enzymes first. Pass the product's name (`"vanillin"`) as `product`.
5. **Check what happens to the product.** `check_product_fate` lists the host reactions that consume it,
   which deletions protect it, and whether its charge at cytosolic pH means it needs a transporter. Treat
   `consumed_by` and `protective_knockouts` as model output; treat the reactive pools and functional-group
   hints as things to check, since the model has no structures.
6. **Simulate, then design.** `simulate_production` first. To compare feedstocks or oxygen settings, use
   `compare_conditions` (comma-separated `substrates`, `oxygen: "both"`) and compare the carbon yields,
   not the mol/mol ones. `flux_distribution` shows the route behind a number: what goes in, what comes out,
   and the busiest reactions. Only run `search_knockouts` when the
   guaranteed yield is low; on genome-scale models it takes minutes, so keep `max_knockouts` at 3 or less.
   For a product the cell has no reason to make (most non-native compounds), FBA finds no deletion that
   helps; use `search_knockouts_moma`, which predicts how a fresh mutant reroutes carbon (about two minutes per
   knockout round on `iJO1366`, so keep `max_knockouts` at 2).
7. **Check before you claim.** Run `check_model` (or `qa_model` for the MEMOTE gate) on the model you used.
   `ok: false` means it can make energy from nothing, and every yield it gives is inflated. Then run
   `check_thermodynamics` with the same product, substrate and pathway: `against_thermodynamics` lists
   reactions the route runs in a direction eQuilibrator rules out. Rows with `method: "group contribution"`
   may be the estimate's fault (fumarate reductase in `e_coli_core` is a known example); say so. To see
   whether a yield depends on such reactions, rerun `simulate_production` with `thermo_bounds: true`.

## Reading results

- Yields are mol product per mol substrate. `production.carbon_yield` is the fraction of the
  substrate's carbon that ends up in product.
- `guaranteed_yield` (least product made while growing as fast as possible) is the realistic number.
  `max_theoretical_yield` ignores growth and is only a ceiling.
- A knockout design is "growth coupled" when its guaranteed yield is above zero: the cell must make
  the product to grow. Report the genes (`summary.genes`), not just reaction ids; those are what a
  lab deletes.

## Caveats to state in answers

FBA predicts what the network can do, not what the cell will do: it ignores gene regulation, enzyme
kinetics and product toxicity. Say so when a result will drive a lab decision. Greedy knockout search
can miss designs a full OptKnock search would find. A `search_knockouts_moma` yield is a prediction for a
freshly made mutant, not a growth-coupled guarantee: evolution can undo it. Say which search a design came from.

A route from `find_precursors` is a graph result, not a validated pathway: the step count ignores whether
the enzyme folds, whether the thermodynamics run in that direction, and whether the intermediate is toxic.
Say so. An acidic product is not "destroyed" by a base the cell makes; it exists as its salt. Real product
losses come from enzymes consuming it, which `check_product_fate` lists, or from reactive pools, which it
flags as heuristics.

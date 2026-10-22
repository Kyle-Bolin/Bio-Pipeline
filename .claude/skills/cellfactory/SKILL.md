---
name: cellfactory
description: Simulate a microbe turning a substrate into a product and design gene knockouts with the cellfactory MCP tools. Use for questions like "how much ethanol can E. coli make from xylose" or "which deletions raise succinate yield".
---

# Using cellfactory

The `cellfactory` MCP server (`cellfactory-mcp`, see `.mcp.json`) runs flux balance analysis on
genome-scale metabolic models. These steps keep answers fast and honest.

## Workflow

1. **Pick the model.** Start on `e_coli_core` (seconds per call) to explore, then confirm on a
   genome-scale model (`iJO1366`, or `iML1515` if BiGG is reachable). `list_models` shows what
   loads offline.
2. **Resolve names first.** If a compound isn't an everyday name (glucose, ethanol, succinate...),
   call `search_metabolites` and pass the id it returns. A product with no `exchange` still works:
   a demand reaction is added, which models the compound accumulating inside the cell.
3. **Set oxygen deliberately.** Fermentation products (ethanol, lactate, succinate, formate) need
   `aerobic: false`; biomass-linked products and most aerobic pathways need `true`. If unsure, run both.
4. **Simulate, then design.** `simulate_production` first. Only run `search_knockouts` when the
   guaranteed yield is low; on genome-scale models it takes minutes, so keep `max_knockouts` at 3 or less.
5. **Check before you claim.** Run `check_model` on the model you used. `ok: false` means it can make
   energy from nothing, and every yield it gives is inflated.

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
can miss designs a full OptKnock search would find.

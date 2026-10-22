---
{module: cellfactory.qa, date: 2026-10-22}
---

# 01 · Model QA gate: models that make something from nothing

A metabolic model is a hand-built list of reactions, each with its stoichiometry (molecules in and out) and allowed direction. FBA exploits any loophole in that list, so before simulating, the gate checks for two that make every later number meaningless.

**Stoichiometric inconsistency.** Every molecule has positive mass. If no set of positive masses balances every reaction, some reaction creates or destroys matter, typically through a typo in a coefficient, and FBA will treat it as free carbon. Unlike the per-reaction mass-balance check, this needs no chemical formulas.

**Energy-generating cycles.** *ATP* is the cell's energy currency; *NADH* carries electrons. Recharging either costs food, as charging a battery costs power. In a model, a loop of reactions can recharge them from nothing when one reaction may run in a thermodynamically impossible direction. Fritzemeier et al. (2017) found such loops in most published models they examined, inflating predicted growth by up to ~25%.

**What this PR found.** MEMOTE, the standard test suite, says iJO1366 makes ATP from nothing, but only after MEMOTE reopens *superoxide dismutase*, which the curators switched off. With it back on, the model runs malate oxidase backwards, turning hydrogen peroxide into oxygen (almost certainly impossible at any useful rate in a real cell) and then "breathes" its home-made oxygen. With the model's own bounds there is no loop, so the gate re-checks every loop MEMOTE reports and only stops on those that survive.

**Misconception.** Passing does not make a model right. It shows only that no conservation law is broken; regulation, enzyme speed and toxicity are untouched.

**What a biologist would ask.** E. coli really makes superoxide dismutase (*sodA*, *sodB*). Why switch that off rather than forbid the backwards oxidase? The model file doesn't say, so this reading is a guess. Item 3 (ΔG'°) is what decides which directions are physically allowed.

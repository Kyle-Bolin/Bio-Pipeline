---
{module: cellfactory.thermodynamics, date: 2026-10-27}
---

# 03 · Reaction thermodynamics: which way can a reaction run?

**ΔG** (Gibbs free energy change, kJ/mol) says which way a reaction can make net progress: negative means forward, positive means backward. Enzymes change how *fast* a reaction goes, never which *way*. **ΔG'°** is the value with every compound at 1 M at a set pH (the prime). A cell never holds 1 M of anything, and the real value is ΔG' = ΔG'° + RT·ln(products/substrates). Each tenfold change in that ratio moves ΔG' by about 5.7 kJ/mol at 25 °C.

**Misconception.** "Positive ΔG'° means the reaction can't run forward." Aldolase is about +24 kJ/mol, yet it runs forward in every glycolysing cell because the next steps keep its products scarce. So this module never judges a reaction by the sign of ΔG'° alone. It asks: across every concentration a cell could plausibly hold (1 µM to 10 mM), plus a margin for the estimate's error, can ΔG' still take both signs? If not, one direction is physically impossible. On iJO1366 that flags malate oxidase, which the model allows to run backwards and make O₂ from H₂O₂. That is the reaction behind the energy cycle in [note 01](01-model-qa-gate.md).

The rule is weak on purpose. Pyruvate kinase and PFK are one-way in living cells, but only because the cell holds their substrates and products far from equilibrium. A concentration *range* can't show that. Measured concentrations can.

**The numbers are estimates.** eQuilibrator's *component contribution* method treats each molecule as a sum of chemical groups and fits it to measured equilibrium constants. For compounds that were never measured it relies on groups alone. That can be badly wrong: it puts succinate reducing ubiquinone at −74 kJ/mol, where redox potentials say a few kJ/mol. Ubiquinone sits in the membrane, and its real potential there is itself debated. Transport reactions are skipped because they also depend on the voltage and pH difference across the membrane.

**What a biologist would ask.** What are the actual concentrations in my strain, under my conditions? Is this step near equilibrium, or held far from it by regulation? Answering either takes metabolomics, not a model.

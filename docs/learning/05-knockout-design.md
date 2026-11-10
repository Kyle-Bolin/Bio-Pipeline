---
{module: cellfactory.strain_design, date: 2026-11-09}
---

# 05 · Knockout design: making the cell want what we want

A cell evolves to grow, not to make our product. Most products are, from the cell's point of view, carbon thrown away. So the wild type makes almost none. The classic fix is **growth coupling**: delete enzymes until the only fast way left to grow also releases the product. Then selection works for us, since the cells that grow fastest are the ones making the product.

**Why the cell would have to.** Growing means balancing books, not just making biomass. Every glucose broken down hands electrons to NAD+, and without oxygen those electrons must go somewhere. *E. coli* has several exits: lactate, ethanol, succinate, and formate and acetate alongside them. Close enough of the others and lactate becomes one of the few left. That is the original OptKnock design (Burgard et al. 2003, the method): delete *adhE* and *pta*, and lactate is forced out; Fong et al. (2005) built and tested those strains.

**Why one-at-a-time search misses designs.** Metabolism is full of parallel routes. Close one and the carbon simply takes another, so no single deletion helps, and the greedy search in `knockouts.py` stops. Aerobic succinate in the core model is the example here: it needs several deletions at once, closing the TCA cycle at succinate dehydrogenase (or fumarase) *and* the acetate route, which overlaps with the deletions Lin et al. (2005) made in real aerobic succinate strains. StrainDesign writes the question as one mixed-integer program (on/off switches for every reaction) and searches every combination at once.

**Three strictness levels.** *OptKnock* assumes the cell, among equally fast flux patterns, picks the one we like: optimistic. *RobustKnock* scores the worst one. *Minimal cut sets* go further and require the product at every growth rate above a floor, not just the fastest ("strong" coupling). The module re-checks every design with the pipeline's own FBA and says which kind it got.

**Reactions are not genes.** A lab deletes genes. One reaction may have two **isozymes** (different proteins doing the same job, so both must go), or be one **complex** of several subunits (one deletion is enough). And one gene can serve several reactions: the core model's transketolase genes run two reactions, so "delete TKT1" in practice also deletes TKT2, and the design stops growing. The gene-level columns show this.

**A model artifact, caught in review.** Left alone, the model's best aerobic succinate design deleted pyruvate dehydrogenase and fed acetyl-CoA through **pyruvate formate-lyase** (PFL) instead. In a real cell PFL is destroyed by oxygen and works essentially only without it, but the core model lets it run in air, and deleting PDH with PFL also closed stops growth. So that design was impossible, and the benchmark now closes PFL first. The lesson generalises: a design is only as good as the model's rules about which enzymes work in which conditions, and a stoichiometric model often doesn't know them.

**Misconception.** "The model says the strain makes product, so a fresh knockout will." Fong et al. (2005) built OptKnock lactate strains; they reached predictions mostly after weeks of **adaptive evolution** (growing them for many generations so faster mutants take over). FBA predicts where evolution can end up, not where a new mutant starts.

**What a biologist would ask.** Will the strain grow at all before evolving? Does regulation fight the design (the aerobic succinate strains of Lin et al. 2005 also deleted a regulator, *iclR*, which these models have no notion of)? Can cells escape through enzymes the model lacks? Is the product toxic? None of these can be read off a stoichiometric model.

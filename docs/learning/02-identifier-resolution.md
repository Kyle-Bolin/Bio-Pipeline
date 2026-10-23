---
{module: cellfactory.identifiers, date: 2026-10-23}
---

# 02 · Identifier resolution: one molecule, many names

Every biology database gives a compound its own id. D-glucose is `glc__D` in BiGG (where the models come from), `C00031` in KEGG, `CHEBI:4167` in ChEBI, `MNXM41` in MetaNetX, whose job is to cross-reference the others. A pathway copied from a KEGG map into a BiGG model only works if each id lands on the right model metabolite. When the match misses, nothing errors: the model just lacks that link and quietly predicts no product.

**Same name, different molecule.** *Stereoisomers* have the same atoms joined in the same order, but arranged differently in 3-D, like left and right hands. Enzymes are shaped to fit one hand, so D-lactate and L-lactate are made by different enzymes (E. coli's *ldhA* makes D-lactate) and a cell can't freely swap one for the other. This PR found two *epimers* (stereoisomers that differ at a single carbon), 2-dehydro-3-deoxy-D-gluconate and -D-galactonate, that share ChEBI ids in iJO1366's annotations. The resolver refuses to guess between them.

**Different name, same molecule.** In water an acid trades protons with its surroundings: acetic acid and acetate are one compound at two charge states, and at a cell's pH (~7) it is almost all acetate. Databases often file the two separately, MetaNetX merges them, and BiGG models use the dominant form. That is why a match made only through MetaNetX is flagged as a caveat. As far as I know MetaNetX's merging can also lump some stereoisomers, but I haven't checked which, so treat such a match as a lead.

**Misconception.** "A common name identifies a compound." It doesn't: this pipeline maps "lactate" to D-lactate, a choice that matters if you wanted L-lactate (the form most PLA plastic is made from).

**What a biologist would ask.** Which stereoisomer does my enzyme actually make, and can the host import, use or tolerate it? An id match says only that the model has a metabolite by that name. Whether the chemistry downstream fits it is a question for the enzyme papers.

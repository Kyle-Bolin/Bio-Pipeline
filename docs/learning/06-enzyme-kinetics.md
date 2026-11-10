---
{module: cellfactory.kinetics, date: 2026-11-10}
---

# 06 · Enzyme speed: how much protein does a flux cost?

Everything so far asks whether a route is *possible*: balanced, pointed downhill, forced by knockouts. None of it asks what it costs. In a cell, every reaction needs an enzyme, and every enzyme is protein the cell had to build. Protein is the cell's biggest expense (around half of an *E. coli* cell's dry weight), so a route that needs a mountain of enzyme is a bad route even if FBA calls it optimal.

**Two numbers.** An enzyme's **kcat** (turnover number) is how many substrate molecules one enzyme molecule converts per second when it has all the substrate it can use. Its **Km** (Michaelis constant) is the substrate concentration at which it runs at half that speed. Flux divided by kcat gives the *least* enzyme a reaction needs; below Km the enzyme spends most of its time waiting, so it needs more. A typical enzyme turns over only about ten molecules per second (Bar-Even et al. 2011), far from the textbook's famous fast ones, which is why the protein bill adds up.

**What this module does.** It looks these numbers up; it doesn't use them yet (the enzyme-constrained FBA that comes next will). Measured values come from BRENDA and SABIO-RK, the two big databases of published enzyme measurements. Where nothing was measured, a machine-learning predictor (DLKcat, TurNuP) can fill in, but those values are always labelled *predicted*: for a single enzyme they are often about tenfold off.

**How a reaction finds its enzyme.** Through its **EC number**, a code for the chemistry done (2.7.1.11 is "put a phosphate from ATP onto fructose 6-phosphate"). An EC number names a reaction type, not a protein, which causes two problems. **Isozymes**, different proteins doing the same chemistry, share one EC: *E. coli* has two phosphofructokinases, PfkA and PfkB, and an EC lookup can't tell their numbers apart. And **multi-enzyme complexes** such as pyruvate dehydrogenase carry several ECs, one per subunit's step, so there is no single kcat to look up; the module says so instead of picking one.

**Whose enzyme.** The same EC in another species is a different protein and can be much faster or slower. The module prefers wild-type measurements from the host, then the same genus, then anything, and the caveats count every step that fell back. It drops mutant enzymes (made in the lab to study the protein, not what the cell has). Genus is a crude stand-in for "closely related"; real evolutionary distance isn't modelled.

**Misconception.** "Km is the enzyme's affinity for its substrate." Only under a special assumption (substrate binding and unbinding much faster than the chemistry) does Km equal the binding constant. In general it mixes binding and catalysis, so a low Km doesn't by itself mean tight binding.

**Test tube versus cell.** Measured kcats come from purified enzyme, often at 25-30 °C rather than 37 °C, with substrate in excess. Inside the cell enzymes are regulated (PFK is slowed by phosphoenolpyruvate, for one) and often not saturated, so most run below kcat. Davidi et al. (2016) found in-cell maximum rates track test-tube kcats reasonably well across *E. coli* enzymes, but with wide scatter for individual ones, so read any one kcat as an order-of-magnitude figure.

**What a biologist would ask.** Which isozyme is the cell actually making under these conditions, and how much of it (that needs proteomics)? Was the kcat measured at a sensible temperature and pH? Is the enzyme switched down by regulation? Is the substrate concentration near or below Km? The table can't answer any of these; it only says how fast the enzyme *can* go, and where that number came from.

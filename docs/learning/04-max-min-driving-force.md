---
{module: cellfactory.mdf, date: 2026-10-28}
---

# 04 · Max-min Driving Force: can the whole route run at once?

[Note 03](03-reaction-thermodynamics.md) judged reactions one at a time. A route is a chain, and neighbouring steps share their intermediates. Step 1 wants its product, intermediate B, scarce. Step 2 wants B plentiful. Each step can be satisfied alone, but B has only one concentration. **Max-min Driving Force** (MDF; Noor et al. 2014) asks the joint question: pick one concentration per metabolite, each between 1 µM and 10 mM, so that the *least* downhill step is as downhill as possible. That best worst step, in kJ/mol, is the MDF.

**Driving force** is −ΔG': how far downhill a step is. It also sets how much of an enzyme's work is wasted. Every enzyme runs both ways, and the ratio of forward to backward events is e^(−ΔG'/RT). At 5.7 kJ/mol that ratio is 10, so one turnover in ten is undone. Near zero, the enzyme mostly spins in place, and the cell has to make far more of it to carry the same net flux. So an MDF of zero or less kills a route. A small positive MDF means a route that works on paper but needs a lot of enzyme.

**Misconception.** "If every step is possible, the route is possible." Glycolysis shows why not. Run backwards (pyruvate to glucose 6-phosphate), no single step is one-way under note 03's rule, yet the route as a whole comes out infeasible here. Cells make glucose with **bypass enzymes** instead, separate enzymes for the pyruvate-kinase and PFK steps that spend extra ATP or release phosphate. Be careful with this example, though. The infeasibility is narrow (about −2 kJ/mol). Textbooks also give regulation as a reason for the bypasses: separate enzymes let the cell switch glycolysis and glucose-making on and off independently. MDF knows nothing about regulation.

**Limits.** Steps without an estimate are left out (transport, for example), and leaving a step out can only make the MDF look better. The 1 µM–10 mM range is a convention, not a measurement, and it applies to cofactors too: it lets ATP outnumber ADP 10,000 to 1, where a living cell keeps them within about 1–10 to 1. So a positive MDF here is optimistic, and an infeasible verdict is the result to trust. Pinning measured cofactor levels (`bounds=`) tightens it.

**What a biologist would ask.** Is the bottleneck enzyme abundant enough to pay for its small driving force? What are the real concentrations? Answering either takes proteomics or metabolomics.

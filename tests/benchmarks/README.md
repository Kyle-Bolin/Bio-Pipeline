# Benchmarks: known answers

These tests check the simulator against results we already know, so a code change (or a new model) that breaks the science shows up as a failing test. Run them with:

```bash
pytest tests/benchmarks -v     # one line per known result
pytest tests/benchmarks -rx    # also print the reason for each known model failure
```

## What's covered

| Group | What it checks | Example |
|---|---|---|
| `GROWTH` | Growth rates published with each model, plus plausibility bands | iJO1366 on glucose = 0.982/h |
| `STOICHIOMETRY` | Theoretical max yields fixed by chemistry; any correct model must match | glucose → 2 ethanol |
| `KNOCKOUTS` | Strain designs from the literature | ΔadhE Δpta forces lactate (OptKnock) |
| `ESSENTIALITY` | Genes whose deletion is/isn't lethal in the lab (Keio collection) | gapA essential, pfkA not |

Cases with `known_failure` are where FBA is known to disagree with the lab, mostly because it has no gene regulation. They run as **strict xfails**: the test passes while the model keeps getting it wrong for the documented reason, and fails loudly if that changes (say a better model or method fixes it), so the note can be updated.

Models that must be downloaded from BiGG (`iML1515`, `iMM904`) are skipped when BiGG can't be reached and run everywhere else (e.g. GitHub Actions).

## Adding a case

Append to the right list in `cases.py`. Every case needs a `source`: a paper with a DOI, or "balanced chemical equation" for stoichiometric limits.

```python
Case("glucose -> 2 lactate", "iJO1366", "max_yield", 2.0, STOICH, product="lactate", aerobic=False),
Case("adhE pta forces lactate", "iJO1366", "guaranteed_yield", (1.5, 2.0), BURGARD,
     product="lactate", aerobic=False, knockouts=["ACALD", "ALCD2x", "PTAr"]),
```

* `metric` is `growth` (1/h), `max_yield` (mol/mol, no growth) or `guaranteed_yield` (mol/mol, minimum at max growth).
* `expected` is a number (checked to `tol`, default 0.01) or a `(low, high)` range. Prefer a range for anything measured in a lab, since FBA predicts an optimum, not the exact measurement.
* `knockouts` are reaction ids. To delete a gene, list every reaction it alone enables (check its GPR rule).
* Substrate and product accept the same names as the CLI (`glucose`, `etoh`, `EX_etoh_e`).

A new kind of check (say, thermodynamic driving force) gets its own list in `cases.py` and its own test function in `test_benchmarks.py`.

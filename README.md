# Bio-Pipeline: cellfactory

Feed **substrate X** to **organism Y** and ask how much **product Z** it can make, then search for gene deletions that make it produce more. Built on [COBRApy](https://opencobra.github.io/cobrapy/) and curated models from [BiGG](http://bigg.ucsd.edu). The teaching notebooks use *E. coli*, glucose and ethanol; the tool works with any model, substrate and product.

This is a starter project for someone who knows code and wants to learn the biology by doing. The README explains the biology in CS terms; the notebooks walk through the analysis step by step.

## Quick start

```bash
# with conda
conda env create -f environment.yml
conda activate cellfactory

# or with plain venv
python -m venv .venv && source .venv/bin/activate
pip install -e ".[checks,dev]"

pytest                 # unit tests + known-result benchmarks, offline (~40 s)
jupyter lab notebooks/ # then open 01_glucose_to_ethanol.ipynb
```

### Command line

```bash
# organism = a BiGG id or a model file; substrate/product = names or ids
cellfactory run --model iML1515 --substrate glucose --product ethanol --anaerobic
cellfactory run --model iYS1720 --substrate glycerol --product succinate --knockouts 3 --checks --json report.json
cellfactory run --model path/to/my_model.xml --substrate xylose --product etoh

# not sure what a compound is called in a model?
cellfactory find --model iML1515 succin
```

Output: growth rate, theoretical max yield, guaranteed and best-case yield at max growth (mol/mol and % of carbon), and optionally a greedy knockout search and model sanity checks. From Python, the same thing is `cellfactory.analyze(model, substrate, product, ...)`.

### For AI agents (MCP)

cellfactory ships an [MCP](https://modelcontextprotocol.io) server, so Claude and other agents can call it as tools and get JSON back:

```bash
pip install -e ".[mcp]"
claude mcp add cellfactory -- cellfactory-mcp   # Claude Code; the repo's .mcp.json does this automatically
```

Tools: `list_models`, `describe_model`, `search_metabolites`, `simulate_production`, `search_knockouts`, `production_envelope`, `check_model`, `knockout_candidates`. Bad input (an unknown compound, a model that can't grow) comes back as an error message saying what to try instead. The skill in `.claude/skills/cellfactory/` teaches Claude the workflow and how to read the results. The same functions are plain Python in `cellfactory.tools`, and `cellfactory run ... --json -` prints the report as JSON.

When you add a public function to the package, expose it in `tools.py` (or list it in `NOT_EXPOSED` with a reason); `tests/test_tools.py` fails until you do.

The notebooks are committed with their outputs, so you can read them on GitHub before installing anything.

| Notebook | What it does |
|---|---|
| [01_glucose_to_ethanol](notebooks/01_glucose_to_ethanol.ipynb) | Load a model, simulate growth with and without oxygen, compute ethanol yields, plot the production envelope, run a greedy knockout search on the small and the genome-scale model |
| [02_sanity_checks](notebooks/02_sanity_checks.ipynb) | Catch model bugs: free-energy loops, mass balance, gene essentiality vs lab data, and a MEMOTE report |
| [03_any_organism](notebooks/03_any_organism.ipynb) | The generic tool: compare organisms and substrates, find compound ids, handle products the organism can't make |

The reusable code lives in `src/cellfactory/`:

| Module | Contents |
|---|---|
| `models.py` | `load_model` (BiGG id or SBML/JSON/MAT file), `resolve_exchange` (names → ids), `set_medium` (swap the carbon source) |
| `simulate.py` | yield at max growth, theoretical yield, production envelope |
| `knockouts.py` | candidate selection and greedy knockout search |
| `checks.py` | gene essentiality, free-energy and mass-balance checks |
| `pipeline.py` | `analyze()`: the whole X → Y → Z run, returning a `Report` |
| `cli.py` | the `cellfactory` command |
| `tools.py` | JSON-in, JSON-out tool functions for agents |
| `mcp_server.py` | the `cellfactory-mcp` MCP server exposing `tools.py` |

Tests live in `tests/`; `tests/benchmarks/` holds known results from the literature (published growth rates, stoichiometric yield limits, OptKnock designs, Keio essential genes) with their sources. See [its README](tests/benchmarks/README.md) to add cases.

### Which model?

Any BiGG id works (browse [bigg.ucsd.edu/models](http://bigg.ucsd.edu/models): yeast `iMM904`, *Bacillus subtilis* `iYO844`, *Pseudomonas putida* `iJN1463`, ...), as does any SBML file. These three are bundled and work offline:

* **`iYS1720`**: *Salmonella* Typhimurium LT2, a close relative of *E. coli*.
* **`e_coli_core`** (95 reactions, 137 genes): central metabolism only. Small enough to read in full; use it to learn and for tests.
* **`iJO1366`** (2,583 reactions, 1,367 genes): full *E. coli* K-12 MG1655 (2011).

The current reference *E. coli* model is **`iML1515`** (2,712 reactions, 1,516 genes, [Monk et al. 2017](https://doi.org/10.1038/nbt.3956)), downloaded from BiGG on first use. The notebooks call `load_model("iML1515", fallback="iJO1366")` so they still run offline.

## The biology, for a CS person

### The cell as a chemical factory

A bacterium takes in raw material (glucose, a 6-carbon sugar) and converts it, through hundreds of chemical **reactions**, into the parts it needs to make a copy of itself: amino acids for proteins, nucleotides for DNA and RNA, lipids for membranes. That network of reactions is its **metabolism**.

The molecules flowing through it are **metabolites**. Each reaction is sped up by an **enzyme**, a protein. Without the enzyme the reaction is effectively off.

### Genes are the source code for enzymes

A **gene** is a stretch of DNA that encodes a protein. Remove (**knock out**) the gene and the cell can't make that enzyme, so the reaction stops. That's the lever metabolic engineers pull.

The mapping isn't one-to-one, so models store it as a boolean **gene-protein-reaction (GPR) rule** per reaction:

* `b0001 and b0002`: the enzyme is a **complex** built from both proteins. Delete either and the reaction stops.
* `b0001 or b0002`: **isozymes**, two different enzymes that do the same job. You must delete both.

*E. coli* genes have names like *ldhA*, but models identify them by **b-numbers** (e.g. `b1380`), stable locus tags from the original genome sequence.

### Energy and electrons: ATP and NADH

Two "currency" molecules show up everywhere:

* **ATP** is the energy currency. Making things costs ATP; breaking down sugar earns it.
* **NADH / NAD⁺** is a rechargeable battery for electrons. Breaking down glucose (**glycolysis**) strips electrons off and stores them as NADH. The battery must be emptied back to NAD⁺, or glycolysis halts once NAD⁺ runs out.

With **oxygen**, the cell empties NADH onto O₂ (**respiration**), which also pays out a lot of ATP, and burns glucose all the way to CO₂.

Without oxygen, the cell has to empty its NADH onto a carbon molecule it makes from glucose and then throws away. That's **fermentation**, and the waste molecules are its products: ethanol, lactate, acetate, succinate. This is why the production problem is a *redox balancing* problem: every route has to come out even on electrons.

### What we want: yield and growth coupling

**Yield** is mol product per mol glucose consumed. For ethanol the ceiling is 2: C₆H₁₂O₆ → 2 C₂H₅OH + 2 CO₂.

Cells evolve to grow fast, not to make our product. So the useful trick is **growth coupling**: delete competing routes until the *only* way for the cell to grow is to make the product. Then the fastest-growing cells are also the best producers, and evolution works for you instead of against you.

### Compartments

Cells have inside and outside. Model metabolite ids carry a suffix: `_c` cytosol (inside), `_p` periplasm (the gap between *E. coli*'s two membranes), `_e` extracellular (the medium). Moving a molecule between compartments is a **transport** reaction, done by a membrane protein.

## The math: flux balance analysis

Let **S** be the stoichiometric matrix: rows are metabolites, columns are reactions, entries are how many of each metabolite a reaction consumes (negative) or produces (positive). Let **v** be the vector of reaction rates (**fluxes**, in mmol per gram dry cell weight per hour).

* **Steady state:** internal metabolites don't pile up, so **S v = 0**.
* **Bounds:** `lb ≤ v ≤ ub`. Irreversible reactions have `lb = 0`. Uptake is capped by what's in the medium (e.g. glucose ≤ 10).
* **Objective:** maximise the **biomass reaction**, a pseudo-reaction that consumes all the building blocks in the proportions found in a real cell. Its flux is the growth rate (1/h).

That's a linear program. A genome-scale one with ~2,700 variables solves in milliseconds, which is what makes brute-force knockout searches feasible.

Related tools used in the notebooks:

* **Flux variability analysis (FVA):** the LP usually has many optimal solutions. FVA finds the min and max of each flux over all of them. We use the *minimum* product flux at max growth as the "guaranteed" yield.
* **Production envelope:** for each growth rate, the min and max product flux. It shows the trade-off between growing and producing.
* **Exchange reactions** (`EX_glc__D_e`, `EX_etoh_e`): the model's I/O ports. Negative flux = uptake, positive = secretion.

### What FBA can't see

FBA knows what the network *can* do, not what the cell *chooses* to do. It ignores:

* **Gene regulation:** cells switch enzymes on and off depending on conditions. A route the model happily uses may be shut off in reality.
* **Kinetics:** how fast an enzyme actually works, and how much of it there is.
* **Toxicity:** high ethanol kills cells; the model doesn't care.

That's why the sanity checks in notebook 02 matter: in `iJO1366`, the knockout search's favourite first move (deleting `ppc`) is one that the gene-essentiality check shows the model gets wrong, because real cells can't use the bypass the model finds. The benchmark suite records these known failures explicitly.

## Next steps

* Try other targets: succinate, lactate, or a non-native compound like lycopene (add the CrtE, CrtB, CrtI reactions yourself; `analyze` then adds a demand reaction for it).
* Replace greedy search with OptKnock or a proper strain-design package ([StrainDesign](https://github.com/klamt-lab/straindesign)).
* Run full gene essentiality on every gene and score it against the complete Keio data set.
* Repeat with yeast (`--model iMM904`), the other workhorse of sugar fermentation.

## Further reading

* Orth, Thiele & Palsson, [What is flux balance analysis?](https://doi.org/10.1038/nbt.1614) *Nature Biotechnology* 2010. Four pages, the best first read.
* Ebrahim et al., [COBRApy](https://doi.org/10.1186/1752-0509-7-74), 2013, and the [COBRApy docs](https://cobrapy.readthedocs.io).
* Lieven et al., [MEMOTE for standardized genome-scale metabolic model testing](https://doi.org/10.1038/s41587-020-0446-y), 2020.
* Burgard, Pharkya & Maranas, [OptKnock](https://doi.org/10.1002/bit.10803), 2003: knockout search as a bilevel optimisation problem.

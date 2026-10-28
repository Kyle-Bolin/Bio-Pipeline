"""Known answers the simulator must reproduce.

Each case is data: what to simulate, what number to expect, and where that
number comes from. To add one, append to a list below; the tests in
test_benchmarks.py pick it up automatically. See README.md in this folder.

``known_failure`` marks cases where FBA is *known* to disagree with the lab.
They are run as strict xfails: they must keep failing for the documented
reason, and if a model change makes one pass, the suite tells you.
"""

from __future__ import annotations

from dataclasses import dataclass, field

ORTH_CORE = "Orth, Fleming & Palsson 2010, EcoSal Plus, doi:10.1128/ecosalplus.10.2.1 (e_coli_core paper)"
ORTH_2011 = "Orth et al. 2011, Mol Syst Biol 7:535, doi:10.1038/msb.2011.65 (iJO1366 paper)"
MONK_2017 = "Monk et al. 2017, Nat Biotechnol 35:904, doi:10.1038/nbt.3956 (iML1515 paper)"
BURGARD = "Burgard, Pharkya & Maranas 2003, Biotechnol Bioeng 84:647, doi:10.1002/bit.10803 (OptKnock)"
FONG = "Fong et al. 2005, Biotechnol Bioeng 91:643, doi:10.1002/bit.20542 (OptKnock lactate strains built and evolved)"
NZN111 = "Bunch, Mat-Jan & Clark 1997, Microbiology 143:187, doi:10.1099/00221287-143-1-187 (pfl ldhA strain NZN111)"
CLARK = "Clark 1989, FEMS Microbiol Rev 63:223, doi:10.1111/j.1574-6968.1989.tb03398.x (E. coli mixed-acid fermentation)"
STOICH = "Balanced chemical equation (no model needed)"
KEIO = "Baba et al. 2006, Mol Syst Biol 2:2006.0008, doi:10.1038/msb4100050 (Keio knockout collection)"


@dataclass
class Case:
    id: str
    model: str
    metric: str  # "growth", "max_yield", "guaranteed_yield"
    expected: float | tuple[float, float]  # a value (checked with ``tol``) or an inclusive (low, high) range
    source: str
    substrate: str = "glucose"
    product: str | None = None
    aerobic: bool = True
    knockouts: list[str] = field(default_factory=list)  # reaction ids
    tol: float = 0.01
    known_failure: str | None = None  # why FBA gets this wrong; the case runs as a strict xfail


# Growth rates (1/h) at 10 mmol glucose/gDW/h on the model's minimal medium.
GROWTH = [
    Case("core aerobic growth", "e_coli_core", "growth", 0.874, ORTH_CORE),
    Case("iJO1366 aerobic growth", "iJO1366", "growth", 0.982, ORTH_2011),
    Case("iML1515 aerobic growth", "iML1515", "growth", 0.877, MONK_2017),
    # Measured MG1655 growth on glucose minimal medium is ~0.6-0.7/h; FBA assumes optimal
    # use of every enzyme, so it should land at or somewhat above that, never far off.
    Case("iJO1366 growth is biologically plausible", "iJO1366", "growth", (0.6, 1.1),
         "Plausibility band: measured E. coli K-12 growth on glucose minimal medium is ~0.6-0.7/h"),
    Case("anaerobic growth is much slower", "iJO1366", "growth", (0.1, 0.4),
         "Plausibility band: fermentation makes ~2 ATP per glucose vs ~20-30 for respiration", aerobic=False),
]

# Maximum theoretical yields (mol/mol): set by chemistry, so every correct model must match.
STOICHIOMETRY = [
    Case("glucose -> 2 ethanol + 2 CO2", "iJO1366", "max_yield", 2.0, STOICH, product="ethanol", aerobic=False),
    Case("glucose -> 2 lactate", "iJO1366", "max_yield", 2.0, STOICH, product="lactate", aerobic=False),
    Case("3 xylose -> 5 ethanol + 5 CO2", "iJO1366", "max_yield", 5 / 3, STOICH, substrate="xylose",
         product="ethanol", aerobic=False),
    # 7 glucose + 6 CO2 -> 12 succinate is the electron-limited maximum, widely quoted as 1.71 mol/mol.
    Case("succinate ceiling 12/7", "iJO1366", "max_yield", 12 / 7, STOICH, product="succinate"),
    Case("core ethanol ceiling", "e_coli_core", "max_yield", 2.0, STOICH, product="ethanol", aerobic=False),
    Case("yeast glucose -> 2 ethanol", "iMM904", "max_yield", 2.0, STOICH, product="ethanol"),
]

# Knockout designs from the literature.
KNOCKOUTS = [
    # Wild type has no reason to make lactate; deleting the ethanol (adhE) and acetate (pta)
    # routes forces it, the first OptKnock design, later built and evolved in the lab.
    Case("wild type makes no lactate", "iJO1366", "guaranteed_yield", 0.0, BURGARD, product="lactate",
         aerobic=False),
    Case("adhE pta forces lactate (OptKnock)", "iJO1366", "guaranteed_yield", (1.5, 2.0), f"{BURGARD}; {FONG}",
         product="lactate", aerobic=False, knockouts=["ACALD", "ALCD2x", "PTAr"]),
    Case("adhE pta forces lactate, core model", "e_coli_core", "guaranteed_yield", (1.5, 2.0), BURGARD,
         product="lactate", aerobic=False, knockouts=["ACALD", "ALCD2x", "PTAr"]),
    # Wild type ferments to a mix that includes ethanol.
    Case("anaerobic wild type secretes ethanol", "iJO1366", "guaranteed_yield", (0.3, 1.2), CLARK,
         product="ethanol", aerobic=False),
    Case("pfl ldhA mutant can't ferment glucose", "iJO1366", "growth", 0.0, NZN111, aerobic=False,
         knockouts=["PFL", "LDH_D"],
         known_failure="The model reoxidises NADH by running pyruvate dehydrogenase anaerobically; in real cells "
                       "PDH is inhibited by the high NADH/NAD+ ratio, which FBA can't represent."),
]

CASES = GROWTH + STOICHIOMETRY + KNOCKOUTS


@dataclass
class EssentialityCase:
    gene: str  # model gene id (E. coli b-number)
    name: str
    essential: bool
    model: str = "iJO1366"
    known_failure: str | None = None


_BYPASS = ("The model routes glucose around the deleted step (Entner-Doudoroff pathway plus the "
           "{route}); real mutants can't grow on glucose because those routes are too slow, toxic, "
           "or switched off by regulation, none of which FBA represents.")

# Aerobic glucose minimal medium (Keio collection + minimal-medium screens).
ESSENTIALITY = [
    EssentialityCase("b1779", "gapA", True, known_failure=_BYPASS.format(route="toxic methylglyoxal bypass")),
    EssentialityCase("b2926", "pgk", True, known_failure=_BYPASS.format(route="toxic methylglyoxal bypass")),
    EssentialityCase("b2779", "eno", True, known_failure=_BYPASS.format(route="serine deaminase and PEP synthase")),
    EssentialityCase("b0720", "gltA", True),
    EssentialityCase("b0004", "thrC", True),
    EssentialityCase("b3960", "argH", True),
    EssentialityCase("b4177", "purA", True),
    EssentialityCase("b1281", "pyrF", True),
    EssentialityCase("b3956", "ppc", True,
                     known_failure="The model refills oxaloacetate via the glyoxylate shunt, which real cells "
                                   "keep repressed on glucose; FBA has no gene regulation."),
    EssentialityCase("b4025", "pgi", False),
    EssentialityCase("b3916", "pfkA", False),
    EssentialityCase("b1852", "zwf", False),
    EssentialityCase("b1380", "ldhA", False),
    EssentialityCase("b2297", "pta", False),
    EssentialityCase("b1241", "adhE", False),
    EssentialityCase("b1101", "ptsG", False),
    # The small core model, for the genes it contains.
    EssentialityCase("b1779", "gapA", True, model="e_coli_core"),
    EssentialityCase("b2779", "eno", True, model="e_coli_core"),
    EssentialityCase("b0720", "gltA", True, model="e_coli_core"),
    EssentialityCase("b3956", "ppc", True, model="e_coli_core",
                     known_failure="The core model refills oxaloacetate via the glyoxylate shunt and PEP "
                                   "carboxykinase; real cells keep these off on glucose."),
    EssentialityCase("b4025", "pgi", False, model="e_coli_core"),
    EssentialityCase("b1380", "ldhA", False, model="e_coli_core"),
]


LIEVEN = "Lieven et al. 2020, Nat Biotechnol 38:272, doi:10.1038/s41587-020-0446-y (MEMOTE)"
FRITZEMEIER = ("Fritzemeier et al. 2017, PLoS Comput Biol 13:e1005494, doi:10.1371/journal.pcbi.1005494 "
               "(energy-generating cycles in published models)")


@dataclass
class QACase:
    id: str
    model: str
    metric: str  # "gate_ok": passes the QA gate; "memote_egc_free": raw MEMOTE finds no energy-generating cycle
    expected: bool
    source: str
    add_reactions: tuple[str, ...] = ()  # reaction strings to add first, to plant a known defect
    expect_hard_fail: str | None = None  # the MEMOTE test that must be among the hard failures
    known_failure: str | None = None


# The model QA gate (qa.py). Needs MEMOTE; skipped where it isn't installed.
QA = [
    QACase("e_coli_core passes the QA gate", "e_coli_core", "gate_ok", True, f"{ORTH_CORE}; {LIEVEN}"),
    QACase("iJO1366 passes the QA gate", "iJO1366", "gate_ok", True, f"{ORTH_2011}; {FRITZEMEIER}"),
    QACase("free ATP regeneration is caught", "e_coli_core", "gate_ok", False, FRITZEMEIER,
           add_reactions=("adp_c + pi_c + h_c --> atp_c + h2o_c",),
           expect_hard_fail="test_detect_energy_generating_cycles"),
    # Two isolated metabolites with a = 2b and b = a: no positive masses satisfy both. They touch no
    # cofactor, so this exercises the consistency check alone, with no energy loop to trip the gate first.
    QACase("mass created from nothing is caught", "e_coli_core", "gate_ok", False, LIEVEN,
           add_reactions=("planted_a_c --> 2 planted_b_c", "planted_b_c --> planted_a_c"),
           expect_hard_fail="test_stoichiometric_consistency"),
    QACase("raw MEMOTE finds no energy cycle in iJO1366", "iJO1366", "memote_egc_free", True, FRITZEMEIER,
           known_failure="MEMOTE resets every reaction to its reversibility, reopening the ones iJO1366 ships "
                         "closed. With superoxide dismutase (SPODM) back on, malate oxidase runs backwards to turn "
                         "H2O2 into O2 and the model respires, making ATP from nothing. The gate re-checks with "
                         "the shipped bounds and passes the model."),
]


MORETTI = ("Moretti et al. 2021, Nucleic Acids Res 49:D570, doi:10.1093/nar/gkaa992 (MetaNetX/MNXref: "
           "cross-references between metabolite namespaces)")
KING_2016 = "King et al. 2016, Nucleic Acids Res 44:D515, doi:10.1093/nar/gkv1049 (BiGG Models database)"


@dataclass
class IdCase:
    id: str
    model: str
    query: str  # what a user types: a database id from KEGG, ChEBI, MetaNetX, ...
    expected: str | tuple[str, ...]  # the model's base metabolite id, or the candidates it must refuse to choose between
    source: str
    known_failure: str | None = None


# Identifier resolution (identifiers.py), offline: only the model's own annotations, MetaNetX never called.
IDENTIFIERS = [
    IdCase("KEGG glucose in the core model", "e_coli_core", "kegg:C00031", "glc__D", KING_2016),
    IdCase("ChEBI ethanol in the core model", "e_coli_core", "CHEBI:16236", "etoh", KING_2016),
    IdCase("MetaNetX glucose in iJO1366", "iJO1366", "MNXM41", "glc__D", MORETTI),
    IdCase("bare KEGG id for succinate in iJO1366", "iJO1366", "C00042", "succ", KING_2016),
    IdCase("old-style HMDB id for D-glucose in iJO1366", "iJO1366", "HMDB00122", "glc__D", KING_2016),
    IdCase("InChIKey for D-glucose in Salmonella", "iYS1720", "WQZGKKKJIJFFOK-GASJEMHNSA-N", "glc__D", KING_2016),
    # Two different sugar acids (epimers: same atoms, one stereocentre flipped) share ChEBI ids in iJO1366's
    # annotations, so the id alone cannot say which is meant. The resolver must refuse, not pick one.
    IdCase("shared ChEBI id is refused", "iJO1366", "CHEBI:57990", ("2ddglcn", "2dh3dgal"), KING_2016),
    IdCase("generic glucose ChEBI id", "iJO1366", "CHEBI:17234", "glc__D", MORETTI,
           known_failure="CHEBI:17234 is 'glucose' with its stereochemistry unspecified; the model annotates only "
                         "the specific D-glucose ids. ChEBI links them by an is-a relation that neither the model's "
                         "annotations nor an exact cross-reference follows, so offline resolution misses it."),
]


LEHNINGER = ("Nelson & Cox, Lehninger Principles of Biochemistry: standard free-energy changes of the glycolytic and "
             "citric acid cycle reactions, and the table of standard reduction potentials (pH 7, 25 °C)")
BEBER_2022 = "Beber et al. 2022, Nucleic Acids Res 50:D603, doi:10.1093/nar/gkab1106 (eQuilibrator 3.0)"
PARK_2016 = ("Park et al. 2016, Nat Chem Biol 12:482, doi:10.1038/nchembio.2077 (measured metabolite concentrations and "
             "in vivo ΔG of central metabolism)")


@dataclass
class ThermoCase:
    id: str
    model: str
    reaction: str
    metric: str  # "dg0_prime" (kJ/mol) or "direction" ("forward", "reverse", "either")
    expected: tuple[float, float] | str
    source: str
    known_failure: str | None = None


# Reaction thermodynamics (thermodynamics.py), offline from eQuilibrator's bundled compound data. Textbook values are
# at pH 7 and zero ionic strength, eQuilibrator's at pH 7.5 and 0.25 M, so the ranges are wide on purpose.
THERMO = [
    ThermoCase("PFK is downhill", "e_coli_core", "PFK", "dg0_prime", (-25.0, -10.0), LEHNINGER),  # -14.2
    ThermoCase("aldolase is uphill at standard conditions", "e_coli_core", "FBA", "dg0_prime", (15.0, 30.0), LEHNINGER),  # +23.8
    ThermoCase("triose phosphate isomerase is mildly uphill", "e_coli_core", "TPI", "dg0_prime", (2.0, 10.0), LEHNINGER),  # +7.5
    ThermoCase("pyruvate kinase is downhill", "e_coli_core", "PYK", "dg0_prime", (-35.0, -20.0), LEHNINGER),  # -31.4
    ThermoCase("malate dehydrogenase is uphill", "e_coli_core", "MDH", "dg0_prime", (20.0, 35.0), LEHNINGER),  # +29.7
    ThermoCase("citrate synthase is downhill", "e_coli_core", "CS", "dg0_prime", (-45.0, -25.0), LEHNINGER),  # -32.2
    ThermoCase("ATP hydrolysis", "e_coli_core", "ATPM", "dg0_prime", (-36.0, -25.0), LEHNINGER),  # -30.5
    # A positive ΔG'° does not make a reaction run backwards: glycolysis pushes aldolase forward by keeping
    # its products scarce. The concentration-range rule has to leave it reversible.
    ThermoCase("aldolase stays reversible despite positive ΔG'°", "e_coli_core", "FBA", "direction", "either",
               LEHNINGER),
    ThermoCase("malate dehydrogenase stays reversible", "e_coli_core", "MDH", "direction", "either", LEHNINGER),
    ThermoCase("ATP hydrolysis is one-way", "e_coli_core", "ATPM", "direction", "forward", LEHNINGER),
    # The reaction behind iJO1366's MEMOTE energy cycle (see QA above): run backwards it turns H2O2 into O2.
    # O2/H2O2 (+0.30 V) against oxaloacetate/malate (-0.17 V) is about -90 kJ/mol forward.
    ThermoCase("malate oxidase can't make O2 from H2O2", "iJO1366", "MOX", "direction", "forward",
               f"{LEHNINGER}; {BEBER_2022}"),
    ThermoCase("pyruvate kinase is one-way in the cell", "e_coli_core", "PYK", "direction", "forward", PARK_2016,
               known_failure="Pyruvate kinase is effectively irreversible in vivo because the cell holds its "
                             "substrates and products far from equilibrium. ΔG'° alone (about -25 kJ/mol) can be "
                             "overturned by a 10^4-fold concentration ratio on each side, so the 1 µM-10 mM rule "
                             "cannot prove it one-way. Measured concentrations would."),
    # Fumarate/succinate (+0.03 V) and ubiquinone/ubiquinol (+0.045 V) are close, so succinate reducing ubiquinone
    # is only a few kJ/mol downhill. Midpoint potentials of membrane quinones are themselves debated.
    ThermoCase("succinate -> ubiquinone is near equilibrium", "e_coli_core", "SUCDi", "dg0_prime", (-20.0, 0.0),
               LEHNINGER,
               known_failure="eQuilibrator estimates ubiquinone-8 by group contribution, for compounds in water, "
                             "though ubiquinone sits in the membrane. It puts this reaction at about -74 kJ/mol. That makes the model's FRD7 (the "
                             "reverse) look impossible and is why apply_directionality skips group-contribution "
                             "estimates by default."),
]


@dataclass
class MDFCase:
    id: str
    model: str
    metric: str  # "mdf" (kJ/mol, point estimate) or "verdict"
    expected: tuple[float, float] | str
    source: str
    product: str | None = None  # screen the max-yield route to this product...
    fluxes: dict[str, float] | None = None  # ...or this hand-written route (reaction id -> flux)
    aerobic: bool = True
    known_failure: str | None = None


# The glycolysis steps of e_coli_core, glucose 6-phosphate to pyruvate. PGK and PGM are written backwards in BiGG.
GLYCOLYSIS = {"PGI": 1, "PFK": 1, "FBA": 1, "TPI": 1, "GAPD": 2, "PGK": -2, "PGM": -2, "ENO": 2, "PYK": 2}

# Max-min Driving Force (mdf.py) over whole routes, from eQuilibrator's bundled data. A route a cell is known to run
# must come out positive; one textbooks say needs bypass enzymes must not.
MDF = [
    MDFCase("glycolysis can run forward", "e_coli_core", "verdict", "feasible", LEHNINGER, fluxes=GLYCOLYSIS),
    # Every step is "either direction" on its own (see THERMO), but not all at once: cells make glucose through
    # bypass enzymes (PEP synthase or PEP carboxykinase, fructose-1,6-bisphosphatase) instead of reversing PYK and PFK.
    MDFCase("glycolysis can't simply run backwards", "e_coli_core", "verdict", "infeasible", LEHNINGER,
            fluxes={r: -v for r, v in GLYCOLYSIS.items()}),
    # The verdict above also needs the optimistic MDF <= 0, which is narrow (about -1 kJ/mol); pin the point value too.
    MDFCase("glycolysis backwards has negative MDF", "e_coli_core", "mdf", (-20.0, 0.0), LEHNINGER,
            fluxes={r: -v for r, v in GLYCOLYSIS.items()}),
    MDFCase("anaerobic ethanol route can run", "e_coli_core", "mdf", (0.0, 100.0), CLARK, product="ethanol",
            aerobic=False),
    MDFCase("anaerobic ethanol route can run, iJO1366", "iJO1366", "mdf", (0.0, 100.0), CLARK, product="ethanol",
            aerobic=False),
    MDFCase("anaerobic lactate route can run", "iJO1366", "mdf", (0.0, 100.0), f"{CLARK}; {FONG}", product="lactate",
            aerobic=False),
    # E. coli makes succinate when fermenting, through fumarate reductase (FRD7 in this model).
    MDFCase("anaerobic succinate route can run", "e_coli_core", "mdf", (0.0, 100.0), CLARK, product="succinate",
            aerobic=False,
            known_failure="The route needs fumarate reductase (FRD7), and eQuilibrator's group-contribution estimate "
                          "for the ubiquinone couple puts it ~74 kJ/mol uphill (see the SUCDi case in THERMO). So MDF "
                          "calls the route 'uncertain' with FRD7 as the bottleneck, not feasible."),
]

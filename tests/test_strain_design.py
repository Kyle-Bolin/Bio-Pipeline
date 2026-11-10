"""Knockout design with StrainDesign: real MILPs on e_coli_core (seconds with SCIP), failure paths mocked."""

import json

import pytest

from cellfactory import (DesignUnavailable, analyze, design_knockouts, greedy_knockout_search, load_model,
                         resolve_exchange, set_medium)
from cellfactory import strain_design
from cellfactory.cli import main
from cellfactory.knockouts import candidate_reactions


@pytest.fixture
def core():
    return load_model("e_coli_core")


def setup(model, product, aerobic):
    substrate = set_medium(model, "glucose", 10.0, aerobic)
    return resolve_exchange(model, product, create_demand=True, online=False), substrate


def test_gene_set_complex_needs_one_gene_isozymes_need_all(core):
    # PDH is one complex (aceE and aceF and lpd): deleting any one subunit is enough, and the pick avoids lpd
    # (b0116), which is also part of 2-oxoglutarate dehydrogenase.
    pdh = strain_design._gene_set(core.reactions.PDH)
    assert len(pdh) == 1 and pdh != ["b0116"]
    # ACALD has two isozymes (mhpF or adhE): both have to go.
    assert strain_design._gene_set(core.reactions.ACALD) == ["b0351", "b1241"]


def test_unknown_method_and_bad_limits(core):
    product, substrate = setup(core, "lactate", False)
    with pytest.raises(ValueError, match="unknown design method"):
        design_knockouts(core, product, substrate, method="optgene")
    with pytest.raises(ValueError, match="at least 1"):
        design_knockouts(core, product, substrate, max_knockouts=0)
    with pytest.raises(ValueError, match="unknown design method"):
        analyze(core, "glucose", "lactate", aerobic=False, qa_gate=False, online_ids=False, design="optgene")


def test_missing_straindesign_is_a_clear_error_and_a_caveat(core, monkeypatch):
    def missing():
        raise DesignUnavailable("StrainDesign is not installed (pip install 'cellfactory[design]')")
    monkeypatch.setattr(strain_design, "_straindesign", missing)
    product, substrate = setup(core, "lactate", False)
    with pytest.raises(DesignUnavailable, match="not installed"):
        design_knockouts(core, product, substrate)
    report = analyze(core, "glucose", "lactate", aerobic=False, qa_gate=False, online_ids=False, design="mcs")
    assert any("Knockout design (Minimal cut sets) not computed" in c for c in report.caveats)
    assert not any(k.startswith("Knockout designs") for k in report.checks)


def test_glpk_is_never_picked_automatically(monkeypatch):
    sd = pytest.importorskip("straindesign")
    monkeypatch.setattr(sd, "avail_solvers", {"glpk"})
    with pytest.raises(DesignUnavailable, match="GLPK"):
        strain_design.milp_solver()
    assert strain_design.milp_solver("glpk") == "glpk"
    with pytest.raises(DesignUnavailable, match="not installed"):
        strain_design.milp_solver("gurobi")
    monkeypatch.setattr(sd, "avail_solvers", {"glpk", "scip", "cplex"})
    assert strain_design.milp_solver() == "cplex"


def test_solver_failure_retries_with_big_m_then_reports_error(core, monkeypatch):
    sd = pytest.importorskip("straindesign")
    product, substrate = setup(core, "lactate", False)
    calls = []

    def broken(model, **kwargs):
        calls.append(kwargs.get("M"))
        raise IndexError("list index out of range")
    monkeypatch.setattr(sd, "compute_strain_designs", broken)
    result = design_knockouts(core, product, substrate, method="robustknock")
    assert calls == [None, strain_design.BIG_M]
    # A malformed problem is not solver trouble: it is raised, not retried.
    calls.clear()

    def bad_problem(model, **kwargs):
        calls.append(kwargs.get("M"))
        raise ValueError("can't parse constraint")
    monkeypatch.setattr(sd, "compute_strain_designs", bad_problem)
    with pytest.raises(ValueError, match="parse"):
        design_knockouts(core, product, substrate, method="robustknock")
    assert calls == [None]
    monkeypatch.setattr(sd, "compute_strain_designs", broken)
    assert result.status == "error" and result.designs.empty
    caveats = strain_design.design_caveats(result)
    assert len(caveats) == 1 and "not computed" in caveats[0] and "IndexError" in caveats[0]


def test_wild_type_too_slow_is_not_sent_to_the_solver(core):
    pytest.importorskip("straindesign")
    product, substrate = setup(core, "lactate", False)
    result = design_knockouts(core, product, substrate, min_growth=5.0)
    assert result.status == "infeasible" and "below" in result.notes[0]


def test_optknock_finds_aerobic_succinate_design_greedy_misses(core):
    # PFL is destroyed by oxygen in the cell but open in e_coli_core, so close it as a real aerobic cell would.
    # Then no single deletion forces succinate out aerobically, so greedy stops at the wild type; the MILP closes
    # succinate dehydrogenase (or fumarase) and the acetate route together.
    pytest.importorskip("straindesign")
    product, substrate = setup(core, "succinate", True)
    core.reactions.PFL.knock_out()
    greedy = greedy_knockout_search(core, product, substrate, max_knockouts=3, candidates=candidate_reactions(core))
    assert greedy.steps[-1].guaranteed_yield < 1e-6
    before = {r.id: r.bounds for r in core.reactions}
    result = design_knockouts(core, product, substrate, method="optknock")
    assert {r.id: r.bounds for r in core.reactions} == before  # the model is left alone
    best = result.best
    assert result.status == "optimal" and result.wild_type_yield < 1e-6
    assert best["guaranteed_yield"] > 0.5 and best["growth"] >= 0.05
    deleted = set(best["knockouts"].split(", "))
    assert deleted & {"SUCDi", "FUM"} and deleted & {"ACKr", "PTAr"}
    caveats = strain_design.design_caveats(result, greedy)
    assert any("beats the greedy search" in c for c in caveats)


def test_reactions_already_closed_are_not_side_effects(core):
    # TKT2 shares its genes with TKT1. Close TKT2 in the "wild type": deleting TKT1's genes then changes nothing more.
    product, substrate = setup(core, "lactate", False)
    row = strain_design._check(core, ["TKT1"], product, substrate, "Biomass_Ecoli_core", 0.05)
    assert row["also_disabled"] == "TKT2"
    core.reactions.TKT2.knock_out()
    row = strain_design._check(core, ["TKT1"], product, substrate, "Biomass_Ecoli_core", 0.05)
    assert row["also_disabled"] == ""


def test_reaction_without_genes_gets_no_gene_level_numbers(core):
    product, substrate = setup(core, "lactate", False)
    core.reactions.ACALD.gene_reaction_rule = ""
    row = strain_design._check(core, ["ACALD"], product, substrate, "Biomass_Ecoli_core", 0.05)
    assert row["genes"] == "" and "no genes: ACALD" in row["also_disabled"]
    assert row["gene_level_guaranteed_yield"] != row["gene_level_guaranteed_yield"]  # NaN


def test_quiet_drops_only_straindesign_records(caplog):
    import logging
    with caplog.at_level(logging.WARNING), strain_design._quiet():
        logging.getLogger("someone.else").warning("keep me")
        record = logging.LogRecord("root", logging.WARNING, "/x/site-packages/straindesign/compute.py", 1, "drop me",
                                   None, None)
        logging.root.handle(record)
    messages = [r.getMessage() for r in caplog.records]
    assert "keep me" in messages and "drop me" not in messages


def test_mcs_gives_strong_coupling_and_gene_level_check(core):
    pytest.importorskip("straindesign")
    product, substrate = setup(core, "lactate", False)
    result = design_knockouts(core, product, substrate, method="mcs", min_yield=1.5)
    designs = result.designs
    assert result.status == "optimal" and len(designs)
    assert (designs["coupling"] == "strong").all()  # the point of cut sets: forced at every growth rate
    assert (designs["guaranteed_yield"] >= 1.5 - 1e-6).all()
    # ACALD's isozyme adhE also runs ALCD2x, so the gene-level design deletes both and must still hold.
    assert "b1241" in designs.iloc[0]["genes"]
    assert designs.iloc[0]["gene_level_guaranteed_yield"] >= 1.5 - 1e-6


def test_gene_level_side_effects_rank_lethal_designs_last(core):
    # TKT1 and TKT2 share their genes (tktA or tktB), so a design deleting TKT1 also loses TKT2 and stops growing.
    pytest.importorskip("straindesign")
    product, substrate = setup(core, "lactate", False)
    result = design_knockouts(core, product, substrate, method="robustknock", max_designs=5)
    d = result.designs
    lethal = d[d["also_disabled"].str.contains("TKT")]
    if lethal.empty:
        pytest.skip("solver returned no design touching transketolase this time")
    assert (lethal["gene_level_growth"] < result.min_growth).all()
    assert lethal["rank"].min() > d[d["also_disabled"] == ""]["rank"].max()
    assert any("stop growth" in c for c in strain_design.design_caveats(result))


def test_no_design_is_reported_honestly(core):
    # Aerobic E. coli has no way to be forced into ethanol at 30% of the ceiling with three deletions.
    pytest.importorskip("straindesign")
    product, substrate = setup(core, "ethanol", True)
    result = design_knockouts(core, product, substrate, method="mcs")
    assert result.status == "infeasible" and result.designs.empty
    assert "no design makes every growing state yield" in strain_design.design_caveats(result)[0]


def test_analyze_and_cli(core, capsys):
    pytest.importorskip("straindesign")
    report = analyze(core, "glucose", "lactate", aerobic=False, qa_gate=False, online_ids=False, design="mcs",
                     max_knockouts=2)
    table = report.checks["Knockout designs (Minimal cut sets)"]
    assert len(table) and "Knockout design (Minimal cut sets), solved with" in " ".join(report.caveats)
    assert any("greedy search" in c for c in report.caveats)  # both ran, so they are compared
    json.dumps(report.to_dict(), default=float)

    assert main(["run", "-p", "lactate", "--anaerobic", "--no-qa", "--offline", "--design", "mcs",
                 "--design-knockouts", "2", "--design-time-limit", "60", "--json", "-"]) == 0
    out = json.loads(capsys.readouterr().out)
    assert out["checks"]["Knockout designs (Minimal cut sets)"][0]["coupling"] == "strong"

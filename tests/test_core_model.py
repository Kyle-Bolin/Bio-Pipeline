"""Fast tests on the 95-reaction core model; they need no internet."""

import pytest

from sugar2compound import (
    free_energy_checks,
    gene_essentiality_check,
    greedy_knockout_search,
    load_model,
    max_theoretical_yield,
    production_at_max_growth,
    set_glucose_medium,
)


@pytest.fixture
def anaerobic_core():
    model = load_model("e_coli_core")
    set_glucose_medium(model, aerobic=False)
    return model


def test_aerobic_growth_is_realistic():
    model = load_model("e_coli_core")
    set_glucose_medium(model, aerobic=True)
    assert model.slim_optimize() == pytest.approx(0.874, abs=0.01)


def test_ethanol_theoretical_yield_is_two(anaerobic_core):
    # glucose (C6) -> 2 ethanol (C2) + 2 CO2
    assert max_theoretical_yield(anaerobic_core, "EX_etoh_e") == pytest.approx(2.0, abs=1e-6)


def test_production_at_max_growth_is_bounded(anaerobic_core):
    res = production_at_max_growth(anaerobic_core, "EX_etoh_e")
    assert 0 < res.growth < 0.874
    assert 0 <= res.target_flux_min <= res.target_flux_max
    assert res.guaranteed_yield <= res.max_theoretical_yield


def test_knockout_search_improves_yield_and_restores_model(anaerobic_core):
    search = greedy_knockout_search(anaerobic_core, "EX_etoh_e", max_knockouts=2)
    yields = [s.guaranteed_yield for s in search.steps]
    assert yields == sorted(yields) and yields[-1] > yields[0]
    assert all(s.growth >= 0.05 for s in search.steps)
    # the search must not leave knockouts behind on the caller's model
    assert all(anaerobic_core.reactions.get_by_id(r).bounds != (0, 0) for r in search.steps[-1].knocked_out)


def test_no_free_lunch():
    assert free_energy_checks(load_model("e_coli_core"))["ok"].all()


def test_glycolysis_genes_are_essential():
    df = gene_essentiality_check(load_model("e_coli_core")).set_index("gene")
    for gene in ("b1779", "b2926", "b2779"):  # gapA, pgk, eno
        assert df.at[gene, "verdict"] == "agree"

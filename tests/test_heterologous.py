"""Non-native pathways (vanillin, lycopene) added to E. coli."""

import pytest

from cellfactory import analyze, load_model
from cellfactory.heterologous import add_pathway


@pytest.fixture(scope="module")
def ijo():
    return load_model("iJO1366")


@pytest.mark.parametrize("product", ["vanillin", "lycopene"])
def test_pathway_is_added_and_mass_balanced(ijo, product):
    model = ijo.copy()
    steps = add_pathway(model, product)
    product_id = steps.pop("product")
    assert set(steps.values()) == {"added"}  # E. coli has none of these enzymes
    for rxn_id in steps:
        assert model.reactions.get_by_id(rxn_id).check_mass_balance() == {}
    assert analyze(model, "glucose", product_id).production.max_theoretical_yield > 0


def test_native_step_is_reused(ijo):
    model = ijo.copy()
    add_pathway(model, "lycopene")
    assert add_pathway(model, "lycopene")["HET_CRTE"] == "native (HET_CRTE)"


def test_lycopene_yield_matches_mep_stoichiometry(ijo):
    # Lycopene is 8 IPP units. The MEP pathway's aerobic ceiling is ~0.8 IPP per glucose,
    # so lycopene tops out near 0.1 mol/mol.
    model = ijo.copy()
    product = add_pathway(model, "lycopene")["product"]
    assert analyze(model, "glucose", product).production.max_theoretical_yield == pytest.approx(0.102, abs=0.005)


def test_unknown_metabolite_raises_without_pathway(ijo):
    with pytest.raises(KeyError):
        analyze(ijo.copy(), "glucose", "vanillin")

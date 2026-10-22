"""Nearest-native-precursor search and the product fate check.

The two pathway tests are the real check on this module. Asked to reach
chorismate and lysine from *E. coli* central metabolism, the search should
rediscover the textbook routes on its own: the seven-step shikimate pathway
(Herrmann & Weaver 1999, Annu. Rev. Plant Physiol. 50:473) and the
diaminopimelate route to L-lysine (Neidhardt, *EcoSal*, ch. 3.6.1.5). It has
no knowledge of either; it only has a reaction graph and a step count. If a
change to the graph, the hub filter or the cost function quietly breaks the
biology, these fail.
"""

import pytest

from cellfactory.models import load_model
from cellfactory.precursors import (
    CURRENCY,
    _base_id,
    _is_transport,
    build_graph,
    fate_to_text,
    nearest_precursors,
    product_fate,
    routes_to_text,
)


@pytest.fixture(scope="module")
def host():
    return load_model("e_coli_core")


@pytest.fixture(scope="module")
def donor():
    return load_model("iJO1366")


# --- the search ------------------------------------------------------------


def test_a_metabolite_the_host_already_makes_needs_no_enzymes(host, donor):
    (route,) = nearest_precursors(host, "succ_c", donor, max_heterologous_steps=2)
    assert route.precursor == "succ_c"
    assert route.heterologous_steps == 0
    assert route.steps == []
    assert route.precursor_max_flux > 0


def test_rediscovers_the_shikimate_pathway_to_chorismate(host, donor):
    """Seven enzymes from PEP to chorismate, in order, with the real aro genes."""
    routes = nearest_precursors(host, "chor_c", donor, max_heterologous_steps=8,
                                max_total_steps=10, top=3)
    best = routes[0]
    assert best.precursor == "pep_c"
    assert best.heterologous_steps == 7
    assert [s.reaction for s in best.steps] == [
        "DDPA", "DHQS", "DHQTi", "SHK3Dr", "SHKK", "PSCVT", "CHORS",
    ]
    # aroA (b0908) and aroC (b2329) are the last two steps; they must come
    # through as a cloning shortlist, not just as reaction ids.
    assert {"b0908", "b2329"} <= set(best.new_genes)
    assert best.precursor_max_flux > 0


def test_rediscovers_the_diaminopimelate_route_to_lysine(host, donor):
    routes = nearest_precursors(host, "lys__L_c", donor, max_heterologous_steps=11,
                                max_total_steps=13, top=3)
    # The textbook precursor is aspartate, which e_coli_core does not contain, so
    # the search correctly backs up one more step to the two halves of the
    # transamination that makes it: glutamate the amino donor, oxaloacetate the
    # carbon skeleton. Both should be offered.
    assert {r.precursor for r in routes} >= {"glu__L_c", "oaa_c"}
    best = routes[0]
    assert best.precursor == "glu__L_c"  # more of it available, and one carbon closer
    steps = [s.reaction for s in best.steps]
    assert steps[0] == "ASPTA"  # oxaloacetate -> aspartate
    assert steps[-1] == "DAPDC"  # the decarboxylation that makes lysine
    assert {"DHDPS", "DAPE"} <= set(steps)
    assert best.heterologous_steps == len(steps)  # none of it is native to e_coli_core


def test_routes_are_chemically_coherent(host, donor):
    """Every route must chain end to end, never reuse a reaction back to back,
    and never depend on an input the host cannot supply."""
    routes = nearest_precursors(host, "chor_c", donor, max_heterologous_steps=8,
                                max_total_steps=10, top=5)
    assert routes
    for route in routes:
        assert route.steps[0].substrate == route.precursor
        assert route.steps[-1].product == route.target
        for before, after in zip(route.steps, route.steps[1:]):
            assert before.product == after.substrate, "route does not chain"
            assert before.reaction != after.reaction, "same reaction used twice in a row"
        for step in route.steps:
            for other in step.also_needs:
                assert other in host.metabolites, f"{step.reaction} needs {other}, which the host lacks"


def test_cofactors_are_never_offered_as_precursors(host, donor):
    """Without the hub filter everything is two steps from everything through ATP."""
    routes = nearest_precursors(host, "chor_c", donor, max_heterologous_steps=8,
                                max_total_steps=10, top=10)
    for route in routes:
        assert _base_id(route.precursor) not in CURRENCY


def test_graph_edges_skip_hub_metabolites(host, donor):
    graph = build_graph(host, donor)
    assert "atp_c" not in graph.producers
    assert not any(e.substrate == "h2o_c" for edges in graph.producers.values() for e in edges)


def test_unreachable_target_returns_nothing_rather_than_a_bad_route(host, donor):
    assert nearest_precursors(host, "chor_c", donor, max_heterologous_steps=2) == []


def test_unknown_compound_suggests_alternatives(host):
    with pytest.raises(KeyError, match="not a metabolite"):
        nearest_precursors(host, "unobtainium", None)


def test_bad_limits_are_rejected(host):
    with pytest.raises(ValueError, match="max_heterologous_steps"):
        nearest_precursors(host, "succ_c", None, max_heterologous_steps=-1)


def test_text_summary_is_readable(host, donor):
    routes = nearest_precursors(host, "chor_c", donor, max_heterologous_steps=8, max_total_steps=10)
    text = routes_to_text(routes)
    assert "Chorismate" in text and "new enzyme" in text and "CHORS" in text
    assert "No native precursor" in routes_to_text([])


# --- product fate ----------------------------------------------------------


def test_succinate_is_eaten_by_succinate_dehydrogenase(host):
    fate = product_fate(host, "succ_c")
    by_id = {row["reaction"]: row for row in fate["consumed_by"]}
    assert by_id["SUCDi"]["active"], "SUCDi carries flux while the cell grows"
    assert by_id["SUCDi"]["safe_to_delete"]
    assert "SUCDi" in fate["protective_knockouts"]
    assert by_id["SUCDi"]["genes"], "a protective knockout is useless without genes"


def test_charge_decides_whether_the_product_can_leave_the_cell(host):
    acid = product_fate(host, "succ_c")
    assert acid["charge"] == -2
    assert "transporter" in acid["charge_note"]

    neutral = product_fate(host, "etoh_c")
    assert neutral["charge"] == 0
    assert "membrane" in neutral["charge_note"]


def test_aldehydes_are_flagged_as_chemically_fragile(host):
    fate = product_fate(host, "acald_c")
    assert "aldehyde" in [g["group"] for g in fate["functional_group_hints"]]
    assert not product_fate(host, "succ_c")["functional_group_hints"]


def test_reactive_pools_are_only_reported_when_the_host_makes_them(host):
    fate = product_fate(host, "succ_c")
    pools = {p["metabolite"] for p in fate["reactive_pools_the_host_makes"]}
    assert "coa_c" in pools, "E. coli core makes free CoA"
    # e_coli_core has no glutathione at all, so it must not be claimed.
    assert "gthrd_c" not in pools
    for pool in fate["reactive_pools_the_host_makes"]:
        assert pool["max_flux"] > 0


def test_transporters_and_placeholder_genes_stay_out_of_knockout_advice(host):
    """Deleting an exporter traps the product; s0001 is not a real gene."""
    fate = product_fate(host, "acald_c")
    by_id = {row["reaction"]: row for row in fate["consumed_by"]}
    assert by_id["ACALDt"]["transport"]
    assert not by_id["ACALDt"]["safe_to_delete"]
    assert "ACALDt" not in fate["protective_knockouts"]
    assert all("s0001" not in row["genes"] for row in fate["consumed_by"])


def test_a_product_the_host_cannot_make_has_nothing_eating_it(host, donor):
    fate = product_fate(host, "chor_c", donor)
    assert fate["in_host_model"] is False
    assert fate["consumed_by"] == []
    assert "no enzyme that touches" in fate_to_text(fate)


def test_is_transport_distinguishes_moving_from_reacting(host):
    assert _is_transport(host.reactions.get_by_id("ACALDt"))
    assert not _is_transport(host.reactions.get_by_id("ACALD"))

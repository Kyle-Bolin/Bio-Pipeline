"""Enzyme kinetics data layer (kcat/Km). Offline: BRENDA is a synthetic fixture file, SABIO-RK is always mocked."""

import io
import json
import urllib.error
from pathlib import Path

import pandas as pd
import pytest

from cellfactory import analyze, enzyme_kinetics, load_kcat_predictions, load_model
from cellfactory import kinetics
from cellfactory.cli import main
from cellfactory.kinetics import KineticValue, KineticsUnavailable, kinetics_caveats, read_brenda, sabiork_lookup

BRENDA = Path(__file__).parent / "data" / "kinetics" / "brenda_sample.txt"  # made-up numbers in BRENDA's layout


@pytest.fixture(scope="module")
def _core():
    return load_model("e_coli_core")


@pytest.fixture
def core(_core):
    return _core.copy()


@pytest.fixture(autouse=True)
def isolated(tmp_path, monkeypatch):
    monkeypatch.setenv("CELLFACTORY_CACHE", str(tmp_path))
    monkeypatch.delenv("CELLFACTORY_OFFLINE", raising=False)
    monkeypatch.delenv("CELLFACTORY_BRENDA", raising=False)


def no_network(ec):
    raise AssertionError("SABIO-RK must not be contacted here")


def row(result, rxn):
    t = result.table
    return t[t["reaction"] == rxn].iloc[0]


# ---------------------------------------------------------------- BRENDA parsing


def test_read_brenda_parses_values_ranges_mutants_and_continuations():
    values = read_brenda(BRENDA)
    pfk = [v for v in values if v.ec == "2.7.1.11" and v.parameter == "kcat"]
    assert sorted(v.value for v in pfk) == [2.0, 100.0, 130.0, 500.0]  # 120-140 -> midpoint; -999 skipped
    assert [v.value for v in pfk if v.mutant] == [2.0]
    assert {v.organism for v in pfk} == {"Escherichia coli", "Homo sapiens"}  # strain and accession dropped
    fba = next(v for v in values if v.ec == "4.1.2.13")
    assert fba.substrate == "D-fructose 1,6-bisphosphate"  # wrapped onto a continuation line
    assert not [v for v in values if v.ec == "5.3.1.1"]  # only a -999 placeholder
    # "#1,3#" lists two E. coli protein entries for one measurement: count it once.
    assert len([v for v in values if v.parameter == "km" and v.substrate.startswith("D-fructose")]) == 1


def test_read_brenda_filters_by_ec():
    assert {v.ec for v in read_brenda(BRENDA, ["5.3.1.9"])} == {"5.3.1.9"}


def test_unreadable_brenda_file_is_a_note_not_a_crash(core, tmp_path):
    result = enzyme_kinetics(core, ["PFK"], brenda=tmp_path / "missing.txt", online=False)
    assert result.sources == [] and "BRENDA not read" in result.notes[0]
    assert row(result, "PFK")["value_type"] == "missing"


# ---------------------------------------------------------------- choosing values


def test_host_wild_type_values_win(core):
    pfk = row(enzyme_kinetics(core, ["PFK"], brenda=BRENDA, online=False), "PFK")
    assert pfk["value_type"] == "measured" and pfk["organism_match"] == "host"
    assert pfk["kcat_per_s"] == 115.0 and (pfk["kcat_min"], pfk["kcat_max"], pfk["n_kcat"]) == (100.0, 130.0, 2)
    assert pfk["substrate_matched"] and pfk["source"] == "BRENDA"
    assert pfk["km_mM"] == 0.03 and pfk["km_substrate"] == "D-fructose 6-phosphate"


def test_organism_tiers(core):
    result = enzyme_kinetics(core, ["PGI", "FBA"], brenda=BRENDA, online=False)
    assert row(result, "PGI")["organism_match"] == "same genus"
    assert row(result, "FBA")["organism_match"] == "other organism"
    assert any("another organism" in c and "PGI" in c and "FBA" in c for c in kinetics_caveats(result))


def test_unknown_host_never_counts_as_host(core):
    core.id = "mystery_bug"
    result = enzyme_kinetics(core, ["PFK"], brenda=BRENDA, online=False, reference="iJO1366")
    assert result.organism is None and "host species" in result.notes[0]
    assert row(result, "PFK")["organism_match"] == "other organism"


def test_explicit_organism_overrides_model_id(core):
    pfk = row(enzyme_kinetics(core, ["PFK"], organism="Homo sapiens", brenda=BRENDA, online=False), "PFK")
    assert pfk["kcat_per_s"] == 500.0 and pfk["organism_match"] == "host"


def test_value_for_another_substrate_is_flagged(core):
    other = [KineticValue("2.7.1.11", "kcat", 7.0, "tagatose 6-phosphate", "Escherichia coli", "SABIO-RK")]
    pfk = row(enzyme_kinetics(core, ["PFK"], lookup=lambda ec: other if ec == "2.7.1.11" else []), "PFK")
    assert pfk["kcat_per_s"] == 7.0 and not pfk["substrate_matched"]
    assert "substrate this reaction doesn't use" in pfk["note"]


def test_only_mutant_data_counts_as_missing(core):
    mutant = [KineticValue("2.7.1.11", "kcat", 2.0, "D-fructose 6-phosphate", "Escherichia coli", "SABIO-RK", True)]
    assert row(enzyme_kinetics(core, ["PFK"], lookup=lambda ec: mutant), "PFK")["value_type"] == "missing"


# ---------------------------------------------------------------- predictions


def test_predictions_fill_gaps_only_and_are_flagged(core):
    preds = pd.DataFrame({"reaction": ["TPI", "PFK"], "kcat": [50.0, 9.0], "predictor": ["DLKcat", "DLKcat"]})
    result = enzyme_kinetics(core, ["PFK", "TPI"], brenda=BRENDA, online=False, predictions=preds)
    assert row(result, "PFK")["value_type"] == "measured" and row(result, "PFK")["kcat_per_s"] == 115.0
    tpi = row(result, "TPI")
    assert tpi["value_type"] == "predicted" and tpi["kcat_per_s"] == 50.0 and tpi["source"] == "DLKcat (predicted)"
    assert "DLKcat predictions" in result.sources
    assert any("PREDICTED" in c and "TPI" in c for c in kinetics_caveats(result))


def test_predictions_by_ec_from_file(core, tmp_path):
    path = tmp_path / "turnup.tsv"
    path.write_text("ec\tkcat\n5.3.1.1\t80\n5.3.1.9\t-1\n")
    preds = load_kcat_predictions(path, predictor="TurNuP")
    assert list(preds["ec"]) == ["5.3.1.1"]  # the negative prediction is dropped
    tpi = row(enzyme_kinetics(core, ["TPI"], online=False, predictions=path), "TPI")
    assert tpi["value_type"] == "predicted" and tpi["kcat_per_s"] == 80.0


def test_predictions_need_the_right_columns(tmp_path):
    with pytest.raises(ValueError, match="'kcat' column"):
        load_kcat_predictions(pd.DataFrame({"gene": ["b3916"], "value": [3]}))


# ---------------------------------------------------------------- EC numbers


def test_core_model_borrows_ecs_from_ijo1366(core):
    pfk = row(enzyme_kinetics(core, ["PFK"], online=False), "PFK")
    assert pfk["ec"] == "2.7.1.11" and "iJO1366" in pfk["ec_source"]


def test_no_borrowing_without_a_reference(core):
    pfk = row(enzyme_kinetics(core, ["PFK"], online=False, reference=None), "PFK")
    assert pfk["ec"] == "" and pfk["note"] == "no EC number annotated"


def test_too_many_or_partial_ecs_are_not_looked_up(core):
    core.reactions.PFK.annotation["ec-code"] = ["2.7.1.11", "2.7.1.90", "2.7.1.144", "2.7.1.56"]
    core.reactions.PGI.annotation["ec-code"] = "5.3.1.-"
    result = enzyme_kinetics(core, ["PFK", "PGI"], brenda=BRENDA, lookup=no_network, online=False)
    assert "too many to pick one enzyme" in row(result, "PFK")["note"]
    assert "partial EC" in row(result, "PGI")["note"]
    assert set(result.table["value_type"]) == {"missing"}


def test_exchanges_biomass_and_geneless_reactions(core):
    result = enzyme_kinetics(core, ["EX_glc__D_e", "Biomass_Ecoli_core", "ATPM", "PFK"], online=False)
    assert list(result.table["reaction"]) == ["ATPM", "PFK"]
    assert row(result, "ATPM")["value_type"] == "no enzyme"


# ---------------------------------------------------------------- SABIO-RK

SABIO_TSV = (
    "EntryID\tOrganism\tECNumber\tEnzymeType\tparameter.type\tparameter.associatedSpecies\tparameter.startValue\t"
    "parameter.unit\n"
    "1\tEscherichia coli\t2.7.1.11\twildtype\tkcat\tD-Fructose 6-phosphate\t6000\tmin^(-1)\n"
    "2\tEscherichia coli\t2.7.1.11\twildtype\tKm\tD-Fructose 6-phosphate\t40\tµM\n"
    "3\tEscherichia coli\t2.7.1.11\tmutant R162S\tkcat\tD-Fructose 6-phosphate\t3\ts^(-1)\n"
    "4\tEscherichia coli\t2.7.1.11\twildtype\tVmax\tD-Fructose 6-phosphate\t9\tumol/min/mg\n"
    "5\tEscherichia coli\t2.7.1.11\twildtype\tkcat\tATP\t\ts^(-1)\n"
)


class FakeResponse(io.BytesIO):
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def test_sabiork_parses_units_and_caches(monkeypatch, tmp_path):
    calls = []

    def fake_urlopen(req, timeout):
        calls.append(req)
        return FakeResponse(SABIO_TSV.encode())

    monkeypatch.setattr(kinetics.urllib.request, "urlopen", fake_urlopen)
    values = sabiork_lookup("2.7.1.11")
    assert [(v.parameter, v.value, v.mutant) for v in values] == [("kcat", 100.0, False), ("km", 0.04, False),
                                                                    ("kcat", 3.0, True)]
    assert b"ECNumber" in calls[0].data
    assert sabiork_lookup("2.7.1.11") == values and len(calls) == 1  # second call from the disk cache
    assert json.loads((tmp_path / "sabiork" / "2.7.1.11.json").read_text())[0]["source"] == "SABIO-RK"


def test_sabiork_404_means_no_data(monkeypatch):
    def not_found(req, timeout):
        raise urllib.error.HTTPError(req.full_url, 404, "Not Found", {}, io.BytesIO(b"No results found for query"))

    monkeypatch.setattr(kinetics.urllib.request, "urlopen", not_found)
    assert sabiork_lookup("5.3.1.1") == []


@pytest.mark.parametrize("error", [urllib.error.URLError("blocked"), TimeoutError("slow"),
                                   urllib.error.HTTPError("u", 503, "down", {}, None)])
def test_sabiork_failures_raise_kineticsunavailable(monkeypatch, error):
    def fail(req, timeout):
        raise error

    monkeypatch.setattr(kinetics.urllib.request, "urlopen", fail)
    with pytest.raises(KineticsUnavailable):
        sabiork_lookup("2.7.1.11")


def test_sabiork_offline_and_bad_ec(monkeypatch):
    with pytest.raises(KineticsUnavailable, match="complete EC"):
        sabiork_lookup("2.7.1.-")
    monkeypatch.setenv("CELLFACTORY_OFFLINE", "1")
    with pytest.raises(KineticsUnavailable, match="OFFLINE"):
        sabiork_lookup("2.7.1.11")


def test_unreachable_sabiork_stops_after_one_try(core):
    calls = []

    def down(ec):
        calls.append(ec)
        raise KineticsUnavailable("connection refused")

    result = enzyme_kinetics(core, ["PFK", "PGI", "FBA", "TPI"], brenda=BRENDA, lookup=down)
    assert len(calls) == 1  # not one timeout per EC number
    assert result.sources == ["BRENDA"] and "SABIO-RK could not be reached" in result.notes[0]
    assert row(result, "PFK")["value_type"] == "measured"  # BRENDA still used


def test_brenda_and_sabiork_combine(core):
    extra = [KineticValue("2.7.1.11", "kcat", 118.0, "D-fructose 6-phosphate", "Escherichia coli", "SABIO-RK")]
    pfk = row(enzyme_kinetics(core, ["PFK"], brenda=BRENDA, lookup=lambda ec: extra if ec == "2.7.1.11" else []),
              "PFK")
    assert pfk["n_kcat"] == 3 and pfk["kcat_per_s"] == 118.0 and pfk["source"] == "BRENDA + SABIO-RK"


def test_brenda_path_from_environment(core, monkeypatch):
    monkeypatch.setenv("CELLFACTORY_BRENDA", str(BRENDA))
    assert row(enzyme_kinetics(core, ["PFK"], online=False), "PFK")["value_type"] == "measured"


# ---------------------------------------------------------------- pipeline and CLI


def test_analyze_reports_route_kinetics(core):
    report = analyze(core, "glucose", "ethanol", aerobic=False, qa_gate=False, online_ids=False, kinetics=True,
                     brenda_file=str(BRENDA))
    table = report.checks["Enzyme kinetics (kcat/Km)"]
    assert {"PFK", "PGI", "FBA", "ALCD2x"} <= set(table["reaction"])
    assert not table["reaction"].str.startswith("EX_").any()
    assert any("no kcat at all" in c for c in report.caveats)
    assert any("test-tube values" in c for c in report.caveats)


def test_analyze_without_any_source_says_so(core):
    report = analyze(core, "glucose", "ethanol", aerobic=False, qa_gate=False, online_ids=False, kinetics=True)
    assert any("not looked up" in c and "every kcat is missing" in c for c in report.caveats)


def test_cli_kinetics(capsys):
    assert main(["run", "-p", "ethanol", "--anaerobic", "--no-qa", "--offline", "--kinetics",
                 "--brenda", str(BRENDA)]) == 0
    out = capsys.readouterr().out
    assert "Enzyme kinetics (kcat/Km)" in out and "PFK" in out

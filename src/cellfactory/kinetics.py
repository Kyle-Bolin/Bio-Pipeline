"""Enzyme speed: look up each reaction's turnover number (kcat) and Michaelis constant (Km).

FBA knows which reactions a cell *could* run, not how much enzyme it would
take to run them. Two measured numbers answer that:

* **kcat** (per second): how many substrate molecules one enzyme molecule turns
  over each second when it is saturated with substrate. Flux divided by kcat is
  the least amount of enzyme a reaction needs.
* **Km** (mM): the substrate concentration at which the enzyme runs at half its
  kcat. Below Km the enzyme mostly waits for substrate and runs well under kcat.

``enzyme_kinetics`` returns one row per reaction with both, saying where each
number came from. This is a data layer only; nothing here constrains the model
(that is the enzyme-constrained FBA that will consume this table).

Where the numbers come from, in order of trust:

1. **Measured, from BRENDA**: read from BRENDA's flat-file download
   (``brenda_YYYY_N.txt``). BRENDA's licence requires each user to register and
   download it themselves, so it is never bundled: pass its path, or set
   ``CELLFACTORY_BRENDA``.
2. **Measured, from SABIO-RK**: its free web service, queried per EC number,
   cached on disk. Best effort: offline, blocked or slow, it becomes a note,
   never a crash or a hang (one failure stops further queries for that call;
   ``CELLFACTORY_OFFLINE=1`` never tries).
3. **Predicted**: a table of kcats from a predictor such as DLKcat or TurNuP,
   run separately (both need PyTorch, protein sequences and substrate
   structures). Used only where no measurement exists, and always tagged
   ``predicted``: a predicted kcat can be ten times off.

Reactions are matched to enzyme data by **EC number** (the Enzyme Commission's
code for the chemistry an enzyme does, such as 2.7.1.11 for
phosphofructokinase), read from the model's ``ec-code`` annotations. A model
without any (the bundled ``e_coli_core``) borrows them from iJO1366 by BiGG
reaction id, which BiGG keeps the same across its models.

Among measurements, wild-type enzymes from the host species beat the same
genus, which beat any other organism; values for one of the reaction's own
substrates beat values for other substrates. The kcat reported is the median
of the best group, with its range and count so the spread stays visible.
"""

from __future__ import annotations

import csv
import http.client
import io
import json
import logging
import os
import re
import statistics
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable

import cobra
import pandas as pd

log = logging.getLogger(__name__)

SABIO_URL = "https://sabiork.h-its.org/sabioRestWebServices/kineticlawsExportTsv"
SABIO_FIELDS = ["EntryID", "Organism", "ECNumber", "EnzymeType", "parameter.type", "parameter.associatedSpecies",
                "parameter.startValue", "parameter.unit"]
TIMEOUT = 10.0  # seconds per SABIO-RK request

# Species each bundled or common BiGG model describes, so callers need not say.
MODEL_ORGANISMS = {
    "e_coli_core": "Escherichia coli", "iJO1366": "Escherichia coli", "iML1515": "Escherichia coli",
    "iAF1260": "Escherichia coli", "iYS1720": "Salmonella enterica", "iMM904": "Saccharomyces cerevisiae",
    "iND750": "Saccharomyces cerevisiae",
}
EC_REFERENCE_MODEL = "iJO1366"  # where ECs are borrowed from for a model with none
MAX_ECS = 3  # more EC numbers than this on one reaction means a lumped or pseudo-reaction: no single enzyme to look up
TIERS = ("host", "same genus", "other organism")

_EC = re.compile(r"^\d+\.\d+\.\d+\.\d+$")

COLUMNS = ["reaction", "ec", "ec_source", "value_type", "kcat_per_s", "kcat_min", "kcat_max", "n_kcat",
           "kcat_organism", "organism_match", "substrate_matched", "source", "km_mM", "km_substrate", "n_km", "note"]


class KineticsUnavailable(RuntimeError):
    """An enzyme-kinetics source could not be read: offline, blocked, timed out, or an unexpected answer."""


@dataclass
class KineticValue:
    """One measured value as a database reports it."""

    ec: str
    parameter: str  # "kcat" (1/s) or "km" (mM)
    value: float
    substrate: str
    organism: str
    source: str  # "BRENDA" or "SABIO-RK"
    mutant: bool = False


@dataclass
class KineticsResult:
    table: pd.DataFrame  # one row per reaction, columns COLUMNS
    organism: str | None  # the host species values were matched against
    sources: list[str] = field(default_factory=list)  # what was consulted: "BRENDA", "SABIO-RK", "<predictor> predictions"
    notes: list[str] = field(default_factory=list)  # sources that failed, and why


# ---------------------------------------------------------------- EC numbers


def _as_list(value) -> list[str]:
    if value is None:
        return []
    return [str(v).strip() for v in (value if isinstance(value, (list, tuple, set)) else [value]) if str(v).strip()]


def _needs_enzyme(rxn: cobra.Reaction) -> bool | None:
    """None: not a reaction at all (exchange, biomass); False: no gene, so no enzyme to look up; True otherwise."""
    if rxn.boundary or rxn.id.lower().startswith("biomass") or rxn.objective_coefficient:
        return None
    return bool(rxn.genes)


def _reaction_ecs(rxn: cobra.Reaction, reference: cobra.Model | None) -> tuple[list[str], str, str]:
    """(complete EC numbers, where they came from, problem or "")."""
    ecs, origin = _as_list(rxn.annotation.get("ec-code")), "model annotation"
    if not ecs and reference is not None and rxn.id in reference.reactions:
        ecs, origin = _as_list(reference.reactions.get_by_id(rxn.id).annotation.get("ec-code")), f"{reference.id} (same reaction id)"
    if not ecs:
        return [], "", "no EC number annotated"
    complete = sorted({e for e in ecs if _EC.match(e)})
    if len(complete) > MAX_ECS:
        return [], origin, (f"{len(complete)} EC numbers annotated: too many to pick one enzyme (a multi-enzyme "
                            "complex, a lumped reaction or a pseudo-reaction)")
    if not complete:
        return [], origin, f"only partial EC numbers ({', '.join(ecs)}): the exact enzyme is unknown"
    return complete, origin, ""


_REFERENCES: dict[str, cobra.Model] = {}  # loaded once per process; only ever read


def _reference(model_id: str) -> cobra.Model:
    if model_id not in _REFERENCES:
        from .models import load_model

        _REFERENCES[model_id] = load_model(model_id)
    return _REFERENCES[model_id]


def _default_reference(model: cobra.Model) -> cobra.Model | None:
    if any("ec-code" in r.annotation for r in model.reactions):
        return None
    try:
        return _reference(EC_REFERENCE_MODEL)
    except Exception as err:  # noqa: BLE001 - a missing reference only means no borrowed ECs
        log.warning("Could not load %s to borrow EC numbers: %s", EC_REFERENCE_MODEL, err)
        return None


# ---------------------------------------------------------------- BRENDA flat file

_ENTRY = re.compile(r"^#([\d,\s]+)#\s*(.*)$")
_VALUE = re.compile(r"^(-?\d+(?:\.\d+)?(?:[eE]-?\d+)?)(?:\s*-\s*(\d+(?:\.\d+)?(?:[eE]-?\d+)?))?\s*\{([^}]*)\}\s*(.*)$")
_REFS = re.compile(r"<[\d,\s]+>")


def _ids(text: str) -> list[str]:
    return [i.strip() for i in text.split(",") if i.strip()]


def _brenda_organism(text: str) -> str:
    # "Escherichia coli P0A796 UniProt <1,2>" -> "Escherichia coli"; strain names and accessions are dropped.
    words = _REFS.sub("", text).split()
    return " ".join(words[:2])


def _mutants(comment: str) -> set[str]:
    """Protein numbers whose commentary describes a mutant enzyme."""
    out = set()
    for part in comment.split(";"):
        m = _ENTRY.match(part.strip())
        if m and "mutant" in m.group(2).lower() and "wild" not in m.group(2).lower():
            out.update(_ids(m.group(1)))
    return out


def read_brenda(path: str | os.PathLike, ecs: Iterable[str] | None = None) -> list[KineticValue]:
    """Turnover numbers and Km values from BRENDA's flat-file download, for ``ecs`` (all if None).

    Reads only the ``PROTEIN`` (``PR``), ``TURNOVER_NUMBER`` (``TN``, 1/s) and
    ``KM_VALUE`` (``KM``, mM) lines. A range ("3.5-4") becomes its midpoint;
    BRENDA's ``-999`` (a value given only in the commentary) is skipped.
    """
    wanted = None if ecs is None else set(ecs)
    out: list[KineticValue] = []
    ec, keep, organisms, entry = None, False, {}, None

    def flush():
        if not entry or not keep:
            return
        tag, text = entry
        m = _ENTRY.match(text)
        if not m:
            return
        proteins, body = _ids(m.group(1)), m.group(2).strip()
        if tag == "PR":
            for p in proteins:
                organisms[p] = _brenda_organism(body)
            return
        v = _VALUE.match(body)
        if not v or v.group(1).startswith("-"):
            return
        lo, hi = float(v.group(1)), float(v.group(2)) if v.group(2) else None
        value = (lo + hi) / 2 if hi is not None else lo
        if value <= 0:
            return
        comment = v.group(4)
        mutant_ids = _mutants(comment[comment.find("(") + 1:comment.rfind(")")] if "(" in comment else "")
        seen = set()  # one line listing two proteins of the same species ("#1,3#") is one measurement, not two
        for p in proteins:
            key = (organisms.get(p, ""), p in mutant_ids)
            if key not in seen:
                seen.add(key)
                out.append(KineticValue(ec, "kcat" if tag == "TN" else "km", value, v.group(3).strip(),
                                        key[0], "BRENDA", mutant=key[1]))

    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            for raw in fh:
                line = raw.rstrip("\n")
                if line.startswith("\t") and entry is not None:
                    entry = (entry[0], entry[1] + " " + line.strip())
                    continue
                flush()
                entry = None
                if line.startswith("ID\t"):
                    ec = line[3:].split()[0] if line[3:].split() else None
                    keep, organisms = bool(ec) and (wanted is None or ec in wanted), {}
                elif line.startswith("///"):
                    ec, keep = None, False
                elif keep:
                    tag, _, rest = line.partition("\t")
                    if tag in ("PR", "TN", "KM"):
                        entry = (tag, rest)
            flush()
    except OSError as err:
        raise KineticsUnavailable(f"could not read the BRENDA file {path} ({err.strerror or err})") from err
    return out


# ---------------------------------------------------------------- SABIO-RK

_KCAT_UNITS = {"s^(-1)": 1.0, "1/s": 1.0, "s^-1": 1.0, "min^(-1)": 1 / 60, "1/min": 1 / 60, "h^(-1)": 1 / 3600}
_KM_UNITS = {"M": 1e3, "mM": 1.0, "µM": 1e-3, "uM": 1e-3, "μM": 1e-3, "nM": 1e-6}


def _cache_dir() -> Path:
    root = os.environ.get("CELLFACTORY_CACHE") or Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache"),
                                                       "cellfactory")
    return Path(root) / "sabiork"


def _parse_sabio_tsv(text: str, ec: str) -> list[KineticValue]:
    out = []
    for row in csv.DictReader(io.StringIO(text), delimiter="\t"):
        kind = (row.get("parameter.type") or "").strip().lower()
        unit = (row.get("parameter.unit") or "").strip()
        try:
            value = float(row.get("parameter.startValue") or "nan")
        except ValueError:
            continue
        if kind == "kcat" and unit in _KCAT_UNITS:
            param, value = "kcat", value * _KCAT_UNITS[unit]
        elif kind == "km" and unit in _KM_UNITS:
            param, value = "km", value * _KM_UNITS[unit]
        else:
            continue
        if not value > 0:  # also drops NaN
            continue
        out.append(KineticValue(ec, param, value, (row.get("parameter.associatedSpecies") or "").strip(),
                                (row.get("Organism") or "").strip(), "SABIO-RK",
                                mutant="mutant" in (row.get("EnzymeType") or "").lower()))
    return out


def sabiork_lookup(ec: str, timeout: float = TIMEOUT, use_cache: bool = True) -> list[KineticValue]:
    """kcat and Km values SABIO-RK holds for EC number ``ec``, from every organism.

    Raises KineticsUnavailable on any network or format problem. Answers,
    including "nothing known", are cached on disk.
    """
    if not _EC.match(ec):
        raise KineticsUnavailable(f"{ec!r} is not a complete EC number")
    if os.environ.get("CELLFACTORY_OFFLINE"):
        raise KineticsUnavailable("CELLFACTORY_OFFLINE is set")
    cache = _cache_dir() / f"{ec}.json"
    if use_cache and cache.exists():
        try:
            return [KineticValue(**d) for d in json.loads(cache.read_text())]
        except (OSError, ValueError, TypeError):
            pass  # unreadable or wrong-shaped cache entry: ask again
    query = [("fields[]", f) for f in SABIO_FIELDS] + [("q", f'ECNumber:"{ec}"')]
    # As SABIO-RK's documented example does: a POST with the query in the URL, not the body.
    req = urllib.request.Request(SABIO_URL + "?" + urllib.parse.urlencode(query), data=b"", method="POST",
                                 headers={"User-Agent": "cellfactory"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            text = resp.read().decode("utf-8", errors="replace")
    except urllib.error.HTTPError as err:
        # SABIO-RK answers 404 "No results found" when nothing matches. Any other 404 (a wrong endpoint, say)
        # is a failure, not "no data", so it is never cached.
        try:
            body = err.read().decode("utf-8", errors="replace") if err.fp else ""
        except OSError:
            body = ""
        if err.code != 404 or "no results" not in body.lower():
            raise KineticsUnavailable(f"SABIO-RK query failed ({err!r})") from err
        text = ""
    except (urllib.error.URLError, http.client.HTTPException, OSError, TimeoutError, ValueError) as err:
        raise KineticsUnavailable(f"SABIO-RK query failed ({err!r})") from err
    if text and "parameter.type" not in text.split("\n", 1)[0]:
        if "no results" in text.lower():
            text = ""
        else:
            raise KineticsUnavailable(f"SABIO-RK sent an unexpected answer ({text[:80]!r})")
    values = _parse_sabio_tsv(text, ec) if text else []
    if use_cache:
        try:
            cache.parent.mkdir(parents=True, exist_ok=True)
            cache.write_text(json.dumps([v.__dict__ for v in values]))
        except OSError:
            pass
    return values


# ---------------------------------------------------------------- predictions


def load_kcat_predictions(source: str | os.PathLike | pd.DataFrame, predictor: str | None = None) -> pd.DataFrame:
    """A predictor's kcats (DLKcat, TurNuP, ...) as a table with ``reaction`` or ``ec``, ``kcat`` and ``predictor``.

    ``source`` is a CSV/TSV file or a DataFrame. Rows keyed by ``reaction``
    match that model reaction id; rows keyed by ``ec`` match any reaction with
    that EC number. ``predictor`` names the tool when the table has no
    ``predictor`` column.
    """
    if isinstance(source, pd.DataFrame):
        df = source.copy()
    else:
        path = Path(source)
        try:
            df = pd.read_csv(path, sep="\t" if path.suffix.lower() in (".tsv", ".tab") else ",")
        except (OSError, ValueError) as err:
            raise ValueError(f"could not read kcat predictions from {path}: {err}") from err
    df.columns = [str(c).strip().lower() for c in df.columns]
    if "kcat" not in df.columns or not ({"reaction", "ec"} & set(df.columns)):
        raise ValueError("kcat predictions need a 'kcat' column (1/s) and a 'reaction' or 'ec' column; "
                         f"got {', '.join(df.columns)}")
    df["kcat"] = pd.to_numeric(df["kcat"], errors="coerce")
    bad = df["kcat"].isna() | (df["kcat"] <= 0)
    if bad.any():
        log.warning("Ignoring %d kcat predictions that are missing, zero or negative", int(bad.sum()))
        df = df[~bad]
    if "predictor" not in df.columns:
        df["predictor"] = predictor or "predictor"
    for col in ("reaction", "ec"):
        if col not in df.columns:
            df[col] = ""
    return df[["reaction", "ec", "kcat", "predictor"]].fillna("").reset_index(drop=True)


# ---------------------------------------------------------------- choosing values


def _match(organism: str, host: str | None) -> str:
    if not host or not organism:
        return "other organism"
    a, b = organism.lower().split(), host.lower().split()
    if a[:2] == b[:2]:
        return "host"
    return "same genus" if a[:1] == b[:1] else "other organism"


def _name(text: str) -> str:
    return re.sub(r"[\s_-]+", " ", text.lower()).strip()


def _best(values: list[KineticValue], substrates: set[str], host: str | None):
    """(values in the most trusted group, its organism tier, whether they name one of the reaction's substrates)."""
    values = [v for v in values if not v.mutant]
    if not values:
        return [], "", False
    matched = [v for v in values if _name(v.substrate) in substrates]
    pool = matched or values
    for tier in TIERS:
        group = [v for v in pool if _match(v.organism, host) == tier]
        if group:
            return group, tier, bool(matched)
    return [], "", False


def _row(rxn: cobra.Reaction, ecs: list[str], ec_source: str, values: list[KineticValue], host: str | None,
         predicted: tuple[float, str] | None, reverse: bool = False) -> dict:
    # A kcat belongs to one direction: match the substrates of the direction the flux actually runs.
    consumed = rxn.products if reverse else rxn.reactants
    substrates = {_name(m.name or "") for m in consumed} | {_name(m.id.rsplit("_", 1)[0]) for m in consumed}
    row = dict.fromkeys(COLUMNS, None)
    row.update(reaction=rxn.id, ec=", ".join(ecs), ec_source=ec_source, n_kcat=0, n_km=0, note="")
    kcats, tier, matched = _best([v for v in values if v.parameter == "kcat"], substrates, host)
    if kcats:
        nums = [v.value for v in kcats]
        row.update(value_type="measured", kcat_per_s=statistics.median(nums), kcat_min=min(nums), kcat_max=max(nums),
                   n_kcat=len(nums), kcat_organism=", ".join(sorted({v.organism for v in kcats})),
                   organism_match=tier, substrate_matched=matched,
                   source=" + ".join(sorted({v.source for v in kcats})))
        if not matched:
            row["note"] = "kcat measured with a substrate this reaction doesn't use"
    elif predicted is not None:
        row.update(value_type="predicted", kcat_per_s=predicted[0], n_kcat=1, source=f"{predicted[1]} (predicted)")
    else:
        row.update(value_type="missing", note="no kcat measured" + (" or predicted" if ecs else ""))
    kms, _, km_matched = _best([v for v in values if v.parameter == "km"], substrates, host)
    if kms:
        by_sub: dict[str, list[float]] = {}
        for v in kms:
            by_sub.setdefault(v.substrate, []).append(v.value)
        sub, nums = max(by_sub.items(), key=lambda kv: (len(kv[1]), kv[0]))
        row.update(km_mM=statistics.median(nums), km_substrate=sub, n_km=len(nums))
    return row


# ---------------------------------------------------------------- main entry


def enzyme_kinetics(
    model: cobra.Model,
    reactions: Iterable[str] | pd.Series | dict[str, float] | None = None,
    organism: str | None = None,
    brenda: str | os.PathLike | None = None,
    online: bool = True,
    predictions: str | os.PathLike | pd.DataFrame | None = None,
    reference: cobra.Model | str | None = "auto",
    lookup: Callable[[str], list[KineticValue]] | None = None,
) -> KineticsResult:
    """kcat and Km for each reaction in ``reactions`` (default: all), with where every number came from.

    ``reactions`` may be reaction ids or fluxes (a Series or dict, such as
    ``mdf.product_route_fluxes``); with fluxes, a reaction running backwards is
    matched against its products, since a kcat is for one direction.

    ``organism`` is the host species ("Escherichia coli"); it defaults from
    the model id for well-known BiGG models, and without one no value can
    count as measured in the host. ``brenda`` is BRENDA's flat-file download
    (default ``$CELLFACTORY_BRENDA``). ``online=False`` never contacts SABIO-RK.
    ``predictions`` is a predictor's output (see ``load_kcat_predictions``),
    used only where nothing was measured. ``reference`` is a model to borrow
    EC numbers from by reaction id; ``"auto"`` uses iJO1366 when the model has
    no EC annotations at all, None borrows nothing. ``lookup`` replaces
    ``sabiork_lookup`` (for tests).

    Exchange and biomass reactions are left out. Reactions with no gene get
    ``value_type`` "no enzyme": spontaneous reactions, diffusion and
    pseudo-reactions need no protein. The model is not changed.
    """
    host = organism or MODEL_ORGANISMS.get(model.id)
    reverse: dict[str, bool] = {}
    if reactions is None:
        ids = [r.id for r in model.reactions]
    elif isinstance(reactions, (pd.Series, dict)):  # fluxes: the sign says which way each reaction runs
        reverse = {rid: float(flux) < 0 for rid, flux in reactions.items()}
        ids = list(reverse)
    else:
        ids = list(reactions)
    if reference == "auto":
        reference = _default_reference(model)
    elif isinstance(reference, str):
        reference = _reference(reference)
    result = KineticsResult(table=pd.DataFrame(columns=COLUMNS), organism=host)
    if not host:
        result.notes.append(f"the host species of {model.id} is unknown (pass organism=), so no value counts as "
                            "measured in the host")

    plan = []  # (reaction, ecs, ec_source, problem)
    for rid in ids:
        rxn = model.reactions.get_by_id(rid)
        needs = _needs_enzyme(rxn)
        if needs is None:
            continue
        if needs is False:
            plan.append((rxn, [], "", None))
            continue
        plan.append((rxn, *_reaction_ecs(rxn, reference)))
    wanted = sorted({ec for _, ecs, _, _ in plan for ec in ecs})

    values: dict[str, list[KineticValue]] = {ec: [] for ec in wanted}
    brenda = brenda or os.environ.get("CELLFACTORY_BRENDA")
    if brenda:
        try:
            for v in read_brenda(brenda, wanted):
                values[v.ec].append(v)
            result.sources.append("BRENDA")
        except KineticsUnavailable as err:
            result.notes.append(f"BRENDA not read: {err}")
    if online and wanted:
        fetch = lookup or sabiork_lookup
        done = 0
        try:
            for ec in wanted:
                values[ec] += fetch(ec)
                done += 1
            result.sources.append("SABIO-RK")
        except KineticsUnavailable as err:
            log.warning("SABIO-RK unavailable: %s", err)
            if done:
                result.sources.append("SABIO-RK (partial)")
                result.notes.append(f"SABIO-RK stopped answering after {done} of {len(wanted)} EC numbers, so the "
                                    f"rest were not looked up there ({err})")
            else:
                result.notes.append(f"SABIO-RK could not be reached, so it was not used ({err})")

    by_rxn: dict[str, tuple[float, str]] = {}
    by_ec: dict[str, tuple[float, str]] = {}
    if predictions is not None:
        pred = load_kcat_predictions(predictions)
        for r in pred.itertuples():
            if r.reaction:
                by_rxn[str(r.reaction)] = (float(r.kcat), str(r.predictor))
            elif r.ec:
                by_ec[str(r.ec)] = (float(r.kcat), str(r.predictor))
        result.sources += sorted({f"{p} predictions" for p in pred["predictor"]})

    rows = []
    for rxn, ecs, ec_source, problem in plan:
        if problem is None:
            row = dict.fromkeys(COLUMNS, None)
            row.update(reaction=rxn.id, ec="", ec_source="", value_type="no enzyme", n_kcat=0, n_km=0,
                       note="no gene: spontaneous, diffusion or a pseudo-reaction")
            rows.append(row)
            continue
        predicted = by_rxn.get(rxn.id) or next((by_ec[e] for e in ecs if e in by_ec), None)
        row = _row(rxn, ecs, ec_source, [v for ec in ecs for v in values.get(ec, [])], host, predicted,
                   reverse=reverse.get(rxn.id, False))
        if problem:
            row["note"] = problem if row["value_type"] == "missing" else f"{row['note']}; {problem}".strip("; ")
        rows.append(row)
    if rows:
        result.table = pd.DataFrame(rows, columns=COLUMNS)
    return result


def kinetics_caveats(result: KineticsResult, label: str = "these reactions") -> list[str]:
    """Plain-words reasons to trust the kcats in ``result`` less, for ``Report.caveats``."""
    t = result.table
    enz = t[t["value_type"] != "no enzyme"]
    n = len(enz)
    out = [f"Enzyme kinetics: {note}." for note in result.notes]
    if not n:
        return out

    def some(df: pd.DataFrame) -> str:
        ids = list(df["reaction"])
        return ", ".join(ids[:8]) + (f" and {len(ids) - 8} more" if len(ids) > 8 else "")

    if not result.sources:
        out.append(f"Enzyme kinetics not looked up for {label}: no BRENDA file, SABIO-RK off or unreachable, and no "
                   "kcat predictions, so every kcat is missing.")
        return out
    pred = enz[enz["value_type"] == "predicted"]
    if len(pred):
        out.append(f"{len(pred)} of {n} enzyme steps in {label} use a PREDICTED kcat, not a measured one "
                   f"({some(pred)}). A predicted kcat can be ten times off, and the error passes straight into any "
                   "yield limited by enzyme capacity.")
    missing = enz[enz["value_type"] == "missing"]
    if len(missing):
        out.append(f"{len(missing)} of {n} enzyme steps in {label} have no kcat at all ({some(missing)}).")
    measured = enz[enz["value_type"] == "measured"]
    foreign = measured[measured["organism_match"] != "host"]
    if len(foreign):
        out.append(f"{len(foreign)} of {n} enzyme steps in {label} use kcats measured in another organism, not "
                   f"{result.organism or 'the host'} ({some(foreign)}); the same reaction's enzyme can differ a lot "
                   "between species.")
    borrowed = enz[enz["ec_source"].fillna("").str.contains("same reaction id")]
    if len(borrowed):
        out.append(f"EC numbers for {len(borrowed)} reactions were borrowed from {EC_REFERENCE_MODEL} by reaction id, "
                   "since the model has none of its own.")
    if len(measured):
        out.append("Measured kcats are test-tube values, at the lab's temperature and pH with substrate in excess; "
                   "inside a cell enzymes usually run below kcat.")
    return out

"""Identifier resolution: map a compound named in any database onto the model's own metabolite.

The same molecule has a different id in every database: D-glucose is
``glc__D`` in BiGG, ``C00031`` in KEGG, ``CHEBI:4167`` in ChEBI, ``GLC`` in
MetaCyc and ``MNXM41`` in MetaNetX, which cross-references all of them. A
substrate or product typed with the wrong namespace's id is the most common
silent failure when stitching a pathway from one source into a model from
another, so ``resolve_metabolite`` tries, in order:

1. the model's own ids and names, plus the ``ALIASES`` of everyday names;
2. the model's own database annotations (``met.annotation``), offline, for
   queries like ``kegg:C00031``, ``CHEBI:4167``, ``MNXM41`` or ``HMDB0000122``;
3. MetaNetX online, only for a database id the annotations don't cover: it
   translates the id to MetaNetX and BiGG ids, which are matched against the model.

If nothing matches, or a database id points at more than one model
metabolite, it says so with near-miss suggestions instead of guessing.

The MetaNetX query is best effort: it needs network access (set
``CELLFACTORY_OFFLINE=1`` to never try), is cached on disk, and any failure
becomes a note in the result, never a crash or a hang.
"""

from __future__ import annotations

import difflib
import http.client
import json
import logging
import os
import re
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import cobra
import pandas as pd

log = logging.getLogger(__name__)

METANETX_SPARQL = "https://rdf.metanetx.org/sparql"
TIMEOUT = 10.0  # seconds; MetaNetX is a nice-to-have, not worth a long stall

# Everyday names -> BiGG metabolite ids, so "glucose" works as well as "glc__D".
ALIASES = {
    "glucose": "glc__D", "fructose": "fru", "xylose": "xyl__D", "arabinose": "arab__L",
    "galactose": "gal", "sucrose": "sucr", "maltose": "malt", "lactose": "lcts", "glycerol": "glyc",
    "ethanol": "etoh", "acetate": "ac", "lactate": "lac__D", "succinate": "succ", "formate": "for",
    "pyruvate": "pyr", "citrate": "cit", "oxygen": "o2", "co2": "co2",
}

# Prefixes people type -> the annotation key BiGG/cobrapy models use (identifiers.org names).
NAMESPACES = {
    "kegg": "kegg.compound", "kegg.compound": "kegg.compound", "keggc": "kegg.compound",
    "chebi": "chebi",
    "metacyc": "biocyc", "metacyc.compound": "biocyc", "metacycm": "biocyc", "biocyc": "biocyc",
    "seed": "seed.compound", "seed.compound": "seed.compound", "seedm": "seed.compound",
    "hmdb": "hmdb",
    "mnx": "metanetx.chemical", "metanetx": "metanetx.chemical", "metanetx.chemical": "metanetx.chemical",
    "bigg": "bigg.metabolite", "bigg.metabolite": "bigg.metabolite", "biggm": "bigg.metabolite",
    "inchikey": "inchi_key", "inchi_key": "inchi_key",
}

# Bare ids whose shape gives away their database.
_BARE = [
    (re.compile(r"^C\d{5}$"), "kegg.compound"),
    (re.compile(r"^MNXM\d+$", re.I), "metanetx.chemical"),
    (re.compile(r"^CHEBI:\d+$", re.I), "chebi"),
    (re.compile(r"^HMDB\d+$", re.I), "hmdb"),
    (re.compile(r"^cpd\d{5}$"), "seed.compound"),
    (re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$"), "inchi_key"),
]


_SAFE_ID = re.compile(r"^[A-Za-z0-9_.:+-]+$")  # what a database id can contain; anything else never reaches SPARQL


class MetaNetXUnavailable(RuntimeError):
    """MetaNetX could not be queried: offline, blocked, timed out, or an unexpected answer."""


def parse_identifier(query: str) -> tuple[str, str] | None:
    """``"kegg:C00031"`` -> ``("kegg.compound", "C00031")``; None if ``query`` doesn't look like a database id.

    The returned id is normalised (see ``_norm``) so it compares equal to the
    model's annotation however it was typed.
    """
    q = query.strip()
    for pattern, ns in _BARE:
        if pattern.match(q):
            return ns, _norm(ns, q)
    if ":" in q:
        prefix, local = q.split(":", 1)
        ns = NAMESPACES.get(prefix.strip().lower())
        if ns and local.strip():
            return ns, _norm(ns, local.strip())
    return None


def _norm(ns: str, value: str) -> str:
    """Spellings of one id -> one form: ``CHEBI:4167``/``4167`` -> ``4167``, ``HMDB00122``/``HMDB0000122`` -> ``122``."""
    v = value.strip()
    if ns == "chebi":
        return v.upper().removeprefix("CHEBI:")
    if ns == "hmdb":
        digits = v.upper().removeprefix("HMDB")
        return str(int(digits)) if digits.isdigit() else digits
    if ns == "biocyc":
        return v.upper().removeprefix("META:")
    if ns in ("metanetx.chemical", "inchi_key"):
        return v.upper()
    return v


def _base(met: cobra.Metabolite) -> str:
    """Metabolite id without its compartment: ``glc__D_e`` -> ``glc__D``."""
    return met.id.rsplit("_", 1)[0] if met.compartment and met.id.endswith(f"_{met.compartment}") else met.id


def annotation_index(model: cobra.Model) -> dict[tuple[str, str], set[str]]:
    """(namespace, normalised id) -> base metabolite ids carrying that annotation in ``model``."""
    index: dict[tuple[str, str], set[str]] = defaultdict(set)
    for met in model.metabolites:
        base = _base(met)
        index[("bigg.metabolite", base)].add(base)
        for ns, values in met.annotation.items():
            for value in values if isinstance(values, list) else [values]:
                index[(ns, _norm(ns, str(value)))].add(base)
    return index


@dataclass
class Resolution:
    query: str
    metabolites: list[cobra.Metabolite] = field(default_factory=list)  # every compartment's copy
    via: str = ""  # "exchange id", "model id", "alias", "name", "annotation <ns>", "MetaNetX <id>", or "" if unresolved
    candidates: list[str] = field(default_factory=list)  # base ids, when an id matched several
    suggestions: list[str] = field(default_factory=list)  # near misses, when nothing matched
    note: str = ""

    @property
    def ok(self) -> bool:
        return bool(self.metabolites)

    @property
    def model_id(self) -> str:
        return _base(self.metabolites[0]) if self.metabolites else ""

    @property
    def caveat(self) -> str | None:
        """A plain-words reason to double-check this match, or None."""
        if self.via.startswith("MetaNetX"):
            return (f"{self.query!r} was matched to {self.model_id} through a MetaNetX cross-reference, not the "
                    "model's own annotations. MetaNetX merges protonation states and some stereoisomers into one "
                    "entry; check this is the molecule you meant.")
        return None

    def error(self, model_id: str) -> str:
        if self.candidates:
            return (f"{self.query!r} matches several metabolites in {model_id}: {', '.join(self.candidates)}. "
                    f"{self.note} Pass the one you mean by its model id.").replace("  ", " ")
        hint = ("Close matches: " + ", ".join(self.suggestions)) if self.suggestions else (
            "Nothing similar is in the model. The organism may lack the pathway; add the reactions with "
            "model.add_reactions() first.")
        return " ".join(p for p in (f"{self.query!r} is not a metabolite in {model_id}.", self.note, hint) if p)


def _mets(model: cobra.Model, bases: set[str]) -> list[cobra.Metabolite]:
    return [m for m in model.metabolites if _base(m) in bases]


def suggest(model: cobra.Model, query: str, n: int = 8) -> list[str]:
    """Near misses for ``query`` among the model's ids and names: substring hits first, then fuzzy ones."""
    q = query.lower()
    labels = {}
    for met in model.metabolites:
        labels.setdefault(_base(met), met.name or "")
    hits = [b for b, name in labels.items() if q in b.lower() or q in name.lower()]
    by_key = {}
    for b, name in labels.items():
        by_key.setdefault(b.lower(), b)
        if name:
            by_key.setdefault(name.lower(), b)
    for key in difflib.get_close_matches(q, list(by_key), n=n, cutoff=0.75):
        if by_key[key] not in hits:
            hits.append(by_key[key])
    return [f"{b}  ({labels[b]})" for b in sorted(hits[:n], key=lambda b: (q not in b.lower(), b))][:n]


def _cache_dir() -> Path:
    root = os.environ.get("CELLFACTORY_CACHE") or Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache"),
                                                       "cellfactory")
    return Path(root) / "metanetx"


def _sparql(ns: str, local: str) -> str:
    # MetaNetX's RDF links each chemical to identifiers.org URIs. Both URI spellings are tried,
    # since identifiers.org has used both. ChEBI's ids carry their own prefix (CHEBI:4167).
    if not _SAFE_ID.match(local):
        raise MetaNetXUnavailable(f"{local!r} is not a valid {ns} id")
    if ns == "chebi":
        forms = (f"https://identifiers.org/CHEBI:{local}", f"http://identifiers.org/chebi/CHEBI:{local}")
    else:
        raw = {"hmdb": f"HMDB{int(local):07d}" if local.isdigit() else local, "biocyc": f"META:{local}"}.get(ns, local)
        forms = (f"https://identifiers.org/{ns}:{raw}", f"http://identifiers.org/{ns}/{raw}")
    uris = " ".join(f"<{u}>" for u in forms)
    if ns == "metanetx.chemical":
        head = f"VALUES ?chem {{ <https://rdf.metanetx.org/chem/{local}> }}"
    else:
        head = f"VALUES ?xref {{ {uris} }} ?chem <https://rdf.metanetx.org/schema/chemXref> ?xref ."
    return ("SELECT DISTINCT ?chem ?bigg WHERE { " + head +
            " OPTIONAL { ?chem <https://rdf.metanetx.org/schema/chemXref> ?bigg . "
            "FILTER(CONTAINS(STR(?bigg), 'bigg.metabolite')) } }")


def metanetx_lookup(ns: str, local: str, timeout: float = TIMEOUT, use_cache: bool = True) -> dict[str, list[str]]:
    """Ask MetaNetX which chemicals carry ``ns:local``. Returns ``{"mnx": [...], "bigg": [...]}``.

    Raises MetaNetXUnavailable on any network or format problem. Answers,
    including "no such id", are cached on disk. ``timeout`` bounds each socket
    operation, not the whole request, so a server trickling bytes can take longer.
    """
    if os.environ.get("CELLFACTORY_OFFLINE"):
        raise MetaNetXUnavailable("CELLFACTORY_OFFLINE is set")
    cache = _cache_dir() / f"{re.sub(r'[^A-Za-z0-9._-]', '_', f'{ns}-{local}')}.json"
    if use_cache and cache.exists():
        try:
            hit = json.loads(cache.read_text())
            if (isinstance(hit, dict) and set(hit) == {"mnx", "bigg"}
                    and all(isinstance(v, list) and all(isinstance(x, str) for x in v) for v in hit.values())):
                return hit
        except (OSError, ValueError):
            pass  # unreadable or wrong-shaped cache entry: ask again
    url = METANETX_SPARQL + "?" + urllib.parse.urlencode({"query": _sparql(ns, local)})
    req = urllib.request.Request(url, headers={"Accept": "application/sparql-results+json",
                                               "User-Agent": "cellfactory"})
    out: dict[str, list[str]] = {"mnx": [], "bigg": []}
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            rows = json.load(resp)["results"]["bindings"]
        for row in rows:
            chem = str(row.get("chem", {}).get("value", "")).rsplit("/", 1)[-1]
            bigg = re.split(r"[:/]", str(row.get("bigg", {}).get("value", "")))[-1]
            for key, value in (("mnx", chem), ("bigg", bigg)):
                if value and value not in out[key]:
                    out[key].append(value)
    except (urllib.error.URLError, http.client.HTTPException, OSError, TimeoutError, ValueError, KeyError,
            TypeError, AttributeError) as err:
        raise MetaNetXUnavailable(f"MetaNetX query failed ({err!r})") from err
    if use_cache:
        try:
            cache.parent.mkdir(parents=True, exist_ok=True)
            cache.write_text(json.dumps(out))
        except OSError:
            pass
    return out


def resolve_metabolite(
    model: cobra.Model,
    query: str,
    online: bool = True,
    lookup: Callable[[str, str], dict[str, list[str]]] | None = None,
) -> Resolution:
    """Find the model metabolite ``query`` refers to. Never raises for a miss; check ``.ok``.

    ``query`` may be a model metabolite id with or without compartment
    (``etoh``, ``etoh_e``), an everyday name from ``ALIASES``, an exact
    metabolite name, or a database id: ``kegg:C00031``, ``C00031``,
    ``CHEBI:4167``, ``MNXM41``, ``metacyc:GLC``, ``seed:cpd00027``,
    ``HMDB0000122``, ``bigg:glc__D`` or an InChIKey. ``online=False`` never
    contacts MetaNetX. ``lookup`` replaces ``metanetx_lookup`` (for tests).
    """
    res = Resolution(query)
    q = query.strip()
    alias = ALIASES.get(q.lower())
    for via, key in (("alias", alias), ("model id", q)):
        if key:
            mets = [m for m in model.metabolites if m.id in (key, f"{key}_e", f"{key}_c")]
            if mets:
                res.metabolites, res.via = mets, via
                return res
    mets = [m for m in model.metabolites if (m.name or "").lower() == q.lower()]
    if mets:
        res.metabolites, res.via = _mets(model, {_base(m) for m in mets}), "name"
        return _ambiguous(res) if len({_base(m) for m in mets}) > 1 else res
    parsed = parse_identifier(q)
    if parsed is None:
        res.suggestions = suggest(model, q)
        return res
    ns, local = parsed
    index = annotation_index(model)
    bases = index.get((ns, local), set())
    if bases:
        res.metabolites, res.via, res.candidates = _mets(model, bases), f"annotation {ns}", sorted(bases)
        return _ambiguous(res) if len(bases) > 1 else _clear(res)
    if not online:
        res.note = f"No metabolite in the model is annotated with {ns}:{local}, and the MetaNetX lookup was off."
        res.suggestions = suggest(model, q)
        return res
    try:
        hit = (lookup or metanetx_lookup)(ns, local)
    except MetaNetXUnavailable as err:
        log.warning("Could not ask MetaNetX about %s: %s", q, err)
        res.note = (f"No metabolite in the model is annotated with {ns}:{local}, and MetaNetX could not be "
                    f"reached to translate it ({err}). Try the model's own id or name instead.")
        res.suggestions = suggest(model, q)
        return res
    bases = set()
    for mnx in hit.get("mnx", []):
        bases |= index.get(("metanetx.chemical", _norm("metanetx.chemical", mnx)), set())
    for bigg in hit.get("bigg", []):
        bases |= index.get(("bigg.metabolite", bigg), set())
    if not bases:
        found = ", ".join(hit.get("bigg", []) + hit.get("mnx", [])) or "nothing"
        res.note = (f"MetaNetX translates {ns}:{local} to {found}, but no metabolite in the model carries those "
                    "ids. The organism may lack it, or the model may be unannotated.")
        res.suggestions = suggest(model, q)
        return res
    res.metabolites, res.candidates = _mets(model, bases), sorted(bases)
    res.via = f"MetaNetX {', '.join(hit.get('mnx', [])) or ns + ':' + local}"
    return _ambiguous(res) if len(bases) > 1 else _clear(res)


def _clear(res: Resolution) -> Resolution:
    res.candidates = []
    return res


def _ambiguous(res: Resolution) -> Resolution:
    res.candidates = sorted({_base(m) for m in res.metabolites})
    res.metabolites = []
    res.note = ("The database treats them as one compound (often stereoisomers or charge forms sharing an id), "
                "but the model keeps them apart.")
    return res


def resolution_table(resolutions: dict[str, Resolution]) -> pd.DataFrame:
    """One row per resolved role (``{"substrate": ..., "product": ...}``) for ``Report.checks``."""
    rows = []
    for role, r in resolutions.items():
        met = r.metabolites[0] if r.metabolites else None
        rows.append((role, r.query, r.model_id, met.name if met else "", met.formula if met else "", r.via,
                     r.caveat or r.note))
    return pd.DataFrame(rows, columns=["role", "query", "model_metabolite", "name", "formula", "matched_via", "note"])

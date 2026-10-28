"""One call that runs the whole analysis: substrate X, organism Y, product Z."""

from __future__ import annotations

import warnings
from dataclasses import asdict, dataclass, field

import cobra

from .checks import free_energy_checks, gene_essentiality_check, mass_balance_problems, predicted_essential_genes
from .knockouts import SearchResult, candidate_reactions, greedy_knockout_search
from .identifiers import Resolution, resolution_table, resolve_metabolite
from .mdf import max_min_driving_force, mdf_caveats, product_route_fluxes
from .models import load_model, resolve_exchange, set_medium
from .qa import memote_qa
from .simulate import ProductionResult, production_at_max_growth
from .thermodynamics import ThermoUnavailable, apply_directionality, reaction_thermodynamics, thermo_caveats


MAX_ROWS = 20  # per table in the text report; the JSON report has everything


@dataclass
class Report:
    model_id: str
    substrate: str
    product: str
    aerobic: bool
    uptake: float
    production: ProductionResult
    carbon_ratio: float  # carbons in product / carbons in substrate, to turn mol/mol into C-mol/C-mol
    knockouts: SearchResult | None = None
    checks: dict = field(default_factory=dict)
    caveats: list[str] = field(default_factory=list)  # reasons to trust these numbers less, in plain words

    @property
    def carbon_yield(self) -> float:
        """Fraction of the substrate's carbon atoms that end up in the product (guaranteed, at max growth)."""
        return self.production.guaranteed_yield * self.carbon_ratio

    def to_text(self) -> str:
        p = self.production
        lines = [
            f"{self.model_id}: {self.substrate} -> {self.product} ({'aerobic' if self.aerobic else 'anaerobic'}, "
            f"uptake {self.uptake} mmol/gDW/h)",
            f"  growth rate                 {p.growth:.3f} /h",
            f"  theoretical max yield       {p.max_theoretical_yield:.3f} mol/mol",
            f"  guaranteed yield @ growth   {p.guaranteed_yield:.3f} mol/mol  ({self.carbon_yield:.0%} of carbon)",
            f"  best-case yield @ growth    {p.best_case_yield:.3f} mol/mol",
        ]
        if self.knockouts is not None:
            lines += ["", "Greedy knockout search:", self.knockouts.to_frame().to_string(index=False)]
        if self.caveats:
            lines += ["", "Caveats:"] + [f"  - {c}" for c in self.caveats]
        for name, df in self.checks.items():
            lines += ["", f"{name}: {len(df)} rows" + (f" (first {MAX_ROWS} shown)" if len(df) > MAX_ROWS else "")]
            lines.append(df.head(MAX_ROWS).to_string(index=False, max_colwidth=60) if len(df) else "  (none)")
        return "\n".join(lines)

    def to_dict(self) -> dict:
        out = {k: v for k, v in asdict(self).items() if k not in ("production", "knockouts", "checks")}
        p = self.production
        out["production"] = {**asdict(p), "guaranteed_yield": p.guaranteed_yield, "best_case_yield": p.best_case_yield,
                             "carbon_yield": self.carbon_yield}
        out["knockouts"] = [asdict(s) for s in self.knockouts.steps] if self.knockouts else None
        out["checks"] = {k: df.to_dict(orient="records") for k, df in self.checks.items()}
        return out


def _carbons(model: cobra.Model, exchange: str) -> int:
    met = next(iter(model.reactions.get_by_id(exchange).metabolites))
    return met.elements.get("C", 0)


def _resolve(model: cobra.Model, query: str, online: bool) -> Resolution:
    if query in model.reactions:
        return Resolution(query, list(model.reactions.get_by_id(query).metabolites), via="exchange id")
    res = resolve_metabolite(model, query, online=online)
    if not res.ok:
        raise KeyError(res.error(model.id))
    return res


def _mdf(model: cobra.Model, product_ex: str, thermo_table, report: Report) -> list[str]:
    # The route at maximum product yield, with no growth: not the route of the growing cell the yields above describe.
    label = f"the zero-growth max-yield route to {product_ex}"
    try:
        result = max_min_driving_force(model, product_route_fluxes(model, product_ex), table=thermo_table)
    except ThermoUnavailable as err:
        return [f"Max-min driving force not computed: {err}."]
    except RuntimeError as err:  # the LP failed to solve
        return [f"Max-min driving force not computed: {err}."]
    except ValueError as err:
        return [f"Max-min driving force not computed: {err}"]
    report.checks["Max-min driving force"] = result.reactions
    return mdf_caveats(result, label)


def analyze(
    model: cobra.Model | str,
    substrate: str,
    product: str,
    uptake: float = 10.0,
    aerobic: bool = True,
    max_knockouts: int = 0,
    min_growth: float = 0.05,
    run_checks: bool = False,
    qa_gate: bool = True,
    online_ids: bool = True,
    thermo: bool = False,
    thermo_bounds: bool = False,
    mdf: bool = False,
) -> Report:
    """Simulate ``model`` eating ``substrate`` and making ``product``.

    ``model`` is a cobra.Model or anything ``load_model`` accepts. Substrate and
    product can be exchange ids, metabolite ids, or names (see ``resolve_exchange``).
    A product without an exchange reaction gets a demand reaction added.
    The model is modified (medium, possibly a demand reaction, and with ``thermo_bounds`` reaction bounds); pass a copy to keep the original.

    The model is first vetted with MEMOTE (see ``qa.memote_qa``; cached per
    model, so it costs seconds once). A hard failure, such as ATP made from
    nothing, raises ``ModelQAError``; warnings land in the report's caveats.
    ``qa_gate=False`` skips the check.

    Substrate and product may also be database ids (``kegg:C00031``,
    ``CHEBI:4167``, ``MNXM41``, ...), matched through the model's annotations
    and, failing that, MetaNetX online; ``online_ids=False`` never goes online.
    How each was matched is in ``checks["Identifier resolution"]``.

    ``thermo=True`` estimates every reaction's ΔG'° with eQuilibrator (the
    ``[thermo]`` extra) and reports, in ``checks["Reaction thermodynamics"]``,
    which directions physics rules out; ``thermo_bounds=True`` also closes
    those directions on reactions the model left reversible, before simulating
    (see ``thermodynamics.apply_directionality``). Without eQuilibrator both
    are skipped with a caveat.

    ``mdf=True`` takes the route the model uses at maximum product yield and
    computes its Max-min Driving Force (see ``mdf.max_min_driving_force``): can
    every step run forward at once at 1 µM-10 mM? The steps are in
    ``checks["Max-min driving force"]`` and the verdict in the caveats. It needs
    eQuilibrator too, and is skipped with a caveat without it.
    """
    if not isinstance(model, cobra.Model):
        model = load_model(model)
    qa = None
    if qa_gate:
        qa = memote_qa(model)
        qa.raise_if_failed()
    ids = {"substrate": _resolve(model, substrate, online_ids), "product": _resolve(model, product, online_ids)}
    substrate_ex = set_medium(model, substrate if substrate in model.reactions else ids["substrate"].model_id,
                              uptake, aerobic)
    product_ex = resolve_exchange(model, product if product in model.reactions else ids["product"].model_id,
                                  create_demand=True, online=False)
    if product_ex == substrate_ex:
        raise ValueError("Substrate and product are the same compound.")
    thermo_table, thermo_notes = None, []
    if thermo or thermo_bounds:
        try:
            thermo_table = reaction_thermodynamics(model)
        except ThermoUnavailable as err:
            thermo_notes = [f"Reaction thermodynamics not checked: {err}."]
        else:
            tightened = apply_directionality(model, thermo_table) if thermo_bounds else None
            thermo_notes = thermo_caveats(thermo_table, tightened=None if tightened is None else len(tightened))
            if tightened is not None and len(tightened) and model.slim_optimize(error_value=0.0) < 1e-6:
                raise ValueError(f"{model.id} can't grow once {', '.join(tightened['reaction'])} are made one-way "
                                 "from their ΔG'° estimates. Run without thermodynamic bounds to compare.")
    if model.slim_optimize(error_value=0.0) < 1e-6:
        raise ValueError(f"{model.id} can't grow on {substrate_ex} under these conditions "
                         f"({'aerobic' if aerobic else 'anaerobic'}). Try the other oxygen setting or another substrate.")
    sub_c, prod_c = _carbons(model, substrate_ex), _carbons(model, product_ex)
    report = Report(
        model_id=model.id, substrate=substrate_ex, product=product_ex, aerobic=aerobic, uptake=uptake,
        production=production_at_max_growth(model, product_ex, substrate_ex),
        carbon_ratio=prod_c / sub_c if sub_c else 0.0,
    )
    report.checks["Identifier resolution"] = resolution_table(ids)
    report.caveats += [r.caveat for r in ids.values() if r.caveat]
    if qa is not None:
        report.checks[f"Model QA gate ({qa.engine})"] = qa.table
        report.caveats += qa.caveats
    if thermo_table is not None:
        report.checks["Reaction thermodynamics"] = thermo_table[thermo_table["dg0_prime"].notna()]
    report.caveats += thermo_notes
    if mdf:
        report.caveats += _mdf(model, product_ex, thermo_table, report)
    if max_knockouts:
        report.knockouts = greedy_knockout_search(model, product_ex, substrate_ex, max_knockouts=max_knockouts,
                                                  min_growth=min_growth, candidates=candidate_reactions(model))
    if run_checks:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="Solver status is 'infeasible'")
            report.checks["Free-energy checks"] = free_energy_checks(model)
            report.checks["Mass-imbalanced reactions"] = mass_balance_problems(model)
            # The built-in lab data is for E. coli on aerobic glucose minimal medium, so compare under that.
            with model:
                set_medium(model, "EX_glc__D_e", 10.0, aerobic=True)
                ess = gene_essentiality_check(model)
            if (ess["verdict"] != "not in model").any():
                report.checks["Gene essentiality vs E. coli lab data (aerobic glucose)"] = ess
            report.checks["Predicted essential genes"] = predicted_essential_genes(model)
    return report

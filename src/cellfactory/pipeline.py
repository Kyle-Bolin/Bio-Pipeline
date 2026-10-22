"""One call that runs the whole analysis: substrate X, organism Y, product Z."""

from __future__ import annotations

import warnings
from dataclasses import asdict, dataclass, field

import cobra

from .checks import free_energy_checks, gene_essentiality_check, mass_balance_problems, predicted_essential_genes
from .knockouts import SearchResult, candidate_reactions, greedy_knockout_search
from .models import load_model, resolve_exchange, set_medium
from .simulate import ProductionResult, production_at_max_growth


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


def analyze(
    model: cobra.Model | str,
    substrate: str,
    product: str,
    uptake: float = 10.0,
    aerobic: bool = True,
    max_knockouts: int = 0,
    min_growth: float = 0.05,
    run_checks: bool = False,
) -> Report:
    """Simulate ``model`` eating ``substrate`` and making ``product``.

    ``model`` is a cobra.Model or anything ``load_model`` accepts. Substrate and
    product can be exchange ids, metabolite ids, or names (see ``resolve_exchange``).
    A product without an exchange reaction gets a demand reaction added.
    The model is modified (medium, possibly a demand reaction); pass a copy to keep the original.
    """
    if not isinstance(model, cobra.Model):
        model = load_model(model)
    substrate_ex = set_medium(model, substrate, uptake, aerobic)
    product_ex = resolve_exchange(model, product, create_demand=True)
    if product_ex == substrate_ex:
        raise ValueError("Substrate and product are the same compound.")
    if model.slim_optimize(error_value=0.0) < 1e-6:
        raise ValueError(f"{model.id} can't grow on {substrate_ex} under these conditions "
                         f"({'aerobic' if aerobic else 'anaerobic'}). Try the other oxygen setting or another substrate.")
    sub_c, prod_c = _carbons(model, substrate_ex), _carbons(model, product_ex)
    report = Report(
        model_id=model.id, substrate=substrate_ex, product=product_ex, aerobic=aerobic, uptake=uptake,
        production=production_at_max_growth(model, product_ex, substrate_ex),
        carbon_ratio=prod_c / sub_c if sub_c else 0.0,
    )
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

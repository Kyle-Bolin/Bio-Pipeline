"""Static PNG charts for the showcase results (the results page draws interactive versions)."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

COLORS = {"ecoli": "#2a78d6", "yeast": "#eb6834"}
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
PRODUCTS = ["vanillin", "lycopene"]


def _style(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=9)
    ax.grid(axis="y", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def _case(results, org, product):
    return next(c for c in results["cases"] if c["organism"] == org and c["product"] == product)


def yields_chart(results, out: Path):
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    width = 0.36
    for i, org in enumerate(COLORS):
        xs = [p + (i - 0.5) * (width + 0.02) for p in range(len(PRODUCTS))]
        vals = [_case(results, org, p)["max_yield_g"] for p in PRODUCTS]
        bars = ax.bar(xs, vals, width, color=COLORS[org], label=results["organisms"][org])
        ax.bar_label(bars, fmt="%.2f", padding=3, color=INK, fontsize=9)
    ax.set_xticks(range(len(PRODUCTS)), [p.capitalize() for p in PRODUCTS], color=INK)
    ax.set_ylabel("max yield (g product / g glucose)", color=MUTED, fontsize=9)
    ax.set_title("Theoretical maximum yield from glucose", loc="left", color=INK, fontsize=11)
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    _style(ax)
    fig.tight_layout()
    fig.savefig(out / "yields.png", dpi=160)
    plt.close(fig)


def envelope_chart(results, out: Path):
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.4))
    for ax, product in zip(axes, PRODUCTS):
        for org in COLORS:
            c = _case(results, org, product)
            xs = [p["growth"] for p in c["envelope"]]
            ys = [p["max_flux"] / c["uptake"] for p in c["envelope"]]
            ax.plot(xs, ys, color=COLORS[org], linewidth=2, label=results["organisms"][org])
        ax.set_title(product.capitalize(), loc="left", color=INK, fontsize=10)
        ax.set_xlabel("growth rate (1/h)", color=MUTED, fontsize=9)
        ax.set_ylim(bottom=0)
        ax.set_xlim(left=0)
        _style(ax)
    axes[0].set_ylabel("max yield (mol / mol glucose)", color=MUTED, fontsize=9)
    axes[1].legend(frameon=False, fontsize=8, loc="upper right")
    fig.suptitle("Growth vs product: every gram of cells costs product", x=0.01, ha="left", color=INK, fontsize=11)
    fig.tight_layout()
    fig.savefig(out / "envelopes.png", dpi=160)
    plt.close(fig)


def knockout_chart(results, out: Path):
    cases = [c for c in results["cases"] if len(c["knockouts"]) > 1]
    if not cases:
        return
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    labels, gains, colors = [], [], []
    for c in cases:
        parent = c["knockouts"][0]["product_yield"]
        for step in c["knockouts"][1:]:
            labels.append(f"{c['product']}, {c['organism']}: −" + ", −".join(g.split(" ")[0] for g in step["genes"]))
            gains.append(100 * (step["product_yield"] / parent - 1))
            colors.append(COLORS[c["organism"]])
    ys = range(len(labels))[::-1]
    bars = ax.barh(list(ys), gains, color=colors, height=0.6)
    ax.bar_label(bars, fmt="+%.0f%%", padding=3, color=INK, fontsize=9)
    ax.set_yticks(list(ys), labels, fontsize=8, color=INK)
    ax.set_xlabel("predicted product yield vs parent strain (MOMA)", color=MUTED, fontsize=9)
    ax.set_title("Gene deletions predicted to raise yield", loc="left", color=INK, fontsize=11)
    _style(ax)
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", color=GRID, linewidth=0.8)
    fig.tight_layout()
    fig.savefig(out / "knockouts.png", dpi=160)
    plt.close(fig)


def make_charts(results: dict, out: Path) -> None:
    out = Path(out)
    yields_chart(results, out)
    envelope_chart(results, out)
    knockout_chart(results, out)
    make_page(results, out)


def make_page(results: dict, out: Path) -> None:
    """Self-contained results page: the template with the results embedded."""
    template = (Path(__file__).parent / "showcase_page.html").read_text()
    page = template.replace("/*DATA*/null", json.dumps(results, default=float))
    (Path(out) / "index.html").write_text(page)

"""Optional figures for the evaluation outputs. Never required for a run to succeed.

Same contract as `m1_genome/figures.py`: matplotlib is not a hard dependency, a missing
PNG is not a reason to fail an analysis, and `save_figures` returns the paths it wrote so
the caller reports what actually landed rather than guessing.

Each figure makes **one claim**, so it can be dropped into a slide on its own:

``seq_cap_loss``
    the seam defect — M1's `--seq-cap` truncates by `m1_priority`, which scores
    hypothetical proteins 0, so it removes the population `triage` exists to find
``score_disagreement``
    `m1_priority` rank against `triage` rank, and how little the two top-50s share
``triage_ablation``
    how many of the top K each triage component moves, with the components whose
    providers never ran marked as a data gap rather than drawn as zero-effect
``outcome_lift``
    precision against the M3 outcome next to that label's base rate, because the two
    being the same height is the finding

Colour: two categorical slots from a validated palette (blue `#2a78d6`, orange `#eb6834`
on light; `#3987e5` / `#d95926` on dark), which pass the lightness, chroma, CVD-separation,
normal-vision and contrast checks in both modes. Identity is never carried by colour alone —
every series is direct-labelled and the never-fired bars are hatched as well as recoloured.
Text stays in ink, never in a series colour.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any, Sequence

#: Validated categorical slots 1 and 2, plus the fixed status colour for a data gap.
SERIES_1 = "#2a78d6"
SERIES_2 = "#eb6834"
STATUS_CRITICAL = "#d03b3b"
INK = "#0b0b0b"
INK_MUTED = "#52514e"
GRID = "#dcdbd6"


def _style(axes) -> None:
    """Recessive axes and grid: the marks carry the chart, not the furniture."""
    for side in ("top", "right"):
        axes.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        axes.spines[side].set_color(GRID)
    axes.tick_params(colors=INK_MUTED, labelsize=9, length=3)
    axes.title.set_color(INK)
    axes.xaxis.label.set_color(INK_MUTED)
    axes.yaxis.label.set_color(INK_MUTED)


def save_figures(
    out_dir: str | Path,
    *,
    stages: Sequence[tuple[str, int, int]] | None = None,
    ranks: tuple[Sequence[int], Sequence[int]] | None = None,
    shared_at_k: int | None = None,
    spearman: float | None = None,
    k: int = 50,
    ablation: dict[str, Any] | None = None,
    outcome: dict[str, Any] | None = None,
) -> list[Path]:
    """Write whatever the available inputs allow, and return the paths written."""
    try:
        import matplotlib
        matplotlib.use("Agg")  # no display in a pipeline; must precede pyplot
        import matplotlib.pyplot as plt
    except ImportError:
        return []

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []

    def emit(figure, name: str, *, bottom: float | None = None) -> None:
        with warnings.catch_warnings():
            # On a very small input the caption can leave no room for tight_layout to
            # satisfy; the figure still renders correctly and a warning here is noise.
            warnings.filterwarnings("ignore", message="Tight layout not applied")
            figure.tight_layout()
        if bottom is not None:
            figure.subplots_adjust(bottom=bottom)
        for suffix in ("png", "svg"):
            path = out_dir / f"{name}.{suffix}"
            figure.savefig(path, facecolor="white")
            written.append(path)
        plt.close(figure)

    # --- 1. the seam defect -------------------------------------------------------------
    # `stages` is [(label, population, gap_like), ...] in pipeline order.
    if stages and len(stages) >= 2:
        stages = list(stages)[:2]
        labels = [row[0] for row in stages]
        totals = [row[1] for row in stages]
        gaps = [row[2] for row in stages]
        retained_all = [100.0 * totals[i] / totals[0] for i in range(len(stages))]
        retained_gap = [100.0 * gaps[i] / gaps[0] if gaps[0] else 0.0
                        for i in range(len(stages))]

        figure, axes = plt.subplots(figsize=(7.2, 3.9), dpi=150)
        positions = range(len(stages))
        width = 0.38
        bars_all = axes.bar([p - width / 2 for p in positions], retained_all, width,
                            color=SERIES_1, label="all proteins", zorder=3)
        bars_gap = axes.bar([p + width / 2 for p in positions], retained_gap, width,
                            color=SERIES_2, label="hypothetical / uncharacterised", zorder=3)
        for bars, counts in ((bars_all, totals), (bars_gap, gaps)):
            for bar, count in zip(bars, counts):
                axes.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1.5,
                          f"{bar.get_height():.0f}%\nn={count:,}", ha="center", va="bottom",
                          fontsize=8.5, color=INK_MUTED, linespacing=1.3)
        axes.set_xticks(list(positions))
        axes.set_xticklabels(labels, fontsize=9)
        axes.set_ylabel("% of the starting population retained")
        axes.set_ylim(0, 128)
        axes.set_yticks([0, 25, 50, 75, 100])
        axes.grid(axis="y", color=GRID, linewidth=0.7, zorder=0)
        axes.set_axisbelow(True)
        axes.set_title(
            "M1's sequence cap removes the proteins M2 is built to find",
            fontsize=11.5, fontweight="600", loc="left", pad=12,
        )
        axes.legend(frameon=False, fontsize=9, loc="upper right",
                    labelcolor=INK_MUTED, ncols=2)
        _style(axes)
        emit(figure, "seq_cap_loss")

    # --- 2. the two rankings disagree ---------------------------------------------------
    if ranks and len(ranks[0]) == len(ranks[1]) and ranks[0]:
        left, right = ranks
        figure, axes = plt.subplots(figsize=(5.6, 5.6), dpi=150)
        axes.scatter(left, right, s=4, color=SERIES_1, alpha=0.22, linewidths=0, zorder=3)
        axes.axvline(k, color=STATUS_CRITICAL, linewidth=1.2, linestyle="--", zorder=4)
        axes.axhline(k, color=STATUS_CRITICAL, linewidth=1.2, linestyle="--", zorder=4)
        top = max(max(left), max(right))
        axes.set_xlim(0, top)
        axes.set_ylim(0, top)
        axes.set_xlabel(f"m1_priority rank  (cut at {k} dashed)")
        axes.set_ylabel(f"triage rank  (cut at {k} dashed)")
        axes.grid(color=GRID, linewidth=0.6, zorder=0)
        axes.set_axisbelow(True)
        headline = "The two scores pick almost disjoint sets"
        axes.set_title(headline, fontsize=11.5, fontweight="600", loc="left", pad=12)
        note = []
        if shared_at_k is not None:
            note.append(f"{shared_at_k} of {k} shared at the cut")
        if spearman is not None:
            note.append(f"Spearman {spearman:+.3f}")
        if note:
            axes.text(0.02, 0.975, "   ·   ".join(note), transform=axes.transAxes,
                      fontsize=9, color=INK_MUTED, va="top", zorder=6,
                      bbox=dict(facecolor="white", edgecolor="none", pad=3))
        _style(axes)
        emit(figure, "score_disagreement")

    # --- 3. which components decide the selection ---------------------------------------
    effects = list((ablation or {}).get("effects") or [])
    if effects:
        effects = sorted(effects, key=lambda e: e.get("churn_at_k") or 0)
        names = [e["component"] for e in effects]
        churn = [e.get("churn_at_k") or 0 for e in effects]
        dead = [not (e.get("nonzero_proteins") or 0) for e in effects]
        cut = (ablation or {}).get("k", k)

        figure, axes = plt.subplots(figsize=(7.4, 3.9), dpi=150)
        for index, (value, is_dead) in enumerate(zip(churn, dead)):
            if is_dead:
                # a ghost spanning the axis, so "not measured" is visibly not "measured 0"
                axes.barh(index, cut, height=0.6, zorder=3, color="none",
                          edgecolor=STATUS_CRITICAL, hatch="///", linewidth=1.0)
            elif value:
                axes.barh(index, value, height=0.6, zorder=3, color=SERIES_1)
        for index, (value, is_dead, effect) in enumerate(zip(churn, dead, effects)):
            if is_dead:
                axes.text(cut * 0.5, index, "provider never ran — not measured",
                          va="center", ha="center", fontsize=8.5, color=STATUS_CRITICAL,
                          zorder=5,
                          bbox=dict(facecolor="white", edgecolor="none", pad=2.5))
            else:
                axes.text(value + cut * 0.015, index,
                          f"{value}   (w {effect['weight']:+.2f})",
                          va="center", fontsize=8.5, color=INK_MUTED, zorder=5)
        axes.set_yticks(range(len(names)))
        axes.set_yticklabels(names, fontsize=9, family="monospace")
        axes.set_xlabel(f"proteins that leave the top {cut} when the component is removed")
        axes.set_xlim(0, cut * 1.28)
        axes.grid(axis="x", color=GRID, linewidth=0.7, zorder=0)
        axes.set_axisbelow(True)
        axes.set_title(
            "annotation_gap, not pdb_evidence, decides the shortlist",
            fontsize=11.5, fontweight="600", loc="left", pad=12,
        )
        _style(axes)
        emit(figure, "triage_ablation")

    # --- 4. the outcome label cannot discriminate ---------------------------------------
    if outcome and outcome.get("base_rate") is not None:
        precision = ((outcome.get("measured") or {}).get("precision"))
        base = outcome["base_rate"]
        if precision is not None:
            figure, axes = plt.subplots(figsize=(6.4, 4.1), dpi=150)
            values = [100.0 * precision, 100.0 * base]
            bars = axes.bar(["precision\n@ the cut", "base rate\n(any selection)"], values,
                            width=0.5, color=[SERIES_1, SERIES_2], zorder=3)
            for bar, value in zip(bars, values):
                axes.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.6,
                          f"{value:.1f}%", ha="center", va="bottom", fontsize=10,
                          color=INK_MUTED)
            axes.set_ylim(90, 103)
            axes.set_yticks([90, 95, 100])
            axes.set_ylabel("% of labelled proteins that are dockable")
            axes.grid(axis="y", color=GRID, linewidth=0.7, zorder=0)
            axes.set_axisbelow(True)
            axes.set_title(
                f"M3's gate passes almost everything — lift {100 * (precision - base):+.1f} pp",
                fontsize=11.5, fontweight="600", loc="left", pad=12,
            )
            _style(axes)
            # Below the axes, not over the marks: red on a blue bar was unreadable.
            figure.text(0.5, 0.015,
                        "Two bars this close mean the label cannot separate a good selection "
                        "from a bad one — it is not a test yet.",
                        ha="center", fontsize=8.5, color=STATUS_CRITICAL)
            emit(figure, "outcome_lift", bottom=0.22)

    return written


__all__ = ["SERIES_1", "SERIES_2", "save_figures"]

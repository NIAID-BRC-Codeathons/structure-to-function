"""Figures for the evaluation outputs. Optional: a missing figure never fails a run.

Same dependency contract as `m1_genome/figures.py` — matplotlib is not required, and
`save_figures` returns the paths it wrote so the caller reports what actually landed.

## Two styles, because the two audiences want opposite things

``publication`` (default)
    No title inside the figure. A journal sets the legend in the caption, and a title baked
    into the artwork has to be edited out in production. Sized to real column widths in
    millimetres, typeset in Arial at 7 pt, written as PDF and SVG with **editable text**
    plus a 600 dpi PNG, and accompanied by `captions.md`.
``slide``
    The headline goes *in* the figure, because a slide has no caption and the audience
    reads the chart for four seconds. Wider, larger type, one claim spelled out on top.

Both styles draw the same marks through the same `_draw_*` functions, so they cannot drift.

## The two settings that decide whether a figure is accepted

matplotlib defaults to `pdf.fonttype = 3` and `svg.fonttype = "path"`. Type 3 fonts are
rejected outright by Nature, Elsevier and IEEE production, and path-outlined SVG text
cannot be corrected by a copy-editor or read by a screen reader. Publication mode sets
`pdf.fonttype = ps.fonttype = 42` (TrueType) and `svg.fonttype = "none"` (live text).

## Colour

Two slots from a validated categorical palette — blue `#2a78d6`, orange `#eb6834` — which
pass the lightness-band, chroma-floor, CVD-separation, normal-vision and surface-contrast
checks in both light and dark. Identity never rests on colour alone: every series is
direct-labelled, and a component whose provider never ran is drawn as a hatched ghost
reading "not measured" rather than as a zero-length bar, because "not measured" and
"measured zero" are different claims and a blank must not leave the reader guessing which.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any, Callable, Sequence

#: Validated categorical slots 1 and 2, plus the fixed status colour for a data gap.
SERIES_1 = "#2a78d6"
SERIES_2 = "#eb6834"
STATUS_CRITICAL = "#d03b3b"
INK = "#0b0b0b"
INK_MUTED = "#52514e"
GRID = "#dcdbd6"

MM = 1.0 / 25.4
#: Journal column widths: 89 mm single, 120 mm one-and-a-half, 183 mm double.
SINGLE_COLUMN = 89 * MM
WIDE_COLUMN = 120 * MM
DOUBLE_COLUMN = 183 * MM

STYLES = ("publication", "slide")

#: PDF first, because it is the format a journal asks for.
PUBLICATION_FORMATS = ("pdf", "svg", "png")
SLIDE_FORMATS = ("png", "svg")


def _rc(style: str) -> dict[str, Any]:
    """matplotlib settings for one style. Applied through `rc_context`, never globally."""
    common = {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        # Without this, mathtext ($w$, $\rho$) silently falls back to DejaVu and the
        # italics sit in a different typeface from every other glyph in the figure.
        "mathtext.fontset": "custom",
        "mathtext.rm": "Arial",
        "mathtext.it": "Arial:italic",
        "mathtext.bf": "Arial:bold",
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "axes.labelcolor": INK_MUTED,
        "text.color": INK,
        "xtick.color": INK_MUTED,
        "ytick.color": INK_MUTED,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    }
    if style == "slide":
        return {**common, "font.size": 9.5, "axes.titlesize": 11.5,
                "axes.labelsize": 9.5, "xtick.labelsize": 9, "ytick.labelsize": 9,
                "savefig.dpi": 200}
    return {
        **common,
        "font.size": 7,
        "axes.titlesize": 8,
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "savefig.dpi": 600,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    }


def _style_axes(axes, *, grid_axis: str | None = "y") -> None:
    """Recessive furniture: the marks carry the chart."""
    for side in ("top", "right"):
        axes.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        axes.spines[side].set_color(GRID)
    if grid_axis:
        axes.grid(axis=grid_axis, color=GRID, linewidth=0.5, zorder=0)
        axes.set_axisbelow(True)


# --------------------------------------------------------------------------- the panels
# Each draws into an Axes it is handed, so a standalone figure and the same panel inside
# the composite are the identical marks.

def _draw_seq_cap_loss(axes, stages: Sequence[tuple[str, int, int]], *, small: bool) -> None:
    stages = list(stages)[:2]
    totals = [row[1] for row in stages]
    gaps = [row[2] for row in stages]
    retained_all = [100.0 * t / totals[0] for t in totals]
    retained_gap = [100.0 * g / gaps[0] if gaps[0] else 0.0 for g in gaps]

    positions, width = range(len(stages)), 0.34
    bars_all = axes.bar([p - width / 2 for p in positions], retained_all, width,
                        color=SERIES_1, label="all CDS", zorder=3)
    bars_gap = axes.bar([p + width / 2 for p in positions], retained_gap, width,
                        color=SERIES_2, label="hypothetical / uncharacterised", zorder=3)
    for bars, counts in ((bars_all, totals), (bars_gap, gaps)):
        for bar, count in zip(bars, counts):
            axes.annotate(f"{bar.get_height():.0f}%\n$n$={count:,}",
                          (bar.get_x() + bar.get_width() / 2, bar.get_height()),
                          textcoords="offset points", xytext=(0, 2.5),
                          ha="center", va="bottom",
                          fontsize=6 if small else 8.5, color=INK_MUTED, linespacing=1.25)
    axes.set_xticks(list(positions))
    axes.set_xticklabels([row[0] for row in stages])
    axes.set_ylabel("Retained (% of starting population)")
    axes.set_ylim(0, 138)
    axes.set_yticks([0, 25, 50, 75, 100])
    axes.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.04),
                ncols=1 if small else 2, handlelength=1.1, handleheight=0.9,
                borderpad=0, labelspacing=0.25, columnspacing=1.2)
    _style_axes(axes)


def _draw_score_disagreement(axes, ranks, *, k: int, shared_at_k, spearman,
                             small: bool) -> None:
    left, right = ranks
    axes.scatter(left, right, s=1.4 if small else 4, color=SERIES_1,
                 alpha=0.2, linewidths=0, rasterized=True, zorder=3)
    axes.axvline(k, color=STATUS_CRITICAL, linewidth=0.8, linestyle=(0, (3, 2)), zorder=4)
    axes.axhline(k, color=STATUS_CRITICAL, linewidth=0.8, linestyle=(0, (3, 2)), zorder=4)
    top = max(max(left), max(right))
    axes.set_xlim(0, top)
    axes.set_ylim(0, top)
    axes.set_xlabel("m1_priority rank")
    axes.set_ylabel("triage rank")
    note = []
    if shared_at_k is not None:
        note.append(f"{shared_at_k}/{k} shared at the cut")
    if spearman is not None:
        note.append(f"Spearman $\\rho$ = {spearman:+.3f}")
    if note:
        axes.text(0.035, 0.975, "\n".join(note), transform=axes.transAxes,
                  fontsize=6 if small else 8.5, color=INK_MUTED, va="top", zorder=6,
                  linespacing=1.4,
                  bbox=dict(facecolor="white", edgecolor="none", pad=2))
    _style_axes(axes, grid_axis=None)
    axes.grid(color=GRID, linewidth=0.5, zorder=0)
    axes.set_axisbelow(True)


def _draw_triage_ablation(axes, ablation: dict[str, Any], *, small: bool) -> None:
    effects = sorted(ablation.get("effects") or [],
                     key=lambda e: e.get("churn_at_k") or 0)
    churn = [e.get("churn_at_k") or 0 for e in effects]
    dead = [not (e.get("nonzero_proteins") or 0) for e in effects]
    cut = ablation.get("k", 50)
    size = 5.8 if small else 8.5

    for index, (value, is_dead) in enumerate(zip(churn, dead)):
        if is_dead:
            axes.barh(index, cut, height=0.62, zorder=3, color="none",
                      edgecolor=STATUS_CRITICAL, hatch="///", linewidth=0.7)
        elif value:
            axes.barh(index, value, height=0.62, zorder=3, color=SERIES_1)
    for index, (value, is_dead, effect) in enumerate(zip(churn, dead, effects)):
        if is_dead:
            axes.text(cut * 0.5, index, "not measured", va="center", ha="center",
                      fontsize=size, color=STATUS_CRITICAL, zorder=5,
                      bbox=dict(facecolor="white", edgecolor="none", pad=1.4))
        else:
            axes.text(value + cut * 0.02, index,
                      f"{value}  ($w$ {effect['weight']:+.2f})",
                      va="center", fontsize=size, color=INK_MUTED, zorder=5)
    axes.set_yticks(range(len(effects)))
    axes.set_yticklabels([e["component"] for e in effects],
                         fontsize=size, family="monospace")
    axes.set_xlabel(f"Proteins leaving the top {cut} when the component is removed")
    axes.set_xlim(0, cut * 1.30)
    axes.set_ylim(-0.7, len(effects) - 0.3)
    _style_axes(axes, grid_axis="x")


def _draw_outcome_lift(axes, outcome: dict[str, Any], *, small: bool) -> None:
    precision = (outcome.get("measured") or {}).get("precision")
    base = outcome["base_rate"]
    values = [100.0 * precision, 100.0 * base]
    bars = axes.bar(["Precision\nat the cut", "Base rate\n(any selection)"], values,
                    width=0.46, color=[SERIES_1, SERIES_2], zorder=3)
    for bar, value in zip(bars, values):
        axes.annotate(f"{value:.1f}%", (bar.get_x() + bar.get_width() / 2, value),
                      textcoords="offset points", xytext=(0, 2.5),
                      ha="center", va="bottom",
                      fontsize=6 if small else 9, color=INK_MUTED)
    # Full baseline, deliberately. Truncating the y-axis would magnify a 0.7-point gap
    # into a visible difference and argue the opposite of what this panel is for: the two
    # bars being the same height *is* the finding.
    axes.set_ylim(0, 118)
    axes.set_yticks([0, 25, 50, 75, 100])
    axes.set_ylabel("Dockable (% of labelled)")
    _style_axes(axes)


# ------------------------------------------------------------------------------ captions

def _captions(*, stages, ranks, shared_at_k, spearman, k, ablation, outcome) -> str:
    """Legends where a journal wants them: outside the artwork, numbered, self-contained."""
    lines = ["# Figure legends", "",
             "Generated by `s2f.evaluation`. Every number is reproducible from the JSON "
             "beside these files in `runs/<id>/eval/`.", ""]
    number = 0
    if stages:
        number += 1
        total_0, gap_0 = stages[0][1], stages[0][2]
        total_1, gap_1 = stages[1][1], stages[1][2]
        overall_loss = 1 - total_1 / total_0
        # A genome smaller than the cap loses nothing, so there is no ratio to quote and
        # claiming one would invent a finding out of a run where the cap never fired.
        differential = (
            f", a {(1 - gap_1 / gap_0) / overall_loss:.1f}-fold higher loss rate"
            if overall_loss > 0 and gap_0 else ""
        )
        headline = (
            "M1's sequence cap depletes the uncharacterised proteins that M2 is designed "
            "to select." if overall_loss > 0 else
            "M1's sequence cap did not fire on this genome."
        )
        lines += [
            f"**Figure {number}. {headline}** Proportion of each population surviving "
            "`--seq-cap` (default 4,000), which truncates by `m1_priority` rank. All CDS: "
            f"{total_1:,}/{total_0:,} retained ({100 * total_1 / total_0:.0f}%). "
            f"Hypothetical or uncharacterised CDS: {gap_1:,}/{gap_0:,} "
            f"({100 * gap_1 / gap_0:.0f}%){differential}. "
            "`m1_priority` scores a hypothetical protein 0, so these sort to the bottom of "
            "the ranking the cap is applied to, while M2's `annotation_gap` component gives "
            "them weight 0.30. Labels give percentage retained and absolute counts.", ""]
    if ranks:
        number += 1
        lines += [
            f"**Figure {number}. The two ranking scores select near-disjoint protein sets.** "
            "Rank under `m1_priority` (host-interaction evidence in the annotation) against "
            "rank under `triage` (structural tractability) for all "
            f"{len(ranks[0]):,} proteins scored by both. Dashed lines mark the top-{k} cut "
            f"used for selection; {shared_at_k} of {k} proteins fall inside both. Spearman "
            f"rho = {spearman:+.3f}. Diagonal streaks are blocks of tied scores resolved by "
            "the shared protein-length tie-break. Points are rasterised for file size; all "
            "other elements are vector.", ""]
    if ablation:
        number += 1
        cut = ablation.get("k", k)
        movers = sum(1 for e in ablation["effects"] if e.get("churn_at_k"))
        lines += [
            f"**Figure {number}. `annotation_gap`, not `pdb_evidence`, determines the "
            "shortlist.** Leave-one-out ablation of the eight `triage` components: each is "
            "set to zero, M2's own `triage_score_from_components` re-applied, and the "
            f"ranking rebuilt. Bars give the number of proteins leaving the top {cut}; "
            f"*w* is the production weight. {movers} of {len(ablation['effects'])} "
            "components alter the selection. Hatched bars mark components whose provider "
            "did not run in this configuration, which is a data gap rather than a measured "
            "null. No production weight was modified.", ""]
    if outcome:
        number += 1
        precision = (outcome.get("measured") or {}).get("precision")
        base = outcome["base_rate"]
        lines += [
            f"**Figure {number}. The M3 structure-quality gate cannot discriminate between "
            f"selections.** Precision of the `triage` top-{k} against M3's "
            f"`usable_for_docking` label ({100 * precision:.1f}%) beside the base rate of "
            f"that label across all labelled proteins ({100 * base:.1f}%): a lift of "
            f"{100 * (precision - base):+.1f} percentage points. The gate passes almost "
            "every structure it is given, so high precision against it is not evidence of "
            "good selection. The comparison is also partly circular, since `pdb_evidence` "
            "contributes 0.40 to the triage score and the presence of a PDB hit is most of "
            "what makes a structure retrievable.", ""]
    return "\n".join(lines)


# --------------------------------------------------------------------------------- entry

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
    style: str = "publication",
    composite: bool = True,
) -> list[Path]:
    """Write whatever the available inputs allow, and return the paths written."""
    if style not in STYLES:
        raise ValueError(f"unknown figure style {style!r}; known: {', '.join(STYLES)}")
    try:
        import matplotlib
        matplotlib.use("Agg")  # no display in a pipeline; must precede pyplot
        import matplotlib.pyplot as plt
    except ImportError:
        return []

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    publication = style == "publication"
    formats = PUBLICATION_FORMATS if publication else SLIDE_FORMATS

    have_stages = bool(stages and len(stages) >= 2 and stages[0][1] and stages[0][2])
    have_ranks = bool(ranks and len(ranks[0]) == len(ranks[1]) and ranks[0])
    have_ablation = bool((ablation or {}).get("effects"))
    have_outcome = bool(outcome and outcome.get("base_rate") is not None
                        and (outcome.get("measured") or {}).get("precision") is not None)

    panels: list[tuple[str, str, float, float, Callable[[Any, bool], None]]] = []
    if have_stages:
        panels.append((
            "seq_cap_loss",
            "M1's sequence cap removes the proteins M2 is built to find",
            SINGLE_COLUMN, 2.5,
            lambda ax, small: _draw_seq_cap_loss(ax, stages, small=small)))
    if have_ranks:
        panels.append((
            "score_disagreement",
            "The two scores pick almost disjoint sets",
            SINGLE_COLUMN, 3.0,
            lambda ax, small: _draw_score_disagreement(
                ax, ranks, k=k, shared_at_k=shared_at_k, spearman=spearman, small=small)))
    if have_ablation:
        effects = ablation["effects"]
        movers = sorted((e for e in effects if e.get("churn_at_k")),
                        key=lambda e: e["churn_at_k"], reverse=True)
        # Derived rather than written: the first version of this headline named the top
        # mover, and went stale the moment M2 gained a ninth component (#60) and the
        # ordering changed. A count cannot.
        headline = (f"{len(movers)} of {len(effects)} triage components move the top {k}"
                    if movers else "No triage component moves the selection")
        panels.append((
            "triage_ablation",
            headline,
            WIDE_COLUMN, 2.6,
            lambda ax, small: _draw_triage_ablation(ax, ablation, small=small)))
    if have_outcome:
        panels.append((
            "outcome_lift",
            "M3's gate passes almost everything",
            SINGLE_COLUMN, 2.5,
            lambda ax, small: _draw_outcome_lift(ax, outcome, small=small)))
    if not panels:
        return []

    def emit(figure, name: str) -> None:
        with warnings.catch_warnings():
            # A dense caption can leave tight_layout nothing to solve; the figure still
            # renders correctly and the warning is noise.
            warnings.filterwarnings("ignore", message="Tight layout not applied")
            figure.tight_layout()
        for suffix in formats:
            path = out_dir / f"{name}.{suffix}"
            figure.savefig(path)
            written.append(path)
        plt.close(figure)

    with plt.rc_context(_rc(style)):
        for name, headline, width, height, draw in panels:
            if publication:
                figure, axes = plt.subplots(figsize=(width, height))
            else:
                figure, axes = plt.subplots(figsize=(width * 1.9, height * 1.35))
                axes.set_title(headline, fontweight="600", loc="left", pad=10)
            draw(axes, publication)
            emit(figure, name)

        if composite and publication and len(panels) > 1:
            columns = 2
            rows = (len(panels) + columns - 1) // columns
            figure, grid = plt.subplots(rows, columns,
                                        figsize=(DOUBLE_COLUMN, 2.6 * rows), squeeze=False)
            flat = [ax for row in grid for ax in row]
            for index, panel in enumerate(panels):
                panel[4](flat[index], True)
                flat[index].text(-0.17, 1.07, "ABCDEFGH"[index],
                                 transform=flat[index].transAxes, fontsize=9,
                                 fontweight="bold", va="top", ha="left", color=INK)
            for spare in flat[len(panels):]:
                spare.set_visible(False)
            emit(figure, "figure_composite")

    caption_path = out_dir / "captions.md"
    caption_path.write_text(
        _captions(stages=list(stages)[:2] if have_stages else None,
                  ranks=ranks if have_ranks else None,
                  shared_at_k=shared_at_k, spearman=spearman, k=k,
                  ablation=ablation if have_ablation else None,
                  outcome=outcome if have_outcome else None),
        encoding="utf-8")
    written.append(caption_path)
    return written


__all__ = [
    "DOUBLE_COLUMN",
    "PUBLICATION_FORMATS",
    "SERIES_1",
    "SERIES_2",
    "SINGLE_COLUMN",
    "SLIDE_FORMATS",
    "STYLES",
    "WIDE_COLUMN",
    "save_figures",
]

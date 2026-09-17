"""Optional publication-quality figures (matplotlib). Never required for a run to succeed.

matplotlib is not in `requirements-m1.txt` as a hard dependency: the HTML report already
carries inline SVG versions of everything here, and these exist only for dropping into a
paper or a slide. If matplotlib is absent, `save_figures` returns an empty list and the run
carries on — a missing PNG is not a reason to fail an analysis.

The functions return the paths they wrote so the caller can report them rather than having
to guess what landed on disk.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any


def save_figures(out_dir: str | Path, *, priorities: list[dict[str, Any]],
                 specialty_counts: dict[str, int],
                 tree: dict[str, Any] | None = None) -> list[Path]:
    """Write the figures that the available libraries allow, and return their paths."""
    try:
        import matplotlib
        matplotlib.use("Agg")  # no display in a pipeline; must precede pyplot
        import matplotlib.pyplot as plt
    except ImportError:
        return []

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []

    # Coerced rather than trusted: these come from `report.json`, whose `genome` section
    # allows unknown properties and does not type these sub-shapes.
    counts = {str(k): v for k, v in (specialty_counts or {}).items()
              if isinstance(v, (int, float)) and not isinstance(v, bool)}
    items = [(k, v) for k, v in sorted(counts.items(), key=lambda kv: kv[1]) if v]
    if items:
        figure, axes = plt.subplots(figsize=(7, 0.5 * len(items) + 1.2), dpi=150)
        axes.barh([k for k, _ in items], [v for _, v in items], color="#0b5cad")
        axes.set_xlabel("gene count")
        axes.set_title("Specialty-gene composition (BV-BRC)")
        for index, (_, value) in enumerate(items):
            axes.text(value, index, f" {value}", va="center", fontsize=9)
        figure.tight_layout()
        path = out_dir / "specialty_composition.png"
        figure.savefig(path)
        plt.close(figure)
        written.append(path)

    scores = [int(p["score"]) for p in priorities or []
              if isinstance(p, dict) and isinstance(p.get("score"), (int, float))
              and not isinstance(p.get("score"), bool)]
    if scores and max(scores) != min(scores):
        figure, axes = plt.subplots(figsize=(6, 3.4), dpi=150)
        axes.hist(scores, bins=range(min(scores), max(scores) + 2),
                  color="#0e7490", edgecolor="white")
        axes.set_xlabel("M1 pathogenesis priority score")
        axes.set_ylabel("proteins")
        axes.set_title("Priority-score distribution")
        figure.tight_layout()
        path = out_dir / "priority_distribution.png"
        figure.savefig(path)
        plt.close(figure)
        written.append(path)

    if tree and tree.get("ok"):
        try:
            from scipy.cluster.hierarchy import dendrogram
        except ImportError:
            return written
        # Redrawn from icoord/dcoord rather than re-clustering: the figure must show the
        # same tree the report does, and re-running linkage could resolve a tie differently.
        figure, axes = plt.subplots(figsize=(7.5, 0.4 * len(tree["ivl"]) + 1.6), dpi=150)
        for xs, ys in zip(tree["dcoord"], tree["icoord"]):
            axes.plot(xs, ys, color="#94a3b8", linewidth=1.4)
        axes.set_yticks([10 * i + 5 for i in range(len(tree["ivl"]))])
        labels = []
        for label in tree["ivl"]:
            meta = tree["leaf_meta"].get(label, {})
            marker = "* " if meta.get("is_query") else ""
            labels.append(f"{marker}{label}")
        axes.set_yticklabels(labels, fontsize=8, fontstyle="italic")
        for tick, label in zip(axes.get_yticklabels(), tree["ivl"]):
            if tree["leaf_meta"].get(label, {}).get("is_pathogen"):
                tick.set_color("#b91c1c")
        axes.set_xlabel("gene-content distance (Jaccard)")
        axes.set_title("Gene-content phylogeny (UPGMA)")
        axes.spines[["top", "right", "left"]].set_visible(False)
        figure.tight_layout()
        path = out_dir / "gene_content_tree.png"
        figure.savefig(path)
        plt.close(figure)
        written.append(path)

    return written

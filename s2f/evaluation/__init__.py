"""Calibration harness: score a ranking against ground truth (issues #49, #12).

Not a pipeline module. It owns no `report.json` key, writes nothing upstream, and reads a
finished run the way a reviewer would. See `docs/08-evaluation.md`.

Four things it measures:

- **calibration** — a ranking against a curated truth set (`truthset`, `metrics`)
- **outcome** — the same ranking against what M3 actually collected, no curation needed
  (`outcome`)
- **ablation** — which of M2's triage components decide the selection (`ablation`)
- **agreement** — how far `m1_priority` and `triage` diverge, which the data contract
  predicts but nobody had measured (`agreement`)
"""

from .ablation import Ablation, ComponentEffect, ablate
from .agreement import Agreement, compare, spearman
from .metrics import Evaluation, Score, evaluate, rank_by_baseline, rank_by_priority
from .outcome import OutcomeSet, outcomes_from_report
from .rankings import Ranking, ablated_triage_score, available_scores, build_ranking
from .truthset import TruthRow, TruthSet, load_truthset, match_proteins, standalone_token

__all__ = [
    "Ablation",
    "Agreement",
    "ComponentEffect",
    "Evaluation",
    "OutcomeSet",
    "Ranking",
    "Score",
    "TruthRow",
    "TruthSet",
    "ablate",
    "ablated_triage_score",
    "available_scores",
    "build_ranking",
    "compare",
    "evaluate",
    "load_truthset",
    "match_proteins",
    "outcomes_from_report",
    "rank_by_baseline",
    "rank_by_priority",
    "spearman",
    "standalone_token",
]

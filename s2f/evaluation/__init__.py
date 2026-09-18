"""Calibration harness: score a ranking against a curated truth set (issue #49).

Not a pipeline module. It owns no `report.json` key, writes nothing upstream, and reads a
finished run the way a reviewer would. See `docs/08-evaluation.md`.
"""

from .metrics import Evaluation, Score, evaluate, rank_by_baseline, rank_by_priority
from .truthset import TruthRow, TruthSet, load_truthset, match_proteins, standalone_token

__all__ = [
    "Evaluation",
    "Score",
    "TruthRow",
    "TruthSet",
    "evaluate",
    "load_truthset",
    "match_proteins",
    "rank_by_baseline",
    "rank_by_priority",
    "standalone_token",
]

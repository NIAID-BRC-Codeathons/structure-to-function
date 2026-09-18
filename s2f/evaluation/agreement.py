"""How far apart two rankings of the same proteins are, and where.

`00a-data-contract.md` states that `m1_priority` and `triage` answer different questions
and are expected to disagree: "a virulence factor with no solved homolog ranks high here
and low there, and that is the pipeline working." That is a design claim, and until it is
measured it stays a claim. If the two agree almost perfectly, one of them is redundant. If
they disagree almost completely, the pipeline's output depends entirely on which one M3
is wired to — and it is wired to `triage`.

Spearman is computed here rather than imported. `scipy` is an optional dependency in this
repo (M1 degrades without it), and a rank correlation over a few thousand integers does not
justify making it required. Ties get average ranks, which is what makes this Spearman
rather than a Pearson correlation on ordinals — and it matters a great deal here, because
both scores are small and heavily tied: 74 HS11286 proteins share `m1_priority` 5.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Sequence

from .rankings import Ranking


def average_ranks(values: Sequence[float]) -> list[float]:
    """Ranks with ties averaged, ascending. The `1, 2.5, 2.5, 4` convention."""
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        shared = (i + j) / 2 + 1
        for position in range(i, j + 1):
            ranks[order[position]] = shared
        i = j + 1
    return ranks


def spearman(left: Sequence[float], right: Sequence[float]) -> float | None:
    """Spearman's rho over paired values, ties averaged.

    Returns None when either side is constant, because a correlation with a zero-variance
    variable is undefined rather than zero — and a score that is constant across the
    population is a finding in itself, not a correlation of 0.
    """
    if len(left) != len(right):
        raise ValueError("paired sequences must be the same length")
    if len(left) < 2:
        return None
    a, b = average_ranks(left), average_ranks(right)
    mean_a, mean_b = sum(a) / len(a), sum(b) / len(b)
    da = [x - mean_a for x in a]
    db = [y - mean_b for y in b]
    numerator = sum(x * y for x, y in zip(da, db))
    denominator = (sum(x * x for x in da) * sum(y * y for y in db)) ** 0.5
    return (numerator / denominator) if denominator else None


@dataclass
class Disagreement:
    """One protein the two rankings place very differently."""

    feature_id: str
    gene: str | None
    product: str | None
    left_rank: int
    right_rank: int
    left_score: float
    right_score: float

    @property
    def gap(self) -> int:
        return abs(self.left_rank - self.right_rank)


@dataclass
class Agreement:
    """Two rankings compared, at the top and over the whole population."""

    left: str
    right: str
    common: int
    spearman: float | None
    overlap_at_k: int
    k: int
    jaccard_at_k: float
    only_left: list[str] = field(default_factory=list)
    only_right: list[str] = field(default_factory=list)
    biggest_disagreements: list[dict[str, Any]] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "left": self.left,
            "right": self.right,
            "common": self.common,
            "spearman": self.spearman,
            "k": self.k,
            "overlap_at_k": self.overlap_at_k,
            "jaccard_at_k": self.jaccard_at_k,
            "only_left": self.only_left,
            "only_right": self.only_right,
            "biggest_disagreements": self.biggest_disagreements,
        }


def compare(
    left: Ranking,
    right: Ranking,
    proteins: Sequence[dict[str, Any]],
    *,
    k: int = 50,
    show: int = 15,
) -> Agreement:
    """Compare two rankings over the proteins both of them scored.

    Restricted to the intersection: a protein one score never reached cannot disagree with
    the other, and padding it in at the bottom would manufacture correlation.
    """
    common = [fid for fid in left.order if fid in right.scores]
    left_position, right_position = left.position(), right.position()
    by_feature = {
        str(p.get("feature_id")): p for p in proteins if p.get("feature_id")
    }

    rho = spearman(
        [left.scores[fid] for fid in common],
        [right.scores[fid] for fid in common],
    )

    left_top, right_top = set(left.top(k)), set(right.top(k))
    union = left_top | right_top
    overlap = len(left_top & right_top)

    disagreements = sorted(
        (
            Disagreement(
                feature_id=fid,
                gene=by_feature.get(fid, {}).get("gene"),
                product=by_feature.get(fid, {}).get("product"),
                left_rank=left_position[fid],
                right_rank=right_position[fid],
                left_score=left.scores[fid],
                right_score=right.scores[fid],
            )
            for fid in common
        ),
        key=lambda d: d.gap,
        reverse=True,
    )[:show]

    return Agreement(
        left=left.name,
        right=right.name,
        common=len(common),
        spearman=rho,
        overlap_at_k=overlap,
        k=k,
        jaccard_at_k=(overlap / len(union)) if union else 0.0,
        only_left=sorted(left_top - right_top),
        only_right=sorted(right_top - left_top),
        biggest_disagreements=[
            {
                "feature_id": d.feature_id,
                "gene": d.gene,
                "product": d.product,
                f"{left.name}_rank": d.left_rank,
                f"{right.name}_rank": d.right_rank,
                f"{left.name}_score": d.left_score,
                f"{right.name}_score": d.right_score,
                "gap": d.gap,
            }
            for d in disagreements
        ],
    )


__all__ = ["Agreement", "Disagreement", "average_ranks", "compare", "spearman"]

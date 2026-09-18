"""Scoring a ranking against a truth set, and auditing where its evidence came from.

Issue #49 asks for "precision/recall of the top 50 against it, recorded before any weight
change". Three things make that harder than it sounds, and each is handled explicitly
rather than averaged away.

**The truth set is sparse.** 5,523 CDS, roughly 200 curated symbols. Classical precision
would divide by 50 and count every unlabelled protein in the top 50 as a false positive,
which is false: most of them are simply proteins nobody has adjudicated. So precision here
is computed over the **labelled subset** of the top K — `positives / (positives +
negatives)` — and `unlabelled` is reported beside it so the coverage is visible. A
precision of 0.8 over 10 labelled proteins is a different claim from 0.8 over 40, and the
output says which one it is.

**Ties make "top 50" ambiguous.** The first live run put 14 proteins at the score ceiling
and had a non-monotonic dip at score 4. If the cut falls inside a tie block, the members
that land above the line did so on the length tie-break, not on evidence. `ceiling_ties`
and `cut_inside_tie_block` report that instead of hiding it.

**The scorer can see the truth set's vocabulary.** `m1_priority` reads
`gene + product + classifications`, so a gene-symbol truth set is not fully independent of
it. :func:`leakage_audit` is the response: for every labelled protein in the top K it
splits the score into evidence a curated database supplied (VFDB, Victors, CARD, NDARO,
drug-target and essentiality tables) and evidence the product string supplied. A positive
that scored only on text-derived components is the regex agreeing with itself, and the
number that matters is how many of them there are.

Nothing here writes to `report.json`. `docs/00-architecture.md` gives every top-level key
an owner and `common.schema.validate` rejects a section it does not know, so inventing an
`evaluation` key would be exactly the improvised schema change the contract forbids. The
module reads the report and writes its own directory; adding the key is proposed as a
separate decision.
"""

from __future__ import annotations

import hashlib
import json
import random
import re
from dataclasses import asdict, dataclass, field
from typing import Any, Callable, Iterable, Sequence

from ..m1_genome.priority import PATHOGENESIS_RE, WEIGHTS
from .truthset import Match, TruthSet

#: Which score components a curated database supplies, and which the product text does.
#: Keyed by the leading words of the `reason` strings `priority.score_protein` writes;
#: if those strings change, `tests/test_eval_metrics.py::test_every_reason_is_classified`
#: fails rather than this silently reclassifying everything as text-derived.
CURATED_REASONS = {
    "Virulence factor": "virulence",
    "Antibiotic-resistance determinant": "amr",
    "Known drug target": "drug_target",
    "Essential-gene homolog": "essential",
    "Transporter": "transporter",
    "Close human homolog": "human_homolog",
}
TEXT_REASONS = {
    "Host-interaction mechanism": "mechanism",
    "Predicted surface or secreted": "surface",
    "Named / characterised gene": "named_gene",
}

#: Ranking strategies the measured ranking is compared against. A score that cannot beat
#: `keyword` is not earning its nine components: that baseline is one regex over the
#: product string, which is the cheapest thing anyone could have shipped instead.
BASELINES = ("keyword", "length", "random")

RANDOM_SEED = 49


def weights_fingerprint(weights: dict[str, int] | None = None) -> str:
    """Short hash of the scoring weights a measurement was taken under.

    Pitfall #12: tuning weights after seeing results invalidates the result. Recording the
    fingerprint is what lets a later run say "these numbers were taken under different
    weights" instead of quietly comparing two different scorers.
    """
    payload = json.dumps(weights if weights is not None else WEIGHTS, sort_keys=True)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:12]


def _score(protein: dict[str, Any]) -> float | None:
    priority = protein.get("m1_priority")
    if not isinstance(priority, dict):
        return None
    value = priority.get("score")
    return float(value) if isinstance(value, (int, float)) else None


def _length(protein: dict[str, Any]) -> int:
    try:
        return int(protein.get("aa_length") or 0)
    except (TypeError, ValueError):
        return 0


def rank_by_priority(proteins: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    """Order proteins by `m1_priority.score`, ties broken by length, as M1 does.

    Re-derived here rather than trusting `m1_priority.rank` so that a report written by an
    older M1, or one where only some proteins carry a priority record, still produces a
    complete and self-consistent ordering.
    """
    scored = [p for p in proteins if _score(p) is not None]
    if not scored:
        raise ValueError(
            "no protein carries m1_priority.score — run `python -m s2f.m1_genome` first, "
            "or pass --dry-run to measure against the committed fixture"
        )
    return sorted(scored, key=lambda p: (_score(p) or 0.0, _length(p)), reverse=True)


def rank_by_baseline(
    proteins: Sequence[dict[str, Any]], baseline: str, *, seed: int = RANDOM_SEED
) -> list[dict[str, Any]]:
    """Order proteins by a trivial strategy, for comparison."""
    candidates = list(proteins)
    if baseline == "length":
        return sorted(candidates, key=_length, reverse=True)
    if baseline == "keyword":
        return sorted(
            candidates,
            key=lambda p: (bool(PATHOGENESIS_RE.search(p.get("product") or "")), _length(p)),
            reverse=True,
        )
    if baseline == "random":
        shuffled = list(candidates)
        random.Random(seed).shuffle(shuffled)
        return shuffled
    raise ValueError(f"unknown baseline {baseline!r}; known: {', '.join(BASELINES)}")


@dataclass
class Score:
    """Precision, recall and their denominators for one ranking at one cut."""

    k: int
    ranked_total: int
    positives: int
    negatives: int
    unlabelled: int
    excluded: int
    recall_denominator: int
    precision: float | None
    recall: float | None
    labelled_coverage: float

    @classmethod
    def compute(
        cls,
        top_k: Sequence[dict[str, Any]],
        by_feature: dict[str, Match],
        *,
        k: int,
        ranked_total: int,
        positives_in_genome: int,
    ) -> "Score":
        hits = [by_feature[fid] for p in top_k if (fid := p.get("feature_id")) in by_feature]
        positives = sum(1 for m in hits if m.label == "positive")
        negatives = sum(1 for m in hits if m.label == "negative")
        excluded = sum(1 for m in hits if m.label == "excluded")
        labelled = positives + negatives
        found = len({m.symbol.lower() for m in hits if m.label == "positive"})
        return cls(
            k=k,
            ranked_total=ranked_total,
            positives=positives,
            negatives=negatives,
            unlabelled=len(top_k) - positives - negatives - excluded,
            excluded=excluded,
            recall_denominator=positives_in_genome,
            precision=(positives / labelled) if labelled else None,
            recall=(found / positives_in_genome) if positives_in_genome else None,
            labelled_coverage=(labelled / len(top_k)) if top_k else 0.0,
        )


@dataclass
class TieReport:
    """What the score ceiling looks like, and whether the cut lands inside it."""

    max_score: float
    at_max: int
    cut_score: float | None
    at_cut_score: int
    cut_inside_tie_block: bool
    distribution: dict[str, int]


@dataclass
class LeakageAudit:
    """Where the top-K evidence came from: a curated database, or the product string."""

    labelled_in_top_k: int
    curated_points: int
    text_points: int
    text_only_positives: list[str] = field(default_factory=list)
    curated_backed_positives: list[str] = field(default_factory=list)
    unclassified_reasons: list[str] = field(default_factory=list)

    @property
    def text_share(self) -> float | None:
        total = self.curated_points + self.text_points
        return (self.text_points / total) if total else None


@dataclass
class Evaluation:
    """Everything one measurement produced, ready to serialise."""

    truthset: str
    truthset_path: str
    truthset_counts: dict[str, int]
    weights: dict[str, int]
    weights_fingerprint: str
    proteins_scored: int
    matched: int
    match_tiers: dict[str, int]
    not_present: list[str]
    positives_in_genome: int
    negatives_in_genome: int
    measured: Score
    baselines: dict[str, Score]
    ties: TieReport
    leakage: LeakageAudit
    by_class: dict[str, dict[str, int]]
    rows: list[dict[str, Any]]

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["leakage"]["text_share"] = self.leakage.text_share
        return payload


def _classify_reason(reason: str) -> tuple[str, str] | None:
    for prefix, component in CURATED_REASONS.items():
        if reason.startswith(prefix):
            return "curated", component
    for prefix, component in TEXT_REASONS.items():
        if reason.startswith(prefix):
            return "text", component
    return None


def leakage_audit(
    top_k: Sequence[dict[str, Any]], by_feature: dict[str, Match]
) -> LeakageAudit:
    """Split top-K score points into curated-database and product-text evidence."""
    audit = LeakageAudit(labelled_in_top_k=0, curated_points=0, text_points=0)
    for protein in top_k:
        match = by_feature.get(protein.get("feature_id") or "")
        if match is None or match.label not in ("positive", "negative"):
            continue
        audit.labelled_in_top_k += 1
        priority = protein.get("m1_priority") or {}
        curated = text = 0
        for entry in priority.get("breakdown") or []:
            reason = str(entry.get("reason") or "")
            points = abs(int(entry.get("points") or 0))
            kind = _classify_reason(reason)
            if kind is None:
                if reason and reason not in audit.unclassified_reasons:
                    audit.unclassified_reasons.append(reason)
                continue
            if kind[0] == "curated":
                curated += points
            else:
                text += points
        audit.curated_points += curated
        audit.text_points += text
        if match.label == "positive":
            if curated:
                audit.curated_backed_positives.append(match.symbol)
            else:
                audit.text_only_positives.append(match.symbol)
    return audit


def tie_report(ranked: Sequence[dict[str, Any]], k: int) -> TieReport:
    """The score distribution, and whether the cut at K splits a block of equal scores."""
    scores = [_score(p) or 0.0 for p in ranked]
    distribution: dict[str, int] = {}
    for value in scores:
        key = f"{value:g}"
        distribution[key] = distribution.get(key, 0) + 1
    max_score = max(scores) if scores else 0.0
    cut_score = scores[k - 1] if 0 < k <= len(scores) else None
    at_cut = sum(1 for value in scores if cut_score is not None and value == cut_score)
    inside = bool(
        cut_score is not None and k < len(scores) and scores[k] == cut_score
    )
    return TieReport(
        max_score=max_score,
        at_max=sum(1 for value in scores if value == max_score),
        cut_score=cut_score,
        at_cut_score=at_cut,
        cut_inside_tie_block=inside,
        distribution=dict(sorted(distribution.items(), key=lambda kv: float(kv[0]), reverse=True)),
    )


def evaluate(
    proteins: Sequence[dict[str, Any]],
    truthset: TruthSet,
    matches: Iterable[Match],
    not_present: Iterable[Any],
    *,
    k: int = 50,
    baselines: Sequence[str] = BASELINES,
    seed: int = RANDOM_SEED,
    ranker: Callable[[Sequence[dict[str, Any]]], list[dict[str, Any]]] = rank_by_priority,
) -> Evaluation:
    """Measure one ranking against one truth set at one cut."""
    by_feature: dict[str, Match] = {m.feature_id: m for m in matches}
    ranked = ranker(proteins)
    k = max(1, min(k, len(ranked)))
    top_k = ranked[:k]

    positives_in_genome = len(
        {m.symbol.lower() for m in by_feature.values() if m.label == "positive"}
    )
    negatives_in_genome = len(
        {m.symbol.lower() for m in by_feature.values() if m.label == "negative"}
    )

    measured = Score.compute(
        top_k, by_feature, k=k, ranked_total=len(ranked),
        positives_in_genome=positives_in_genome,
    )
    baseline_scores = {}
    for name in baselines:
        alt = rank_by_baseline(ranked, name, seed=seed)[:k]
        baseline_scores[name] = Score.compute(
            alt, by_feature, k=k, ranked_total=len(ranked),
            positives_in_genome=positives_in_genome,
        )

    by_class: dict[str, dict[str, int]] = {}
    in_top_k = {p.get("feature_id") for p in top_k}
    for match in by_feature.values():
        bucket = by_class.setdefault(match.klass, {"in_genome": 0, "in_top_k": 0})
        bucket["in_genome"] += 1
        if match.feature_id in in_top_k:
            bucket["in_top_k"] += 1

    tiers: dict[str, int] = {}
    for match in by_feature.values():
        tiers[match.tier] = tiers.get(match.tier, 0) + 1

    position = {p.get("feature_id"): i + 1 for i, p in enumerate(ranked)}
    ranked_by_id = {p.get("feature_id"): p for p in ranked}
    rows = sorted(
        (
            {
                "rank": position.get(match.feature_id),
                "score": _score(ranked_by_id.get(match.feature_id) or {}),
                "feature_id": match.feature_id,
                "symbol": match.symbol,
                "label": match.label,
                "class": match.klass,
                "match_tier": match.tier,
                "in_top_k": match.feature_id in in_top_k,
                "gene": match.gene,
                "product": match.product,
            }
            for match in by_feature.values()
        ),
        key=lambda row: row["rank"] or 10**9,
    )

    return Evaluation(
        truthset=truthset.name,
        truthset_path=truthset.path,
        truthset_counts=truthset.counts(),
        weights=dict(WEIGHTS),
        weights_fingerprint=weights_fingerprint(),
        proteins_scored=len(ranked),
        matched=len(by_feature),
        match_tiers=dict(sorted(tiers.items())),
        not_present=sorted(row.symbol for row in not_present),
        positives_in_genome=positives_in_genome,
        negatives_in_genome=negatives_in_genome,
        measured=measured,
        baselines=baseline_scores,
        ties=tie_report(ranked, k),
        leakage=leakage_audit(top_k, by_feature),
        by_class=dict(sorted(by_class.items())),
        rows=rows,
    )


__all__ = [
    "BASELINES",
    "CURATED_REASONS",
    "Evaluation",
    "LeakageAudit",
    "RANDOM_SEED",
    "Score",
    "TEXT_REASONS",
    "TieReport",
    "evaluate",
    "leakage_audit",
    "rank_by_baseline",
    "rank_by_priority",
    "tie_report",
    "weights_fingerprint",
]

"""Leave-one-out ablation of M2's triage score.

The question the weights cannot answer on their own: which components actually decide
the selection? `pdb_evidence` carries weight 0.40 and `annotation_gap`
0.30, but a weight is an intention. A component only matters if zeroing it moves proteins
across the cut, and whether that movement is good or bad needs a label.

For each component this reports, against the full score as the reference:

- **selection churn** — how many of the top K change, which is the only measure that does
  not need any ground truth at all
- **precision and recall** against the curated host-interaction truth set
- **precision** against M3's outcome labels, where a run has them

The component list comes from `WEIGHTS`, never a literal, so a component added upstream is
swept automatically. `ligandable_homolog` arrived in #60 after this was written and needed
no change here.

And separately, whether the component ever fired. A component that is zero for every
protein cannot change anything, and reading that as "inert" would be wrong: it usually
means the provider behind it never ran. `OPTIONAL_COMPONENTS` is M2's own name for that
distinction and it is carried through here rather than collapsed.

## This changes no weights

`triage_score_from_components` is re-applied with one input zeroed. The production weights
are untouched, nothing is refit, and no result here licenses a weight change on its own:
the numbers come from one genome and one truth set, and pitfall #12's rule that tuning
after seeing results invalidates the result applies to anyone who acts on them. What an
ablation buys is knowing which components are worth arguing about.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Sequence

from .metrics import Score
from .rankings import (
    OPTIONAL_COMPONENTS,
    TRIAGE_COMPONENTS,
    TRIAGE_WEIGHTS,
    Ranking,
    ablated_triage_score,
    build_ranking,
    component_availability,
    triage_components,
)
from .truthset import Match


@dataclass
class ComponentEffect:
    """What removing one component does."""

    component: str
    weight: float
    optional: bool

    #: Population behaviour of the component itself, before any ranking is built.
    nonzero_proteins: int
    providers_available: int | None
    mean_abs_contribution: float

    #: Selection churn against the full score.
    retained_at_k: int
    churn_at_k: int

    #: Against the curated truth set, and against M3's outcomes where present.
    truth: Score | None = None
    truth_delta_precision: float | None = None
    outcome: Score | None = None
    outcome_delta_precision: float | None = None

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["truth"] = asdict(self.truth) if self.truth else None
        payload["outcome"] = asdict(self.outcome) if self.outcome else None
        return payload


@dataclass
class Ablation:
    """The full leave-one-out sweep."""

    score: str
    k: int
    weights: dict[str, float]
    reference_truth: Score | None
    reference_outcome: Score | None
    effects: list[ComponentEffect] = field(default_factory=list)
    never_fired: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "score": self.score,
            "k": self.k,
            "weights": self.weights,
            "reference_truth": asdict(self.reference_truth) if self.reference_truth else None,
            "reference_outcome": (
                asdict(self.reference_outcome) if self.reference_outcome else None
            ),
            "effects": [effect.to_dict() for effect in self.effects],
            "never_fired": self.never_fired,
            "notes": self.notes,
        }

    def load_bearing(self) -> list[ComponentEffect]:
        """Components whose removal moves at least one protein across the cut."""
        return [effect for effect in self.effects if effect.churn_at_k > 0]


def _score_against(
    ranked: Ranking,
    proteins_by_id: dict[str, dict[str, Any]],
    by_feature: dict[str, Match] | None,
    *,
    k: int,
) -> Score | None:
    if not by_feature:
        return None
    top = [proteins_by_id[fid] for fid in ranked.top(k) if fid in proteins_by_id]
    positives = len({
        m.symbol.lower() for m in by_feature.values() if m.label == "positive"
    })
    return Score.compute(
        top, by_feature, k=k, ranked_total=len(ranked.order),
        positives_in_genome=positives,
    )


def ablate(
    proteins: Sequence[dict[str, Any]],
    *,
    truth_matches: dict[str, Match] | None = None,
    outcome_matches: dict[str, Match] | None = None,
    k: int = 50,
) -> Ablation:
    """Leave-one-out over every triage component."""
    scored = [p for p in proteins if triage_components(p) is not None]
    if not scored:
        raise ValueError(
            "no protein carries triage.components — run `python -m s2f.m2_triage` first"
        )
    proteins_by_id = {str(p["feature_id"]): p for p in scored if p.get("feature_id")}

    reference = build_ranking(scored, "triage")
    reference_top = set(reference.top(k))
    reference_truth = _score_against(reference, proteins_by_id, truth_matches, k=k)
    reference_outcome = _score_against(reference, proteins_by_id, outcome_matches, k=k)

    effects: list[ComponentEffect] = []
    never_fired: list[str] = []

    for component in TRIAGE_COMPONENTS:
        values = [
            (triage_components(p) or {}).get(component, 0.0) for p in scored
        ]
        nonzero = sum(1 for value in values if value)
        mean_abs = (sum(abs(value) for value in values) / len(values)) if values else 0.0
        optional = component in OPTIONAL_COMPONENTS
        available = (
            sum(1 for p in scored if component_availability(p).get(component) is True)
            if optional else None
        )
        if nonzero == 0:
            never_fired.append(component)

        ablated = build_ranking(
            scored, f"triage-{component}",
            score_of=lambda p, c=component: ablated_triage_score(p, c),
        )
        ablated_top = set(ablated.top(k))
        truth = _score_against(ablated, proteins_by_id, truth_matches, k=k)
        outcome = _score_against(ablated, proteins_by_id, outcome_matches, k=k)

        effects.append(
            ComponentEffect(
                component=component,
                weight=TRIAGE_WEIGHTS[component],
                optional=optional,
                nonzero_proteins=nonzero,
                providers_available=available,
                mean_abs_contribution=round(mean_abs * abs(TRIAGE_WEIGHTS[component]), 6),
                retained_at_k=len(reference_top & ablated_top),
                churn_at_k=len(reference_top - ablated_top),
                truth=truth,
                truth_delta_precision=_delta(reference_truth, truth),
                outcome=outcome,
                outcome_delta_precision=_delta(reference_outcome, outcome),
            )
        )

    notes: list[str] = []
    if never_fired:
        notes.append(
            "never fired on this run, so their ablation is uninformative rather than a "
            "finding that they do not matter: " + ", ".join(never_fired)
        )
    inert = [e.component for e in effects if e.churn_at_k == 0 and e.nonzero_proteins]
    if inert:
        notes.append(
            "fired but changed nothing in the top "
            f"{k}: " + ", ".join(inert)
        )
    return Ablation(
        score="triage",
        k=k,
        weights=dict(TRIAGE_WEIGHTS),
        reference_truth=reference_truth,
        reference_outcome=reference_outcome,
        effects=sorted(effects, key=lambda e: e.churn_at_k, reverse=True),
        never_fired=never_fired,
        notes=notes,
    )


def _delta(reference: Score | None, ablated: Score | None) -> float | None:
    """Ablated minus reference precision. Negative means the component was helping."""
    if reference is None or ablated is None:
        return None
    if reference.precision is None or ablated.precision is None:
        return None
    return round(ablated.precision - reference.precision, 6)


__all__ = ["Ablation", "ComponentEffect", "ablate"]

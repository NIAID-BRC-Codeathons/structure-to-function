"""The rankings under test, and the counterfactual re-scores an ablation needs.

Two scores in this pipeline order proteins, and `00a-data-contract.md` is explicit that
they answer different questions and are expected to disagree:

``m1_priority``
    host-interaction evidence in the annotation, before any structure has been looked at.
``triage``
    structural tractability — can we fold and dock this — with PDB evidence at weight 0.40.

Both are read here through one accessor so every measurement in this package can be
pointed at either. Nothing in this module re-runs a search or touches the network.

## Counterfactual re-scoring

`m2_triage.score.triage_score_from_components` is a pure function of the component
values M2 writes out per protein, precisely so "the selection is reproducible from the recorded
numbers alone". That makes leave-one-out ablation free: zero a component, re-apply the
same function, re-rank, and compare. No refit, no network, no second run.

This is measurement of counterfactual rankings, not tuning. Pitfall #12 forbids changing
weights after seeing results; it does not forbid asking what a different weight vector
*would* have selected, and answering that is how you find out which components are load
bearing before anyone proposes a change.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any, Callable, Iterable, Sequence

from ..m2_triage.score import WEIGHTS as TRIAGE_WEIGHTS
from ..m2_triage.score import OPTIONAL_COMPONENTS, triage_score_from_components

#: Scores this package can measure.
SCORES = ("m1_priority", "triage")

#: Triage components, in the order M2 declares them.
TRIAGE_COMPONENTS = tuple(TRIAGE_WEIGHTS)

#: Which triage components are evidence a curated database supplied, and which the
#: product string supplied. `annotation_gap` rewards a protein for reading "hypothetical",
#: which is the one component computed purely from annotation text. `surface_bonus` and
#: `membrane_penalty` are `OPTIONAL_COMPONENTS`: a real provider (SignalP, DeepTMHMM,
#: PSORTb) when one ran, a product-text heuristic otherwise, so they are attributed per
#: protein from `components_available` rather than assumed either way.
TRIAGE_CURATED = ("pdb_evidence", "virulence_amr", "essential", "drug_target",
                  "human_homolog_penalty")
TRIAGE_TEXT = ("annotation_gap",)


def _as_float(value: Any) -> float | None:
    return float(value) if isinstance(value, (int, float)) and not isinstance(value, bool) else None


def m1_priority_score(protein: dict[str, Any]) -> float | None:
    record = protein.get("m1_priority")
    return _as_float(record.get("score")) if isinstance(record, dict) else None


def triage_score(protein: dict[str, Any]) -> float | None:
    record = protein.get("triage")
    return _as_float(record.get("score")) if isinstance(record, dict) else None


def score_accessor(name: str) -> Callable[[dict[str, Any]], float | None]:
    """The score function for a named ranking."""
    if name == "m1_priority":
        return m1_priority_score
    if name == "triage":
        return triage_score
    raise ValueError(f"unknown score {name!r}; known: {', '.join(SCORES)}")


def triage_components(protein: dict[str, Any]) -> dict[str, float] | None:
    """The eight recorded component values for one protein, or None if M2 never scored it."""
    record = protein.get("triage")
    if not isinstance(record, dict):
        return None
    components = record.get("components")
    if not isinstance(components, dict):
        return None
    return {name: _as_float(components.get(name)) or 0.0 for name in TRIAGE_COMPONENTS}


def component_availability(protein: dict[str, Any]) -> dict[str, bool]:
    """Which optional components had a provider behind them for this protein.

    M2 records this because a component whose provider never ran contributes 0 — the same
    number as "we checked and it is false". An ablation that ignored the distinction would
    report `membrane_penalty` as inert on a run where DeepTMHMM simply was not installed.
    """
    record = protein.get("triage")
    available = record.get("components_available") if isinstance(record, dict) else None
    if not isinstance(available, dict):
        return {}
    return {name: bool(value) for name, value in available.items()}


def ablated_triage_score(
    protein: dict[str, Any], drop: str | Iterable[str]
) -> float | None:
    """Recompute `triage.score` with one or more components zeroed.

    Uses M2's own `triage_score_from_components`, so an ablated score is produced by the
    same arithmetic as the real one and cannot drift from it.
    """
    components = triage_components(protein)
    if components is None:
        return None
    dropped = {drop} if isinstance(drop, str) else set(drop)
    unknown = dropped - set(TRIAGE_COMPONENTS)
    if unknown:
        raise ValueError(
            f"unknown triage component(s): {', '.join(sorted(unknown))}; "
            f"known: {', '.join(TRIAGE_COMPONENTS)}"
        )
    kept = {name: (0.0 if name in dropped else value) for name, value in components.items()}
    return triage_score_from_components(_Components(kept))


class _Components:
    """Adapter so a plain dict can be handed to M2's scoring function."""

    __slots__ = ("_values",)

    def __init__(self, values: dict[str, float]) -> None:
        self._values = values

    def as_dict(self) -> dict[str, float]:
        return self._values


def _length(protein: dict[str, Any]) -> int:
    try:
        return int(protein.get("aa_length") or 0)
    except (TypeError, ValueError):
        return 0


@dataclass
class Ranking:
    """One ordered list of proteins, and the score that ordered it."""

    name: str
    order: list[str] = field(default_factory=list)
    scores: dict[str, float] = field(default_factory=dict)
    scored: list[dict[str, Any]] = field(default_factory=list)

    def top(self, k: int) -> list[str]:
        return self.order[: max(0, k)]

    def position(self) -> dict[str, int]:
        """`feature_id -> 1-based rank`."""
        return {feature_id: i + 1 for i, feature_id in enumerate(self.order)}


def build_ranking(
    proteins: Sequence[dict[str, Any]],
    name: str,
    score_of: Callable[[dict[str, Any]], float | None] | None = None,
) -> Ranking:
    """Order proteins by a score, ties broken by length descending.

    The tie-break matches what M1 and M2 both do. It is not evidence, and where it decides
    membership of a cut `metrics.tie_report` says so.
    """
    accessor = score_of or score_accessor(name)
    scored = [(p, accessor(p)) for p in proteins]
    keep = [(p, s) for p, s in scored if s is not None and p.get("feature_id")]
    if not keep:
        raise ValueError(
            f"no protein carries a {name!r} score — run the module that writes it first, "
            "or pass --dry-run to measure against the committed fixture"
        )
    keep.sort(key=lambda item: (item[1], _length(item[0])), reverse=True)
    return Ranking(
        name=name,
        order=[p["feature_id"] for p, _ in keep],
        scores={p["feature_id"]: s for p, s in keep},
        scored=[p for p, _ in keep],
    )


#: A product string that says the annotation does not know what this protein is. Used only
#: to describe how the *population* moves through the pipeline, never to score anything —
#: M2's own `annotation_gap` component is the scoring definition and this must not be
#: confused with it. Kept deliberately broad, because the question it answers is "how many
#: of the proteins this project exists to characterise survive each stage".
GAP_LIKE = re.compile(r"hypothetical|uncharacteri[sz]ed|putative|\bDUF\d", re.I)


def is_gap_like(protein: dict[str, Any]) -> bool:
    """Does the annotation admit it does not know what this protein is?"""
    return bool(GAP_LIKE.search(protein.get("product") or ""))


def population_stages(
    proteins: Sequence[dict[str, Any]]
) -> list[tuple[str, int, int]]:
    """`(stage, population, gap_like)` as proteins pass through M1's cap and M2's selection.

    The seam this exposes: M1's `--seq-cap` truncates by `m1_priority`, which scores a
    hypothetical protein 0 because it has no gene name, no mechanism keyword and no
    specialty hit. M2's `annotation_gap` then pays 0.30 for exactly that property. The two
    policies are each locally sensible and together they cancel, which is only visible
    from outside both modules.
    """
    scored = [p for p in proteins if p.get("triage")]
    selected = [p for p in scored if (p.get("triage") or {}).get("selected") is True]
    stages = [("in genome", proteins), ("survived M1 seq-cap", scored)]
    if selected:
        stages.append(("selected by triage", selected))
    return [
        (label, len(pop), sum(1 for p in pop if is_gap_like(p)))
        for label, pop in stages
        if pop
    ]


def available_scores(proteins: Sequence[dict[str, Any]]) -> list[str]:
    """Which of `SCORES` this report actually carries."""
    return [
        name for name in SCORES
        if any(score_accessor(name)(p) is not None for p in proteins)
    ]


__all__ = [
    "OPTIONAL_COMPONENTS",
    "Ranking",
    "SCORES",
    "TRIAGE_COMPONENTS",
    "TRIAGE_CURATED",
    "TRIAGE_TEXT",
    "TRIAGE_WEIGHTS",
    "ablated_triage_score",
    "GAP_LIKE",
    "available_scores",
    "build_ranking",
    "is_gap_like",
    "population_stages",
    "component_availability",
    "m1_priority_score",
    "score_accessor",
    "triage_components",
    "triage_score",
]

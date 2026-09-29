"""Objective outcome labels, taken from what the pipeline actually produced.

The truth set in `truthset.py` is curated: someone read the literature and made a call.
This module needs no curation at all. M2's `triage` score is a **prediction** — "this
protein is worth folding and docking" — and M3's quality gate is the **outcome**: a
structure was collected and passed, or it did not. Scoring the prediction against the
outcome asks whether the pipeline's central selection step allocates the downstream budget
well, and nobody has to agree about biology first.

Four states, from `structures[]`:

===================  ==========================================================
`dockable`           collected and `usable_for_docking is True` — the positive
`unusable`           collected but the quality gate rejected it
`no_structure`       `prediction_required`: nothing existing to reuse
`failed`             the collection itself errored
===================  ==========================================================

`failed` is **excluded**, not counted as a negative. Chen's workflow contract is explicit
that a technical failure is never recorded as a negative finding — a download that timed
out is not evidence that a protein has no usable structure — and the same rule that keeps
a failed sequence search from becoming "no homolog" applies here. It is the same role the
`excluded` label plays in the curated set, arrived at from the other direction.

`no_structure` **is** a negative, and deliberately so. M2 selected the protein on a claim
of structural tractability; M3 finding nothing to reuse is that claim failing. It is
recorded as its own class rather than merged with `unusable`, because the two suggest
different fixes: more sources against a prediction step.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable, Sequence

from .truthset import Match

#: Outcome class -> the label the metrics machinery scores it as.
OUTCOME_LABELS = {
    "dockable": "positive",
    "unusable": "negative",
    "no_structure": "negative",
    "failed": "excluded",
}

TIER = "m3_outcome"


@dataclass
class OutcomeSet:
    """Outcome labels for the proteins M3 reached, plus why the rest are missing."""

    matches: list[Match]
    counts: dict[str, int]
    policy_version: str | None
    covered: int
    #: Proteins M2 selected that M3 never produced a record for. Not a negative: M3 may
    #: have been run with `--limit`, or not run at all after the last M2 pass.
    selected_but_uncollected: list[str]

    def by_feature(self) -> dict[str, Match]:
        return {match.feature_id: match for match in self.matches}

    @property
    def base_rate(self) -> float | None:
        """Share of the labelled population that is `dockable`.

        The number a precision@K has to be read against. On the first HS11286 run M3
        passed 298 of 300, so precision@50 came out at 100% against a base rate of 99.3%:
        a lift of half a point, which is not evidence that the ranking is selecting well.
        It is evidence that the label cannot tell good selections from bad ones, because
        M2 selects largely on PDB evidence and M3 then succeeds because PDB evidence is
        exactly what makes a structure retrievable. Reporting precision without the base
        rate beside it would turn that into a false headline.
        """
        labelled = self.counts.get("dockable", 0) + self.counts.get("unusable", 0) \
            + self.counts.get("no_structure", 0)
        return (self.counts.get("dockable", 0) / labelled) if labelled else None

    def lift(self, precision: float | None) -> float | None:
        """Precision minus base rate. Near zero means the label is not discriminating."""
        base = self.base_rate
        if base is None or precision is None:
            return None
        return round(precision - base, 6)

    def discriminates(self, *, floor: float = 0.05) -> bool:
        """Is there enough variance in the label for a precision figure to mean anything?

        A label that is 99% one class cannot separate a good ranking from a bad one. The
        threshold is a judgement call, stated here rather than left implicit.
        """
        base = self.base_rate
        return base is not None and floor <= base <= (1 - floor)


def classify(record: dict[str, Any]) -> str:
    """The outcome class of one `structures[]` record."""
    status = str(record.get("collection_status") or "")
    if status == "failed":
        return "failed"
    if status == "prediction_required":
        return "no_structure"
    usable = record.get("usable_for_docking")
    if usable is True:
        return "dockable"
    if usable is False:
        return "unusable"
    # Collected, but the gate never ran — M3 records `None` for that. Treat it the way a
    # missing provider is treated elsewhere: unmeasured, not failed.
    return "failed"


def outcomes_from_report(
    report: dict[str, Any], proteins: Sequence[dict[str, Any]]
) -> OutcomeSet:
    """Build outcome labels from `structures[]`, joined to `proteins[]` on feature_id."""
    structures = report.get("structures") or []
    if not isinstance(structures, list):
        raise ValueError("report structures section is not a list")

    by_feature = {
        str(p.get("feature_id")): p for p in proteins if p.get("feature_id")
    }
    matches: list[Match] = []
    counts: dict[str, int] = {name: 0 for name in OUTCOME_LABELS}
    seen: set[str] = set()
    policy_version: str | None = None

    for record in structures:
        if not isinstance(record, dict):
            continue
        feature_id = str(record.get("feature_id") or "")
        if not feature_id or feature_id in seen:
            continue
        seen.add(feature_id)
        quality = record.get("quality")
        if isinstance(quality, dict) and quality.get("policy_version"):
            policy_version = str(quality["policy_version"])
        klass = classify(record)
        counts[klass] += 1
        protein = by_feature.get(feature_id, {})
        matches.append(
            Match(
                feature_id=feature_id,
                symbol=str(record.get("accession") or record.get("source") or feature_id),
                label=OUTCOME_LABELS[klass],
                klass=klass,
                tier=TIER,
                gene=protein.get("gene"),
                product=protein.get("product"),
            )
        )

    selected = [
        str(p.get("feature_id"))
        for p in proteins
        if (p.get("triage") or {}).get("selected") is True and p.get("feature_id")
    ]
    return OutcomeSet(
        matches=matches,
        counts=counts,
        policy_version=policy_version,
        covered=len(matches),
        selected_but_uncollected=sorted(set(selected) - seen),
    )


def coverage_warning(outcomes: OutcomeSet, k: int, ranked_top_k: Iterable[str]) -> str | None:
    """Why an outcome-based figure may not mean what it looks like.

    M3 only collects for proteins M2 selected, so outcome labels exist for the selection
    and nowhere else. Measuring precision@K against them when K reaches past the selected
    set silently changes the denominator, and the run has to say so.
    """
    labelled = outcomes.by_feature()
    missing = [fid for fid in ranked_top_k if fid not in labelled]
    if not missing:
        return None
    return (
        f"{len(missing)} of the top {k} have no M3 outcome record, because M3 collects only "
        "for proteins M2 selected. The outcome figures describe the selected set, not the "
        "whole ranking."
    )


__all__ = [
    "OUTCOME_LABELS",
    "TIER",
    "OutcomeSet",
    "classify",
    "coverage_warning",
    "outcomes_from_report",
]

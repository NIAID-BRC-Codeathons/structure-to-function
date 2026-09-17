"""Interim triage scoring for M2 (issue #12).

Weights live in ``WEIGHTS`` and are mirrored in ``docs/02a-m2-pdb-evidence.md``. Pitfall #12:
change them in both places in one commit and add a dated row to that doc's change log, saying
whether the change came before or after seeing a run's results.

Every component is written out per protein, so the selection is reproducible from the recorded
numbers alone (``triage_score_from_components``).
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field

from .bvbrc_input import Protein
from .pdb_evidence import SequenceHit

# Scoring gates for a hit to count as evidence at all.
MIN_EVALUE = 1e-5
MIN_QUERY_COVERAGE = 0.5

#: Bump with any weight change, alongside a dated row in docs/02a-m2-pdb-evidence.md.
WEIGHTS_VERSION = "2026-09-16"

WEIGHTS = {
    "pdb_evidence": 0.40,
    "virulence_amr": 0.20,
    "essential": 0.15,
    "drug_target": 0.15,
    "annotation_gap": 0.30,
    "human_homolog_penalty": -0.25,
}

LIGAND_BONUS = 0.10
PARTNER_BONUS = 0.05
ANNOTATION_BONUS = 0.05
# Best possible raw hit score: identity 1.0 x coverage 1.0 x quality 1.0, plus every bonus.
# pdb_evidence is normalized by this instead of clipped at 1.0, which would collapse the
# distinction between the strongest hits (see the change log in docs/02a-m2-pdb-evidence.md).
MAX_HIT_SCORE = 1.0 + LIGAND_BONUS + PARTNER_BONUS + ANNOTATION_BONUS

HUMAN_TAXONOMY_ID = 9606

#: Identity at or above which a human homolog is a selectivity risk worth reporting (pitfall #4).
#: The penalty itself is continuous; this only drives the flag M6 states and M4 reads.
CLOSE_HUMAN_HOMOLOG_IDENTITY = 40.0


def hit_qualifies(hit: SequenceHit) -> bool:
    """A hit counts only if it is significant and covers half the query."""
    return hit.evalue <= MIN_EVALUE and hit.query_coverage >= MIN_QUERY_COVERAGE


def structure_quality(hit: SequenceHit) -> float:
    """Quality factor from experimental method and resolution."""
    method = (hit.experimental_method or "").upper()
    if "NMR" in method:
        return 0.7
    if "X-RAY" in method or "ELECTRON MICROSCOPY" in method or method in {"EM", "X-RAY"}:
        if hit.resolution is None:
            return 0.6
        if hit.resolution <= 2.5:
            return 1.0
        if hit.resolution <= 3.5:
            return 0.8
        return 0.6
    return 0.6


def hit_score(hit: SequenceHit) -> float:
    """seq x quality + bonuses, per docs/02a-m2-pdb-evidence.md."""
    if not hit_qualifies(hit):
        return 0.0
    base = hit.identity * hit.query_coverage * structure_quality(hit)
    bonus = 0.0
    if hit.has_ligand:
        bonus += LIGAND_BONUS
    if hit.has_partner:
        bonus += PARTNER_BONUS
    if hit.has_ec or hit.has_go:
        bonus += ANNOTATION_BONUS
    return base + bonus


def best_hit(hits: list[SequenceHit]) -> SequenceHit | None:
    """Best qualifying hit: score, then distinct-entry support is handled by the caller."""
    qualifying = [hit for hit in hits if hit_qualifies(hit)]
    if not qualifying:
        return None
    return max(qualifying, key=lambda hit: (hit_score(hit), hit.bitscore))


@dataclass
class TriageComponents:
    """The scored components, written to proteins.tsv so scoring is reproducible."""

    pdb_evidence: float = 0.0
    virulence_amr: float = 0.0
    essential: float = 0.0
    drug_target: float = 0.0
    annotation_gap: float = 0.0
    human_homolog_penalty: float = 0.0

    def as_dict(self) -> dict[str, float]:
        return asdict(self)


@dataclass
class ProteinScore:
    """One protein's triage outcome plus the flags M3/M4 need."""

    feature_id: str
    product: str
    components: TriageComponents
    score: float
    retrieval_status: str
    best_hit: SequenceHit | None = None
    distinct_entries: int = 0
    qualifying_hits: int = 0
    flags: dict[str, object] = field(default_factory=dict)
    rank: int = 0
    selected: bool = False
    reason: str = ""
    error: str = ""


def triage_score_from_components(components: TriageComponents) -> float:
    """The single definition of the score; used by the pipeline and by the tests."""
    values = components.as_dict()
    # human_homolog_penalty carries a negative weight, so one loop covers every component.
    return round(sum(weight * values[name] for name, weight in WEIGHTS.items()), 6)


def score_protein(
    protein: Protein,
    hits: list[SequenceHit],
    retrieval_status: str,
    error: str = "",
    *,
    human_identity: float | None = None,
    human_homolog_source: str = "",
) -> ProteinScore:
    """Score one protein from its PDB hits and M1 specialty rows.

    ``human_identity`` overrides the specialty table with evidence we computed ourselves
    (issue #29). BV-BRC precomputes `Human Homolog` rows for public genomes only, so on a blinded
    genome the specialty value is absent and the penalty would silently never fire.
    """
    qualifying = [hit for hit in hits if hit_qualifies(hit)]
    top = best_hit(hits)
    pdb_evidence = (hit_score(top) / MAX_HIT_SCORE) if top is not None else 0.0

    if human_identity is None:
        human_identity = protein.human_homolog_identity
        human_homolog_source = human_homolog_source or "bvbrc_specialty"
    components = TriageComponents(
        pdb_evidence=round(pdb_evidence, 6),
        virulence_amr=1.0 if protein.has_virulence_or_amr else 0.0,
        essential=1.0 if protein.is_essential_ortholog else 0.0,
        drug_target=1.0 if protein.is_drug_target else 0.0,
        annotation_gap=1.0 if protein.is_uncharacterized else 0.0,
        human_homolog_penalty=round(human_identity / 100.0, 6) if human_identity else 0.0,
    )

    flags = {
        "has_ligand_in_entry": bool(top and top.has_ligand),
        "has_partner_in_entry": bool(top and top.has_partner),
        "holo_homolog": bool(top and top.has_ligand),
        "metals_in_entry": ";".join(top.metals) if top and top.metals else "",
        "transporter": protein.is_transporter,
        "metal_resistance": protein.has_metal_resistance,
        "human_homolog_identity": human_identity if human_identity is not None else "",
        "close_human_homolog": bool(
            human_identity is not None and human_identity >= CLOSE_HUMAN_HOMOLOG_IDENTITY
        ),
        "human_homolog_source": human_homolog_source or "",
        "no_pdb_hit": not qualifying,
        "pdb_hit_organism": top.organism if top else "",
        "human_pdb_hit": bool(top and top.taxonomy_id == HUMAN_TAXONOMY_ID),
        "uniprot_of_best_hit": ";".join(top.uniprot_ids) if top else "",
    }

    return ProteinScore(
        feature_id=protein.feature_id,
        product=protein.product,
        components=components,
        score=triage_score_from_components(components),
        retrieval_status=retrieval_status,
        best_hit=top,
        distinct_entries=len({hit.entry_id for hit in qualifying}),
        qualifying_hits=len(qualifying),
        flags=flags,
        error=error,
    )


def _reason(score: ProteinScore) -> str:
    """Short, human-readable justification recorded for every protein, kept or discarded."""
    parts: list[str] = []
    top = score.best_hit
    if top is not None:
        parts.append(
            f"PDB {top.entity_id} identity {top.identity:.0%}, coverage {top.query_coverage:.0%}"
            + (f", {top.resolution:.2f} A" if top.resolution is not None else "")
        )
        if top.ligands:
            parts.append("ligand-bound homolog: " + ",".join(top.ligands[:3]))
    elif score.retrieval_status == "query-failed":
        parts.append("PDB search failed; no structural evidence recorded")
    else:
        parts.append("no qualifying PDB hit")

    if score.components.virulence_amr:
        parts.append("virulence/AMR")
    if score.components.essential:
        parts.append("essential-gene ortholog")
    if score.components.drug_target:
        parts.append("known drug target")
    if score.components.annotation_gap:
        parts.append("uncharacterized product")
    if score.components.human_homolog_penalty:
        parts.append(f"human homolog {score.components.human_homolog_penalty:.0%} identity (penalty)")
    return "; ".join(parts)


def rank_and_select(scores: list[ProteinScore], *, top_n: int = 50) -> list[ProteinScore]:
    """Rank all proteins, mark the top N selected, and give every protein a reason."""
    ordered = sorted(
        scores,
        key=lambda s: (
            -s.score,
            -s.components.pdb_evidence,
            -s.distinct_entries,
            s.feature_id,
        ),
    )
    for index, score in enumerate(ordered, start=1):
        score.rank = index
        score.selected = index <= top_n
        score.reason = _reason(score)
    return ordered

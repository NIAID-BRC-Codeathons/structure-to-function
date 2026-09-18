"""Interim triage scoring for M2 (issue #12).

Weights live in ``WEIGHTS`` and are mirrored in ``docs/02a-m2-pdb-evidence.md``. Pitfall #12:
change them in both places in one commit and add a dated row to that doc's change log, saying
whether the change came before or after seeing a run's results.

Every component is written out per protein, so the selection is reproducible from the recorded
numbers alone (``triage_score_from_components``).
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Sequence

from .bvbrc_input import Protein, same_genus, same_species
from .pdb_evidence import SequenceHit

# Scoring gates for a hit to count as evidence at all.
MIN_EVALUE = 1e-5
MIN_QUERY_COVERAGE = 0.5

#: Bump with any weight change, alongside a dated row in docs/02a-m2-pdb-evidence.md.
WEIGHTS_VERSION = "2026-09-18"

WEIGHTS = {
    "pdb_evidence": 0.40,
    "virulence_amr": 0.20,
    "essential": 0.15,
    "drug_target": 0.15,
    "annotation_gap": 0.30,
    "surface_bonus": 0.10,
    "ligandable_homolog": 0.10,
    "membrane_penalty": -0.15,
    "human_homolog_penalty": -0.25,
}

#: Components that need a provider that may not have run. A component whose provider was absent
#: contributes 0 — the same number as "we checked and it is false" — so the two are recorded
#: separately per protein and summarised per run. Otherwise a missing DeepTMHMM run reads as
#: "no membrane proteins in this genome" (issue #10's flag_row makes the same distinction).
OPTIONAL_COMPONENTS = ("surface_bonus", "membrane_penalty", "ligandable_homolog")

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

#: A human homolog at or above this identity disqualifies a protein from selection outright.
#: Selectivity is a requirement, not a preference: a compound that hits our target and its human
#: counterpart is a toxicity problem, not a slightly weaker candidate. The continuous
#: ``human_homolog_penalty`` still ranks, but it can be outvoted — in the 2026-09-18 lambda0 run
#: ``essential`` (+0.15) beat a 50.6%-identity penalty (-0.127) and returned ribosomal protein
#: S12p to the top 50, whose human mitochondrial counterpart MRPS12 is an established
#: aminoglycoside liability. Coverage is gated upstream by ``--human-min-coverage`` (default
#: 0.5), so an identity reaching this check has already passed a coverage floor.
DISQUALIFYING_HUMAN_HOMOLOG_IDENTITY = 40.0

#: Term kinds that amount to a functional assignment, versus a domain match alone.
FUNCTION_TERM_KINDS = ("ec", "ko", "cog", "go")
DOMAIN_TERM_KINDS = ("pfam", "interpro", "signature")
PARTIAL_ANNOTATION_GAP = 0.5


def annotation_gap_value(protein: Protein, annotation: Any = None) -> float:
    """How much of this protein's function is unknown, judged on evidence rather than wording.

    ``is_uncharacterized`` reads the BV-BRC product string, which is written before any
    annotation provider runs. Scoring the gap from that alone ignores the annotation layer
    entirely: in the 2026-09-18 lambda0 run, 91 of the 178 proteins collecting a full
    ``annotation_gap`` had InterProScan terms, and 17 of the 34 uncharacterized proteins in the
    top 50 did. A protein whose product says "hypothetical" but which InterProScan assigns an EC
    number is not an annotation gap.

    1.0  uncharacterized product, no functional terms at all
    0.5  uncharacterized product with a domain match but no EC/KO/COG/GO assignment
    0.0  a characterized product, or terms that name the function
    """
    if not protein.is_uncharacterized:
        return 0.0
    terms_of = getattr(annotation, "terms_of", None)
    if terms_of is None:
        # No annotation layer ran for this protein, so the product string is all we have.
        return 1.0
    if any(terms_of(kind) for kind in FUNCTION_TERM_KINDS):
        return 0.0
    if any(terms_of(kind) for kind in DOMAIN_TERM_KINDS):
        return PARTIAL_ANNOTATION_GAP
    return 1.0


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
    surface_bonus: float = 0.0
    ligandable_homolog: float = 0.0
    membrane_penalty: float = 0.0
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
    #: Which optional components had a provider behind them for this protein. A component that
    #: is False here contributed 0 because nothing measured it, not because it is absent.
    components_available: dict[str, bool] = field(default_factory=dict)
    rank: int = 0
    selected: bool = False
    reason: str = ""
    #: Why this protein may never be selected, whatever it scored. Empty when it is eligible.
    disqualified: str = ""
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
    essential: bool | None = None,
    essential_source: str = "",
    annotation: Any = None,
    query_organism: str = "",
    query_taxonomy: Any = None,
    chembl_targets: Sequence[str] = (),
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
    if essential is None:
        essential = protein.is_essential_ortholog
        essential_source = essential_source or "bvbrc_specialty"

    # Localization and membrane state come from #10. A flag with no source behind it is
    # unknown, not false: `annotate` leaves it empty when no provider covered the protein.
    surface_known = membrane_known = False
    surface_value = membrane_value = 0.0
    surface_source = membrane_source = ""
    if annotation is not None:
        membrane_source = getattr(annotation, "membrane_source", "") or ""
        if membrane_source:
            membrane_known = True
            membrane_value = 1.0 if getattr(annotation, "membrane", False) else 0.0
        surface_source = (
            getattr(annotation, "signal_source", "") or getattr(annotation, "localization_source", "") or ""
        )
        if surface_source:
            surface_known = True
            surface_value = 1.0 if (
                getattr(annotation, "surface_exposed", False) or getattr(annotation, "secreted", False)
            ) else 0.0

    # Ligandability: a small molecule has been observed bound to the homolog, or compounds have
    # been assayed against it. This is evidence the site is tractable, **not** evidence of
    # potency — a ligand in a PDB entry is not a binding measurement, which is why additives,
    # glycans and metals are already excluded from `ligands` upstream. M4 needs a holo template
    # (pitfall #2); M3 prefers one to align to.
    ligandable = 0.0
    ligandable_basis = ""
    ligandable_known = False
    if top is not None:
        ligandable_known = True
        if top.ligands and top.has_partner:
            # Entry-level ligand presence is not evidence that *our* chain is ligandable. A
            # ribosomal protein's best hit is a ribosome, and ribosome entries routinely carry
            # antibiotics that survive the additive/glycan/metal filter, so every subunit would
            # inherit credit for a ligand bound to the rRNA or to another chain. Measured on
            # 2026-09-18: this promoted L33p and L16p to ranks 13-14. Attributing the ligand
            # properly needs per-entity contacts we do not fetch, so withhold rather than guess.
            ligandable_basis = (
                f"ligand in {top.entry_id} not credited: entry has partner entities, so the "
                "ligand cannot be attributed to this chain"
            )
        elif top.ligands:
            ligandable = 1.0
            ligandable_basis = f"ligand bound in {top.entry_id}: {','.join(top.ligands[:3])}"
    if chembl_targets:
        ligandable_known = True
        if not ligandable:
            ligandable = 1.0
            ligandable_basis = "compounds assayed against the mapped target: " + ";".join(
                list(chembl_targets)[:3]
            )

    amr_value, amr_basis = protein.amr_evidence
    virulence_amr = max(1.0 if protein.is_virulence_factor else 0.0, amr_value)
    components = TriageComponents(
        pdb_evidence=round(pdb_evidence, 6),
        virulence_amr=virulence_amr,
        essential=1.0 if essential else 0.0,
        drug_target=1.0 if protein.is_drug_target else 0.0,
        annotation_gap=annotation_gap_value(protein, annotation),
        surface_bonus=surface_value,
        ligandable_homolog=ligandable,
        membrane_penalty=membrane_value,
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
        "essential_source": essential_source or "",
        "amr_basis": amr_basis,
        "antibiotic_target_not_resistance": protein.is_antibiotic_target,
        "ligandable_basis": ligandable_basis,
        "surface_exposed_source": surface_source,
        "membrane_source": membrane_source,
        # "no qualifying hit" must not absorb "the search failed" — a network failure would
        # otherwise be laundered into a structural claim downstream (found by @Ashita2619
        # while building #43 on top of this file).
        "no_pdb_hit": not qualifying and retrieval_status != "query-failed",
        "pdb_search_failed": retrieval_status == "query-failed",
        "pdb_hit_organism": top.organism if top else "",
        "human_pdb_hit": bool(top and top.taxonomy_id == HUMAN_TAXONOMY_ID),
        "uniprot_of_best_hit": ";".join(top.uniprot_ids) if top else "",
        # A hit against our own organism is not a cross-organism transfer: the test genome has
        # its own structures in the PDB, so a 100% identity "discovery" can be a self-match
        # (raised by @cmmann21 on #12). Recorded, never scored differently — but M6 must be able
        # to say why an annotation was easy.
        "same_species_hit": _same_species(query_taxonomy, query_organism, top),
        "same_genus_hit": _same_genus(query_taxonomy, query_organism, top),
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
        components_available={
            "surface_bonus": surface_known,
            "ligandable_homolog": ligandable_known,
            "membrane_penalty": membrane_known,
        },
        error=error,
    )


def _same_species(taxonomy: Any, organism_name: str, top: SequenceHit | None) -> bool:
    """Prefer M1's lineage (issue #55); fall back to name comparison when standalone."""
    if top is None:
        return False
    if taxonomy is not None and getattr(taxonomy, "determined", False):
        return bool(taxonomy.same_species_as(top.organism, getattr(top, "taxonomy_id", None)))
    return bool(organism_name and same_species(organism_name, top.organism))


def _same_genus(taxonomy: Any, organism_name: str, top: SequenceHit | None) -> bool:
    if top is None:
        return False
    if taxonomy is not None and getattr(taxonomy, "determined", False):
        return bool(taxonomy.same_genus_as(top.organism, getattr(top, "taxonomy_id", None)))
    return bool(organism_name and same_genus(organism_name, top.organism))


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
    if score.components.ligandable_homolog:
        parts.append("ligandable homolog" + (f" ({score.flags.get('ligandable_basis','')})"
                                             if score.flags.get("ligandable_basis") else ""))
    if score.components.surface_bonus:
        parts.append("surface-exposed or secreted")
    if score.components.membrane_penalty:
        parts.append("predicted membrane protein (penalty)")
    if score.components.human_homolog_penalty:
        parts.append(f"human homolog {score.components.human_homolog_penalty:.0%} identity (penalty)")
    if score.flags.get("same_species_hit"):
        parts.append("best structural hit is the same species — not a cross-organism transfer")
    elif score.flags.get("same_genus_hit"):
        parts.append("best structural hit is the same genus")
    if score.flags.get("antibiotic_target_not_resistance"):
        parts.append("antibiotic target in a susceptible species, not a resistance gene")
    if score.disqualified:
        parts.append(f"DISQUALIFIED: {score.disqualified}")
    return "; ".join(parts)


def human_homolog_disqualification(score: ProteinScore) -> str:
    """Empty when the protein may be selected; otherwise the reason it may not be."""
    try:
        identity = float(score.flags.get("human_homolog_identity"))
    except (TypeError, ValueError):
        return ""
    if identity < DISQUALIFYING_HUMAN_HOMOLOG_IDENTITY:
        return ""
    source = score.flags.get("human_homolog_source") or "unknown source"
    return f"human homolog {identity:.1f}% identity ({source})"


def rank_and_select(scores: list[ProteinScore], *, top_n: int = 50) -> list[ProteinScore]:
    """Rank every protein, mark the top N *eligible* ones, and give each one a reason.

    A disqualified protein keeps its rank and its score — the ranking still records what the
    evidence was worth, and the report can show what was removed and why — but it can never be
    selected, and its slot passes to the next eligible protein so the caller still receives
    ``top_n`` candidates.
    """
    ordered = sorted(
        scores,
        key=lambda s: (
            -s.score,
            -s.components.pdb_evidence,
            -s.distinct_entries,
            s.feature_id,
        ),
    )
    selected = 0
    for index, score in enumerate(ordered, start=1):
        score.rank = index
        score.disqualified = human_homolog_disqualification(score)
        if score.disqualified:
            score.selected = False
        elif selected < top_n:
            score.selected = True
            selected += 1
        else:
            score.selected = False
        score.reason = _reason(score)
    return ordered

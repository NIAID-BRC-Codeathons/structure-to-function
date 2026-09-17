from s2f.m2_triage.bvbrc_input import Protein, SpecialtyHit
from s2f.m2_triage.pdb_evidence import SequenceHit
from s2f.m2_triage.score import (
    MAX_HIT_SCORE,
    TriageComponents,
    hit_qualifies,
    hit_score,
    rank_and_select,
    score_protein,
    structure_quality,
    triage_score_from_components,
)


def make_hit(**kwargs) -> SequenceHit:
    defaults = dict(
        entity_id="1ABC_1",
        identity=1.0,
        evalue=1e-40,
        bitscore=400.0,
        query_beg=1,
        query_end=100,
        query_length=100,
        subject_length=100,
        experimental_method="X-ray",
        resolution=2.0,
    )
    defaults.update(kwargs)
    return SequenceHit(**defaults)


def make_protein(product: str = "Flavoprotein MioC", specialty=()) -> Protein:
    return Protein(
        feature_id="fig|1.1.peg.1",
        sequence="MKALIV",
        product=product,
        specialty=[
            SpecialtyHit(property_name=name, source=source, product="", identity=identity, query_coverage=None)
            for name, source, identity in specialty
        ],
    )


def test_gates_reject_weak_or_partial_hits() -> None:
    assert hit_qualifies(make_hit())
    assert not hit_qualifies(make_hit(evalue=1e-4))
    # 49 of 100 query residues aligned is below the 0.5 coverage gate.
    assert not hit_qualifies(make_hit(query_beg=1, query_end=49))
    assert hit_score(make_hit(evalue=1e-4)) == 0.0


def test_structure_quality_tiers() -> None:
    assert structure_quality(make_hit(resolution=2.5)) == 1.0
    assert structure_quality(make_hit(resolution=3.0)) == 0.8
    assert structure_quality(make_hit(resolution=4.0)) == 0.6
    assert structure_quality(make_hit(experimental_method="NMR", resolution=None)) == 0.7
    assert structure_quality(make_hit(experimental_method="", resolution=None)) == 0.6
    assert structure_quality(make_hit(experimental_method="Electron Microscopy", resolution=2.0)) == 1.0


def test_hit_score_adds_documented_bonuses() -> None:
    plain = make_hit()
    assert hit_score(plain) == 1.0

    with_ligand = make_hit(ligands=["MK1"])
    assert round(hit_score(with_ligand), 6) == 1.10

    everything = make_hit(ligands=["MK1"], partner_entities=1, has_go=True)
    assert round(hit_score(everything), 6) == round(MAX_HIT_SCORE, 6)


def test_metals_and_additives_do_not_earn_the_ligand_bonus() -> None:
    # classify_components strips them upstream, so a metals-only entry carries no ligands.
    metals_only = make_hit(ligands=[], metals=["ZN"])
    assert hit_score(metals_only) == 1.0


def test_pdb_evidence_is_normalized_not_clipped() -> None:
    strong = score_protein(make_protein(), [make_hit(ligands=["MK1"], partner_entities=1, has_go=True)], "found")
    weaker = score_protein(make_protein(), [make_hit(ligands=["MK1"])], "found")

    assert strong.components.pdb_evidence == 1.0
    # Both exceeded 1.0 before normalization; they must stay distinguishable.
    assert weaker.components.pdb_evidence < strong.components.pdb_evidence


def test_score_is_reproducible_from_recorded_components() -> None:
    protein = make_protein(
        product="hypothetical protein",
        specialty=[("Virulence Factor", "VFDB", 84.0), ("Essential Gene", "PATRIC", None)],
    )
    scored = score_protein(protein, [make_hit()], "found")

    recomputed = triage_score_from_components(
        TriageComponents(**scored.components.as_dict())
    )
    assert recomputed == scored.score
    assert scored.components.virulence_amr == 1.0
    assert scored.components.essential == 1.0
    assert scored.components.annotation_gap == 1.0
    assert scored.components.drug_target == 0.0


def test_human_homolog_is_a_penalty_scaled_by_identity() -> None:
    plain = score_protein(make_protein(), [make_hit()], "found")
    human = score_protein(
        make_protein(specialty=[("Human Homolog", "Human", 80.0)]), [make_hit()], "found"
    )

    assert human.components.human_homolog_penalty == 0.8
    assert round(plain.score - human.score, 6) == 0.2  # 0.25 weight x 0.80 identity
    assert human.score < plain.score


def test_no_hit_and_failed_query_are_distinguished() -> None:
    no_hit = score_protein(make_protein(), [], "no-hit")
    failed = score_protein(make_protein(), [], "query-failed", error="HTTP 500")

    assert no_hit.components.pdb_evidence == 0.0
    assert failed.components.pdb_evidence == 0.0
    assert no_hit.flags["no_pdb_hit"] and failed.flags["no_pdb_hit"]
    assert no_hit.retrieval_status != failed.retrieval_status
    assert failed.error == "HTTP 500"

    ranked = rank_and_select([no_hit, failed], top_n=1)
    reasons = {score.retrieval_status: score.reason for score in ranked}
    assert "no qualifying PDB hit" in reasons["no-hit"]
    assert "search failed" in reasons["query-failed"]


def test_rank_and_select_keeps_every_protein_with_a_reason() -> None:
    strong = score_protein(make_protein("Aspartate--ammonia ligase"), [make_hit(ligands=["ASN"])], "found")
    weak = score_protein(make_protein("hypothetical protein"), [], "no-hit")

    ranked = rank_and_select([weak, strong], top_n=1)

    assert [s.rank for s in ranked] == [1, 2]
    assert ranked[0].feature_id == strong.feature_id
    assert ranked[0].selected is True
    assert ranked[1].selected is False
    assert all(score.reason for score in ranked)


def test_best_hit_prefers_higher_scoring_structure() -> None:
    low_res = make_hit(entity_id="1LOW_1", resolution=4.0)
    high_res = make_hit(entity_id="1HIG_1", resolution=1.5)
    scored = score_protein(make_protein(), [low_res, high_res], "found")

    assert scored.best_hit is not None
    assert scored.best_hit.entity_id == "1HIG_1"
    assert scored.qualifying_hits == 2
    assert scored.distinct_entries == 2

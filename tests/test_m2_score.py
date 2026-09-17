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
    assert no_hit.retrieval_status != failed.retrieval_status
    assert failed.error == "HTTP 500"

    # `no_pdb_hit` means "we searched and found nothing", so a failed search must not set it:
    # downstream consumers treat that file as "no homolog exists" and would launder a network
    # failure into a structural claim.
    assert no_hit.flags["no_pdb_hit"] is True
    assert no_hit.flags["pdb_search_failed"] is False
    assert failed.flags["no_pdb_hit"] is False
    assert failed.flags["pdb_search_failed"] is True

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


# --- issue #12: the components #10 unlocked, and the AMR correction --------------------


class FakeAnnotation:
    """Stands in for function.FunctionalAnnotation; only the fields scoring reads."""

    def __init__(self, *, membrane=False, membrane_source="", secreted=False,
                 surface_exposed=False, signal_source="", localization_source=""):
        self.membrane = membrane
        self.membrane_source = membrane_source
        self.secreted = secreted
        self.surface_exposed = surface_exposed
        self.signal_source = signal_source
        self.localization_source = localization_source


def specialty(property_name, classification="", identity=None, source="PATRIC"):
    return SpecialtyHit(
        property_name=property_name, source=source, product="", identity=identity,
        query_coverage=None, classification=classification,
    )


def test_a_drug_target_in_a_susceptible_species_is_not_resistance_evidence() -> None:
    """gyrA and rpoB are what fluoroquinolones and rifamycins hit, not resistance genes (#30)."""
    target = Protein(
        feature_id="fig|1.1.peg.1", sequence="MKALIV", product="DNA gyrase subunit A",
        specialty=[specialty("Antibiotic Resistance", "antibiotic target in susceptible species")],
    )
    scored = score_protein(target, [make_hit()], "found")

    assert scored.components.virulence_amr == 0.0       # not resistance
    assert scored.components.drug_target == 1.0          # but it *is* a drug target
    assert scored.flags["antibiotic_target_not_resistance"] is True

    # rank_and_select writes the reason, and a reader must see why it is not resistance.
    ranked = rank_and_select([scored], top_n=1)
    assert "antibiotic target in a susceptible species" in ranked[0].reason
    assert "not a resistance gene" in ranked[0].reason


def test_a_real_resistance_mechanism_still_counts() -> None:
    efflux = Protein(
        feature_id="fig|1.1.peg.2", sequence="MKALIV", product="efflux pump",
        specialty=[specialty("Antibiotic Resistance", "['efflux pump conferring antibiotic resistance']")],
    )
    scored = score_protein(efflux, [make_hit()], "found")

    assert scored.components.virulence_amr == 1.0
    assert "efflux pump" in scored.flags["amr_basis"]


def test_an_unclassified_amr_row_is_half_not_full_credit() -> None:
    """The export sometimes omits classification; that is uncertainty, not confirmation."""
    unknown = Protein(
        feature_id="fig|1.1.peg.3", sequence="MKALIV", product="something",
        specialty=[specialty("Antibiotic Resistance", "")],
    )
    scored = score_protein(unknown, [make_hit()], "found")

    assert scored.components.virulence_amr == 0.5
    assert "no classification" in scored.flags["amr_basis"]


def test_a_virulence_factor_still_scores_in_full() -> None:
    virulent = Protein(
        feature_id="fig|1.1.peg.4", sequence="MKALIV", product="adhesin",
        specialty=[specialty("Virulence Factor", identity=84.0)],
    )
    assert score_protein(virulent, [make_hit()], "found").components.virulence_amr == 1.0


def test_surface_exposure_adds_and_membrane_subtracts() -> None:
    plain = score_protein(make_protein(), [make_hit()], "found")
    surface = score_protein(
        make_protein(), [make_hit()], "found",
        annotation=FakeAnnotation(secreted=True, surface_exposed=True, signal_source="signalp6"),
    )
    membrane = score_protein(
        make_protein(), [make_hit()], "found",
        annotation=FakeAnnotation(membrane=True, membrane_source="deeptmhmm"),
    )

    assert round(surface.score - plain.score, 6) == 0.10
    assert round(membrane.score - plain.score, 6) == -0.15   # deprioritized, not excluded
    assert membrane.score > 0  # pitfall #3 says deprioritize, not drop


def test_an_unmeasured_flag_is_recorded_as_unmeasured_not_false() -> None:
    """A missing DeepTMHMM run must not read as "no membrane proteins in this genome"."""
    unmeasured = score_protein(make_protein(), [make_hit()], "found", annotation=FakeAnnotation())
    measured_false = score_protein(
        make_protein(), [make_hit()], "found",
        annotation=FakeAnnotation(membrane=False, membrane_source="deeptmhmm"),
    )

    # Both contribute 0 to the score...
    assert unmeasured.components.membrane_penalty == 0.0
    assert measured_false.components.membrane_penalty == 0.0
    # ...but only one of them was actually looked at.
    assert unmeasured.components_available["membrane_penalty"] is False
    assert measured_false.components_available["membrane_penalty"] is True


def test_no_annotation_at_all_leaves_both_components_unmeasured() -> None:
    scored = score_protein(make_protein(), [make_hit()], "found")

    assert scored.components_available == {"surface_bonus": False, "membrane_penalty": False}
    assert scored.components.surface_bonus == 0.0


def test_score_is_still_reproducible_from_the_recorded_components() -> None:
    scored = score_protein(
        make_protein(product="hypothetical protein"), [make_hit()], "found",
        annotation=FakeAnnotation(membrane=True, membrane_source="deeptmhmm"),
    )
    assert triage_score_from_components(TriageComponents(**scored.components.as_dict())) == scored.score

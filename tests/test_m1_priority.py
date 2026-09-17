"""M1 pathogenesis priority: the score, its reasons, and its separation from M2's triage."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from s2f.m1_genome.priority import (CATEGORIES, WEIGHTS, classify_mechanisms,
                                    index_specialty, mechanism_hypothesis,
                                    pathogenesis_proteins, rank_proteins, score_protein,
                                    specialty_type)

FIXTURE = Path(__file__).resolve().parents[1] / "fixtures/m1/bvbrc_api/klebsiella_hs11286.json"
ADHESIN = "fig|1125630.4.peg.998"
SIDEROPHORE = "fig|1125630.4.peg.3382"
REGULATOR = "fig|1125630.4.peg.1591"
TET = "fig|1125630.4.peg.5465"
HYPOTHETICAL = "fig|1125630.4.peg.5521"


@pytest.fixture(scope="module")
def payload() -> dict:
    return json.loads(FIXTURE.read_text())


@pytest.fixture
def ranked(payload) -> dict:
    return rank_proteins(payload["features"], payload["specialty"])


# --- mechanism classification ----------------------------------------------
def test_mechanisms_come_back_in_canonical_order_not_match_order() -> None:
    # The first category is what the hypothesis sentence leads with, so a protein must
    # not change its apparent primary mechanism because a regex ran in a different order.
    found = classify_mechanisms("fimbrial adhesin and haemolysin toxin")
    assert found == ["adhesion", "toxin"]
    assert [c for c, *_ in CATEGORIES].index("adhesion") < \
           [c for c, *_ in CATEGORIES].index("toxin")


def test_a_curated_classification_can_add_a_mechanism_the_text_missed() -> None:
    # "EcpR" says nothing; VFDB's Adherence classification does.
    assert classify_mechanisms("Transcriptional regulator EcpR", ["Adherence"]) == ["adhesion"]


def test_an_amr_classification_string_is_not_read_as_a_mechanism() -> None:
    # This real BV-BRC string contains "regulator ... resistance genes" and must not be
    # mistaken for a virulence-regulation mechanism claim.
    assert classify_mechanisms(
        "", ["regulator modulating expression of antibiotic resistance genes"]) == []


def test_a_generic_two_component_system_is_not_virulence_regulation() -> None:
    """A typical enterobacterium encodes ~30 two-component systems.

    Most regulate osmolarity, phosphate or nitrogen. Treating the generic phrase as
    host-interaction machinery scored ordinary metabolic regulators as virulence
    regulators — this is the guard against that.
    """
    assert classify_mechanisms("Two-component response regulator OmpR") == []
    assert classify_mechanisms("transcriptional regulator, LysR family") == []


def test_a_generic_regulator_counts_once_a_database_says_it_is_a_virulence_factor() -> None:
    # The curated flag supplies the context the product string does not.
    assert classify_mechanisms(
        "Two-component response regulator", is_virulence=True) == ["regulation"]


def test_a_named_virulence_regulator_is_recognised_from_text_alone() -> None:
    for product in ["PhoP two-component regulator", "quorum sensing regulator LuxR",
                    "ToxR transcriptional activator", "virulence regulator VirF"]:
        assert classify_mechanisms(product) == ["regulation"], product


def test_specialty_property_names_map_to_scoring_types() -> None:
    assert specialty_type("Virulence Factor") == "virulence"
    assert specialty_type("Antibiotic Resistance") == "amr"
    assert specialty_type("Human Homolog") == "human_homolog"
    assert specialty_type(None) == "other"


# --- the hypothesis sentence ------------------------------------------------
def test_the_hypothesis_is_always_phrased_as_one() -> None:
    text = mechanism_hypothesis("fimbrial adhesin", ["adhesion"], True)
    assert "likely" in text
    assert "Adhesion & attachment" in text


def test_an_unevidenced_protein_says_so_instead_of_guessing() -> None:
    assert "Uncharacterised" in mechanism_hypothesis("hypothetical protein", [], False)
    assert "not resolvable" in mechanism_hypothesis("Some protein", [], True)
    assert "Housekeeping" in mechanism_hypothesis("DNA gyrase subunit A", [], False)


# --- scoring ---------------------------------------------------------------
def test_every_point_awarded_has_a_reason_attached(payload) -> None:
    index = index_specialty(payload["specialty"])
    feature = next(f for f in payload["features"] if f["patric_id"] == SIDEROPHORE)
    priority = score_protein(feature, index[SIDEROPHORE])
    assert priority["score"] == sum(item["points"] for item in priority["breakdown"])
    assert all(item["reason"] for item in priority["breakdown"])


def test_the_weights_used_are_recorded_with_the_score(payload) -> None:
    # So a report can print the weights it actually ran with rather than prose that drifts.
    priority = score_protein({"patric_id": "x", "product": "toxin"}, None)
    assert priority["weights"] == WEIGHTS


def test_a_human_homolog_is_penalised(payload) -> None:
    feature = {"patric_id": "x", "product": "fimbrial adhesin", "gene": "fimA"}
    without = score_protein(feature, {"types": {"virulence"}, "classes": set(), "hits": []})
    with_homolog = score_protein(
        feature, {"types": {"virulence", "human_homolog"}, "classes": set(), "hits": []})
    # A target that looks like a host protein is a selectivity risk, so it must rank lower.
    assert with_homolog["score"] == without["score"] + WEIGHTS["human_homolog"]
    assert with_homolog["score"] < without["score"]


def test_selection_requires_evidence_not_merely_a_positive_score() -> None:
    # A named surface protein scores 2 but has no virulence, AMR, drug-target or mechanism
    # evidence, so it must not be handed to M2 as a candidate.
    priority = score_protein(
        {"patric_id": "x", "product": "outer membrane lipoprotein", "gene": "ompX"}, None)
    assert priority["score"] > 0
    assert priority["selected"] is False


def test_an_unevidenced_hypothetical_scores_zero_and_is_not_selected(ranked) -> None:
    assert ranked[HYPOTHETICAL]["score"] == 0
    assert ranked[HYPOTHETICAL]["selected"] is False
    assert ranked[HYPOTHETICAL]["categories"] == []


def test_real_virulence_and_amr_proteins_outrank_hypotheticals(ranked) -> None:
    assert ranked[ADHESIN]["rank"] < ranked[HYPOTHETICAL]["rank"]
    assert ranked[SIDEROPHORE]["rank"] < ranked[HYPOTHETICAL]["rank"]
    assert ranked[ADHESIN]["selected"] and ranked[SIDEROPHORE]["selected"]


def test_the_adhesin_is_scored_on_both_text_and_curated_evidence(ranked) -> None:
    priority = ranked[ADHESIN]
    assert "adhesion" in priority["categories"]
    assert "virulence" in priority["specialty_types"]
    reasons = " ".join(item["reason"] for item in priority["breakdown"])
    assert "Virulence factor" in reasons
    assert "Adhesion" in reasons


def test_the_amr_regulator_gets_no_mechanism_category(ranked) -> None:
    # BasR/PmrA is a real AMR determinant but its classification is about resistance
    # regulation, not host interaction.
    assert ranked[REGULATOR]["categories"] == []
    assert "amr" in ranked[REGULATOR]["specialty_types"]


def test_ranks_are_dense_and_start_at_one(ranked, payload) -> None:
    ranks = sorted(p["rank"] for p in ranked.values())
    assert ranks == list(range(1, len(payload["features"]) + 1))


def test_ties_break_on_length_so_the_order_is_deterministic() -> None:
    features = [
        {"patric_id": "short", "product": "toxin", "aa_length": 100},
        {"patric_id": "long", "product": "toxin", "aa_length": 900},
    ]
    ranked = rank_proteins(features, [])
    assert ranked["long"]["rank"] < ranked["short"]["rank"]
    assert ranked["long"]["score"] == ranked["short"]["score"]


def test_a_non_integer_length_does_not_break_ranking() -> None:
    # BV-BRC writes coordinates and lengths as strings, and sometimes as ''.
    ranked = rank_proteins(
        [{"patric_id": "a", "product": "toxin", "aa_length": "300"},
         {"patric_id": "b", "product": "toxin", "aa_length": ""}], [])
    assert ranked["a"]["rank"] == 1


def test_features_without_an_id_are_skipped_not_ranked_as_none() -> None:
    ranked = rank_proteins([{"product": "toxin"}, {"patric_id": "a", "product": "x"}], [])
    assert list(ranked) == ["a"]


# --- separation from M2 -----------------------------------------------------
def test_the_priority_record_never_carries_m2s_triage_fields(ranked) -> None:
    """`triage` is M2's key with different required fields; M1 must not write it.

    Overwriting it would destroy the structural evidence M3 gates on
    (`docs/00a-data-contract.md`).
    """
    for priority in ranked.values():
        assert "triage" not in priority
        assert set(priority) == {"score", "breakdown", "categories", "specialty_types",
                                 "mechanism_hypothesis", "selected", "weights", "rank"}


def test_priority_and_triage_disagreeing_is_not_a_bug(ranked) -> None:
    # M1 ranks on annotation evidence, M2 on structural tractability. The adhesin is a
    # top M1 candidate precisely because nothing structural has been looked at yet.
    assert ranked[ADHESIN]["score"] > 0
    assert "components" not in ranked[ADHESIN]  # M2's triage shape, not M1's


# --- pathogenesis candidates ------------------------------------------------
def test_the_two_kinds_of_evidence_stay_distinguishable(payload) -> None:
    rows = pathogenesis_proteins(payload["features"], payload["specialty"][:4])
    kinds = {row["evidence"] for row in rows}
    # A curated VFDB hit and a product-keyword guess are very different claims.
    assert kinds == {"virulence specialty gene", "product keyword match"}


def test_a_protein_is_listed_once_with_its_strongest_evidence(payload) -> None:
    rows = pathogenesis_proteins(payload["features"], payload["specialty"][:4])
    ids = [row["patric_id"] for row in rows]
    assert len(ids) == len(set(ids))
    adhesin = next(r for r in rows if r["patric_id"] == ADHESIN)
    assert adhesin["evidence"] == "virulence specialty gene"

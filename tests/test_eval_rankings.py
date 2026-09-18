"""Reading the two scores, and counterfactual re-scoring for the ablation."""

from __future__ import annotations

import pytest

from s2f.evaluation.rankings import (
    SCORES,
    TRIAGE_COMPONENTS,
    TRIAGE_WEIGHTS,
    ablated_triage_score,
    available_scores,
    build_ranking,
    component_availability,
    m1_priority_score,
    score_accessor,
    triage_components,
    triage_score,
)
from s2f.m2_triage.score import triage_score_from_components


def with_triage(fid, components, *, length=100, selected=False, gene=None, product=""):
    full = {name: 0.0 for name in TRIAGE_COMPONENTS}
    full.update(components)

    class _C:
        def as_dict(self):
            return full

    return {
        "feature_id": fid,
        "gene": gene,
        "product": product,
        "aa_length": length,
        "triage": {
            "score": triage_score_from_components(_C()),
            "components": full,
            "selected": selected,
            "reason": "",
        },
    }


def with_m1(fid, score, *, length=100):
    return {"feature_id": fid, "aa_length": length, "m1_priority": {"score": score}}


# --------------------------------------------------------------------------- accessors

def test_accessors_read_their_own_key_and_nothing_else():
    protein = {**with_m1("f", 4), **with_triage("f", {"pdb_evidence": 1.0})}
    assert m1_priority_score(protein) == 4
    assert triage_score(protein) == pytest.approx(0.4)
    assert m1_priority_score({"feature_id": "x"}) is None
    assert triage_score({"feature_id": "x"}) is None


def test_booleans_are_not_scores():
    """`True` is an int in Python; a score of True would sort above everything."""
    assert m1_priority_score({"m1_priority": {"score": True}}) is None


def test_unknown_score_is_refused():
    with pytest.raises(ValueError, match="unknown score"):
        score_accessor("vibes")


def test_available_scores_reports_what_the_report_carries():
    assert available_scores([with_m1("a", 1)]) == ["m1_priority"]
    assert available_scores([with_triage("a", {})]) == ["triage"]
    both = {**with_m1("a", 1), **with_triage("a", {})}
    assert available_scores([both]) == list(SCORES)
    assert available_scores([{"feature_id": "a"}]) == []


# ----------------------------------------------------------------------------- ranking

def test_build_ranking_orders_by_score_then_length():
    proteins = [with_m1("a", 3, length=10), with_m1("b", 5), with_m1("c", 3, length=900)]
    assert build_ranking(proteins, "m1_priority").order == ["b", "c", "a"]


def test_build_ranking_refuses_a_score_nothing_carries():
    with pytest.raises(ValueError, match="no protein carries a 'triage' score"):
        build_ranking([with_m1("a", 1)], "triage")


def test_proteins_without_a_feature_id_are_dropped():
    assert build_ranking([with_m1("a", 1), {"m1_priority": {"score": 9}}],
                         "m1_priority").order == ["a"]


def test_position_is_one_based():
    ranking = build_ranking([with_m1("a", 2), with_m1("b", 1)], "m1_priority")
    assert ranking.position() == {"a": 1, "b": 2}
    assert ranking.top(1) == ["a"]


# --------------------------------------------------------------------- counterfactuals

def test_ablated_score_uses_m2s_own_arithmetic():
    protein = with_triage("f", {"pdb_evidence": 1.0, "annotation_gap": 1.0})
    assert protein["triage"]["score"] == pytest.approx(0.4 + 0.3)
    assert ablated_triage_score(protein, "annotation_gap") == pytest.approx(0.4)
    assert ablated_triage_score(protein, "pdb_evidence") == pytest.approx(0.3)


def test_ablating_several_components_at_once():
    protein = with_triage("f", {"pdb_evidence": 1.0, "annotation_gap": 1.0, "essential": 1.0})
    assert ablated_triage_score(
        protein, ["pdb_evidence", "annotation_gap"]
    ) == pytest.approx(TRIAGE_WEIGHTS["essential"])


def test_ablating_a_penalty_raises_the_score():
    """human_homolog_penalty carries a negative weight, so removing it must help."""
    protein = with_triage("f", {"pdb_evidence": 1.0, "human_homolog_penalty": 1.0})
    assert ablated_triage_score(protein, "human_homolog_penalty") > protein["triage"]["score"]


def test_ablating_everything_gives_zero():
    protein = with_triage("f", {name: 1.0 for name in TRIAGE_COMPONENTS})
    assert ablated_triage_score(protein, TRIAGE_COMPONENTS) == pytest.approx(0.0)


def test_unknown_component_is_refused_by_name():
    with pytest.raises(ValueError, match="unknown triage component"):
        ablated_triage_score(with_triage("f", {}), "enthusiasm")


def test_ablation_of_an_unscored_protein_is_none():
    assert ablated_triage_score({"feature_id": "f"}, "pdb_evidence") is None


def test_components_default_missing_entries_to_zero():
    protein = {"feature_id": "f", "triage": {"score": 0, "components": {"pdb_evidence": 0.5}}}
    assert triage_components(protein) == {
        **{name: 0.0 for name in TRIAGE_COMPONENTS}, "pdb_evidence": 0.5
    }


def test_component_availability_is_read_not_assumed():
    protein = with_triage("f", {})
    protein["triage"]["components_available"] = {"surface_bonus": True,
                                                 "membrane_penalty": False}
    assert component_availability(protein) == {"surface_bonus": True,
                                               "membrane_penalty": False}
    assert component_availability({"feature_id": "f"}) == {}

"""M3 outcome labels, and rank agreement between the two scores."""

from __future__ import annotations

import pytest

from s2f.evaluation.agreement import average_ranks, compare, spearman
from s2f.evaluation.outcome import OUTCOME_LABELS, classify, coverage_warning, outcomes_from_report
from s2f.evaluation.rankings import build_ranking
from tests.test_eval_rankings import with_m1, with_triage


# ------------------------------------------------------------------------ outcome classes

@pytest.mark.parametrize(
    "record,expected",
    [
        ({"collection_status": "collected", "usable_for_docking": True}, "dockable"),
        ({"collection_status": "collected", "usable_for_docking": False}, "unusable"),
        ({"collection_status": "prediction_required", "usable_for_docking": None}, "no_structure"),
        ({"collection_status": "failed", "usable_for_docking": None}, "failed"),
        ({"collection_status": "collected", "usable_for_docking": None}, "failed"),
        ({}, "failed"),
    ],
)
def test_classify(record, expected):
    assert classify(record) == expected


def test_a_technical_failure_is_excluded_not_a_negative():
    """Chen's workflow contract: a technical failure is never a negative finding.

    A download that 503'd is not evidence that the protein has no usable structure, so it
    must not land in the denominator of a precision figure.
    """
    assert OUTCOME_LABELS["failed"] == "excluded"
    assert OUTCOME_LABELS["dockable"] == "positive"
    assert OUTCOME_LABELS["unusable"] == "negative"


def test_no_structure_is_a_negative_and_kept_distinct_from_unusable():
    """M2 selected it on a tractability claim; M3 finding nothing is that claim failing."""
    assert OUTCOME_LABELS["no_structure"] == "negative"
    assert "no_structure" != "unusable"


def build_report():
    proteins = [
        with_triage("a", {"pdb_evidence": 0.9}, selected=True, product="A"),
        with_triage("b", {"pdb_evidence": 0.8}, selected=True, product="B"),
        with_triage("c", {"pdb_evidence": 0.7}, selected=True, product="C"),
        with_triage("d", {"pdb_evidence": 0.6}, selected=True, product="D"),
        with_triage("e", {"pdb_evidence": 0.1}, selected=False, product="E"),
    ]
    report = {
        "proteins": proteins,
        "structures": [
            {"feature_id": "a", "source": "pdb", "collection_status": "collected",
             "usable_for_docking": True, "accession": "1ABC",
             "quality": {"policy_version": "provisional-2026-09-17"}},
            {"feature_id": "b", "source": "pdb", "collection_status": "collected",
             "usable_for_docking": False, "accession": "2DEF"},
            {"feature_id": "c", "source": "predicted",
             "collection_status": "prediction_required", "usable_for_docking": None},
            {"feature_id": "d", "source": "pdb", "collection_status": "failed",
             "usable_for_docking": None},
        ],
    }
    return report, proteins


def test_outcomes_are_read_from_structures():
    report, proteins = build_report()
    outcomes = outcomes_from_report(report, proteins)
    assert outcomes.counts == {"dockable": 1, "unusable": 1, "no_structure": 1, "failed": 1}
    assert outcomes.policy_version == "provisional-2026-09-17"
    assert outcomes.covered == 4
    labels = {m.feature_id: m.label for m in outcomes.matches}
    assert labels == {"a": "positive", "b": "negative", "c": "negative", "d": "excluded"}


def test_outcomes_carry_the_product_for_the_report():
    report, proteins = build_report()
    outcomes = outcomes_from_report(report, proteins)
    assert outcomes.by_feature()["a"].product == "A"


def test_selected_proteins_m3_never_reached_are_reported_not_scored():
    report, proteins = build_report()
    report["structures"] = [s for s in report["structures"] if s["feature_id"] != "d"]
    outcomes = outcomes_from_report(report, proteins)
    assert outcomes.selected_but_uncollected == ["d"]
    assert "d" not in outcomes.by_feature()


def test_duplicate_structure_records_are_collapsed():
    report, proteins = build_report()
    report["structures"].append(dict(report["structures"][0]))
    assert outcomes_from_report(report, proteins).counts["dockable"] == 1


def test_a_report_with_no_structures_yields_nothing_rather_than_failing():
    report, proteins = build_report()
    report["structures"] = []
    outcomes = outcomes_from_report(report, proteins)
    assert outcomes.matches == [] and outcomes.covered == 0


def test_coverage_warning_fires_when_the_cut_reaches_past_the_selected_set():
    report, proteins = build_report()
    outcomes = outcomes_from_report(report, proteins)
    assert coverage_warning(outcomes, 5, ["a", "b", "c", "d", "e"]) is not None
    assert coverage_warning(outcomes, 2, ["a", "b"]) is None


# ------------------------------------------------------------------------------ agreement

def test_average_ranks_share_ties():
    assert average_ranks([10, 20, 20, 30]) == [1.0, 2.5, 2.5, 4.0]
    assert average_ranks([5, 5, 5]) == [2.0, 2.0, 2.0]


def test_spearman_endpoints():
    assert spearman([1, 2, 3, 4], [1, 2, 3, 4]) == pytest.approx(1.0)
    assert spearman([1, 2, 3, 4], [4, 3, 2, 1]) == pytest.approx(-1.0)


def test_spearman_is_none_when_a_side_is_constant():
    """Undefined, not zero — and a score constant across the population is a finding."""
    assert spearman([1, 2, 3], [7, 7, 7]) is None
    assert spearman([1], [2]) is None


def test_spearman_refuses_unequal_lengths():
    with pytest.raises(ValueError, match="same length"):
        spearman([1, 2], [1])


def test_spearman_handles_heavy_ties():
    """Both real scores are small integers with large tie blocks."""
    assert spearman([5, 5, 5, 1], [5, 5, 5, 1]) == pytest.approx(1.0)


def build_pair():
    proteins = [
        {**with_m1("a", 6, length=100), **with_triage("a", {"pdb_evidence": 0.1})},
        {**with_m1("b", 5, length=100), **with_triage("b", {"pdb_evidence": 0.9})},
        {**with_m1("c", 4, length=100), **with_triage("c", {"pdb_evidence": 0.5})},
        {**with_m1("d", 1, length=100), **with_triage("d", {"pdb_evidence": 0.7})},
    ]
    return proteins


def test_compare_reports_overlap_and_the_proteins_each_side_picked_alone():
    proteins = build_pair()
    left = build_ranking(proteins, "m1_priority")
    right = build_ranking(proteins, "triage")
    result = compare(left, right, proteins, k=2)
    assert result.common == 4
    assert result.overlap_at_k == 1                      # only 'b' is in both top 2
    assert result.only_left == ["a"]
    assert result.only_right == ["d"]
    assert result.jaccard_at_k == pytest.approx(1 / 3)


def test_compare_names_the_biggest_disagreement():
    proteins = build_pair()
    result = compare(build_ranking(proteins, "m1_priority"),
                     build_ranking(proteins, "triage"), proteins, k=2)
    worst = result.biggest_disagreements[0]
    assert worst["feature_id"] in {"a", "d"}
    assert worst["gap"] == 3
    assert "m1_priority_rank" in worst and "triage_rank" in worst


def test_compare_uses_only_proteins_both_scores_reached():
    """Padding an unscored protein in at the bottom would manufacture correlation."""
    proteins = build_pair() + [with_m1("e", 3)]
    result = compare(build_ranking(proteins, "m1_priority"),
                     build_ranking([p for p in proteins if "triage" in p], "triage"),
                     proteins, k=2)
    assert result.common == 4


def test_agreement_serialises():
    proteins = build_pair()
    payload = compare(build_ranking(proteins, "m1_priority"),
                      build_ranking(proteins, "triage"), proteins, k=2).to_dict()
    assert payload["left"] == "m1_priority" and payload["right"] == "triage"
    assert set(payload) >= {"spearman", "overlap_at_k", "jaccard_at_k", "only_left"}


# ------------------------------------------------------------------ base rate and lift
# A precision figure read without its base rate is the easiest false headline in this
# package. On the first real run M3 passed 298 of 300, so precision@50 was 100% against a
# base rate of 99.3% — a lift of half a point.

def test_base_rate_is_over_the_labelled_population_only():
    report, proteins = build_report()
    outcomes = outcomes_from_report(report, proteins)
    # 1 dockable, 1 unusable, 1 no_structure labelled; the failed one is excluded
    assert outcomes.base_rate == pytest.approx(1 / 3)


def test_lift_is_precision_minus_base_rate():
    report, proteins = build_report()
    outcomes = outcomes_from_report(report, proteins)
    assert outcomes.lift(1.0) == pytest.approx(1 - 1 / 3)
    assert outcomes.lift(outcomes.base_rate) == pytest.approx(0.0)
    assert outcomes.lift(None) is None


def test_base_rate_and_lift_are_none_without_labels():
    report, proteins = build_report()
    report["structures"] = []
    outcomes = outcomes_from_report(report, proteins)
    assert outcomes.base_rate is None
    assert outcomes.lift(1.0) is None


def test_a_near_constant_label_is_reported_as_non_discriminating():
    """298 of 300 passing is what the first HS11286 run actually produced."""
    report, proteins = build_report()
    report["structures"] = [
        {"feature_id": f"p{i}", "source": "pdb", "collection_status": "collected",
         "usable_for_docking": i >= 2}
        for i in range(300)
    ]
    outcomes = outcomes_from_report(report, proteins)
    assert outcomes.base_rate == pytest.approx(298 / 300)
    assert outcomes.discriminates() is False


def test_a_balanced_label_discriminates():
    report, proteins = build_report()
    report["structures"] = [
        {"feature_id": f"p{i}", "source": "pdb", "collection_status": "collected",
         "usable_for_docking": i % 2 == 0}
        for i in range(100)
    ]
    assert outcomes_from_report(report, proteins).discriminates() is True


def test_failures_stay_out_of_the_base_rate():
    report, proteins = build_report()
    report["structures"] = [
        {"feature_id": "a", "source": "pdb", "collection_status": "collected",
         "usable_for_docking": True},
        {"feature_id": "b", "source": "pdb", "collection_status": "failed",
         "usable_for_docking": None},
    ]
    outcomes = outcomes_from_report(report, proteins)
    assert outcomes.counts["failed"] == 1
    assert outcomes.base_rate == pytest.approx(1.0)   # not 0.5

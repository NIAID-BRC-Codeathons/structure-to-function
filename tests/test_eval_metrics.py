"""Precision, recall, ties, baselines and the leakage audit (issue #49)."""

from __future__ import annotations

import pytest

from s2f.evaluation.metrics import (
    BASELINES,
    CURATED_REASONS,
    TEXT_REASONS,
    _classify_reason,
    evaluate,
    leakage_audit,
    rank_by_baseline,
    rank_by_priority,
    tie_report,
    weights_fingerprint,
)
from s2f.evaluation.truthset import load_truthset, match_proteins
from s2f.m1_genome.priority import WEIGHTS, score_protein

HEADER = "symbol\tlabel\tclass\tbasis\tnote\n"


def truth(tmp_path, body: str):
    path = tmp_path / "t.truth.tsv"
    path.write_text(HEADER + body, encoding="utf-8")
    return load_truthset(path)


def protein(fid: str, score: float, *, gene=None, product="", length=100, breakdown=None):
    return {
        "feature_id": fid,
        "gene": gene,
        "product": product,
        "aa_length": length,
        "m1_priority": {"score": score, "breakdown": breakdown or []},
    }


# ------------------------------------------------------------------ reason attribution

def test_every_reason_priority_can_emit_is_classified():
    """If `priority.score_protein` renames a reason, this fails instead of the audit
    silently reclassifying that evidence as product-text derived."""
    cases = [
        ({"product": "x"}, {"types": {"virulence"}}),
        ({"product": "x"}, {"types": {"amr"}}),
        ({"product": "x"}, {"types": {"drug_target"}}),
        ({"product": "x"}, {"types": {"essential"}}),
        ({"product": "x"}, {"types": {"transporter"}}),
        ({"product": "x"}, {"types": {"human_homolog"}}),
        ({"product": "type III secretion system effector"}, None),
        ({"product": "outer membrane lipoprotein"}, None),
        ({"product": "x", "gene": "abcD"}, None),
    ]
    seen = set()
    for feature, specialty in cases:
        record = score_protein(feature, specialty)
        for entry in record["breakdown"]:
            reason = entry["reason"]
            assert _classify_reason(reason) is not None, f"unattributed reason: {reason!r}"
            seen.add(_classify_reason(reason)[1])
    assert seen == set(CURATED_REASONS.values()) | set(TEXT_REASONS.values())


def test_curated_and_text_components_together_cover_every_weight():
    """No scoring component may be left out of the split."""
    components = set(CURATED_REASONS.values()) | set(TEXT_REASONS.values())
    assert components == set(WEIGHTS)


def test_unknown_reason_is_surfaced_not_swallowed(tmp_path):
    """A reason the audit cannot attribute is named, so the split is never silently wrong."""
    proteins = [
        protein("f1", 5, gene="wza", breakdown=[{"points": 3, "reason": "Something new"}])
    ]
    truthset = truth(tmp_path, "wza\tpositive\tcapsule\tbecause\t\n")
    matches, _ = match_proteins(proteins, truthset)
    audit = leakage_audit(proteins, {m.feature_id: m for m in matches})
    assert audit.unclassified_reasons == ["Something new"]
    assert audit.curated_points == 0 and audit.text_points == 0


# ------------------------------------------------------------------------- fingerprint

def test_weights_fingerprint_is_stable_and_changes_with_the_weights():
    assert weights_fingerprint() == weights_fingerprint(dict(WEIGHTS))
    nudged = dict(WEIGHTS)
    nudged["mechanism"] = nudged["mechanism"] + 1
    assert weights_fingerprint(nudged) != weights_fingerprint()


# ------------------------------------------------------------------------------ ranking

def test_rank_by_priority_orders_by_score_then_length():
    proteins = [
        protein("a", 3, length=100),
        protein("b", 5, length=10),
        protein("c", 3, length=900),
    ]
    assert [p["feature_id"] for p in rank_by_priority(proteins)] == ["b", "c", "a"]


def test_rank_by_priority_refuses_a_report_m1_never_scored():
    with pytest.raises(ValueError, match="no protein carries m1_priority"):
        rank_by_priority([{"feature_id": "a"}])


def test_unscored_proteins_are_dropped_not_ranked_last():
    ranked = rank_by_priority([protein("a", 1), {"feature_id": "b"}])
    assert [p["feature_id"] for p in ranked] == ["a"]


@pytest.mark.parametrize("baseline", BASELINES)
def test_baselines_are_deterministic(baseline):
    proteins = [protein(f"f{i}", i, length=i * 10, product="toxin") for i in range(8)]
    first = [p["feature_id"] for p in rank_by_baseline(proteins, baseline, seed=1)]
    second = [p["feature_id"] for p in rank_by_baseline(proteins, baseline, seed=1)]
    assert first == second


def test_keyword_baseline_puts_pathogenesis_vocabulary_first():
    proteins = [
        protein("plain", 0, product="hypothetical protein", length=900),
        protein("viru", 0, product="haemolysin toxin", length=10),
    ]
    assert [p["feature_id"] for p in rank_by_baseline(proteins, "keyword")] == ["viru", "plain"]


def test_unknown_baseline_is_refused():
    with pytest.raises(ValueError, match="unknown baseline"):
        rank_by_baseline([protein("a", 1)], "vibes")


# -------------------------------------------------------------------------------- ties

def test_tie_report_finds_the_ceiling_and_a_split_block():
    ranked = [protein(f"f{i}", 7) for i in range(4)] + [protein(f"g{i}", 3) for i in range(4)]
    report = tie_report(ranked, k=2)
    assert report.max_score == 7
    assert report.at_max == 4
    assert report.cut_score == 7
    assert report.cut_inside_tie_block is True


def test_tie_report_knows_when_the_cut_is_clean():
    ranked = [protein("a", 9), protein("b", 5), protein("c", 1)]
    report = tie_report(ranked, k=1)
    assert report.cut_inside_tie_block is False
    assert report.distribution == {"9": 1, "5": 1, "1": 1}


# ---------------------------------------------------------------------------- scoring

def build(tmp_path):
    """Two positives and one negative in the top 3, one positive out of reach, one absent."""
    truthset = truth(
        tmp_path,
        "wza\tpositive\tcapsule\tone\t\n"
        "mrkD\tpositive\tfimbriae\ttwo\t\n"
        "ureC\tpositive\turease\tthree\t\n"
        "skp\tnegative\tomp_biogenesis\tfour\t\n"
        "rmpA\tpositive\tcapsule\tfive\t\n",
    )
    proteins = [
        protein("p1", 9, gene="wza",
                breakdown=[{"points": 2, "reason": "Host-interaction mechanism: x"}]),
        protein("p2", 8, gene="mrkD",
                breakdown=[{"points": 3, "reason": "Virulence factor (VFDB / Victors)"}]),
        protein("p3", 7, gene="skp",
                breakdown=[{"points": 1, "reason": "Predicted surface or secreted (x)"}]),
        protein("p4", 1, gene="ureC"),
        protein("p5", 0, gene=None, product="hypothetical protein"),
    ]
    matches, absent = match_proteins(proteins, truthset)
    return truthset, proteins, matches, absent


def test_precision_is_over_the_labelled_subset_not_over_k(tmp_path):
    truthset, proteins, matches, absent = build(tmp_path)
    result = evaluate(proteins, truthset, matches, absent, k=4, baselines=())
    # top 4 = wza, mrkD, skp, ureC -> 3 positive, 1 negative, 0 unlabelled
    assert (result.measured.positives, result.measured.negatives) == (3, 1)
    assert result.measured.precision == pytest.approx(0.75)
    assert result.measured.labelled_coverage == pytest.approx(1.0)


def test_unlabelled_proteins_are_not_counted_against_precision(tmp_path):
    truthset, proteins, matches, absent = build(tmp_path)
    result = evaluate(proteins, truthset, matches, absent, k=5, baselines=())
    assert result.measured.unlabelled == 1
    assert result.measured.precision == pytest.approx(0.75)  # unchanged by the hypothetical
    assert result.measured.labelled_coverage == pytest.approx(0.8)


def test_recall_denominator_excludes_symbols_absent_from_the_genome(tmp_path):
    truthset, proteins, matches, absent = build(tmp_path)
    result = evaluate(proteins, truthset, matches, absent, k=2, baselines=())
    # 4 curated positives, but rmpA is not in this genome, so the denominator is 3.
    assert result.positives_in_genome == 3
    assert result.measured.recall_denominator == 3
    assert result.measured.recall == pytest.approx(2 / 3)
    assert result.not_present == ["rmpA"]


def test_k_is_clamped_to_the_number_of_scored_proteins(tmp_path):
    truthset, proteins, matches, absent = build(tmp_path)
    result = evaluate(proteins, truthset, matches, absent, k=10_000, baselines=())
    assert result.measured.k == len(proteins)


def test_baselines_are_measured_at_the_same_cut(tmp_path):
    truthset, proteins, matches, absent = build(tmp_path)
    result = evaluate(proteins, truthset, matches, absent, k=3)
    assert set(result.baselines) == set(BASELINES)
    assert all(score.k == 3 for score in result.baselines.values())


def test_per_class_counts_separate_genome_from_top_k(tmp_path):
    truthset, proteins, matches, absent = build(tmp_path)
    result = evaluate(proteins, truthset, matches, absent, k=2, baselines=())
    assert result.by_class["urease"] == {"in_genome": 1, "in_top_k": 0}
    assert result.by_class["capsule"] == {"in_genome": 1, "in_top_k": 1}


# ------------------------------------------------------------------------ leakage audit

def test_leakage_audit_splits_curated_from_text_evidence(tmp_path):
    truthset, proteins, matches, absent = build(tmp_path)
    result = evaluate(proteins, truthset, matches, absent, k=3, baselines=())
    leak = result.leakage
    assert leak.labelled_in_top_k == 3
    assert leak.curated_points == 3      # mrkD's VFDB hit
    assert leak.text_points == 3         # wza mechanism 2 + skp surface 1
    assert leak.curated_backed_positives == ["mrkD"]
    assert leak.text_only_positives == ["wza"]
    assert leak.text_share == pytest.approx(0.5)


def test_leakage_text_share_is_none_when_nothing_scored(tmp_path):
    truthset, proteins, matches, absent = build(tmp_path)
    result = evaluate(proteins, truthset, matches, absent, k=1, baselines=())
    # wza alone: text only, so curated is zero but the share is defined.
    assert result.leakage.text_share == pytest.approx(1.0)


def test_evaluation_serialises(tmp_path):
    truthset, proteins, matches, absent = build(tmp_path)
    payload = evaluate(proteins, truthset, matches, absent, k=3).to_dict()
    assert payload["weights_fingerprint"] == weights_fingerprint()
    assert payload["leakage"]["text_share"] == pytest.approx(0.5)
    assert {row["symbol"] for row in payload["rows"]} == {"wza", "mrkD", "skp", "ureC"}

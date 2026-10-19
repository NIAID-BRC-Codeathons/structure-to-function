"""Leave-one-out ablation of the triage score."""

from __future__ import annotations

import pytest

from s2f.evaluation.ablation import ablate
from s2f.evaluation.rankings import TRIAGE_COMPONENTS
from s2f.evaluation.truthset import Match
from tests.test_eval_rankings import with_triage


def match(fid, label, symbol=None, klass="test"):
    return Match(feature_id=fid, symbol=symbol or fid, label=label, klass=klass,
                 tier="gene", gene=None, product=None)


def population():
    """Scores, and what each ablation does to the top 2 — worked out by hand:

        full                 a .56  c .38  d .34  b .32  e .02   -> {a, c}
        drop pdb_evidence    c .30  d .30  a .20  b .00  e .00   -> {c, d}   churn 2
        drop virulence_amr   c .38  a .36  d .34  b .32  e .02   -> {c, a}   churn 0
        drop annotation_gap  a .56  b .32  c .08  d .04  e .02   -> {a, b}   churn 1

    So `annotation_gap` at weight 0.30 lifts two hypothetical proteins over a real PDB hit,
    which is the design question the ablation exists to ask.
    """
    return [
        with_triage("a", {"pdb_evidence": 0.9, "virulence_amr": 1.0}, length=400),
        with_triage("b", {"pdb_evidence": 0.8}, length=300),
        with_triage("c", {"pdb_evidence": 0.2, "annotation_gap": 1.0}, length=200),
        with_triage("d", {"pdb_evidence": 0.1, "annotation_gap": 1.0}, length=100),
        with_triage("e", {"pdb_evidence": 0.05}, length=50),
    ]


def test_it_refuses_a_report_m2_never_scored():
    with pytest.raises(ValueError, match="no protein carries triage.components"):
        ablate([{"feature_id": "a"}], k=2)


def test_every_component_is_swept():
    result = ablate(population(), k=2)
    assert {effect.component for effect in result.effects} == set(TRIAGE_COMPONENTS)
    assert result.score == "triage" and result.k == 2


def test_effects_are_ordered_by_how_much_they_move_the_selection():
    result = ablate(population(), k=2)
    churn = [effect.churn_at_k for effect in result.effects]
    assert churn == sorted(churn, reverse=True)


def test_dropping_the_dominant_component_churns_the_selection():
    result = ablate(population(), k=2)
    pdb = next(e for e in result.effects if e.component == "pdb_evidence")
    assert pdb.churn_at_k > 0
    assert pdb.weight == pytest.approx(0.4)
    assert pdb.nonzero_proteins == 5


def test_a_component_that_never_fires_is_flagged_rather_than_called_inert():
    """Zero everywhere usually means the provider never ran, not that it does not matter."""
    result = ablate(population(), k=2)
    assert "membrane_penalty" in result.never_fired
    membrane = next(e for e in result.effects if e.component == "membrane_penalty")
    assert membrane.nonzero_proteins == 0 and membrane.churn_at_k == 0
    assert any("never fired" in note for note in result.notes)


def test_a_component_that_fires_but_changes_nothing_is_reported_separately():
    result = ablate(population(), k=2)
    viru = next(e for e in result.effects if e.component == "virulence_amr")
    assert viru.nonzero_proteins == 1          # it fires
    assert viru.churn_at_k == 0                # but only reorders inside the cut
    assert any("changed nothing" in note for note in result.notes)


def test_annotation_gap_lifts_hypothetical_proteins_over_a_real_pdb_hit():
    """Weight 0.30 for being uncharacterised is the ablation's most interesting question."""
    result = ablate(population(), k=2)
    gap = next(e for e in result.effects if e.component == "annotation_gap")
    assert gap.nonzero_proteins == 2
    assert gap.churn_at_k == 1                 # removing it drops 'c' for 'b'


def test_optional_components_carry_their_provider_count():
    proteins = population()
    for protein in proteins:
        protein["triage"]["components_available"] = {"surface_bonus": True,
                                                     "membrane_penalty": False}
    result = ablate(proteins, k=2)
    surface = next(e for e in result.effects if e.component == "surface_bonus")
    membrane = next(e for e in result.effects if e.component == "membrane_penalty")
    assert surface.optional is True and surface.providers_available == 5
    assert membrane.optional is True and membrane.providers_available == 0
    pdb = next(e for e in result.effects if e.component == "pdb_evidence")
    assert pdb.optional is False and pdb.providers_available is None


def test_precision_deltas_are_measured_against_the_full_score():
    truth = {
        "a": match("a", "positive"),
        "b": match("b", "positive"),
        "c": match("c", "negative"),
        "d": match("d", "negative"),
    }
    result = ablate(population(), truth_matches=truth, k=2)
    assert result.reference_truth is not None
    assert result.reference_truth.precision == pytest.approx(0.5)   # top 2 = a(+), c(-)
    pdb = next(e for e in result.effects if e.component == "pdb_evidence")
    assert pdb.truth_delta_precision is not None
    membrane = next(e for e in result.effects if e.component == "membrane_penalty")
    assert membrane.truth_delta_precision == 0      # it changed nothing, so nor did precision


def test_outcome_deltas_are_absent_when_m3_has_not_run():
    result = ablate(population(), truth_matches={"a": match("a", "positive")}, k=2)
    assert result.reference_outcome is None
    assert all(effect.outcome_delta_precision is None for effect in result.effects)


def test_outcome_deltas_appear_when_m3_labels_exist():
    outcomes = {
        "a": match("a", "positive", klass="dockable"),
        "c": match("c", "negative", klass="no_structure"),
    }
    result = ablate(population(), outcome_matches=outcomes, k=2)
    assert result.reference_outcome is not None
    assert result.reference_outcome.precision == pytest.approx(0.5)


def test_load_bearing_lists_only_components_that_move_the_cut():
    result = ablate(population(), k=2)
    movers = result.load_bearing()
    assert movers and all(effect.churn_at_k > 0 for effect in movers)
    assert "membrane_penalty" not in {effect.component for effect in movers}


def test_ablation_serialises():
    payload = ablate(population(), k=2).to_dict()
    assert payload["score"] == "triage"
    assert len(payload["effects"]) == len(TRIAGE_COMPONENTS)
    assert payload["weights"]["pdb_evidence"] == pytest.approx(0.4)
    assert payload["effects"][0]["truth"] is None


def test_the_reference_ranking_is_never_mutated():
    proteins = population()
    before = [p["triage"]["score"] for p in proteins]
    ablate(proteins, k=2)
    assert [p["triage"]["score"] for p in proteins] == before

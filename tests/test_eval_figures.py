"""Optional evaluation figures.

Same contract as M1's: matplotlib is not a hard dependency and a missing PNG never fails a
run. These tests assert the contract, not the pixels.
"""

from __future__ import annotations

import builtins
import json

import pytest

from s2f.evaluation.__main__ import main
from s2f.evaluation.figures import save_figures
from s2f.evaluation.rankings import is_gap_like, population_stages
from tests.test_eval_rankings import with_triage

matplotlib = pytest.importorskip("matplotlib", reason="figures are optional")

STAGES = [("in genome", 100, 40), ("survived M1 seq-cap", 70, 15)]
ABLATION = {
    "k": 10,
    "effects": [
        {"component": "pdb_evidence", "weight": 0.4, "churn_at_k": 6, "nonzero_proteins": 90},
        {"component": "annotation_gap", "weight": 0.3, "churn_at_k": 4, "nonzero_proteins": 20},
        {"component": "essential", "weight": 0.15, "churn_at_k": 0, "nonzero_proteins": 5},
        {"component": "surface_bonus", "weight": 0.1, "churn_at_k": 0, "nonzero_proteins": 0},
    ],
}
OUTCOME = {"base_rate": 0.993, "measured": {"precision": 1.0}}


def test_every_figure_is_written_when_its_input_is_present(tmp_path):
    written = save_figures(
        tmp_path, stages=STAGES,
        ranks=([1, 2, 3, 4], [4, 3, 1, 2]), shared_at_k=1, spearman=0.37, k=2,
        ablation=ABLATION, outcome=OUTCOME,
    )
    names = {path.stem for path in written}
    assert names == {"seq_cap_loss", "score_disagreement", "triage_ablation", "outcome_lift"}
    assert {path.suffix for path in written} == {".png", ".svg"}
    assert all(path.stat().st_size > 0 for path in written)


def test_each_input_is_independent(tmp_path):
    """A run with only M1 still gets the figure M1 supports, and no others."""
    written = save_figures(tmp_path, stages=STAGES)
    assert {path.stem for path in written} == {"seq_cap_loss"}


def test_nothing_is_written_when_there_is_nothing_to_draw(tmp_path):
    assert save_figures(tmp_path) == []


def test_mismatched_rank_vectors_are_skipped_not_plotted(tmp_path):
    written = save_figures(tmp_path, ranks=([1, 2, 3], [1, 2]))
    assert written == []


def test_a_missing_matplotlib_returns_empty_rather_than_raising(tmp_path, monkeypatch):
    real_import = builtins.__import__

    def refuse(name, *args, **kwargs):
        if name == "matplotlib" or name.startswith("matplotlib."):
            raise ImportError("no matplotlib")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", refuse)
    assert save_figures(tmp_path, stages=STAGES, ablation=ABLATION) == []


# ------------------------------------------------------------------- the seam definition

def test_gap_like_matches_the_vocabulary_an_annotation_uses_to_admit_ignorance():
    assert is_gap_like({"product": "hypothetical protein"})
    assert is_gap_like({"product": "Uncharacterized ferredoxin-like protein YfhL"})
    assert is_gap_like({"product": "Putative oxidoreductase YncB"})
    assert is_gap_like({"product": "DUF1234 domain protein"})
    assert not is_gap_like({"product": "Type 1 fimbrial adhesin FimH"})
    assert not is_gap_like({"product": ""})
    assert not is_gap_like({})


def test_population_stages_tracks_the_gap_like_share_through_the_cap():
    proteins = (
        [with_triage(f"s{i}", {}, product="hypothetical protein") for i in range(3)]
        + [with_triage(f"k{i}", {}, product="Type 1 fimbrial adhesin FimH") for i in range(2)]
        # never reached M2: no triage record at all
        + [{"feature_id": f"d{i}", "product": "hypothetical protein"} for i in range(5)]
    )
    stages = population_stages(proteins)
    assert stages[0] == ("in genome", 10, 8)
    assert stages[1] == ("survived M1 seq-cap", 5, 3)


def test_population_stages_adds_the_selection_when_m2_made_one():
    proteins = [
        with_triage("a", {}, selected=True, product="hypothetical protein"),
        with_triage("b", {}, selected=False, product="FimH"),
    ]
    stages = population_stages(proteins)
    assert [row[0] for row in stages][-1] == "selected by triage"
    assert stages[-1] == ("selected by triage", 1, 1)


# --------------------------------------------------------------------------------- CLI

def test_the_cli_writes_figures_only_when_asked(tmp_path):
    plain, drawn = tmp_path / "plain", tmp_path / "drawn"
    assert main(["--dry-run", "--run", str(plain), "--score", "triage"]) == 0
    assert not (plain / "eval" / "figures").exists()
    assert main(["--dry-run", "--run", str(drawn), "--score", "triage", "--figures"]) == 0
    assert (drawn / "eval" / "figures" / "seq_cap_loss.png").is_file()
    assert (drawn / "eval" / "figures" / "triage_ablation.svg").is_file()

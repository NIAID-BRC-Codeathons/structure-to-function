"""Optional evaluation figures.

Same contract as M1's: matplotlib is not a hard dependency and a missing PNG never fails a
run. These tests assert the contract, not the pixels.
"""

from __future__ import annotations

import builtins
import json

import pytest

from s2f.evaluation.__main__ import main
from s2f.evaluation.figures import (
    DOUBLE_COLUMN, PUBLICATION_FORMATS, SINGLE_COLUMN, save_figures,
)
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


ALL_INPUTS = dict(
    stages=STAGES, ranks=([1, 2, 3, 4], [4, 3, 1, 2]), shared_at_k=1, spearman=0.37,
    k=2, ablation=ABLATION, outcome=OUTCOME,
)


def test_every_panel_plus_a_composite_and_captions(tmp_path):
    written = save_figures(tmp_path, **ALL_INPUTS)
    names = {path.stem for path in written}
    assert names == {"seq_cap_loss", "score_disagreement", "triage_ablation",
                     "outcome_lift", "figure_composite", "captions"}
    assert {p.suffix for p in written if p.suffix != ".md"} == {".pdf", ".svg", ".png"}
    assert all(path.stat().st_size > 0 for path in written)


def test_each_input_is_independent(tmp_path):
    """A run with only M1 gets the figure M1 supports, and no composite for one panel."""
    written = save_figures(tmp_path, stages=STAGES)
    assert {path.stem for path in written} == {"seq_cap_loss", "captions"}


def test_slide_style_puts_the_headline_in_the_figure_and_skips_pdf(tmp_path):
    written = save_figures(tmp_path, style="slide", **ALL_INPUTS)
    assert {p.suffix for p in written if p.suffix != ".md"} == {".png", ".svg"}
    assert "figure_composite" not in {p.stem for p in written}


def test_unknown_style_is_refused(tmp_path):
    with pytest.raises(ValueError, match="unknown figure style"):
        save_figures(tmp_path, style="poster", stages=STAGES)


# ------------------------------------------------- what makes a figure submittable
# matplotlib's defaults fail both of these, and a journal rejects the file for it.

def test_pdf_text_is_embedded_truetype_never_type3(tmp_path):
    save_figures(tmp_path, **ALL_INPUTS)
    for path in tmp_path.glob("*.pdf"):
        raw = path.read_bytes()
        assert b"/Type3" not in raw, f"{path.name} embeds Type 3 fonts"
        assert b"FontFile2" in raw, f"{path.name} does not embed a TrueType face"


def test_svg_text_stays_live_rather_than_outlined(tmp_path):
    save_figures(tmp_path, **ALL_INPUTS)
    for path in tmp_path.glob("*.svg"):
        assert b"<text" in path.read_bytes(), f"{path.name} has no live text"


def test_panels_are_sized_to_journal_columns(tmp_path):
    import re

    save_figures(tmp_path, **ALL_INPUTS)
    def width_mm(name):
        raw = (tmp_path / name).read_bytes()
        box = re.search(rb"/MediaBox\s*\[\s*0\s+0\s+([\d.]+)", raw)
        return float(box.group(1)) / 72 * 25.4
    # tight bbox trims a little, so assert the band rather than the exact millimetre
    assert width_mm("seq_cap_loss.pdf") <= SINGLE_COLUMN * 25.4 + 1
    assert width_mm("figure_composite.pdf") <= DOUBLE_COLUMN * 25.4 + 1
    assert width_mm("figure_composite.pdf") > width_mm("seq_cap_loss.pdf")


def test_captions_are_numbered_and_carry_the_numbers(tmp_path):
    """A legend belongs outside the artwork, and must stand on its own."""
    save_figures(tmp_path, **ALL_INPUTS)
    captions = (tmp_path / "captions.md").read_text(encoding="utf-8")
    assert captions.count("**Figure ") == 4
    assert "Figure 1." in captions and "Figure 4." in captions
    assert "99.3%" in captions          # the base rate
    assert "rho = +0.370" in captions   # the correlation, as passed in
    assert "data gap rather than a measured null" in captions


def test_nothing_is_written_when_there_is_nothing_to_draw(tmp_path):
    assert save_figures(tmp_path) == []
    assert not (tmp_path / "captions.md").exists()


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
    figures = drawn / "eval" / "figures"
    assert (figures / "seq_cap_loss.pdf").is_file()
    assert (figures / "triage_ablation.svg").is_file()
    assert (figures / "captions.md").is_file()


def test_the_cli_can_still_produce_slide_figures(tmp_path):
    run_dir = tmp_path / "slides"
    assert main(["--dry-run", "--run", str(run_dir), "--score", "triage",
                 "--figures", "--figure-style", "slide"]) == 0
    figures = run_dir / "eval" / "figures"
    assert (figures / "seq_cap_loss.png").is_file()
    assert not (figures / "seq_cap_loss.pdf").exists()


def test_a_genome_under_the_cap_does_not_get_a_fabricated_loss_ratio(tmp_path):
    """Nothing dropped means no differential to quote; the caption must not invent one."""
    written = save_figures(tmp_path, stages=[("in genome", 900, 300),
                                             ("survived M1 seq-cap", 900, 300)])
    assert written
    captions = (tmp_path / "captions.md").read_text(encoding="utf-8")
    assert "did not fire on this genome" in captions
    assert "fold higher loss rate" not in captions

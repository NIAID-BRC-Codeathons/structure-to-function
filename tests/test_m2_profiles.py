"""Scoring profiles: the triage policy is a choice, recorded with the ranking (#93)."""

from __future__ import annotations

import csv
import json

import pytest

from s2f.m2_triage.__main__ import main as m2_main
from s2f.m2_triage.profiles import DEFAULT_PROFILE, PROFILES, load_profile
from s2f.m2_triage.rescore import main as rescore_main
from s2f.m2_triage.score import WEIGHTS, WEIGHTS_VERSION


def _rows(path):
    with open(path, encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_the_default_profile_is_exactly_the_documented_weights() -> None:
    profile = PROFILES[DEFAULT_PROFILE]
    assert profile.weights == WEIGHTS
    assert profile.version == WEIGHTS_VERSION
    assert profile.disqualify_human_homologs


def test_every_profile_names_every_component() -> None:
    for profile in PROFILES.values():
        assert set(profile.weights) == set(WEIGHTS), profile.name


def test_a_weights_file_must_name_every_component(tmp_path) -> None:
    partial = tmp_path / "w.json"
    partial.write_text(json.dumps({"weights": {"pdb_evidence": 1.0}}))
    with pytest.raises(ValueError, match="missing"):
        load_profile(weights_file=str(partial))
    typo = tmp_path / "t.json"
    typo.write_text(json.dumps({"weights": {**WEIGHTS, "pdb_evidance": 0.1}}))
    with pytest.raises(ValueError, match="unknown"):
        load_profile(weights_file=str(typo))


def test_a_valid_weights_file_loads_with_its_source(tmp_path) -> None:
    path = tmp_path / "mine.json"
    path.write_text(json.dumps({"weights": {**WEIGHTS, "annotation_gap": 0.5},
                                "disqualify_human_homologs": False}))
    profile = load_profile(weights_file=str(path))
    assert profile.name == "mine" and profile.source == str(path)
    assert profile.weights["annotation_gap"] == 0.5 and not profile.disqualify_human_homologs


def test_run_json_records_the_profile_that_produced_the_ranking(tmp_path) -> None:
    assert m2_main(["--run", str(tmp_path / "a"), "--dry-run"]) == 0
    assert m2_main(["--run", str(tmp_path / "b"), "--dry-run", "--profile", "evidence-only"]) == 0
    a = json.loads((tmp_path / "a/m2_pdb/run.json").read_text())
    b = json.loads((tmp_path / "b/m2_pdb/run.json").read_text())
    assert a["profile"]["name"] == DEFAULT_PROFILE and a["weights"] == WEIGHTS
    assert b["profile"]["name"] == "evidence-only"
    assert b["weights"] == PROFILES["evidence-only"].weights


def test_a_bad_profile_fails_before_any_work(tmp_path) -> None:
    bad = tmp_path / "bad.json"
    bad.write_text("{}")
    assert m2_main(["--run", str(tmp_path / "r"), "--dry-run", "--weights", str(bad)]) == 2
    assert not (tmp_path / "r/m2_pdb/run.json").exists()


def test_rescoring_under_the_default_profile_reproduces_m2_exactly(tmp_path) -> None:
    assert m2_main(["--run", str(tmp_path), "--dry-run"]) == 0
    assert rescore_main(["--run", str(tmp_path)]) == 0
    original = {r["feature_id"]: r for r in _rows(tmp_path / "m2_pdb/proteins.tsv")}
    rescored = _rows(tmp_path / f"m2_pdb/rescored/{DEFAULT_PROFILE}.tsv")
    assert len(rescored) == len(original)
    for row in rescored:
        assert row["rank"] == original[row["feature_id"]]["rank"]
        assert float(row["score"]) == pytest.approx(float(original[row["feature_id"]]["triage_score"]))
        assert row["selected"] == original[row["feature_id"]]["selected"]


def test_rescoring_never_touches_the_original_ranking(tmp_path) -> None:
    assert m2_main(["--run", str(tmp_path), "--dry-run"]) == 0
    before = (tmp_path / "m2_pdb/proteins.tsv").read_text()
    assert rescore_main(["--run", str(tmp_path), "--profile", "evidence-only"]) == 0
    assert (tmp_path / "m2_pdb/proteins.tsv").read_text() == before
    manifest = json.loads((tmp_path / "m2_pdb/rescored/evidence-only.json").read_text())
    assert manifest["profile"]["name"] == "evidence-only"


def test_evidence_only_does_not_disqualify_human_homologs(tmp_path) -> None:
    assert m2_main(["--run", str(tmp_path), "--dry-run"]) == 0
    assert rescore_main(["--run", str(tmp_path), "--profile", "evidence-only"]) == 0
    assert not any(r["disqualified"] for r in
                   _rows(tmp_path / "m2_pdb/rescored/evidence-only.tsv"))


def test_old_output_without_every_component_is_refused(tmp_path) -> None:
    (tmp_path / "m2_pdb").mkdir()
    (tmp_path / "m2_pdb/proteins.tsv").write_text("rank\tfeature_id\tpdb_evidence\n1\tx\t0.5\n")
    assert rescore_main(["--run", str(tmp_path)]) == 2

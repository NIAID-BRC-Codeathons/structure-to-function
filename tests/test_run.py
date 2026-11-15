"""The pipeline runner forwards what it claims to, and says when evidence is missing (#81)."""

from __future__ import annotations

import json

from s2f.m2_triage.__main__ import build_parser as m2_parser
from s2f.run import build_parser, missing_evidence, stage_args

M2_FLAGS = ["--full-evidence", "--human-proteome", "data/h.fasta", "--human-min-coverage", "0.5",
            "--essentiality-keyword", "Mycoplasma", "--essentiality-limit", "7",
            "--essentiality-min-identity", "35", "--essentiality-min-coverage", "0.6"]


def _m2_args(tmp_path, argv):
    args = build_parser().parse_args(["--run", str(tmp_path), *argv])
    return m2_parser().parse_args(stage_args("m2", args))


def test_every_m2_flag_the_runner_supports_arrives(tmp_path) -> None:
    got = _m2_args(tmp_path, M2_FLAGS)
    assert got.human_homology and got.essentiality
    assert got.human_proteome == "data/h.fasta"
    assert got.human_min_coverage == 0.5
    assert got.essentiality_keyword == "Mycoplasma"
    assert got.essentiality_limit == 7
    assert got.essentiality_min_identity == 35.0
    assert got.essentiality_min_coverage == 0.6


def test_each_evidence_switch_works_on_its_own(tmp_path) -> None:
    got = _m2_args(tmp_path, ["--human-homology"])
    assert got.human_homology and not got.essentiality
    got = _m2_args(tmp_path, ["--essentiality"])
    assert got.essentiality and not got.human_homology


def test_unset_tuning_flags_leave_m2_defaults_alone(tmp_path) -> None:
    assert _m2_args(tmp_path, []) == m2_parser().parse_args(
        ["--run", str(tmp_path), "--report"])


def _manifest(tmp_path, **sections):
    out = tmp_path / "m2_pdb"
    out.mkdir()
    manifest = {"weights": {"essential": 0.15, "human_homolog_penalty": -0.25,
                            "membrane_penalty": -0.15},
                "human_homology": {"enabled": False}, "essentiality": {"enabled": False},
                "components_available": {"membrane_penalty": {"measured_for": 0, "of": 5}},
                **sections}
    (out / "run.json").write_text(json.dumps(manifest))


def test_a_run_without_external_evidence_says_so(tmp_path) -> None:
    _manifest(tmp_path)
    lines = missing_evidence(str(tmp_path))
    assert any(line.startswith("human_homolog_penalty") and "--human-homology" in line
               for line in lines)
    assert any(line.startswith("essential ") and "--essentiality" in line for line in lines)
    assert any(line.startswith("membrane_penalty") for line in lines)


def test_a_failed_provider_reports_its_error_and_a_full_run_is_quiet(tmp_path) -> None:
    _manifest(tmp_path, human_homology={"enabled": False, "error": "diamond not found"},
              essentiality={"enabled": True},
              components_available={"membrane_penalty": {"measured_for": 3, "of": 5}})
    assert missing_evidence(str(tmp_path)) == [
        "human_homolog_penalty (weight -0.25) failed: diamond not found"]


def test_no_m2_manifest_means_no_warning(tmp_path) -> None:
    assert missing_evidence(str(tmp_path)) == []

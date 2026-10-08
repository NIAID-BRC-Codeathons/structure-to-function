"""Sequence-only runs record what nobody looked at as not available, never as zero (#90)."""

from __future__ import annotations

import json
from pathlib import Path

from s2f.full_report.render import build_page
from s2f.m1_genome.__main__ import main as m1_main
from s2f.m2_triage.availability import not_available
from s2f.run import missing_evidence

ROOT = Path(__file__).resolve().parents[1]
QUERY = ROOT / "fixtures/m1/fasta_query/query.faa"
GENOME_REPORT = ROOT / "fixtures/report.fixture.json"


def test_a_genome_route_has_nothing_unavailable() -> None:
    assert not_available("api", {}) == {}
    assert not_available("cga", {}) == {}
    assert not_available("", {}) == {}  # standalone M2: M1 did not run here


def test_the_fasta_route_marks_every_specialty_component() -> None:
    got = not_available("fasta", {})
    assert set(got) == {"virulence_amr", "drug_target", "essential", "human_homolog_penalty"}
    assert "--human-homology" in got["human_homolog_penalty"]


def test_a_provider_that_ran_makes_its_component_available() -> None:
    got = not_available("fasta", {"human_homology": True, "essentiality": True})
    assert set(got) == {"virulence_amr", "drug_target"}


def _m2_manifest(run_dir: Path, unavailable: dict) -> None:
    (run_dir / "m2_pdb").mkdir(parents=True)
    (run_dir / "m2_pdb/run.json").write_text(json.dumps({
        "weights": {"virulence_amr": 0.2, "drug_target": 0.15, "essential": 0.15,
                    "human_homolog_penalty": -0.25},
        "human_homology": {"enabled": False}, "essentiality": {"enabled": False},
        "not_available": unavailable}))


def test_the_runner_says_not_available_instead_of_bvbrc_rows(tmp_path) -> None:
    _m2_manifest(tmp_path, not_available("fasta", {}))
    lines = missing_evidence(str(tmp_path))
    assert len(lines) == 4
    assert all("not available" in line for line in lines)
    assert not any("BV-BRC's precomputed rows" in line for line in lines)


def test_a_genome_run_and_a_fasta_run_are_distinguishable_from_report_json(tmp_path) -> None:
    assert m1_main(["--run", str(tmp_path / "q"), "--from-fasta", str(QUERY)]) == 0
    fasta = json.loads((tmp_path / "q/report.json").read_text())
    genome = json.loads(GENOME_REPORT.read_text())
    assert fasta["genome"]["annotation_route"] == "fasta"
    assert genome["genome"].get("annotation_route") != "fasta"
    assert fasta["genome"]["not_available"] and not genome["genome"].get("not_available")


def test_the_full_report_shows_not_available_only_when_something_was(tmp_path) -> None:
    assert m1_main(["--run", str(tmp_path / "q"), "--from-fasta", str(QUERY)]) == 0
    fasta = json.loads((tmp_path / "q/report.json").read_text())
    page = build_page(fasta, {"not_available": not_available("fasta", {})})
    assert "Not available on this run" in page and "virulence_amr" in page
    assert "Not available on this run" not in build_page(json.loads(GENOME_REPORT.read_text()))

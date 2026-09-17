"""End-to-end check of the acceptance criterion: --dry-run works from fixtures, no network."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from s2f.m2_triage.__main__ import main
from s2f.m2_triage.score import TriageComponents, triage_score_from_components


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_dry_run_produces_a_ranked_package(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    assert main(["--dry-run", "--run", str(run_dir), "--top", "2"]) == 0

    out = run_dir / "m2_pdb"
    proteins = _read_tsv(out / "proteins.tsv")
    manifest = json.loads((out / "run.json").read_text())

    assert manifest["offline"] is True
    assert manifest["failures"]["search"] == []
    # Four fixture proteins, two of which share a sequence: three searches.
    assert manifest["counts"]["proteins_scored"] == 4
    assert manifest["counts"]["unique_sequences"] == 3
    assert [p["rank"] for p in proteins] == ["1", "2", "3", "4"]
    assert sum(1 for p in proteins if p["selected"] == "True") == 2
    assert all(p["reason"] for p in proteins)


def test_dry_run_ranks_the_sars_cov_2_rbd_first_with_a_6m0j_hit(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    main(["--dry-run", "--run", str(run_dir), "--top", "2"])
    proteins = _read_tsv(run_dir / "m2_pdb" / "proteins.tsv")

    top = proteins[0]
    assert top["feature_id"] == "fig|999999.1.peg.1"
    assert float(top["best_identity"]) == 1.0
    assert top["best_entity_id"].startswith("6M0J")
    assert top["retrieval_status"] == "found"


def test_dry_run_scores_are_reproducible_from_the_written_components(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    main(["--dry-run", "--run", str(run_dir)])

    for row in _read_tsv(run_dir / "m2_pdb" / "proteins.tsv"):
        components = TriageComponents(
            pdb_evidence=float(row["pdb_evidence"]),
            virulence_amr=float(row["virulence_amr"]),
            essential=float(row["essential"]),
            drug_target=float(row["drug_target"]),
            annotation_gap=float(row["annotation_gap"]),
            human_homolog_penalty=float(row["human_homolog_penalty"]),
        )
        assert triage_score_from_components(components) == float(row["triage_score"])


def test_dry_run_separates_no_hit_proteins_and_penalizes_the_human_homolog(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    main(["--dry-run", "--run", str(run_dir)])
    out = run_dir / "m2_pdb"

    no_hits = _read_tsv(out / "no_pdb_hit.tsv")
    assert [row["feature_id"] for row in no_hits] == ["fig|999999.1.peg.3"]
    assert no_hits[0]["retrieval_status"] == "no-hit"

    by_id = {row["feature_id"]: row for row in _read_tsv(out / "proteins.tsv")}
    lysozyme = by_id["fig|999999.1.peg.2"]
    assert float(lysozyme["human_homolog_penalty"]) == 0.6
    assert lysozyme["holo_homolog"] == "True"


def test_limit_restricts_work_without_breaking_the_manifest(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    main(["--dry-run", "--run", str(run_dir), "--limit", "1"])
    manifest = json.loads((run_dir / "m2_pdb" / "run.json").read_text())

    assert manifest["counts"]["proteins_loaded"] == 4
    assert manifest["counts"]["proteins_scored"] == 1
    assert manifest["parameters"]["limit"] == 1

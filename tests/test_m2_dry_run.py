"""End-to-end check of the acceptance criterion: --dry-run works from fixtures, no network."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from s2f.m2_triage.__main__ import main
from s2f.m2_triage.score import WEIGHTS, TriageComponents, triage_score_from_components


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


def test_dry_run_finds_the_sars_cov_2_rbd_against_6m0j(tmp_path: Path) -> None:
    """End-to-end sanity: the RBD fixture must recover its own structure at identity 1.0.

    It is no longer rank 1. Since `ligandable_homolog` landed, the lysozyme fixture outranks it:
    lysozyme's best hit has a real ligand bound, while 6M0J's only non-polymer components are a
    glycan, a metal and an additive — correctly excluded from `ligands`. The pipeline check here
    is the hit, not the position.
    """
    run_dir = tmp_path / "run"
    main(["--dry-run", "--run", str(run_dir), "--top", "2"])
    proteins = {row["feature_id"]: row for row in _read_tsv(run_dir / "m2_pdb" / "proteins.tsv")}

    rbd = proteins["fig|999999.1.peg.1"]
    assert float(rbd["best_identity"]) == 1.0
    assert rbd["best_entity_id"].startswith("6M0J")
    assert rbd["retrieval_status"] == "found"
    assert rbd["selected"] == "True"
    assert float(rbd["ligandable_homolog"]) == 0.0  # NAG/ZN/CL are not functional ligands


def test_dry_run_scores_are_reproducible_from_the_written_components(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    main(["--dry-run", "--run", str(run_dir)])

    for row in _read_tsv(run_dir / "m2_pdb" / "proteins.tsv"):
        # Every weighted component, read back from the TSV. Building this from WEIGHTS rather
        # than a hand-written list means a new component cannot be added without the
        # reproducibility claim covering it.
        components = TriageComponents(**{name: float(row[name]) for name in WEIGHTS})
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

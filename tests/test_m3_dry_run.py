"""Offline end-to-end test for the M3 existing-structure collector."""

from __future__ import annotations

import json
from pathlib import Path

from s2f.common.io import read_report
from s2f.common.schema import validate
from s2f.m3_fold.__main__ import main


def test_dry_run_collects_fixture_structures_without_network(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"

    assert main(["--dry-run", "--run", str(run_dir)]) == 0

    report = read_report(run_dir)
    structures = report["structures"]
    assert len(structures) == 3
    assert {record["source"] for record in structures} == {"pdb"}
    assert all(record["collection_status"] == "collected" for record in structures)
    assert all(record["usable_for_docking"] is None for record in structures)
    assert all((run_dir / record["path"]).exists() for record in structures)
    assert all("elapsed_seconds" in record for record in structures)
    assert all(record["method"] == "download_existing_structure" for record in structures)
    assert all("params" in record for record in structures)
    validate("structures", structures)

    manifest = json.loads((run_dir / "m3_fold" / "run.json").read_text())
    assert manifest["module"] == "m3_fold"
    assert manifest["command"][-2:] == ["--run", str(run_dir)]
    assert manifest["counts"]["selected"] == 3
    assert manifest["counts"]["collected"] == 3
    assert manifest["counts"]["failed"] == 0
    assert len(manifest["input_report_sha256"]) == 64
    assert len(manifest["output_report_sha256"]) == 64
    assert manifest["input_report_sha256"] != manifest["output_report_sha256"]
    assert manifest["environment"]["python"]
    assert "commit" in manifest["git"]

    # M3 replaces only its own section; the fixture's other module sections survive.
    assert "proteins" in report
    assert "kg" in report
    assert "ligands" in report

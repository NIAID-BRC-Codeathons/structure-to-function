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
    assert len(structures) == 4
    assert {record["source"] for record in structures} == {"pdb"}
    assert all(record["collection_status"] == "collected" for record in structures)
    assert all(record["usable_for_docking"] is True for record in structures)
    assert all(record["quality_gate"]["provisional"] is True for record in structures)
    assert all((run_dir / record["path"]).exists() for record in structures)
    assert all("elapsed_seconds" in record for record in structures)
    assert all(record["method"] == "download_existing_structure" for record in structures)
    assert all("params" in record for record in structures)
    validate("structures", structures)

    manifest = json.loads((run_dir / "m3_fold" / "run.json").read_text())
    assert manifest["module"] == "m3_fold"
    assert manifest["command"][-2:] == ["--run", str(run_dir)]
    assert manifest["counts"]["selected"] == 4
    assert manifest["counts"]["collected"] == 4
    assert manifest["counts"]["usable_for_docking"] == 4
    assert manifest["counts"]["not_usable_for_docking"] == 0
    assert manifest["counts"]["quality_pending"] == 0
    assert manifest["parameters"]["confidence_gate_applied"] is True
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


def test_the_dry_run_defines_a_site_from_the_holo_fixture(tmp_path: Path) -> None:
    """#14, offline: 6MU0 binds ribose-5-phosphate, and its site reaches report.json."""
    run_dir = tmp_path / "run"
    assert main(["--dry-run", "--run", str(run_dir)]) == 0

    structures = read_report(run_dir)["structures"]
    holo = [record for record in structures if record["accession"] == "6MU0_1"]
    assert len(holo) == 1
    pockets = holo[0]["pockets"]
    assert len(pockets) == 1
    assert pockets[0]["component"] == "5RP"
    assert pockets[0]["aligned_to_known_site"] is True
    assert len(pockets[0]["residues"]) >= 5
    validate("structures", structures)

    # The three apo proteins get no pocket, and each says why rather than looking unexamined.
    apo = [record for record in structures if not record["pockets"]]
    assert len(apo) == 3
    assert all(record["pocket_search"]["status"] == "no_pocket" for record in apo)
    assert all("P2Rank" in record["pocket_search"]["reason"] for record in apo)

    counts = json.loads((run_dir / "m3_fold" / "run.json").read_text())["counts"]
    assert counts["with_pocket"] == 1
    assert counts["usable_pockets"] == 1
    assert counts["cofactor_required_pockets"] == 0


def test_the_dry_run_prepares_a_receptor_and_a_box(tmp_path: Path) -> None:
    """#14 part B, offline and with no preparation tool required."""
    run_dir = tmp_path / "run"
    assert main(["--dry-run", "--run", str(run_dir), "--no-receptors"]) == 0
    assert not (run_dir / "receptors").exists()

    run_dir = tmp_path / "with"
    assert main(["--dry-run", "--run", str(run_dir)]) == 0

    holo = [r for r in read_report(run_dir)["structures"] if r["accession"] == "6MU0_1"][0]
    receptor = holo["receptor"]
    assert receptor["chains"] == ["A"]
    assert receptor["atoms"] > 100
    # The feature id ends in a number, which must survive the filename.
    assert receptor["pdb_path"].endswith("peg.463.pdb")
    assert (run_dir / receptor["pdb_path"]).exists()
    # Fidelity depends on what is installed, so assert the honesty rule, not the tool.
    assert receptor["fidelity"] in {"full", "pdb_only", "reduced"}
    if receptor["fidelity"] == "pdb_only":
        assert receptor["pdbqt_path"] is None

    box = holo["pockets"][0]["box"]
    assert {"center_x", "size_x", "atoms"} <= set(box)

    counts = json.loads((run_dir / "m3_fold" / "run.json").read_text())["counts"]
    assert counts["receptors_prepared"] == 1
    assert counts["search_boxes"] == 1

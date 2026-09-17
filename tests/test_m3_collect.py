"""Existing-structure collection tests for M3 issue 13."""

from __future__ import annotations

from pathlib import Path

from s2f.common.schema import validate
from s2f.m3_fold.collect import collect_existing, choose_existing


def _protein(
    feature_id: str,
    *,
    selected: bool = True,
    experimental: dict | None = None,
    predicted: dict | None = None,
) -> dict:
    return {
        "feature_id": feature_id,
        "flags": {
            "experimental_homolog": experimental,
            "predicted_model": predicted,
        },
        "triage": {"selected": selected},
    }


EXPERIMENTAL = {
    "pdb_entity": "1abc_2",
    "identity": 0.91,
    "coverage": 0.95,
    "resolution": 2.1,
    "method": "X-ray",
    "holo": True,
    "ligands": ["ATP"],
}

PREDICTED = {
    "source": "afdb",
    "entry_id": "AF-P12345-F1",
    "mean_plddt": 88.4,
    "coverage": 1.0,
    "url": "https://alphafold.ebi.ac.uk/files/AF-P12345-F1-model_v6.cif",
}


def test_experimental_structure_is_preferred_over_afdb() -> None:
    candidate = choose_existing(
        _protein("fig|1.1.peg.1", experimental=EXPERIMENTAL, predicted=PREDICTED)
    )

    assert candidate is not None
    assert candidate.source == "pdb"
    assert candidate.accession == "1abc_2"
    assert candidate.url.endswith("/1ABC.cif")
    assert candidate.confidence == 2.1
    assert candidate.experimental_homolog["holo"] is True


def test_afdb_is_used_when_there_is_no_experimental_structure() -> None:
    candidate = choose_existing(_protein("fig|1.1.peg.2", predicted=PREDICTED))

    assert candidate is not None
    assert candidate.source == "afdb"
    assert candidate.accession == "AF-P12345-F1"
    assert candidate.confidence == 88.4
    assert candidate.suffix == ".cif"


def test_collection_writes_selected_structures_and_marks_the_gap(tmp_path: Path) -> None:
    report = {
        "proteins": [
            _protein("fig|1.1.peg.1", experimental=EXPERIMENTAL, predicted=PREDICTED),
            _protein("fig|1.1.peg.2", predicted=PREDICTED),
            _protein("fig|1.1.peg.3"),
            _protein("fig|1.1.peg.4", selected=False, experimental=EXPERIMENTAL),
        ]
    }
    fetched: list[str] = []

    def fetch(candidate) -> bytes:
        fetched.append(candidate.accession)
        return f"data_{candidate.accession}\n".encode()

    summary = collect_existing(report, tmp_path, fetch)

    assert summary.selected == 3
    assert summary.collected == 2
    assert summary.prediction_required == 1
    assert summary.failed == 0
    assert fetched == ["1abc_2", "AF-P12345-F1"]
    assert all(record["feature_id"] != "fig|1.1.peg.4" for record in summary.structures)

    by_id = {record["feature_id"]: record for record in summary.structures}
    pdb = by_id["fig|1.1.peg.1"]
    assert pdb["collection_status"] == "collected"
    assert pdb["usable_for_docking"] is None
    assert pdb["pockets"] == []
    assert (tmp_path / pdb["path"]).read_bytes() == b"data_1abc_2\n"

    missing = by_id["fig|1.1.peg.3"]
    assert missing["collection_status"] == "prediction_required"
    assert missing["path"] is None
    assert "prediction required" in missing["reason"]
    validate("structures", summary.structures)


def test_download_failure_is_explicit_and_does_not_abort_other_proteins(tmp_path: Path) -> None:
    report = {
        "proteins": [
            _protein("fig|1.1.peg.1", experimental=EXPERIMENTAL),
            _protein("fig|1.1.peg.2", predicted=PREDICTED),
        ]
    }

    def fetch(candidate) -> bytes:
        if candidate.source == "pdb":
            raise OSError("fixture says network unavailable")
        return b"data_AFDB\n"

    summary = collect_existing(report, tmp_path, fetch)

    assert summary.collected == 1
    assert summary.failed == 1
    failed = summary.structures[0]
    assert failed["collection_status"] == "failed"
    assert failed["usable_for_docking"] is False
    assert "network unavailable" in failed["reason"]
    assert summary.structures[1]["collection_status"] == "collected"


def test_empty_download_is_recorded_as_a_failure(tmp_path: Path) -> None:
    report = {"proteins": [_protein("fig|1.1.peg.1", experimental=EXPERIMENTAL)]}

    summary = collect_existing(report, tmp_path, lambda _candidate: b"")

    assert summary.failed == 1
    assert "was empty" in summary.structures[0]["reason"]


def test_limit_applies_to_selected_proteins_only(tmp_path: Path) -> None:
    report = {
        "proteins": [
            _protein("not-selected", selected=False, experimental=EXPERIMENTAL),
            _protein("first-selected", experimental=EXPERIMENTAL),
            _protein("second-selected", predicted=PREDICTED),
        ]
    }

    summary = collect_existing(report, tmp_path, lambda _candidate: b"structure", limit=1)

    assert summary.selected == 1
    assert [record["feature_id"] for record in summary.structures] == ["first-selected"]

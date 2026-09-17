"""Human-homology search tests (issue #29).

Parsing and gating run everywhere. The one test that shells out to DIAMOND skips when it is not
installed, so CI stays tool-free.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from s2f.m2_triage.human_homology import (
    DEFAULT_SENSITIVITY,
    STATUS_HIT,
    STATUS_NO_HIT,
    DiamondMissing,
    HumanHit,
    diamond_version,
    parse_diamond_hits,
    search,
    write_fasta,
)

# qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore
ROWS = [
    ["fig|1.1.peg.1", "sp|P06576|ATPB_HUMAN", "63.7", "460", "470", "529", "5", "464", "60", "519", "5.11e-206", "640"],
    ["fig|1.1.peg.1", "sp|P25705|ATPA_HUMAN", "24.0", "300", "470", "553", "10", "310", "20", "320", "1e-20", "90"],
    ["fig|1.1.peg.2", "sp|O15235|RT12_HUMAN", "50.6", "83", "139", "138", "41", "123", "45", "127", "2.83e-22", "85.5"],
    # 50% identity over 20 residues of a 200-residue protein: real alignment, not a homolog.
    ["fig|1.1.peg.3", "sp|Q00000|SHORT_HUMAN", "50.0", "20", "200", "300", "1", "20", "1", "20", "1e-6", "40"],
]


def test_best_hit_per_protein_is_chosen_by_bitscore() -> None:
    hits = parse_diamond_hits(ROWS)

    assert hits["fig|1.1.peg.1"].accession == "P06576"
    assert hits["fig|1.1.peg.1"].entry_name == "ATPB_HUMAN"
    assert hits["fig|1.1.peg.1"].identity == 63.7
    assert hits["fig|1.1.peg.1"].bitscore == 640


def test_coverage_is_computed_from_the_query_span() -> None:
    hits = parse_diamond_hits(ROWS)

    # rpsL: residues 41-123 of 139 = 83/139
    assert hits["fig|1.1.peg.2"].coverage == pytest.approx(83 / 139, abs=1e-4)
    assert hits["fig|1.1.peg.2"].counted is True


def test_a_short_alignment_does_not_count_as_a_homolog() -> None:
    """BV-BRC records an identity with no coverage at all; 50% over 20 residues is not homology."""
    hits = parse_diamond_hits(ROWS, min_coverage=0.5)

    short = hits["fig|1.1.peg.3"]
    assert short.status == STATUS_HIT
    assert short.identity == 50.0
    assert short.counted is False
    assert "below the 50% gate" in short.note


def test_the_coverage_gate_is_adjustable() -> None:
    assert parse_diamond_hits(ROWS, min_coverage=0.05)["fig|1.1.peg.3"].counted is True
    assert parse_diamond_hits(ROWS, min_coverage=0.9)["fig|1.1.peg.2"].counted is False


def test_malformed_rows_are_skipped() -> None:
    hits = parse_diamond_hits([["too", "few"], ["fig|1.1.peg.9", "sp|X|Y", "nope"] , *ROWS])
    assert "fig|1.1.peg.9" not in hits
    assert len(hits) == 3


def test_annotation_row_carries_source_and_reference() -> None:
    """report.json requires `source` on every annotation (issue #2)."""
    hit = parse_diamond_hits(ROWS)["fig|1.1.peg.2"]

    row = hit.as_annotation(retrieved_at="2026-09-17T12:00:00Z", source_version="2026_03")

    assert row["source"] == "human_homology_diamond"
    assert row["hit"] == "O15235"
    assert row["identity"] == pytest.approx(0.506, abs=1e-3)  # fraction, per the schema
    assert row["coverage"] == pytest.approx(0.597, abs=1e-3)
    assert row["counted_as_homolog"] is True
    assert "UP000005640" in row["reference"]


def test_a_protein_with_no_hit_says_so_explicitly() -> None:
    """"We looked and found nothing" must be distinguishable from "we never looked"."""
    missing = HumanHit(feature_id="fig|1.1.peg.7", status=STATUS_NO_HIT, note="no human hit")
    row = missing.as_annotation(retrieved_at="2026-09-17T12:00:00Z")

    assert row["hit"] is None
    assert row["description"] == "no human homolog"
    assert row["counted_as_homolog"] is False
    assert missing.identity is None  # not zero


def test_write_fasta_round_trips(tmp_path: Path) -> None:
    path = tmp_path / "q.faa"
    written = write_fasta(path, [("a", "MKALIV" * 20), ("b", "MSTKNK"), ("", "SKIPPED"), ("c", "")])

    assert written == 2
    text = path.read_text()
    assert text.startswith(">a\n")
    assert ">b\n" in text and "SKIPPED" not in text
    assert max(len(line) for line in text.splitlines() if not line.startswith(">")) <= 60


def test_missing_diamond_is_a_clear_error(monkeypatch) -> None:
    monkeypatch.setattr(shutil, "which", lambda _name: None)

    with pytest.raises(DiamondMissing, match="brew install diamond"):
        diamond_version()


@pytest.mark.skipif(shutil.which("diamond") is None, reason="diamond not installed")
def test_search_runs_diamond_end_to_end(tmp_path: Path) -> None:
    """A real DIAMOND run on a two-sequence reference, including the no-hit record."""
    reference = tmp_path / "ref.fasta"
    target = (
        "MAAQASPSPLSKLWSPSLLFSRQAAPAQLRCVSAAVFRHSGPRLVATHRVSSTRSLVSSHSHAVAVRSSGGSAHHVR"
        "GTRLCEQFVSSMQSRQVAALQLEHLARPTLRSFSTSLLSTRSFSTSTLSRSSRSSQSSQSSHTHSSQSSHSHSHSHS"
    )
    reference.write_text(f">sp|P00001|TEST_HUMAN Test protein\n{target}\n>sp|P00002|OTHER_HUMAN Other\nMKALIVMKALIVMKALIVMKALIVMKALIVMKALIV\n")

    result = search(
        [("fig|1.1.peg.1", target), ("fig|1.1.peg.2", "WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW")],
        reference,
        min_coverage=0.5,
        threads=1,
        workdir=tmp_path,
    )

    assert result.sensitivity == DEFAULT_SENSITIVITY
    assert result.reference_sequences == 2
    assert result.diamond_version.startswith("diamond")

    hit = result.hits["fig|1.1.peg.1"]
    assert hit.status == STATUS_HIT
    assert hit.accession == "P00001"
    assert hit.identity == pytest.approx(100.0, abs=0.1)
    assert hit.counted is True

    # Every protein gets a record, including the one with nothing.
    assert result.hits["fig|1.1.peg.2"].status == STATUS_NO_HIT
    assert result.summary()["proteins_with_hit"] == 1

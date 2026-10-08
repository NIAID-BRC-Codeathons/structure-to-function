"""M1's --from-fasta route: protein sequences in, no genome, no BV-BRC (#89)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from s2f.common.schema import validate_report
from s2f.m1_genome.__main__ import main
from s2f.m1_genome.fasta_route import query_proteins, read_protein_fasta
from s2f.m2_triage.bvbrc_input import load_input
from s2f.run import build_parser, stage_args

QUERY = Path(__file__).resolve().parents[1] / "fixtures/m1/fasta_query/query.faa"


def _run(tmp_path, *extra) -> dict:
    run_dir = tmp_path / "q1"
    assert main(["--run", str(run_dir), "--from-fasta", str(QUERY), *extra]) == 0
    return json.loads((run_dir / "report.json").read_text())


def test_a_fasta_run_writes_a_valid_report_that_says_what_it_is(tmp_path) -> None:
    report = _run(tmp_path)
    validate_report(report)
    assert report["genome"]["annotation_route"] == "fasta"
    assert report["run"]["modules"]["m1"]["bvbrc_used"] is False
    assert [p["feature_id"] for p in report["proteins"]] == [
        "gyrA", "tr|TEST00001|TEST1_KLEPN", "sp|P0DTC2|SPIKE_SARS2"]


def test_genome_only_evidence_is_not_available_rather_than_empty(tmp_path) -> None:
    report = _run(tmp_path)
    assert {"specialty", "subsystems", "contigs", "tree"} <= set(report["genome"]["not_available"])
    for protein in report["proteins"]:
        assert "specialty" not in protein          # absent, not an empty list
        assert "specialty" in protein["not_available"]


def test_the_organism_comes_from_the_flag_never_from_headers(tmp_path) -> None:
    assert _run(tmp_path)["genome"]["taxonomy"] is None
    report = _run(tmp_path, "--organism", "Klebsiella pneumoniae")
    assert report["genome"]["taxonomy"] == {"scientific_name": "Klebsiella pneumoniae",
                                            "source": "user"}
    # The headers carried two different organisms in brackets and OS= fields; none leak
    # into proteins.faa, where M2 would otherwise infer an organism from them.
    faa = (tmp_path / "q1/m1/proteins.faa").read_text()
    assert "[" not in faa and "OS=" not in faa


def test_m2_reads_the_fasta_route_output(tmp_path) -> None:
    _run(tmp_path)
    bundle = load_input(tmp_path / "q1/m1")
    by_id = {p.feature_id: p for p in bundle.proteins}
    assert len(by_id) == 3
    assert by_id["tr|TEST00001|TEST1_KLEPN"].product == "Uncharacterized protein"
    assert by_id["tr|TEST00001|TEST1_KLEPN"].is_uncharacterized
    assert not by_id["tr|TEST00001|TEST1_KLEPN"].sequence.endswith("*")
    assert bundle.organism == ""


def test_duplicate_ids_are_kept_apart_and_empty_records_dropped() -> None:
    proteins = query_proteins([("x first", "MKV"), ("x second", "MAL"), ("y empty", "")])
    assert [p["feature_id"] for p in proteins] == ["x", "x_2"]


def test_a_nucleotide_or_malformed_file_is_refused(tmp_path) -> None:
    bad = tmp_path / "bad.faa"
    bad.write_text("ACGT\n>late header\nMKV\n")
    with pytest.raises(ValueError):
        read_protein_fasta(bad)
    assert main(["--run", str(tmp_path / "r"), "--from-fasta", str(bad)]) == 2
    dna = tmp_path / "dna.faa"
    dna.write_text(">contig_1\n" + "ACGTTGCA" * 10 + "\n")
    with pytest.raises(ValueError, match="nucleotide"):
        read_protein_fasta(dna)


def test_routes_are_mutually_exclusive(tmp_path) -> None:
    assert main(["--run", str(tmp_path / "r"), "--from-fasta", str(QUERY),
                 "--from-bvbrc-api"]) == 2


def test_the_runner_forwards_the_route_to_m1(tmp_path) -> None:
    args = build_parser().parse_args(["--run", str(tmp_path), "--from-fasta", str(QUERY),
                                      "--organism", "Klebsiella pneumoniae"])
    m1 = stage_args("m1", args)
    assert m1[m1.index("--from-fasta") + 1] == str(QUERY)
    assert m1[m1.index("--organism") + 1] == "Klebsiella pneumoniae"

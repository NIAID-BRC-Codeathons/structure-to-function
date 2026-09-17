"""M1 CGA parsing (issue #5), against the trimmed real output in fixtures/m1."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from s2f.m1_genome.parse import (
    CgaLayoutError,
    genome_section,
    load_cga,
    proteins_section,
    run_fields,
    write_m1_dir,
)

FIXTURE = Path(__file__).resolve().parents[1] / "fixtures" / "m1" / "cga_sample"
GYRA = "fig|243273.147.peg.4"
GPMI = "fig|243273.147.peg.485"
HYPOTHETICAL = "fig|243273.147.peg.11"


@pytest.fixture(scope="module")
def run():
    return load_cga(FIXTURE)


def test_rejects_a_directory_that_is_not_cga_output(tmp_path):
    with pytest.raises(CgaLayoutError):
        load_cga(tmp_path)


def test_proteins_are_keyed_by_feature_id_with_integer_coordinates(run):
    proteins = proteins_section(run)
    assert len(proteins) == 5
    gyra = next(p for p in proteins if p["feature_id"] == GYRA)
    # BV-BRC writes coordinates as strings; the schema requires integers or null.
    assert isinstance(gyra["start"], int) and isinstance(gyra["end"], int)
    assert gyra["aa_length"] == 836
    assert gyra["contig"] == "243273.147.con.0001"


def test_specialty_rows_keep_classification_and_null_identity(run):
    """Most AMR rows are 'antibiotic target in susceptible species', not resistance,
    and k-mer rows have no identity. Both facts have to survive parsing."""
    proteins = {p["feature_id"]: p for p in proteins_section(run)}
    kmer = next(s for s in proteins[GYRA]["specialty"] if s["evidence"] == "K-mer Search")
    assert kmer["classification"] == "antibiotic target in susceptible species"
    assert kmer["identity"] is None and kmer["coverage"] is None
    assert kmer["antibiotics"]

    diamond = next(s for s in proteins[GPMI]["specialty"] if s["evidence"] == "DIAMOND")
    assert diamond["identity"] is not None and diamond["coverage"] is not None
    assert diamond["database"] == "TTD"


def test_kmer_confidence_is_carried_for_annotated_features(run):
    proteins = {p["feature_id"]: p for p in proteins_section(run)}
    assert proteins[GYRA]["kmer_confidence"]["hit_count"]


def test_protein_without_specialty_is_empty_not_missing(run):
    proteins = {p["feature_id"]: p for p in proteins_section(run)}
    assert proteins[HYPOTHETICAL]["specialty"] == []


def test_genome_section_reports_quality_and_tree(run):
    genome = genome_section(run)
    assert genome["genome_id"] == "243273.147"
    assert genome["taxon_id"] == 243273
    assert genome["quality"]["genome_quality"] == "Good"
    assert genome["quality"]["genome_quality_flags"] == []
    assert genome["tree_newick"].startswith("(")
    assert len(genome["tree_ingroup"]) == 10
    # CGA's own close_genomes is empty; closest_genomes comes from Minhash.
    assert genome["closest_genomes"] == []


def test_genome_section_carries_a_minhash_call_when_given_one(run):
    call = {
        "called_by": "minhash",
        "called_rank": "species",
        "top_hit": {"distance": 0.0},
        "hits": [{"genome_id": "243273.25", "genome_name": "M. genitalium G37",
                  "distance": 0.0, "pvalue": 0.0, "counts": "1000/1000"}],
    }
    genome = genome_section(run, call)
    assert genome["taxonomy"]["called_by"] == "minhash"
    assert genome["closest_genomes"][0]["mash_distance"] == 0.0
    assert genome["closest_genomes"][0]["ani"] is None  # issue #7 still owes us this


def test_run_fields_report_the_job_and_its_tools(run):
    fields = run_fields(run)
    assert fields["cga"]["job_id"] == "23587516"
    assert fields["cga"]["elapsed_seconds"] == pytest.approx(164.5, abs=0.5)
    assert fields["cga"]["genetic_code"] == 4
    assert fields["tool_versions"]["cga_analysis_events"]


def test_write_m1_dir_produces_what_m2_reads(run, tmp_path):
    counts = write_m1_dir(run, tmp_path / "m1")
    assert counts == {"proteins_faa": 5, "feature_rows": 5, "specialty_rows": 4}

    faa = (tmp_path / "m1" / "proteins.faa").read_text()
    assert faa.count(">") == 5
    assert faa.startswith(">fig|")

    from s2f.m2_triage.bvbrc_input import load_input

    bundle = load_input(tmp_path / "m1")
    assert len(bundle.proteins) == 5
    loaded = {p.feature_id: p for p in bundle.proteins}
    assert loaded[GYRA].product.startswith("DNA gyrase subunit A")
    assert loaded[GYRA].specialty, "M2 must see the specialty rows M1 wrote"
    assert loaded[HYPOTHETICAL].is_uncharacterized


def test_poor_quality_is_detected(run, tmp_path):
    assert run.is_poor is False
    poor = load_cga(FIXTURE)
    poor.genome["genome_quality"] = "Poor"
    assert poor.is_poor is True
    assert "Poor" in poor.quality_reason()


class TestUploadSafety:
    """p3-cp exits 0 without overwriting, so a reused filename silently feeds Minhash
    the previous run's genome. That happened once; these pin the guard."""

    def test_ls_size_is_parsed_from_a_p3_ls_row(self):
        from s2f.m1_genome.cga import LS_SIZE

        row = "-rw-- someone@bvbrc 20866 Sep 17 09:32 mgen_G37_cga"
        assert LS_SIZE.match(row).group(1) == "20866"

    def test_upload_rejects_a_workspace_copy_of_the_wrong_size(self, tmp_path, monkeypatch):
        from s2f.m1_genome import cga

        local = tmp_path / "run_a.contigs.fna"
        local.write_text(">contig_1\nACGT\n")

        monkeypatch.setattr(cga, "_run", lambda *a, **k: "")
        monkeypatch.setattr(cga, "remote_size", lambda path: 999999)
        with pytest.raises(cga.CgaError, match="not this run's file|is not "):
            cga.upload(local, "/someone@bvbrc/home/s2f")

    def test_upload_returns_the_path_when_sizes_agree(self, tmp_path, monkeypatch):
        from s2f.m1_genome import cga

        local = tmp_path / "run_a.contigs.fna"
        local.write_text(">contig_1\nACGT\n")

        monkeypatch.setattr(cga, "_run", lambda *a, **k: "")
        monkeypatch.setattr(cga, "remote_size", lambda path: local.stat().st_size)
        assert cga.upload(local, "/someone@bvbrc/home/s2f") == \
            "/someone@bvbrc/home/s2f/uploads/run_a.contigs.fna"

    def test_upload_rejects_a_missing_workspace_copy(self, tmp_path, monkeypatch):
        from s2f.m1_genome import cga

        local = tmp_path / "run_a.contigs.fna"
        local.write_text(">contig_1\nACGT\n")

        monkeypatch.setattr(cga, "_run", lambda *a, **k: "")
        monkeypatch.setattr(cga, "remote_size", lambda path: None)
        with pytest.raises(cga.CgaError, match="left nothing"):
            cga.upload(local, "/someone@bvbrc/home/s2f")

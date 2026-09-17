"""The shared genome fixture on main (PR #25) must load as M1 input.

Layout: ``fixtures/genomes/<name>/<name>.faa`` plus ``<name>.features.json`` and
``<name>.sp_gene.json`` — JSON, and FASTA headers that differ from the p3-CLI export.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from s2f.m2_triage.bvbrc_input import _feature_id_from_header, load_input

FIXTURE = Path(__file__).resolve().parents[1] / "fixtures" / "genomes" / "mgen_G37"


@pytest.mark.skipif(not FIXTURE.exists(), reason="shared genome fixture not present")
def test_loads_the_shared_g37_fixture() -> None:
    bundle = load_input(FIXTURE)

    assert len(bundle.proteins) == 542
    # Nothing falls out of the join in either direction.
    assert bundle.specialty_without_protein == []
    assert bundle.table_rows_without_sequence == []
    assert bundle.sequences_without_table_row == []

    protein = bundle.proteins[0]
    assert protein.feature_id == "fig|243273.25.peg.308"
    assert protein.product == "hypothetical protein"
    assert protein.locus_tag == "MG_267"
    assert protein.sequence.startswith("MTLLFKLVKIAIL")


@pytest.mark.skipif(not FIXTURE.exists(), reason="shared genome fixture not present")
def test_g37_specialty_properties_reach_the_scoring_flags() -> None:
    bundle = load_input(FIXTURE)

    assert sum(1 for p in bundle.proteins if p.is_essential_ortholog) == 148
    assert sum(1 for p in bundle.proteins if p.human_homolog_identity) == 5
    assert sum(1 for p in bundle.proteins if p.is_uncharacterized) == 190

    atp_alpha = next(p for p in bundle.proteins if p.feature_id == "fig|243273.25.peg.468")
    assert atp_alpha.human_homolog_identity == 52.0


@pytest.mark.parametrize(
    ("header", "feature_id", "product"),
    [
        # p3-CLI / web export
        (
            "fig|1125630.4.peg.5196  DNA replication protein",
            "fig|1125630.4.peg.5196",
            "DNA replication protein",
        ),
        # The fixture layout: extra pipe-separated IDs, organism in trailing brackets.
        (
            "fig|243273.25.peg.308|MG_267|VBIMycGen98045_0308| hypothetical protein "
            "[Mycoplasma genitalium G37 | 243273.25]",
            "fig|243273.25.peg.308",
            "hypothetical protein",
        ),
        ("fig|1.1.peg.1", "fig|1.1.peg.1", ""),
        ("some_other_id  a product", "some_other_id", "a product"),
    ],
)
def test_feature_id_from_header_handles_both_layouts(header, feature_id, product) -> None:
    assert _feature_id_from_header(header) == (feature_id, product)


def test_json_tables_are_read_like_csv(tmp_path: Path) -> None:
    m1 = tmp_path / "m1"
    m1.mkdir()
    (m1 / "genome.faa").write_text(">fig|9.9.peg.1|LT_1| Test protein [Org | 9.9]\nMKALIV\n")
    (m1 / "genome.features.json").write_text(
        json.dumps([{"patric_id": "fig|9.9.peg.1", "product": "Test protein", "refseq_locus_tag": "LT_1"}])
    )
    (m1 / "genome.sp_gene.json").write_text(
        json.dumps(
            [
                {"patric_id": "fig|9.9.peg.1", "property": "Essential Gene", "source": "PATRIC"},
                {"patric_id": "fig|9.9.peg.1", "property": "Human Homolog", "source": "Human", "identity": 52},
            ]
        )
    )

    bundle = load_input(m1)

    assert len(bundle.proteins) == 1
    protein = bundle.proteins[0]
    assert protein.locus_tag == "LT_1"
    assert protein.is_essential_ortholog
    assert protein.human_homolog_identity == 52.0


def test_nucleotide_fna_is_not_mistaken_for_protein_input(tmp_path: Path) -> None:
    m1 = tmp_path / "m1"
    m1.mkdir()
    (m1 / "genome.fna").write_text(">contig_1\nACGTACGTACGT\n")
    (m1 / "genome.faa").write_text(">fig|9.9.peg.1  Test protein\nMKALIV\n")

    bundle = load_input(m1)

    assert [p.sequence for p in bundle.proteins] == ["MKALIV"]

"""Query-taxonomy resolution tests (issue #55).

The point of the change: M1 determines the organism, so M2 should read it rather than infer it
from FASTA header brackets. The lineage also settles two things spelling cannot — the current
genus name, and identity at a named rank.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from s2f.m2_triage.taxonomy import (
    SOURCE_HEADERS,
    SOURCE_M1,
    SOURCE_OVERRIDE,
    SOURCE_UNDETERMINED,
    QueryTaxonomy,
    from_report,
    resolve,
)

# Real values: BV-BRC calls the genome "Mycoplasma genitalium G37", NCBI's lineage says
# "Mycoplasmoides". The species is 2097, the strain 243273.
G37_GENOME = {
    "genome_id": "243273.25",
    "taxon_id": 243273,
    "taxonomy": {
        "scientific_name": "Mycoplasma genitalium G37",
        "lineage_names": [
            "Bacteria", "Mycoplasmatota", "Mollicutes", "Mycoplasmoidales",
            "Mycoplasmoidaceae", "Mycoplasmoides", "Mycoplasmoides genitalium",
            "Mycoplasmoides genitalium G37",
        ],
        "lineage_ids": [2, 544448, 31969, 2790995, 2790998, 2995234, 2097, 243273],
    },
}


def write_report(run_dir: Path, genome: dict | None) -> Path:
    run_dir.mkdir(parents=True, exist_ok=True)
    payload: dict = {"schema_version": "0.1.0"}
    if genome is not None:
        payload["genome"] = genome
    (run_dir / "report.json").write_text(json.dumps(payload), encoding="utf-8")
    return run_dir


def test_the_genus_comes_from_the_lineage_not_the_name(tmp_path: Path) -> None:
    taxonomy = from_report(write_report(tmp_path, G37_GENOME))

    assert taxonomy.source == SOURCE_M1
    assert taxonomy.scientific_name == "Mycoplasma genitalium G37"
    # The name says Mycoplasma; NCBI and the PDB say Mycoplasmoides. No stem heuristic needed.
    assert taxonomy.genus == "mycoplasmoides"
    assert taxonomy.genus_taxon_id() == 2995234


def test_a_species_level_hit_matches_a_strain_level_query(tmp_path: Path) -> None:
    """@cmmann21's case on #12: hit annotated 2097, query genome 243273.25. IDs never match."""
    taxonomy = from_report(write_report(tmp_path, G37_GENOME))

    assert taxonomy.species_taxon_ids() == {2097, 243273}
    assert taxonomy.same_species_as("Mycoplasmoides genitalium", 2097) is True
    assert taxonomy.same_species_as("Mycoplasmoides genitalium G37", 243273) is True


def test_same_genus_without_same_species(tmp_path: Path) -> None:
    taxonomy = from_report(write_report(tmp_path, G37_GENOME))

    assert taxonomy.same_species_as("Mycoplasmoides pneumoniae M129", 2104) is False
    assert taxonomy.same_genus_as("Mycoplasmoides pneumoniae M129", 2104) is True


def test_an_unrelated_organism_matches_nothing(tmp_path: Path) -> None:
    taxonomy = from_report(write_report(tmp_path, G37_GENOME))

    assert taxonomy.same_species_as("Escherichia coli", 562) is False
    assert taxonomy.same_genus_as("Escherichia coli", 562) is False


def test_an_unknown_hit_organism_is_not_called_foreign(tmp_path: Path) -> None:
    """An unresolvable hit is unflagged, not flagged as a cross-organism transfer."""
    taxonomy = from_report(write_report(tmp_path, G37_GENOME))

    assert taxonomy.same_species_as("", None) is False
    assert taxonomy.same_genus_as("", None) is False


def test_precedence_override_then_m1_then_headers(tmp_path: Path) -> None:
    run_dir = write_report(tmp_path / "run", G37_GENOME)

    explicit = resolve(run_dir, override="Klebsiella pneumoniae", header_organism="from headers")
    assert explicit.source == SOURCE_OVERRIDE
    assert explicit.scientific_name == "Klebsiella pneumoniae"

    from_m1 = resolve(run_dir, header_organism="from headers")
    assert from_m1.source == SOURCE_M1

    standalone = resolve(tmp_path / "no-report", header_organism="Mycoplasma genitalium G37")
    assert standalone.source == SOURCE_HEADERS
    assert standalone.determined is True


def test_no_report_and_no_headers_is_undetermined(tmp_path: Path) -> None:
    """Zero same-species hits and zero measurements are not the same claim."""
    taxonomy = resolve(tmp_path / "nothing")

    assert taxonomy.source == SOURCE_UNDETERMINED
    assert taxonomy.determined is False
    assert taxonomy.as_dict()["scientific_name"] is None


@pytest.mark.parametrize("genome", [None, {}, {"taxonomy": {}}, {"taxonomy": None}])
def test_a_report_without_a_usable_genome_section_is_undetermined(tmp_path, genome) -> None:
    taxonomy = from_report(write_report(tmp_path / str(abs(hash(str(genome)))), genome))
    assert taxonomy.determined is False


def test_a_corrupt_report_does_not_raise(tmp_path: Path) -> None:
    tmp_path.mkdir(parents=True, exist_ok=True)
    (tmp_path / "report.json").write_text("{not json", encoding="utf-8")

    assert from_report(tmp_path).determined is False


def test_unpaired_lineage_falls_back_to_names(tmp_path: Path) -> None:
    """Names and IDs must line up to be used as pairs; otherwise the IDs are dropped."""
    genome = {
        "taxon_id": 243273,
        "taxonomy": {
            "scientific_name": "Mycoplasma genitalium G37",
            "lineage_names": ["Bacteria", "Mycoplasmoides", "Mycoplasmoides genitalium"],
            "lineage_ids": [2, 2995234],  # one short
        },
    }
    taxonomy = from_report(write_report(tmp_path, genome))

    assert taxonomy.lineage_ids == []
    assert taxonomy.genus == "mycoplasmoides"      # names still usable
    assert taxonomy.species_taxon_ids() == {243273}  # only the genome's own taxon


def test_without_a_lineage_it_falls_back_to_name_comparison() -> None:
    """The standalone path: no M1 report, so the string heuristics still have to work."""
    taxonomy = QueryTaxonomy(scientific_name="Mycoplasma genitalium G37", source=SOURCE_HEADERS)

    assert taxonomy.has_lineage is False
    assert taxonomy.same_species_as("Mycoplasmoides genitalium G37") is True   # renaming tolerated
    assert taxonomy.same_genus_as("Mycoplasmoides pneumoniae M129") is True
    assert taxonomy.same_genus_as("Streptomyces coelicolor") is False

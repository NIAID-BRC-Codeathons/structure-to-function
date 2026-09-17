"""Quality-gate tests for M3 issue 13."""

from __future__ import annotations

from s2f.m3_fold.quality import QualityThresholds, assess_structure


PDB_EVIDENCE = {
    "identity": 0.91,
    "coverage": 0.95,
    "resolution": 2.1,
    "method": "X-ray",
}

AFDB_EVIDENCE = {
    "mean_plddt": 82.5,
    "coverage": 0.95,
}

AFDB_CIF = b"""data_AFDB
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.label_atom_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.B_iso_or_equiv
_atom_site.pdbx_PDB_model_num
ATOM 1 CA A 1 92.0 1
ATOM 2 CA A 2 65.0 1
ATOM 3 CA A 3 45.0 1
ATOM 4 CA A 4 88.0 1
#
"""


def test_good_experimental_homolog_passes_with_a_readable_reason() -> None:
    decision = assess_structure("pdb", PDB_EVIDENCE, b"coordinates", QualityThresholds())

    assert decision.usable_for_docking is True
    assert "identity 91.0%" in decision.reason
    assert all(check["passed"] for check in decision.details["checks"])
    assert decision.details["pae"]["status"] == "not_applicable"


def test_experimental_homolog_fails_and_reports_each_failed_check() -> None:
    evidence = {**PDB_EVIDENCE, "coverage": 0.40, "resolution": 4.0}

    decision = assess_structure("pdb", evidence, b"coordinates", QualityThresholds())

    assert decision.usable_for_docking is False
    assert "sequence_coverage" in decision.reason
    assert "resolution_angstrom" in decision.reason


def test_nmr_structure_does_not_require_a_resolution_value() -> None:
    evidence = {**PDB_EVIDENCE, "method": "Solution NMR", "resolution": None}

    decision = assess_structure("pdb", evidence, b"coordinates", QualityThresholds())

    assert decision.usable_for_docking is True
    resolution_check = decision.details["checks"][2]
    assert resolution_check["operator"] == "not_applicable"


def test_afdb_gate_masks_low_confidence_regions_without_renumbering() -> None:
    decision = assess_structure("afdb", AFDB_EVIDENCE, AFDB_CIF, QualityThresholds())

    assert decision.usable_for_docking is True
    assert decision.details["observed"]["assessed_residues"] == 4
    assert decision.details["observed"]["masked_residues"] == 2
    assert decision.details["masked_regions"] == [
        {
            "chain": "A",
            "start": 2,
            "end": 3,
            "residue_count": 2,
            "mean_plddt": 55.0,
            "threshold": 70.0,
            "action": "exclude_from_pocket_search",
        }
    ]
    assert decision.details["pae"]["status"] == "not_available"


def test_afdb_gate_rejects_low_global_confidence() -> None:
    evidence = {**AFDB_EVIDENCE, "mean_plddt": 65.0}

    decision = assess_structure("afdb", evidence, AFDB_CIF, QualityThresholds())

    assert decision.usable_for_docking is False
    assert "mean_plddt 65 is below minimum 70" in decision.reason


def test_quality_thresholds_can_be_overridden() -> None:
    strict = QualityThresholds(min_pdb_identity=0.95)

    decision = assess_structure("pdb", PDB_EVIDENCE, b"coordinates", strict)

    assert decision.usable_for_docking is False
    assert decision.details["thresholds"]["min_pdb_identity"] == 0.95

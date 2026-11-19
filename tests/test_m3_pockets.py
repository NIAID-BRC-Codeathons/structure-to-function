"""Binding sites transferred from what the homolog entry binds (issue #14).

The geometry cases use a small hand-built mmCIF so the expected answer can be read off the
coordinates. The last tests use `fixtures/m3/6MU0.cif`, a real entry, so the module is also
exercised against a file nobody wrote for it.
"""

from __future__ import annotations

from pathlib import Path

from s2f.m3_fold.pockets import (CONTACT_DISTANCE, find_holo_sites, find_pockets,
                                 reject_low_confidence)

FIXTURE = Path(__file__).resolve().parents[1] / "fixtures/m3/6MU0.cif"

HEADER = """data_T
#
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
B 2
D 3
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.B_iso_or_equiv
_atom_site.pdbx_PDB_model_num
"""


def _cif(rows: list[str]) -> bytes:
    return (HEADER + "\n".join(rows) + "\n#\n").encode("utf-8")


def _atom(serial, atom, comp, chain, seq, x, y, z, b=90.0, het=False, model=1):
    group = "HETATM" if het else "ATOM"
    return (f"{group} {serial} {atom} {comp} {chain} {seq} {seq} {chain} "
            f"{x} {y} {z} {b} {model}")


#: Three residues line the ligand, one sits far away, one belongs to another entity.
CONTACT_CASE = _cif([
    _atom(1, "CA", "ASP", "A", 1, 0.0, 0.0, 0.0),
    _atom(2, "CA", "HIS", "A", 2, 3.0, 0.0, 0.0),
    _atom(3, "CA", "THR", "A", 3, 0.0, 3.0, 0.0),
    _atom(4, "CA", "LEU", "A", 4, 30.0, 30.0, 30.0),
    _atom(5, "CA", "VAL", "B", 1, 1.0, 1.0, 1.0),
    _atom(6, "MG", "MG", "D", 101, 1.0, 1.0, 0.0, b=20.0, het=True),
])


def test_a_site_is_the_residues_that_line_the_bound_component() -> None:
    pockets = find_holo_sites(CONTACT_CASE, metals=["MG"], entity_id="1")

    assert len(pockets) == 1
    pocket = pockets[0]
    assert [residue.seq_id for residue in pocket.residues] == [1, 2, 3]
    assert pocket.component == "MG"
    assert pocket.component_kind == "metal"
    # Read from the entry that binds it, so no superposition was involved.
    assert pocket.aligned_to_known_site is True
    assert pocket.center == (1.0, 1.0, 0.0)


def test_only_our_entity_lines_the_site() -> None:
    """Chain B is a different entity: its residue is near the metal but not our pocket."""
    pockets = find_holo_sites(CONTACT_CASE, metals=["MG"], entity_id="1")
    assert {residue.chain for residue in pockets[0].residues} == {"A"}

    both = find_holo_sites(CONTACT_CASE, metals=["MG"], entity_id="")
    assert {residue.chain for residue in both[0].residues} == {"A", "B"}


def test_a_metal_site_says_a_cofactor_is_required() -> None:
    pocket = find_holo_sites(CONTACT_CASE, metals=["MG"], entity_id="1")[0]

    assert pocket.cofactor_required is True
    assert pocket.cofactor_present is True
    assert pocket.params["metals_in_entry"] == ["MG"]


def test_a_component_nobody_reported_is_not_a_site() -> None:
    """The ligand list comes from M2, which already dropped waters and additives."""
    assert find_holo_sites(CONTACT_CASE, ligands=[], metals=[], entity_id="1") == []
    assert find_holo_sites(CONTACT_CASE, ligands=["ATP"], entity_id="1") == []


def test_a_brush_of_one_residue_is_not_a_pocket() -> None:
    sparse = _cif([
        _atom(1, "CA", "ASP", "A", 1, 0.0, 0.0, 0.0),
        _atom(2, "CA", "LEU", "A", 2, 30.0, 0.0, 0.0),
        _atom(3, "MG", "MG", "D", 101, 1.0, 0.0, 0.0, het=True),
    ])
    assert find_holo_sites(sparse, metals=["MG"], entity_id="1") == []


def test_a_predicted_model_carries_no_site_to_read() -> None:
    """The site is in the holo entry, not in a model of our protein (#14, pitfalls #2)."""
    assert find_holo_sites(CONTACT_CASE, metals=["MG"], structure_is_template=False) == []


def test_a_pocket_in_a_low_confidence_region_is_rejected() -> None:
    poor = _cif([
        _atom(1, "CA", "ASP", "A", 1, 0.0, 0.0, 0.0, b=40.0),
        _atom(2, "CA", "HIS", "A", 2, 3.0, 0.0, 0.0, b=35.0),
        _atom(3, "CA", "THR", "A", 3, 0.0, 3.0, 0.0, b=30.0),
        _atom(4, "MG", "MG", "D", 101, 1.0, 1.0, 0.0, het=True),
    ])
    pockets = reject_low_confidence(
        find_holo_sites(poor, metals=["MG"], entity_id="1"), source="afdb")

    assert pockets[0].usable is False
    assert "below 70" in pockets[0].reason

    # The same B-factor column in an experimental entry is not a confidence claim.
    experimental = reject_low_confidence(
        find_holo_sites(poor, metals=["MG"], entity_id="1"), source="pdb")
    assert experimental[0].usable is True
    assert experimental[0].confidence is None


def test_a_drug_like_ligand_outranks_a_metal() -> None:
    both = _cif([
        _atom(1, "CA", "ASP", "A", 1, 0.0, 0.0, 0.0),
        _atom(2, "CA", "HIS", "A", 2, 3.0, 0.0, 0.0),
        _atom(3, "CA", "THR", "A", 3, 0.0, 3.0, 0.0),
        _atom(4, "MG", "MG", "D", 101, 1.0, 1.0, 0.0, het=True),
        _atom(5, "C1", "5RP", "B", 201, 1.0, 1.0, 2.0, het=True),
    ])
    pockets = find_holo_sites(both, ligands=["5RP"], metals=["MG"], entity_id="1")

    assert [pocket.component for pocket in pockets] == ["5RP", "MG"]
    # The entry holds a metal, so even the ligand site may be incomplete without it.
    assert all(pocket.cofactor_required for pocket in pockets)


# --- against a real entry ---------------------------------------------------------------


def test_the_real_entry_gives_the_known_active_site() -> None:
    """6MU0 is ribose-5-phosphate isomerase B with its substrate bound."""
    payload = FIXTURE.read_bytes()

    pockets, summary = find_pockets(
        payload,
        source="pdb",
        holo_template={"pdb_entity": "6MU0_1", "ligands": ["5RP"], "metals": []},
        structure_accession="6MU0_1",
    )

    assert summary["status"] == "ok"
    assert summary["template_is_collected_structure"] is True
    assert len(pockets) == 1
    pocket = pockets[0]
    assert pocket["component"] == "5RP"
    assert pocket["cofactor_required"] is False
    # The catalytic aspartate and histidine of RpiB line the substrate.
    lining = {(residue["auth_seq_id"], residue["comp_id"]) for residue in pocket["residues"]}
    assert {("10", "ASP"), ("11", "HIS")} <= lining
    assert all(residue["min_distance"] <= CONTACT_DISTANCE for residue in pocket["residues"])
    assert pocket["residues"] == sorted(
        pocket["residues"], key=lambda residue: int(residue["auth_seq_id"]))


def test_a_protein_with_no_holo_homolog_says_why_it_has_no_pocket() -> None:
    """Uncovered is not the same claim as "no pocket exists" (the #90 distinction)."""
    pockets, summary = find_pockets(
        FIXTURE.read_bytes(), source="pdb", holo_template=None, structure_accession="6MU0_1")

    assert pockets == []
    assert summary["status"] == "no_pocket"
    assert "P2Rank" in summary["reason"]


def test_a_site_in_another_entry_is_not_silently_transferred() -> None:
    pockets, summary = find_pockets(
        FIXTURE.read_bytes(),
        source="afdb",
        holo_template={"pdb_entity": "6MU0_1", "ligands": ["5RP"], "metals": []},
        structure_accession="AF-P47543-F1",
    )

    assert pockets == []
    assert summary["template_is_collected_structure"] is False
    assert "superposition" in summary["reason"]

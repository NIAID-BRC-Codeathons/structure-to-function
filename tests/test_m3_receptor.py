"""Receptors and search boxes for M4 (issue #14, part B).

No test here depends on a preparation tool being installed: `find_tool` is injected, so the
suite behaves the same on a node with OpenBabel and on CI without it.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from s2f.m3_fold.receptor import (autodock_type, docking_box, prepare_receptor,
                                  receptor_atoms, vina_config)
from s2f.m3_fold.pockets import find_pockets

FIXTURE = Path(__file__).resolve().parents[1] / "fixtures/m3/6MU0.cif"
NO_TOOL = lambda: None  # noqa: E731 - a stand-in for a node with nothing installed
FAKE_TOOL = lambda: "obabel"  # noqa: E731


@pytest.fixture
def payload() -> bytes:
    return FIXTURE.read_bytes()


@pytest.fixture
def pocket(payload: bytes) -> dict:
    pockets, _ = find_pockets(
        payload,
        source="pdb",
        holo_template={"pdb_entity": "6MU0_1", "ligands": ["5RP"], "metals": []},
        structure_accession="6MU0_1",
    )
    return pockets[0]


def test_the_receptor_is_our_entity_without_solvent_or_ligand(payload: bytes) -> None:
    atoms = receptor_atoms(payload, "1")

    assert atoms
    assert {atom.auth_chain for atom in atoms} == {"A"}
    # Docking into a site that still holds its ligand finds the ligand.
    assert not any(atom.comp_id == "5RP" for atom in atoms)
    assert not any(atom.comp_id == "HOH" for atom in atoms)
    assert not any(atom.atom_id.startswith("H") for atom in atoms)


def test_without_a_preparation_tool_no_pdbqt_is_invented(tmp_path: Path, payload: bytes) -> None:
    prepared = prepare_receptor(
        payload, entity_id="1", destination_stem=tmp_path / "p", run_dir=tmp_path,
        find_tool=NO_TOOL)

    assert prepared.pdb_path == "p.pdb"
    assert prepared.pdbqt_path is None
    assert prepared.fidelity == "pdb_only"
    assert "skipped rather than approximated" in prepared.reason
    assert (tmp_path / "p.pdb").exists()


def test_the_fallback_converter_says_what_it_is(tmp_path: Path, payload: bytes) -> None:
    prepared = prepare_receptor(
        payload, entity_id="1", destination_stem=tmp_path / "p", run_dir=tmp_path,
        find_tool=NO_TOOL, allow_fallback=True)

    assert prepared.fidelity == "reduced"
    assert prepared.prepared_by == "builtin_fallback"
    assert "no added" in prepared.reason and "hydrogens" in prepared.reason
    assert "not for a reported docking result" in prepared.reason
    line = (tmp_path / "p.pdbqt").read_text().splitlines()[0]
    assert line.startswith("ATOM")
    assert "+0.000" in line


def test_a_feature_id_ending_in_a_number_keeps_it(tmp_path: Path, payload: bytes) -> None:
    """`Path.with_suffix` treats the `.463` of fig|...peg.463 as an extension and eats it."""
    prepared = prepare_receptor(
        payload, entity_id="1", destination_stem=tmp_path / "fig_243273.25.peg.463",
        run_dir=tmp_path, find_tool=NO_TOOL)

    assert prepared.pdb_path == "fig_243273.25.peg.463.pdb"
    assert (tmp_path / "fig_243273.25.peg.463.pdb").exists()


def test_a_structure_with_no_readable_atoms_prepares_nothing(tmp_path: Path) -> None:
    prepared = prepare_receptor(
        b"data_X\n#\n", entity_id="1", destination_stem=tmp_path / "p", run_dir=tmp_path,
        find_tool=NO_TOOL)

    assert prepared.pdb_path is None
    assert prepared.fidelity == "none"
    assert not list(tmp_path.iterdir())


def test_a_tool_that_writes_nothing_is_not_reported_as_success(
    tmp_path: Path, payload: bytes, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(
        "s2f.m3_fold.receptor.run_preparation_tool", lambda *_, **__: (False, "obabel broke"))

    prepared = prepare_receptor(
        payload, entity_id="1", destination_stem=tmp_path / "p", run_dir=tmp_path,
        find_tool=FAKE_TOOL)

    assert prepared.pdbqt_path is None
    assert prepared.fidelity == "pdb_only"
    assert "obabel broke" in prepared.reason


# --- search boxes -------------------------------------------------------------------------


def test_the_box_is_drawn_round_the_site_not_a_cube_at_its_centre(
    payload: bytes, pocket: dict
) -> None:
    box = docking_box(payload, pocket)

    assert box is not None
    # A real site is not a cube: the axes differ, which a fixed-size box cannot express.
    assert len({box["size_x"], box["size_y"], box["size_z"]}) > 1
    assert all(box[f"size_{axis}"] >= 16.0 for axis in "xyz")
    # The ligand sits inside the box it defines.
    centre = pocket["center"]
    for value, axis in zip(centre, "xyz"):
        assert abs(value - box[f"center_{axis}"]) <= box[f"size_{axis}"] / 2


def test_a_tight_site_still_gets_a_searchable_box(payload: bytes) -> None:
    single = {"component": "", "residues": [
        {"chain": "A", "auth_seq_id": "11", "comp_id": "HIS"}]}

    box = docking_box(payload, single, padding=0.0)

    assert box is not None
    assert box["size_x"] == box["size_y"] == box["size_z"] == 16.0


def test_a_pocket_with_no_residues_has_no_box(payload: bytes) -> None:
    assert docking_box(payload, {"residues": []}) is None


def test_the_vina_config_states_its_own_parameters(payload: bytes, pocket: dict) -> None:
    config = vina_config(docking_box(payload, pocket), "receptor.pdbqt")

    assert "receptor = receptor.pdbqt" in config
    assert "exhaustiveness = 8" in config
    assert "num_modes = 9" in config


def test_aromatic_carbons_are_typed_apart_from_aliphatic_ones() -> None:
    assert autodock_type("CG", "PHE", "C") == "A"
    assert autodock_type("CB", "PHE", "C") == "C"
    assert autodock_type("ND1", "HIS", "N") == "NA"
    assert autodock_type("N", "ALA", "N") == "N"
    assert autodock_type("OD1", "ASP", "O") == "OA"

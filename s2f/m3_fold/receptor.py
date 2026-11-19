"""Docking-ready receptors and search boxes for M4 (issue #14, part B).

Two of the three things M4 needs can be produced exactly, with no external tool: a cleaned
receptor and a search box drawn around the site M3 defined. The third, a PDBQT, cannot —
assigning atom types and adding the polar hydrogens Vina needs is what Meeko, OpenBabel and
PDB2PQR exist for, and a homemade substitute that merely *parses* is a receptor that docks
and lies about it.

So the default is: always write the receptor PDB and the box; write a PDBQT when a real
preparation tool is on PATH; otherwise say so and write none. `--pdbqt-fallback` turns on the
built-in converter adapted from Ping Hu's `bacteria_structure_docking.py`, which produces a
geometrically valid PDBQT with AutoDock atom types, no added hydrogens and zero charges. That
is enough for a smoke test and not enough for a result; every record it writes carries
`fidelity: "reduced"` and the reason, and M4 is expected to refuse it unless asked.
"""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Iterable

from .mmcif import Atom, entity_chains, read_atoms

#: Padding added around the site's own atoms to give a ligand room to move. Vina's own
#: guidance is a box a little larger than the ligand; 4 A each side is the usual starting
#: point and is recorded in every box record rather than assumed.
BOX_PADDING = 4.0

#: No box smaller than this on any axis: a tight site would otherwise leave no room to search.
MIN_BOX_EDGE = 16.0

#: Preparation tools, best first. Meeko is AutoDock's own; OpenBabel is the common fallback.
PDBQT_TOOLS = ("mk_prepare_receptor.py", "obabel")


@dataclass(frozen=True)
class PreparedReceptor:
    """What was written for one structure, and what it is fit for."""

    pdb_path: str | None
    pdbqt_path: str | None
    prepared_by: str
    fidelity: str
    reason: str
    chains: list[str]
    atoms: int

    def as_dict(self) -> dict[str, Any]:
        return {
            "pdb_path": self.pdb_path,
            "pdbqt_path": self.pdbqt_path,
            "prepared_by": self.prepared_by,
            "fidelity": self.fidelity,
            "reason": self.reason,
            "chains": self.chains,
            "atoms": self.atoms,
        }


def available_tool() -> str | None:
    """The first receptor-preparation tool on PATH, or None.

    Discovery only: nothing here is installed, and no path is hard-coded, the same contract
    M2's annotation providers keep.
    """
    for tool in PDBQT_TOOLS:
        if shutil.which(tool):
            return tool
    return None


def receptor_atoms(payload: bytes, entity_id: str = "") -> list[Atom]:
    """The polymer atoms a receptor is built from: our entity, first model, no solvent.

    Ligands and waters are dropped because docking into a site that still contains its
    ligand finds the ligand. The bound components stay in the collected structure file,
    which is never modified, so nothing is lost.
    """
    chains = entity_chains(payload, entity_id) if entity_id else set()
    return [
        atom
        for atom in read_atoms(payload)
        if not atom.is_hetatm
        and atom.has_coordinates
        and (not chains or atom.chain in chains)
        and not _is_hydrogen(atom)
    ]


def _is_hydrogen(atom: Atom) -> bool:
    name = atom.atom_id.strip()
    return name.startswith("H") or (len(name) > 1 and name[0].isdigit() and name[1] == "H")


def write_receptor_pdb(atoms: Iterable[Atom], destination: Path) -> int:
    """Write the receptor as PDB, the format every preparation tool reads. Returns atoms."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    lines = []
    for serial, atom in enumerate(atoms, 1):
        name = atom.atom_id[:4]
        # PDB column rules: a one- or two-letter element starts in column 14, not 13.
        padded = f" {name:<3}" if len(name) < 4 else name
        element = _element(atom)
        lines.append(
            f"ATOM  {serial:>5} {padded}{'':1}{atom.comp_id[:3]:>3} "
            f"{(atom.auth_chain or atom.chain)[:1]:>1}"
            f"{(atom.auth_seq_id or atom.seq_id or 0):>4}{'':1}   "
            f"{atom.x:>8.3f}{atom.y:>8.3f}{atom.z:>8.3f}"
            f"{1.00:>6.2f}{(atom.b_iso_value()):>6.2f}{'':10}{element:>2}"
        )
    lines.append("END")
    destination.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return len(lines) - 1


def _with_extension(stem: Path, extension: str) -> Path:
    """Append an extension. `Path.with_suffix` would eat the `.463` of a feature id."""
    return stem.with_name(stem.name + extension)


def _element(atom: Atom) -> str:
    name = atom.atom_id.strip()
    for character in name:
        if character.isalpha():
            return character.upper()
    return "C"


# --- AutoDock typing, adapted from Ping Hu's bacteria_structure_docking.py ----------------

#: Carbons in aromatic rings are type A, not C; AutoDock scores them differently.
_AROMATIC_CARBONS = {
    "PHE": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TYR": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "HIS": {"CG", "CD2", "CE1"},
    "TRP": {"CG", "CD1", "CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
}


def autodock_type(atom_name: str, residue: str, element: str) -> str:
    if element == "C":
        return "A" if atom_name in _AROMATIC_CARBONS.get(residue, ()) else "C"
    if element == "N":
        return "NA" if residue == "HIS" and atom_name in {"ND1", "NE2"} else "N"
    return {"O": "OA", "S": "SA", "P": "P"}.get(element, element or "C")


def write_pdbqt_fallback(atoms: Iterable[Atom], destination: Path) -> int:
    """A PDBQT with AutoDock types, no hydrogens and zero charges. Reduced fidelity, stated."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    lines = []
    for serial, atom in enumerate(atoms, 1):
        element = _element(atom)
        name = atom.atom_id.strip()[:4]
        padded = f" {name:<3}" if len(name) < 4 else name
        lines.append(
            f"ATOM  {serial:>5} {padded}{'':1}{atom.comp_id[:3]:>3} "
            f"{(atom.auth_chain or atom.chain)[:1]:>1}"
            f"{(atom.auth_seq_id or atom.seq_id or 0):>4}{'':1}   "
            f"{atom.x:>8.3f}{atom.y:>8.3f}{atom.z:>8.3f}"
            f"{1.00:>6.2f}{atom.b_iso_value():>6.2f}    +0.000 "
            f"{autodock_type(name, atom.comp_id[:3], element):<2}"
        )
    lines.append("END")
    destination.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return len(lines) - 1


def run_preparation_tool(tool: str, pdb_path: Path, destination: Path) -> tuple[bool, str]:
    """Convert with a real preparation tool. Returns (wrote_a_file, message)."""
    if tool == "obabel":
        command = ["obabel", str(pdb_path), "-O", str(destination), "-xr"]
    else:
        base = str(destination)[: -len(".pdbqt")] if str(destination).endswith(".pdbqt") else str(destination)
        command = [tool, "--read_pdb", str(pdb_path), "-o", base, "-p"]
    try:
        finished = subprocess.run(command, capture_output=True, text=True, timeout=300)
    except (OSError, subprocess.SubprocessError) as exc:
        return False, f"{tool} failed: {exc}"
    stem = destination.name[: -len(".pdbqt")] if destination.name.endswith(".pdbqt") else destination.name
    candidates = [destination, destination.with_name(stem + "_rigid.pdbqt")]
    for candidate in candidates:
        if candidate.exists() and candidate.stat().st_size:
            if candidate != destination:
                candidate.replace(destination)
            return True, f"prepared with {tool}"
    detail = (finished.stderr or finished.stdout or "").strip().splitlines()
    return False, f"{tool} produced no PDBQT" + (f": {detail[-1]}" if detail else "")


def prepare_receptor(
    payload: bytes,
    *,
    entity_id: str,
    destination_stem: Path,
    run_dir: Path,
    allow_fallback: bool = False,
    find_tool: Callable[[], str | None] = available_tool,
) -> PreparedReceptor:
    """Write the receptor for one structure, and a PDBQT where that can be done honestly."""
    atoms = receptor_atoms(payload, entity_id)
    if not atoms:
        return PreparedReceptor(
            pdb_path=None,
            pdbqt_path=None,
            prepared_by="none",
            fidelity="none",
            reason="no polymer atoms could be read for this entity",
            chains=[],
            atoms=0,
        )

    pdb_path = _with_extension(destination_stem, ".pdb")
    written = write_receptor_pdb(atoms, pdb_path)
    chains = sorted({atom.auth_chain or atom.chain for atom in atoms})
    relative = str(pdb_path.relative_to(run_dir))

    pdbqt_path = _with_extension(destination_stem, ".pdbqt")
    tool = find_tool()
    if tool:
        ok, message = run_preparation_tool(tool, pdb_path, pdbqt_path)
        if ok:
            return PreparedReceptor(
                pdb_path=relative,
                pdbqt_path=str(pdbqt_path.relative_to(run_dir)),
                prepared_by=tool,
                fidelity="full",
                reason=message,
                chains=chains,
                atoms=written,
            )
        reason = message
    else:
        reason = (
            "no receptor preparation tool on PATH (" + ", ".join(PDBQT_TOOLS) + ")"
        )

    if not allow_fallback:
        return PreparedReceptor(
            pdb_path=relative,
            pdbqt_path=None,
            prepared_by="none",
            fidelity="pdb_only",
            reason=reason + "; receptor PDB written, PDBQT skipped rather than approximated",
            chains=chains,
            atoms=written,
        )

    write_pdbqt_fallback(atoms, pdbqt_path)
    return PreparedReceptor(
        pdb_path=relative,
        pdbqt_path=str(pdbqt_path.relative_to(run_dir)),
        prepared_by="builtin_fallback",
        fidelity="reduced",
        reason=(
            reason + "; built-in converter used: AutoDock atom types, but no added "
            "hydrogens and zero partial charges. Fit for a pipeline smoke test, not for a "
            "reported docking result"
        ),
        chains=chains,
        atoms=written,
    )


# --- search boxes -------------------------------------------------------------------------


def docking_box(
    payload: bytes,
    pocket: dict[str, Any],
    *,
    padding: float = BOX_PADDING,
    minimum_edge: float = MIN_BOX_EDGE,
) -> dict[str, Any] | None:
    """A search box around one pocket's own atoms, not a fixed cube at its centre.

    Drawn from the lining residues and the bound component together, so the box follows the
    shape of the real site. A cube centred on the ligand is the usual shortcut and it either
    clips an elongated site or wastes search volume on a compact one.
    """
    residues = pocket.get("residues") or []
    if not residues:
        return None
    wanted = {
        (str(residue.get("chain")), str(residue.get("auth_seq_id")))
        for residue in residues
        if isinstance(residue, dict)
    }
    if not wanted:
        return None

    component = str(pocket.get("component") or "")
    points = [
        atom
        for atom in read_atoms(payload)
        if atom.has_coordinates
        and (
            ((atom.auth_chain or atom.chain), str(atom.auth_seq_id)) in wanted
            or (atom.is_hetatm and component and atom.comp_id == component)
        )
    ]
    if not points:
        return None

    box: dict[str, Any] = {"padding": padding, "minimum_edge": minimum_edge}
    for axis, values in (
        ("x", [atom.x for atom in points]),
        ("y", [atom.y for atom in points]),
        ("z", [atom.z for atom in points]),
    ):
        low, high = min(values), max(values)
        size = max((high - low) + 2 * padding, minimum_edge)
        box[f"center_{axis}"] = round((low + high) / 2, 3)
        box[f"size_{axis}"] = round(size, 3)
    box["atoms"] = len(points)
    return box


def vina_config(box: dict[str, Any], receptor: str, *, exhaustiveness: int = 8,
                num_modes: int = 9) -> str:
    """The box as an AutoDock Vina config file.

    Written for M4 to use or ignore; M3 runs no docking. Exhaustiveness and the mode count
    are Vina's own defaults, recorded here so a run states them rather than inheriting them.
    """
    return (
        f"receptor = {receptor}\n"
        f"center_x = {box['center_x']}\ncenter_y = {box['center_y']}\n"
        f"center_z = {box['center_z']}\n"
        f"size_x = {box['size_x']}\nsize_y = {box['size_y']}\nsize_z = {box['size_z']}\n"
        f"exhaustiveness = {exhaustiveness}\nnum_modes = {num_modes}\n"
    )

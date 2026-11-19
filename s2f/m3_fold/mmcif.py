"""A small read-only mmCIF loop reader, shared by M3's quality gate and pocket finder.

Deliberately not a structure library. It reads the two loops M3 needs — ``_atom_site`` for
coordinates and confidence, ``_struct_asym`` for the chains an entity owns — and returns
rows of strings. A file it cannot parse yields no rows, so the caller fails an auditable
check rather than silently passing on a structure nobody read.
"""

from __future__ import annotations

import shlex
from dataclasses import dataclass


@dataclass(frozen=True)
class Loop:
    """One mmCIF ``loop_`` as column names plus rows of raw string values."""

    headers: list[str]
    rows: list[list[str]]

    def column(self, name: str) -> int | None:
        try:
            return self.headers.index(name)
        except ValueError:
            return None

    def value(self, row: list[str], name: str) -> str | None:
        index = self.column(name)
        if index is None or index >= len(row):
            return None
        token = row[index]
        return None if token in {".", "?"} else token


def read_loop(payload: bytes, required_column: str) -> Loop:
    """The first ``loop_`` carrying ``required_column``, or an empty loop.

    ``required_column`` is a full mmCIF tag, for example ``_atom_site.Cartn_x``. A tag ending
    in ``.`` matches any column in that category, which is how the atom-site loop is found
    whether or not the file carries coordinates. Values are returned exactly as written;
    quoting is resolved and numbers are not converted, because which columns are numeric
    differs between the loops this reads.
    """
    lines = payload.decode("utf-8", errors="replace").splitlines()
    index = 0
    while index < len(lines):
        if lines[index].strip() != "loop_":
            index += 1
            continue

        index += 1
        headers: list[str] = []
        while index < len(lines) and lines[index].lstrip().startswith("_"):
            headers.append(lines[index].split()[0])
            index += 1
        if required_column.endswith("."):
            if not any(header.startswith(required_column) for header in headers):
                continue
        elif required_column not in headers:
            continue

        tokens: list[str] = []
        while index < len(lines):
            line = lines[index].strip()
            if line in {"#", "loop_"} or line.startswith("_") or line.startswith("data_"):
                break
            if line:
                try:
                    tokens.extend(shlex.split(line, comments=False, posix=True))
                except ValueError:
                    return Loop([], [])
            index += 1

        width = len(headers)
        if not width or len(tokens) % width:
            return Loop([], [])
        rows = [tokens[start : start + width] for start in range(0, len(tokens), width)]
        return Loop(headers, rows)
    return Loop([], [])


@dataclass(frozen=True)
class Atom:
    """One atom-site row, reduced to what M3 asks of it."""

    is_hetatm: bool
    chain: str
    auth_chain: str
    comp_id: str
    seq_id: int | None
    auth_seq_id: str | None
    atom_id: str
    #: None when the file carries no coordinate columns; geometry skips such atoms.
    x: float | None
    y: float | None
    z: float | None
    b_factor: float | None

    @property
    def has_coordinates(self) -> bool:
        return None not in (self.x, self.y, self.z)


def read_atoms(payload: bytes, *, first_model_only: bool = True) -> list[Atom]:
    """Every atom of the first model, or an empty list if the file cannot be read."""
    loop = read_loop(payload, "_atom_site.")
    if not loop.rows:
        return []

    atoms: list[Atom] = []
    for row in loop.rows:
        if first_model_only:
            model = loop.value(row, "_atom_site.pdbx_PDB_model_num")
            if model is not None and model != "1":
                continue
        coordinates: list[float | None] = []
        for axis in ("Cartn_x", "Cartn_y", "Cartn_z"):
            token = loop.value(row, f"_atom_site.{axis}")
            try:
                coordinates.append(float(token) if token is not None else None)
            except ValueError:
                coordinates.append(None)
        x, y, z = coordinates
        chain = loop.value(row, "_atom_site.label_asym_id") or ""
        auth_chain = loop.value(row, "_atom_site.auth_asym_id") or chain
        seq_token = loop.value(row, "_atom_site.label_seq_id")
        try:
            seq_id = int(seq_token) if seq_token is not None else None
        except ValueError:
            seq_id = None
        b_token = loop.value(row, "_atom_site.B_iso_or_equiv")
        try:
            b_factor = float(b_token) if b_token is not None else None
        except ValueError:
            b_factor = None
        atoms.append(
            Atom(
                is_hetatm=(loop.value(row, "_atom_site.group_PDB") or "ATOM").upper() == "HETATM",
                chain=chain,
                auth_chain=auth_chain,
                comp_id=(loop.value(row, "_atom_site.label_comp_id") or "").upper(),
                seq_id=seq_id,
                auth_seq_id=loop.value(row, "_atom_site.auth_seq_id"),
                atom_id=loop.value(row, "_atom_site.label_atom_id") or "",
                x=x,
                y=y,
                z=z,
                b_factor=b_factor,
            )
        )
    return atoms


def entity_chains(payload: bytes, entity_id: str) -> set[str]:
    """``label_asym_id`` values belonging to one entity, from ``_struct_asym``.

    An RCSB entry often holds many copies and many different chains; a ligand bound to a
    neighbouring entity is not our homolog's site. An empty set means the file does not say,
    and the caller must decide whether to fall back to every chain.
    """
    loop = read_loop(payload, "_struct_asym.entity_id")
    if not loop.rows or not entity_id:
        return set()
    return {
        asym
        for row in loop.rows
        if (loop.value(row, "_struct_asym.entity_id") or "") == str(entity_id)
        and (asym := loop.value(row, "_struct_asym.id"))
    }

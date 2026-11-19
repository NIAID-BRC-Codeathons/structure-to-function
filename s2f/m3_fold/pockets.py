"""Binding sites for M4, transferred from what the homolog entry actually binds (issue #14).

The primary method needs no external tool and no network. M2 names the ligands and metals
bound in each qualifying PDB entry, and M3 has already downloaded that entry's coordinates,
so the site is read straight out of the file: the residues of *our* entity's chains within a
cutoff of the named component's atoms.

What this deliberately does not do:

- It does not invent a pocket. A protein whose homolog binds nothing gets no pocket from
  here and is reported as uncovered, never as "no pocket exists" (the distinction issue #90
  draws). P2Rank fills that gap when it is installed; `find_pockets` is where it plugs in.
- It does not transfer a site between different coordinate frames. A site is only read from
  the structure M3 collected, so no superposition is involved and none is implied. The
  `aligned_to_known_site` flag the contract asks for is therefore true exactly when the
  collected structure *is* the entry that binds the component.
- It does not treat a metal as a drug site. A metal marks a real site, and a predicted model
  will not contain the metal, which is the failure pitfalls.md #2 is about: those pockets
  carry `cofactor_required` so M4 and M6 can say so rather than dock into a hole.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import sqrt
from typing import Any, Iterable

from .mmcif import Atom, entity_chains, read_atoms

#: Residues with any atom within this distance of a bound component line its site. 4.5 A is
#: the usual first-shell contact distance; it is recorded in every pocket's `params`.
CONTACT_DISTANCE = 4.5

#: Fewer contacts than this is a surface brush, not a pocket.
MIN_POCKET_RESIDUES = 3

#: Below this mean pLDDT over the pocket's own residues the site is not trustworthy, whatever
#: the model's global score says: pitfalls.md says a pocket in a low-confidence region is not
#: a pocket. Matches QualityThresholds.min_local_plddt.
MIN_POCKET_PLDDT = 70.0


@dataclass(frozen=True)
class PocketResidue:
    chain: str
    seq_id: int | None
    auth_seq_id: str | None
    comp_id: str
    min_distance: float
    #: The mmCIF B_iso_or_equiv column as written: a real B-factor in an experimental entry,
    #: pLDDT in an AlphaFold model. Named for the column, not for either meaning.
    b_iso: float | None

    def as_dict(self) -> dict[str, Any]:
        record = {
            "chain": self.chain,
            "seq_id": self.seq_id,
            "auth_seq_id": self.auth_seq_id,
            "comp_id": self.comp_id,
            "min_distance": round(self.min_distance, 2),
        }
        if self.b_iso is not None:
            record["b_iso"] = round(self.b_iso, 2)
        return record


@dataclass
class Pocket:
    """One binding site, with the evidence that produced it."""

    pocket_id: str
    method: str
    component: str
    component_kind: str
    residues: list[PocketResidue]
    center: tuple[float, float, float] | None
    aligned_to_known_site: bool
    cofactor_required: bool
    cofactor_present: bool
    usable: bool = True
    reason: str = ""
    #: A detector's own score, where the method has one. Site transfer does not: the site is
    #: observed in the entry, not predicted, so there is nothing to score.
    score: float | None = None
    confidence: float | None = None
    params: dict[str, Any] = field(default_factory=dict)

    def as_dict(self) -> dict[str, Any]:
        return {
            # `id`, `method`, `residues`, `score` and `aligned_to_known_site` are the keys
            # docs/00-architecture.md fixed for this record; the rest is added evidence.
            "id": self.pocket_id,
            "method": self.method,
            "score": self.score,
            "component": self.component,
            "component_kind": self.component_kind,
            "residues": [residue.as_dict() for residue in self.residues],
            "residue_count": len(self.residues),
            "center": (
                [round(value, 3) for value in self.center] if self.center is not None else None
            ),
            "aligned_to_known_site": self.aligned_to_known_site,
            "cofactor_required": self.cofactor_required,
            "cofactor_present": self.cofactor_present,
            "usable": self.usable,
            "reason": self.reason,
            "confidence": round(self.confidence, 2) if self.confidence is not None else None,
            "params": dict(self.params),
        }


def _residue_key(atom: Atom) -> tuple[str, int | None, str | None, str]:
    return (atom.chain, atom.seq_id, atom.auth_seq_id, atom.comp_id)


def _distance(left: Atom, right: Atom) -> float:
    return sqrt(
        (left.x - right.x) ** 2 + (left.y - right.y) ** 2 + (left.z - right.z) ** 2
    )


def _component_groups(
    atoms: Iterable[Atom], names: set[str]
) -> dict[tuple[str, int | None, str | None, str], list[Atom]]:
    """Group the atoms of each named component copy, keyed by its residue identity."""
    groups: dict[tuple[str, int | None, str | None, str], list[Atom]] = {}
    for atom in atoms:
        if atom.is_hetatm and atom.comp_id in names and atom.has_coordinates:
            groups.setdefault(_residue_key(atom), []).append(atom)
    return groups


def find_holo_sites(
    payload: bytes,
    *,
    ligands: Iterable[str] = (),
    metals: Iterable[str] = (),
    entity_id: str = "",
    contact_distance: float = CONTACT_DISTANCE,
    min_residues: int = MIN_POCKET_RESIDUES,
    structure_is_template: bool = True,
) -> list[Pocket]:
    """Sites lining the components M2 named, strongest first.

    `ligands` and `metals` come from `report.json`, where M2 already dropped waters,
    cryoprotectants and glycans, so nothing here re-filters chemistry by name.

    `entity_id` restricts the lining residues to the chains of our homolog's entity. Without
    it every chain is eligible, which on a large assembly can return a neighbour's site.

    `structure_is_template` is false when the coordinates are a predicted model rather than
    the entry that binds the component. The site is then not present in the file at all, so
    no pocket can be read from it; the caller gets an empty list rather than a guess.
    """
    ligand_names = {name.upper() for name in ligands if name}
    metal_names = {name.upper() for name in metals if name}
    if not structure_is_template or not (ligand_names | metal_names):
        return []

    atoms = read_atoms(payload)
    if not atoms:
        return []

    chains = entity_chains(payload, entity_id) if entity_id else set()
    polymer = [
        atom
        for atom in atoms
        if not atom.is_hetatm and atom.has_coordinates and (not chains or atom.chain in chains)
    ]
    if not polymer:
        return []

    groups = _component_groups(atoms, ligand_names | metal_names)
    metals_in_file = {
        key[3] for key in _component_groups(atoms, metal_names)
    }

    pockets: list[Pocket] = []
    for key, component_atoms in groups.items():
        chain, _, auth_seq_id, comp_id = key
        kind = "ligand" if comp_id in ligand_names else "metal"

        nearest: dict[tuple[str, int | None, str | None, str], tuple[float, Atom]] = {}
        for atom in polymer:
            best = None
            for other in component_atoms:
                distance = _distance(atom, other)
                if best is None or distance < best:
                    best = distance
            if best is None or best > contact_distance:
                continue
            residue_key = _residue_key(atom)
            if residue_key not in nearest or best < nearest[residue_key][0]:
                nearest[residue_key] = (best, atom)

        if len(nearest) < min_residues:
            continue

        residues = sorted(
            (
                PocketResidue(
                    chain=atom.auth_chain or atom.chain,
                    seq_id=atom.seq_id,
                    auth_seq_id=atom.auth_seq_id,
                    comp_id=atom.comp_id,
                    min_distance=distance,
                    b_iso=atom.b_factor,
                )
                for distance, atom in nearest.values()
            ),
            key=lambda residue: (residue.chain, residue.seq_id or 0, residue.auth_seq_id or ""),
        )
        center = (
            sum(atom.x for atom in component_atoms) / len(component_atoms),
            sum(atom.y for atom in component_atoms) / len(component_atoms),
            sum(atom.z for atom in component_atoms) / len(component_atoms),
        )
        # A metal site needs its metal; a ligand site whose entry also holds a metal may need
        # that too. Either way the question is whether a model would carry it (pitfalls #2).
        needs_cofactor = kind == "metal" or bool(metals_in_file)
        pockets.append(
            Pocket(
                pocket_id=f"{comp_id}_{chain}{auth_seq_id or ''}",
                method="holo_site_transfer",
                component=comp_id,
                component_kind=kind,
                residues=residues,
                center=center,
                # True because the site was read from the entry that binds the component,
                # not aligned from somewhere else.
                aligned_to_known_site=True,
                cofactor_required=needs_cofactor,
                cofactor_present=True,
                params={
                    "contact_distance": contact_distance,
                    "min_residues": min_residues,
                    "entity_id": entity_id or None,
                    "entity_chains": sorted(chains) or None,
                    "metals_in_entry": sorted(metals_in_file) or None,
                },
            )
        )

    pockets.sort(
        key=lambda pocket: (pocket.component_kind != "ligand", -len(pocket.residues))
    )
    return pockets


def reject_low_confidence(
    pockets: list[Pocket],
    *,
    source: str,
    min_plddt: float = MIN_POCKET_PLDDT,
) -> list[Pocket]:
    """Mark pockets whose own residues are poorly modelled as unusable.

    Only meaningful where the B-factor column holds pLDDT, which is AlphaFold's convention;
    for an experimental entry the column is a real B-factor and no confidence claim is made.
    """
    for pocket in pockets:
        scored = [
            residue.b_iso for residue in pocket.residues if residue.b_iso is not None
        ]
        if source != "afdb" or not scored:
            continue
        mean_plddt = sum(scored) / len(scored)
        pocket.confidence = mean_plddt
        pocket.params["min_pocket_plddt"] = min_plddt
        if mean_plddt < min_plddt:
            pocket.usable = False
            pocket.reason = (
                f"mean pLDDT {mean_plddt:.1f} over the pocket's own residues is below "
                f"{min_plddt:.0f}; a pocket in a low-confidence region is not a pocket"
            )
    return pockets


def find_pockets(
    payload: bytes,
    *,
    source: str,
    holo_template: dict[str, Any] | None,
    structure_accession: str = "",
    contact_distance: float = CONTACT_DISTANCE,
    min_plddt: float = MIN_POCKET_PLDDT,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Every pocket M3 can define for one collected structure, plus how it got there.

    Returns the pocket records for `structures[].pockets` and a summary naming the method
    used and why, so a protein with no pocket says which case it fell into rather than
    looking like one nobody looked at.
    """
    template = holo_template or {}
    entity = str(template.get("pdb_entity") or "")
    # Only the entry that binds the component holds the site. M3 collects one structure per
    # protein, so this is true when that structure is the holo entry itself.
    is_template = bool(entity) and _same_entry(entity, structure_accession)

    pockets = find_holo_sites(
        payload,
        ligands=template.get("ligands") or (),
        metals=template.get("metals") or (),
        entity_id=entity.split("_", 1)[1] if "_" in entity else "",
        contact_distance=contact_distance,
        structure_is_template=is_template,
    )
    reject_low_confidence(pockets, source=source, min_plddt=min_plddt)

    usable = [pocket for pocket in pockets if pocket.usable]
    summary: dict[str, Any] = {
        "method": "holo_site_transfer",
        "pockets": len(pockets),
        "usable": len(usable),
        "holo_template": entity or None,
        "template_is_collected_structure": is_template,
    }
    if not pockets:
        summary["status"] = "no_pocket"
        summary["reason"] = _no_pocket_reason(template, entity, is_template)
    else:
        summary["status"] = "ok" if usable else "all_rejected"
    return [pocket.as_dict() for pocket in pockets], summary


def _same_entry(entity_id: str, accession: str) -> bool:
    """Whether an entity such as ``1YDX_1`` names the same entry as a collected accession."""
    if not accession:
        return False
    return entity_id.split("_", 1)[0].upper() == accession.split("_", 1)[0].upper()


def _no_pocket_reason(template: dict[str, Any], entity: str, is_template: bool) -> str:
    if not entity:
        return (
            "no qualifying homolog binds a ligand or metal; a pocket detector such as "
            "P2Rank is required for this protein and none ran"
        )
    if not is_template:
        return (
            f"the binding site is in {entity}, which is not the structure collected here; "
            "transferring it needs a superposition M3 does not yet do"
        )
    bound = ", ".join([*(template.get("ligands") or []), *(template.get("metals") or [])])
    return f"no residue of our entity lines {bound or 'the bound components'} in this entry"

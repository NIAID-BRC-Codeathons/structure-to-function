"""The query genome's taxonomy, taken from M1 rather than guessed (issue #55).

M2's same-species and same-genus flags need to know what organism the query genome is. M1
already determines that and writes it to ``report.json``; the string heuristics in
``bvbrc_input`` exist only for the standalone entry point, where there is no ``genome`` section
to read (``docs/02a-m2-pdb-evidence.md``, "Input contract").

Two things the lineage gives that spelling cannot:

- **The authoritative genus.** NCBI lineage names carry the current name, so a renamed genus
  needs no stem-similarity test: BV-BRC's ``Mycoplasma genitalium`` has lineage
  ``… Mycoplasmoides, Mycoplasmoides genitalium …`` — the same spelling the PDB uses.
- **Rank-resolved identity.** ``taxon_lineage_ids`` carries the species-level taxon alongside
  the strain, so a hit annotated with the species (2097) matches a strain-level query
  (243273.25) without ID equality ever holding. That mismatch is exactly the case @cmmann21
  raised on #12.

A hit whose taxonomy we cannot resolve is left unflagged rather than guessed: an unknown
organism is not evidence of a foreign one.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .bvbrc_input import organism_key, same_genus as same_genus_by_name
from .bvbrc_input import same_species as same_species_by_name

SOURCE_M1 = "m1_genome_section"
SOURCE_OVERRIDE = "--organism"
SOURCE_HEADERS = "fasta_headers"
SOURCE_UNDETERMINED = "undetermined"


@dataclass
class QueryTaxonomy:
    """What we know about the query genome, and where it came from."""

    scientific_name: str = ""
    taxon_id: int | None = None
    lineage_names: list[str] = field(default_factory=list)
    lineage_ids: list[int] = field(default_factory=list)
    source: str = SOURCE_UNDETERMINED

    @property
    def determined(self) -> bool:
        return bool(self.scientific_name)

    @property
    def has_lineage(self) -> bool:
        return bool(self.lineage_names)

    @property
    def genus(self) -> str:
        """The genus as NCBI spells it, falling back to the first token of the name."""
        if self.lineage_names:
            # The lineage runs kingdom -> strain; the genus is the last entry that is a single
            # word, since species and strain names carry a space.
            for name in reversed(self.lineage_names):
                if name and " " not in name.strip():
                    return name.strip().lower()
        return organism_key(self.scientific_name)[0]

    @property
    def species_name(self) -> str:
        genus, epithet = organism_key(self.scientific_name)
        return f"{genus} {epithet}".strip() if epithet else ""

    def species_taxon_ids(self) -> set[int]:
        """Taxon IDs at species rank and below — what a hit may legitimately carry."""
        ids: set[int] = set()
        if self.lineage_names and self.lineage_ids:
            for name, taxon in zip(self.lineage_names, self.lineage_ids):
                if name and " " in name.strip():  # species or below
                    ids.add(int(taxon))
        if self.taxon_id:
            ids.add(int(self.taxon_id))
        return ids

    def genus_taxon_id(self) -> int | None:
        if not (self.lineage_names and self.lineage_ids):
            return None
        for name, taxon in zip(reversed(self.lineage_names), reversed(self.lineage_ids)):
            if name and " " not in name.strip():
                return int(taxon)
        return None

    # --- the comparisons the flags are built on ---------------------------------------

    def same_species_as(self, hit_organism: str, hit_taxon_id: int | None = None) -> bool:
        if hit_taxon_id and hit_taxon_id in self.species_taxon_ids():
            return True
        if not self.determined or not hit_organism:
            return False
        if self.has_lineage:
            # Authoritative genus spelling, so this is a plain comparison rather than a
            # stem-similarity guess.
            hit_genus, hit_epithet = organism_key(hit_organism)
            query_epithet = organism_key(self.scientific_name)[1]
            if not hit_epithet or not query_epithet:
                return False
            return hit_epithet == query_epithet and hit_genus == self.genus
        return same_species_by_name(self.scientific_name, hit_organism)

    def same_genus_as(self, hit_organism: str, hit_taxon_id: int | None = None) -> bool:
        if self.same_species_as(hit_organism, hit_taxon_id):
            return True
        if not self.determined or not hit_organism:
            return False
        if self.has_lineage:
            return organism_key(hit_organism)[0] == self.genus
        return same_genus_by_name(self.scientific_name, hit_organism)

    def as_dict(self) -> dict[str, Any]:
        return {
            "scientific_name": self.scientific_name or None,
            "taxon_id": self.taxon_id,
            "genus": self.genus or None,
            "source": self.source,
            "lineage_depth": len(self.lineage_names),
        }


def from_report(run_dir: str | Path) -> QueryTaxonomy:
    """Read M1's `genome` section. Returns an undetermined record when there is none."""
    path = Path(run_dir) / "report.json"
    if not path.exists():
        return QueryTaxonomy()
    try:
        report = json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return QueryTaxonomy()

    genome = report.get("genome") or {}
    taxonomy = genome.get("taxonomy") or {}
    name = (taxonomy.get("scientific_name") or "").strip()
    lineage_names = [str(n) for n in (taxonomy.get("lineage_names") or []) if n]
    lineage_ids_raw = taxonomy.get("lineage_ids") or genome.get("taxon_lineage_ids") or []
    lineage_ids: list[int] = []
    for value in lineage_ids_raw:
        try:
            lineage_ids.append(int(value))
        except (TypeError, ValueError):
            continue
    if len(lineage_ids) != len(lineage_names):
        lineage_ids = []  # unpaired lineage is unusable; fall back to names alone

    taxon_id = genome.get("taxon_id")
    try:
        taxon_id = int(taxon_id) if taxon_id else None
    except (TypeError, ValueError):
        taxon_id = None

    if not name and not taxon_id:
        return QueryTaxonomy()
    return QueryTaxonomy(
        scientific_name=name,
        taxon_id=taxon_id,
        lineage_names=lineage_names,
        lineage_ids=lineage_ids,
        source=SOURCE_M1,
    )


def resolve(run_dir: str | Path, *, override: str = "", header_organism: str = "") -> QueryTaxonomy:
    """M1's genome section, then --organism, then the FASTA headers, then undetermined."""
    if override:
        return QueryTaxonomy(scientific_name=override.strip(), source=SOURCE_OVERRIDE)
    from_m1 = from_report(run_dir)
    if from_m1.determined:
        return from_m1
    if header_organism:
        return QueryTaxonomy(scientific_name=header_organism.strip(), source=SOURCE_HEADERS)
    return QueryTaxonomy()

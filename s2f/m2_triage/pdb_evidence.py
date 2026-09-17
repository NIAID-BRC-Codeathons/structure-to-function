"""RCSB PDB sequence search and entity metadata for M2 (issue #8).

Parameters and cutoffs are documented in ``docs/02a-m2-pdb-evidence.md``; change them there in
the same commit. Retrieval status per protein is ``found``, ``no-hit`` or ``query-failed`` — a
failed query is never recorded as an absence of hits.
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from typing import Any, Iterable

from ..common.http import CachedJsonClient, HttpError

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL_URL = "https://data.rcsb.org/graphql"

EVALUE_CUTOFF = 1e-3
IDENTITY_CUTOFF = 0.25
ROWS_PER_QUERY = 25
GRAPHQL_BATCH = 100

# Crystallization additives, buffers, cryoprotectants and ions: present in the entry but not
# evidence of a functional ligand. Metals are excluded from the ligand bonus but still reported,
# since a metal in a homolog's site matters to M3 (pitfalls.md #2).
ADDITIVES = {
    "HOH", "DOD", "GOL", "EDO", "PEG", "PG4", "PGE", "PGO", "P6G", "1PE", "2PE", "MPD", "MRD",
    "SO4", "PO4", "ACT", "ACY", "FMT", "CIT", "FLC", "TRS", "TAR", "MES", "EPE", "IMD", "DMS",
    "DMF", "URE", "NH4", "AZI", "BME", "OXL", "SCN", "NO3", "CO3", "BCT", "IOD", "BR", "CL",
    "FLU", "UNX", "UNL",
}
METALS = {
    "NA", "K", "MG", "CA", "ZN", "MN", "FE", "FE2", "NI", "CU", "CU1", "CO", "CD", "HG", "PT",
    "AU", "AG", "PB", "SR", "BA", "CS", "RB", "LI", "AL", "W", "MO", "V",
}
# N-/O-linked glycans and common sugar modifications: post-translational decoration, not binders.
GLYCANS = {
    "NAG", "NDG", "BMA", "MAN", "BGC", "GLC", "GAL", "GLA", "FUC", "FUL", "XYS", "XYP", "SIA",
    "NGA", "A2G", "RIB", "SUC",
}


@dataclass
class SequenceHit:
    """One RCSB polymer-entity hit with its alignment and entity metadata."""

    entity_id: str
    identity: float
    evalue: float
    bitscore: float
    query_beg: int
    query_end: int
    query_length: int
    subject_length: int
    description: str = ""
    experimental_method: str = ""
    resolution: float | None = None
    organism: str = ""
    taxonomy_id: int | None = None
    uniprot_ids: list[str] = field(default_factory=list)
    ligands: list[str] = field(default_factory=list)
    metals: list[str] = field(default_factory=list)
    partner_entities: int = 0
    has_ec: bool = False
    has_go: bool = False
    annotation_types: list[str] = field(default_factory=list)

    @property
    def entry_id(self) -> str:
        return self.entity_id.split("_")[0]

    @property
    def query_coverage(self) -> float:
        if self.query_length <= 0:
            return 0.0
        covered = self.query_end - self.query_beg + 1
        return max(0.0, min(1.0, covered / self.query_length))

    @property
    def has_ligand(self) -> bool:
        return bool(self.ligands)

    @property
    def has_partner(self) -> bool:
        return self.partner_entities > 0


@dataclass
class SearchResult:
    """Search outcome for one distinct sequence."""

    status: str  # found | no-hit | query-failed
    hits: list[SequenceHit] = field(default_factory=list)
    error: str = ""


def search_payload(sequence: str, *, rows: int = ROWS_PER_QUERY) -> dict[str, Any]:
    """Build the RCSB sequence-search request documented in 02a-m2-pdb-evidence.md."""
    return {
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "sequence_type": "protein",
                "value": sequence,
                "evalue_cutoff": EVALUE_CUTOFF,
                "identity_cutoff": IDENTITY_CUTOFF,
            },
        },
        "return_type": "polymer_entity",
        "request_options": {
            "paginate": {"start": 0, "rows": rows},
            "scoring_strategy": "sequence",
            "results_verbosity": "verbose",
        },
    }


def parse_search_response(payload: dict[str, Any]) -> list[SequenceHit]:
    """Turn one search response into hits, keeping the best match_context per entity."""
    hits: list[SequenceHit] = []
    for result in payload.get("result_set") or []:
        entity_id = str(result.get("identifier") or "")
        if not entity_id:
            continue
        contexts = [
            context
            for service in result.get("services") or []
            for node in service.get("nodes") or []
            for context in node.get("match_context") or []
        ]
        if not contexts:
            continue
        best = max(contexts, key=lambda c: float(c.get("bitscore") or 0.0))
        hits.append(
            SequenceHit(
                entity_id=entity_id,
                identity=float(best.get("sequence_identity") or 0.0),
                evalue=float(best.get("evalue") if best.get("evalue") is not None else 1.0),
                bitscore=float(best.get("bitscore") or 0.0),
                query_beg=int(best.get("query_beg") or 1),
                query_end=int(best.get("query_end") or 0),
                query_length=int(best.get("query_length") or 0),
                subject_length=int(best.get("subject_length") or 0),
            )
        )
    return hits


ENTITY_QUERY = """
query($ids: [String!]!) {
  polymer_entities(entity_ids: $ids) {
    rcsb_id
    rcsb_polymer_entity { pdbx_description rcsb_ec_lineage { id } }
    rcsb_polymer_entity_container_identifiers { uniprot_ids }
    rcsb_entity_source_organism { ncbi_scientific_name ncbi_taxonomy_id }
    rcsb_polymer_entity_annotation { type }
    entry {
      rcsb_entry_info { resolution_combined experimental_method }
      nonpolymer_entities { pdbx_entity_nonpoly { comp_id } }
      polymer_entities { rcsb_id }
    }
  }
}
"""


@dataclass
class EntityMetadata:
    description: str = ""
    experimental_method: str = ""
    resolution: float | None = None
    organism: str = ""
    taxonomy_id: int | None = None
    uniprot_ids: list[str] = field(default_factory=list)
    ligands: list[str] = field(default_factory=list)
    metals: list[str] = field(default_factory=list)
    partner_entities: int = 0
    has_ec: bool = False
    has_go: bool = False
    annotation_types: list[str] = field(default_factory=list)


def classify_components(comp_ids: Iterable[str]) -> tuple[list[str], list[str]]:
    """Split entry components into (functional ligands, metals), dropping additives/glycans."""
    ligands: list[str] = []
    metals: list[str] = []
    for raw in comp_ids:
        comp = (raw or "").strip().upper()
        if not comp:
            continue
        if comp in METALS:
            if comp not in metals:
                metals.append(comp)
            continue
        if comp in ADDITIVES or comp in GLYCANS:
            continue
        if comp not in ligands:
            ligands.append(comp)
    return ligands, metals


def parse_entity_metadata(payload: dict[str, Any]) -> dict[str, EntityMetadata]:
    """Map GraphQL results by rcsb_id; unknown IDs are dropped by the API, not positional."""
    out: dict[str, EntityMetadata] = {}
    entities = ((payload.get("data") or {}).get("polymer_entities")) or []
    for entity in entities:
        if not entity:
            continue
        entity_id = str(entity.get("rcsb_id") or "")
        if not entity_id:
            continue
        core = entity.get("rcsb_polymer_entity") or {}
        identifiers = entity.get("rcsb_polymer_entity_container_identifiers") or {}
        organisms = entity.get("rcsb_entity_source_organism") or []
        annotations = entity.get("rcsb_polymer_entity_annotation") or []
        entry = entity.get("entry") or {}
        info = entry.get("rcsb_entry_info") or {}
        resolutions = info.get("resolution_combined") or []
        comp_ids = [
            ((nonpoly.get("pdbx_entity_nonpoly") or {}).get("comp_id") or "")
            for nonpoly in (entry.get("nonpolymer_entities") or [])
        ]
        ligands, metals = classify_components(comp_ids)
        annotation_types = sorted({str(a.get("type")) for a in annotations if a.get("type")})
        out[entity_id.upper()] = EntityMetadata(
            description=core.get("pdbx_description") or "",
            experimental_method=info.get("experimental_method") or "",
            resolution=float(resolutions[0]) if resolutions else None,
            organism=(organisms[0].get("ncbi_scientific_name") or "") if organisms else "",
            taxonomy_id=(organisms[0].get("ncbi_taxonomy_id") if organisms else None),
            uniprot_ids=[str(u) for u in (identifiers.get("uniprot_ids") or [])],
            ligands=ligands,
            metals=metals,
            partner_entities=max(0, len(entry.get("polymer_entities") or []) - 1),
            has_ec=bool(core.get("rcsb_ec_lineage")),
            has_go="GO" in annotation_types,
            annotation_types=annotation_types,
        )
    return out


class PdbEvidenceClient:
    """Sequence search plus batched entity metadata, cached and retried."""

    def __init__(self, client: CachedJsonClient, *, rows: int = ROWS_PER_QUERY) -> None:
        self.client = client
        self.rows = rows
        self.metadata_failures: list[dict[str, Any]] = []

    def search_sequence(self, sequence: str) -> SearchResult:
        """Search one sequence. HTTP 204 (empty body) is RCSB's no-hit response."""
        try:
            payload = self.client.post_json("rcsb_sequence", SEARCH_URL, search_payload(sequence, rows=self.rows))
        except HttpError as exc:
            return SearchResult(status="query-failed", error=str(exc))
        hits = parse_search_response(payload or {})
        if not hits:
            return SearchResult(status="no-hit")
        return SearchResult(status="found", hits=hits)

    def search_sequences(
        self, sequences: dict[str, str], *, workers: int = 4
    ) -> dict[str, SearchResult]:
        """Search distinct sequences concurrently, keyed by the caller's sequence hash."""
        keys = list(sequences)
        if not keys:
            return {}
        if workers <= 1:
            return {key: self.search_sequence(sequences[key]) for key in keys}
        with ThreadPoolExecutor(max_workers=workers) as pool:
            results = list(pool.map(lambda key: self.search_sequence(sequences[key]), keys))
        return dict(zip(keys, results))

    def fetch_entity_metadata(self, entity_ids: Iterable[str]) -> dict[str, EntityMetadata]:
        """Fetch metadata for hit entities in batches of GRAPHQL_BATCH."""
        unique = sorted({str(entity_id).upper() for entity_id in entity_ids if entity_id})
        metadata: dict[str, EntityMetadata] = {}
        for start in range(0, len(unique), GRAPHQL_BATCH):
            batch = unique[start : start + GRAPHQL_BATCH]
            payload = {"query": ENTITY_QUERY, "variables": {"ids": batch}}
            try:
                response = self.client.post_json("rcsb_graphql", GRAPHQL_URL, payload)
            except HttpError as exc:
                self.metadata_failures.append({"entity_ids": batch, "error": str(exc)})
                continue
            if (response or {}).get("errors"):
                self.metadata_failures.append(
                    {"entity_ids": batch, "error": str(response["errors"])[:500]}
                )
            metadata.update(parse_entity_metadata(response or {}))
        return metadata


def apply_metadata(hits: Iterable[SequenceHit], metadata: dict[str, EntityMetadata]) -> None:
    """Copy entity metadata onto hits in place; hits with no metadata keep search-only fields."""
    for hit in hits:
        meta = metadata.get(hit.entity_id.upper())
        if meta is None:
            continue
        hit.description = meta.description
        hit.experimental_method = meta.experimental_method
        hit.resolution = meta.resolution
        hit.organism = meta.organism
        hit.taxonomy_id = meta.taxonomy_id
        hit.uniprot_ids = list(meta.uniprot_ids)
        hit.ligands = list(meta.ligands)
        hit.metals = list(meta.metals)
        hit.partner_entities = meta.partner_entities
        hit.has_ec = meta.has_ec
        hit.has_go = meta.has_go
        hit.annotation_types = list(meta.annotation_types)

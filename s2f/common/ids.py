"""Identifier mapping: BV-BRC feature_id <-> UniProt <-> PDB <-> ChEMBL <-> gene name (issue #3).

This is the only place in the pipeline that maps identifiers (pitfall #11). Modules import it;
nobody writes their own.

Routes are tried in order and the one that produced a mapping is recorded on every result, so a
weak mapping is never mistaken for a strong one:

0. ``proteome_checksum`` / ``proteome_locus`` — one bulk fetch of the genome's taxon gives every
   UniProt entry with its CRC64 checksum and ordered locus names, so most proteins are matched
   locally with no per-protein request. This is the fast path: a 542-protein genome resolves in
   one call instead of ~1,800.
1. ``sequence_crc64`` — the protein's exact sequence, as a UniProt CRC64 checksum, looked up in
   UniParc. Same exactness as ``proteome_checksum``, for proteins the taxon index does not cover.
2. ``patric_crossref`` — UniParc indexes BV-BRC/PATRIC feature IDs directly, so the feature ID
   itself can be looked up when the sequence has drifted.
3. ``locus_tag`` — UniProt search by locus tag within the genome's taxon, verified against the
   entry's ordered locus names. UniProt writes ``MG401`` where BV-BRC writes ``MG_401``, so the
   comparison ignores punctuation.
4. ``pdb_hit`` — UniProt accessions the caller already has from PDB sequence hits (M2). A homolog
   structure's accession is **not** our protein, so this is recorded as a related accession, never
   as the protein's own identity.

When no route succeeds the result is an explicit unmapped record with ``uniprot = None``. A
plausible-looking wrong mapping is worse than none.
"""

from __future__ import annotations

import re
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass, field
from typing import Any, Iterable, Sequence

from .http import CachedJsonClient, HttpError

UNIPARC_SEARCH = "https://rest.uniprot.org/uniparc/search"
UNIPARC_ENTRY = "https://rest.uniprot.org/uniparc/{upi}.json"
UNIPROTKB_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
UNIPROTKB_STREAM = "https://rest.uniprot.org/uniprotkb/stream"

UNIPROTKB_FIELDS = "accession,id,protein_name,gene_names,organism_id,xref_pdb,xref_chembl,xref_refseq"

# CRC64-ISO, reflected, as UniProt stores in `sequence.crc64`. Verified against P0DTC2.
_CRC64_POLY = 0xD800000000000000
_CRC64_TABLE: list[int] = []
for _byte in range(256):
    _part = _byte
    for _ in range(8):
        _part = (_part >> 1) ^ (_CRC64_POLY if _part & 1 else 0)
    _CRC64_TABLE.append(_part)

ROUTE_PROTEOME_SEQUENCE = "proteome_checksum"
ROUTE_PROTEOME_LOCUS = "proteome_locus"
ROUTE_SEQUENCE = "sequence_crc64"
ROUTE_PATRIC = "patric_crossref"
ROUTE_LOCUS = "locus_tag"
ROUTE_PDB_HIT = "pdb_hit"
ROUTE_UNMAPPED = "unmapped"

_PUNCTUATION = re.compile(r"[^A-Z0-9]")


def crc64(sequence: str) -> str:
    """UniProt's CRC64 checksum for a sequence, uppercase hex."""
    crc = 0
    for char in sequence.strip().upper().encode("ascii", errors="ignore"):
        crc = _CRC64_TABLE[(crc ^ char) & 0xFF] ^ (crc >> 8)
    return f"{crc:016X}"


def _normalize_locus(value: str) -> str:
    """``MG_401`` and ``MG401`` are the same locus written two ways."""
    return _PUNCTUATION.sub("", (value or "").upper())


@dataclass
class Xrefs:
    """Stable cross-reference record written onto a protein.

    ``uniprot`` is the protein's own accession, or ``None`` when nothing mapped.
    ``related_uniprot`` holds accessions of homologs (e.g. from PDB hits) and is never promoted
    into ``uniprot``.
    """

    feature_id: str
    uniprot: str | None = None
    uniprot_entry_name: str | None = None
    uniparc: str | None = None
    gene_name: str | None = None
    locus_tag: str | None = None
    taxon_id: int | None = None
    protein_name: str | None = None
    refseq: list[str] = field(default_factory=list)
    pdb: list[str] = field(default_factory=list)
    chembl: list[str] = field(default_factory=list)
    related_uniprot: list[str] = field(default_factory=list)
    route: str = ROUTE_UNMAPPED
    note: str = ""

    @property
    def mapped(self) -> bool:
        return self.uniprot is not None

    def as_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["mapped"] = self.mapped
        return payload


def _select_uniprot_crossref(
    crossrefs: Sequence[dict[str, Any]], *, feature_id: str = "", taxon_id: int | None = None
) -> tuple[str | None, str]:
    """Pick one active UniProtKB accession from a UniParc entry, and say why.

    Preference: the entry that also lists our PATRIC feature ID, then a taxon match, then
    Swiss-Prot over TrEMBL. Ambiguity is reported in the note rather than hidden.
    """
    candidates = [
        ref
        for ref in crossrefs
        if str(ref.get("database", "")).startswith("UniProtKB") and ref.get("active")
    ]
    if not candidates:
        return None, "UniParc entry has no active UniProtKB cross-reference"

    patric_ids = {
        str(ref.get("id", ""))
        for ref in crossrefs
        if str(ref.get("database", "")).upper() in {"PATRIC", "SEED"} and ref.get("active")
    }
    feature_match = bool(feature_id) and feature_id in patric_ids

    def rank(ref: dict[str, Any]) -> tuple[int, int]:
        same_taxon = (ref.get("organism") or {}).get("taxonId") == taxon_id if taxon_id else False
        swissprot = str(ref.get("database")) == "UniProtKB/Swiss-Prot"
        return (0 if same_taxon else 1, 0 if swissprot else 1)

    ordered = sorted(candidates, key=rank)
    chosen = ordered[0]
    notes = []
    if feature_match:
        notes.append("UniParc entry lists this feature ID")
    if len(ordered) > 1:
        others = ", ".join(str(ref.get("id")) for ref in ordered[1:4])
        notes.append(f"{len(ordered)} active accessions; also {others}")
    return str(chosen.get("id")), "; ".join(notes)


@dataclass
class TaxonIndex:
    """Every UniProt entry for one taxon, indexed by sequence checksum and locus tag."""

    taxon_id: int
    entries: int = 0
    by_checksum: dict[str, dict[str, Any]] = field(default_factory=dict)
    by_locus: dict[str, dict[str, Any]] = field(default_factory=dict)


def build_taxon_index(taxon_id: int, results: Sequence[dict[str, Any]]) -> TaxonIndex:
    """Index a taxon's UniProt entries. First entry wins on a duplicate key."""
    index = TaxonIndex(taxon_id=taxon_id, entries=len(results))
    for entry in results:
        checksum = ((entry.get("sequence") or {}).get("crc64") or "").upper()
        if checksum:
            index.by_checksum.setdefault(checksum, entry)
        for gene in entry.get("genes") or []:
            for oln in gene.get("orderedLocusNames") or []:
                key = _normalize_locus(oln.get("value", ""))
                if key:
                    index.by_locus.setdefault(key, entry)
    return index


class IdMapper:
    """Map BV-BRC proteins to UniProt and onward, with every response cached."""

    def __init__(self, client: CachedJsonClient, *, taxon_id: int | None = None) -> None:
        self.client = client
        self.taxon_id = taxon_id
        self.failures: list[dict[str, str]] = []
        self._indexes: dict[int, TaxonIndex] = {}

    def taxon_index(self, taxon_id: int | None = None) -> TaxonIndex | None:
        """Fetch (once, then cached) every UniProt entry for the genome's taxon."""
        taxon = taxon_id if taxon_id is not None else self.taxon_id
        if not taxon:
            return None
        if taxon in self._indexes:
            return self._indexes[taxon]
        try:
            payload = self.client.get_json(
                "uniprotkb_stream",
                UNIPROTKB_STREAM,
                {
                    "query": f"taxonomy_id:{taxon}",
                    "format": "json",
                    "fields": UNIPROTKB_FIELDS + ",sequence",
                },
            )
        except HttpError as exc:
            self.failures.append({"feature_id": f"taxon:{taxon}", "error": str(exc)})
            self._indexes[taxon] = TaxonIndex(taxon_id=taxon)
            return self._indexes[taxon]
        index = build_taxon_index(taxon, payload.get("results") or [])
        self._indexes[taxon] = index
        return index

    # --- individual routes -------------------------------------------------

    def _uniparc_entry(self, upi: str) -> dict[str, Any]:
        return self.client.get_json("uniparc_entry", UNIPARC_ENTRY.format(upi=upi))

    def _uniparc_search(self, query: str) -> str | None:
        payload = self.client.get_json(
            "uniparc_search", UNIPARC_SEARCH, {"query": query, "format": "json", "size": 1}
        )
        results = payload.get("results") or []
        return str(results[0].get("uniParcId")) if results else None

    def _by_sequence(self, sequence: str) -> str | None:
        if not sequence:
            return None
        return self._uniparc_search(f"checksum:{crc64(sequence)}")

    def _by_feature_id(self, feature_id: str) -> str | None:
        # Quoted and database-scoped: an unquoted feature ID matches loosely.
        return self._uniparc_search(f'database:PATRIC AND dbid:"{feature_id}"')

    def _by_locus_tag(self, locus_tag: str, taxon_id: int | None) -> tuple[str | None, str]:
        if not locus_tag or not taxon_id:
            return None, ""
        payload = self.client.get_json(
            "uniprotkb_search",
            UNIPROTKB_SEARCH,
            {
                "query": f"{locus_tag} AND taxonomy_id:{taxon_id}",
                "format": "json",
                "fields": UNIPROTKB_FIELDS,
                "size": 5,
            },
        )
        wanted = _normalize_locus(locus_tag)
        for entry in payload.get("results") or []:
            names = {
                _normalize_locus(oln.get("value", ""))
                for gene in entry.get("genes") or []
                for oln in gene.get("orderedLocusNames") or []
            }
            if wanted in names:
                return str(entry.get("primaryAccession")), "locus tag verified against entry"
        return None, ""

    def _fetch_entry(self, accession: str) -> dict[str, Any] | None:
        payload = self.client.get_json(
            "uniprotkb_search",
            UNIPROTKB_SEARCH,
            {
                "query": f"accession:{accession}",
                "format": "json",
                "fields": UNIPROTKB_FIELDS,
                "size": 1,
            },
        )
        results = payload.get("results") or []
        return results[0] if results else None

    # --- assembly ----------------------------------------------------------

    def _enrich(self, xrefs: Xrefs) -> Xrefs:
        """Fill gene name, PDB, ChEMBL and RefSeq from the resolved accession."""
        if not xrefs.uniprot:
            return xrefs
        entry = self._fetch_entry(xrefs.uniprot)
        if entry is None:
            xrefs.note = "; ".join(filter(None, [xrefs.note, "accession not found in UniProtKB"]))
            return xrefs
        return self._apply_entry(xrefs, entry)

    def _apply_entry(self, xrefs: Xrefs, entry: dict[str, Any]) -> Xrefs:
        """Copy the fields we expose out of one UniProtKB entry."""
        xrefs.uniprot = xrefs.uniprot or str(entry.get("primaryAccession") or "") or None
        xrefs.uniprot_entry_name = entry.get("uniProtkbId")
        description = entry.get("proteinDescription") or {}
        recommended = (description.get("recommendedName") or {}).get("fullName") or {}
        xrefs.protein_name = recommended.get("value")
        genes = entry.get("genes") or []
        if genes:
            xrefs.gene_name = (genes[0].get("geneName") or {}).get("value")
        organism = entry.get("organism") or {}
        if organism.get("taxonId"):
            xrefs.taxon_id = int(organism["taxonId"])

        for crossref in entry.get("uniProtKBCrossReferences") or []:
            database = crossref.get("database")
            identifier = str(crossref.get("id", ""))
            if not identifier:
                continue
            if database == "PDB" and identifier not in xrefs.pdb:
                xrefs.pdb.append(identifier)
            elif database == "ChEMBL" and identifier not in xrefs.chembl:
                xrefs.chembl.append(identifier)
            elif database == "RefSeq" and identifier not in xrefs.refseq:
                xrefs.refseq.append(identifier)
        return xrefs

    def map_protein(
        self,
        feature_id: str,
        sequence: str = "",
        *,
        locus_tag: str | None = None,
        taxon_id: int | None = None,
        pdb_hit_uniprot_ids: Iterable[str] = (),
    ) -> Xrefs:
        """Map one protein. Never guesses: an unmapped protein comes back with uniprot=None."""
        taxon = taxon_id if taxon_id is not None else self.taxon_id
        xrefs = Xrefs(feature_id=feature_id, locus_tag=locus_tag, taxon_id=taxon)

        # Fast path: the taxon's whole proteome, fetched once, matched locally.
        index = self.taxon_index(taxon)
        if index is not None:
            entry = index.by_checksum.get(crc64(sequence)) if sequence else None
            route = ROUTE_PROTEOME_SEQUENCE
            if entry is None and locus_tag:
                entry = index.by_locus.get(_normalize_locus(locus_tag))
                route = ROUTE_PROTEOME_LOCUS
            if entry is not None:
                xrefs.route = route
                xrefs.note = f"matched in the taxon {taxon} proteome ({index.entries} entries)"
                return self._apply_entry(xrefs, entry)

        try:
            for route, upi in (
                (ROUTE_SEQUENCE, self._by_sequence(sequence)),
                (ROUTE_PATRIC, self._by_feature_id(feature_id)),
            ):
                if not upi:
                    continue
                entry = self._uniparc_entry(upi)
                accession, note = _select_uniprot_crossref(
                    entry.get("uniParcCrossReferences") or [],
                    feature_id=feature_id,
                    taxon_id=taxon,
                )
                xrefs.uniparc = upi
                if accession:
                    xrefs.uniprot = accession
                    xrefs.route = route
                    xrefs.note = note
                    return self._enrich(xrefs)
                xrefs.note = note

            accession, note = self._by_locus_tag(locus_tag or "", taxon)
            if accession:
                xrefs.uniprot = accession
                xrefs.route = ROUTE_LOCUS
                xrefs.note = note
                return self._enrich(xrefs)
        except HttpError as exc:
            self.failures.append({"feature_id": feature_id, "error": str(exc)})
            xrefs.route = ROUTE_UNMAPPED
            xrefs.note = f"lookup failed: {exc}"
            return xrefs

        related = [str(accession) for accession in pdb_hit_uniprot_ids if accession]
        if related:
            xrefs.related_uniprot = sorted(set(related))
            xrefs.route = ROUTE_PDB_HIT
            xrefs.note = (
                "no mapping for this protein; accessions listed are the PDB homolog's, not ours"
            )
            return xrefs

        xrefs.route = ROUTE_UNMAPPED
        xrefs.note = xrefs.note or "no UniParc, locus-tag or PDB-homolog mapping found"
        return xrefs

    def map_many(
        self, requests_: Sequence[dict[str, Any]], *, workers: int = 4
    ) -> dict[str, Xrefs]:
        """Map many proteins concurrently. Each item is the kwargs of map_protein()."""
        if not requests_:
            return {}
        if workers <= 1:
            results = [self.map_protein(**item) for item in requests_]
        else:
            with ThreadPoolExecutor(max_workers=workers) as pool:
                results = list(pool.map(lambda item: self.map_protein(**item), requests_))
        return {xrefs.feature_id: xrefs for xrefs in results}

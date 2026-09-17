"""AlphaFold DB lookups by UniProt accession (issue #8).

M3 must not re-fold what AlphaFold DB already has (`03-m3-fold.md`), so M2 records what exists
and how good it is. Two numbers decide whether a model is usable, and both are recorded rather
than assumed:

- **Confidence.** ``globalMetricValue`` is the mean pLDDT. AlphaFold's own bands: >90 very high,
  70–90 confident, 50–70 low, <50 very low. A pocket in a low-confidence region is not a pocket
  (pitfall in `03-m3-fold.md`), so the band travels with the model.
- **Coverage.** A prediction may cover only a fragment of the accession. SARS-CoV-2 orf1ab
  (7,096 aa) returns a model for residues 1368–1493 — 2% of the protein. Recording only "a model
  exists" would be misleading, so the residue range and the covered fraction are recorded too.

A predicted model also has no cofactors, metals or ligands (pitfall #2), which is why
``holo_homolog`` from the PDB side stays the stronger signal for docking.

API behaviour, verified: 200 with a list = found; 404 = no model (e.g. titin); 400 = malformed
accession; the list may hold several fragments.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Iterable

from ..common.http import CachedJsonClient, HttpError

AFDB_PREDICTION = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"

STATUS_FOUND = "found"
STATUS_NO_MODEL = "no-model"
STATUS_INVALID = "invalid-accession"
STATUS_FAILED = "query-failed"
STATUS_NOT_QUERIED = "not-queried"

# AlphaFold's published pLDDT bands.
BAND_VERY_HIGH = 90.0
BAND_CONFIDENT = 70.0
BAND_LOW = 50.0


def confidence_band(mean_plddt: float | None) -> str:
    if mean_plddt is None:
        return ""
    if mean_plddt >= BAND_VERY_HIGH:
        return "very high"
    if mean_plddt >= BAND_CONFIDENT:
        return "confident"
    if mean_plddt >= BAND_LOW:
        return "low"
    return "very low"


@dataclass
class AlphaFoldModel:
    """One AlphaFold DB prediction, or an explicit reason there is none."""

    accession: str = ""
    status: str = STATUS_NOT_QUERIED
    entry_id: str = ""
    version: str = ""
    mean_plddt: float | None = None
    uniprot_start: int | None = None
    uniprot_end: int | None = None
    sequence_length: int | None = None
    fragments: int = 0
    cif_url: str = ""
    pdb_url: str = ""
    pae_url: str = ""
    created: str = ""
    error: str = ""

    @property
    def found(self) -> bool:
        return self.status == STATUS_FOUND

    @property
    def confidence_band(self) -> str:
        return confidence_band(self.mean_plddt)

    @property
    def covered_residues(self) -> int | None:
        if self.uniprot_start is None or self.uniprot_end is None:
            return None
        return max(0, self.uniprot_end - self.uniprot_start + 1)

    @property
    def coverage(self) -> float | None:
        """Fraction of the accession's sequence the model covers."""
        covered = self.covered_residues
        if covered is None or not self.sequence_length:
            return None
        return round(min(1.0, covered / self.sequence_length), 4)

    def usable_for_docking(self, *, min_plddt: float = BAND_CONFIDENT, min_coverage: float = 0.8) -> tuple[bool, str]:
        """M3's gate, with the reason it needs to record either way."""
        if not self.found:
            return False, f"no AlphaFold model ({self.status})"
        if self.mean_plddt is not None and self.mean_plddt < min_plddt:
            return False, f"mean pLDDT {self.mean_plddt} is {self.confidence_band}, below {min_plddt}"
        coverage = self.coverage
        if coverage is not None and coverage < min_coverage:
            return (
                False,
                f"model covers residues {self.uniprot_start}-{self.uniprot_end} "
                f"({coverage:.0%} of the accession), below {min_coverage:.0%}",
            )
        return True, f"mean pLDDT {self.mean_plddt} ({self.confidence_band})"

    def as_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["confidence_band"] = self.confidence_band
        payload["coverage"] = self.coverage
        return payload


def parse_prediction(accession: str, payload: Any) -> AlphaFoldModel:
    """Turn one API response into a model record. Picks the highest-confidence fragment."""
    if isinstance(payload, dict) and payload.get("_http_status"):
        status = int(payload["_http_status"])
        if status == 404:
            return AlphaFoldModel(accession=accession, status=STATUS_NO_MODEL)
        if status == 400:
            return AlphaFoldModel(
                accession=accession, status=STATUS_INVALID, error="AlphaFold rejected the identifier format"
            )
        return AlphaFoldModel(accession=accession, status=STATUS_FAILED, error=f"HTTP {status}")

    entries = payload if isinstance(payload, list) else []
    if not entries:
        return AlphaFoldModel(accession=accession, status=STATUS_NO_MODEL)

    def confidence(entry: dict[str, Any]) -> float:
        try:
            return float(entry.get("globalMetricValue") or 0.0)
        except (TypeError, ValueError):
            return 0.0

    best = max(entries, key=confidence)
    sequence = best.get("uniprotSequence") or ""
    return AlphaFoldModel(
        accession=accession,
        status=STATUS_FOUND,
        entry_id=str(best.get("entryId") or ""),
        version=str(best.get("latestVersion") or ""),
        mean_plddt=confidence(best) or None,
        uniprot_start=_int_or_none(best.get("uniprotStart")),
        uniprot_end=_int_or_none(best.get("uniprotEnd")),
        sequence_length=len(sequence) if sequence else None,
        fragments=len(entries),
        cif_url=str(best.get("cifUrl") or ""),
        pdb_url=str(best.get("pdbUrl") or ""),
        pae_url=str(best.get("paeDocUrl") or ""),
        created=str(best.get("modelCreatedDate") or ""),
    )


def _int_or_none(value: Any) -> int | None:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


class AlphaFoldClient:
    """Look up AlphaFold DB predictions, cached like every other source."""

    def __init__(self, client: CachedJsonClient) -> None:
        self.client = client
        self.failures: list[dict[str, str]] = []

    def lookup(self, accession: str) -> AlphaFoldModel:
        if not accession:
            return AlphaFoldModel(status=STATUS_NOT_QUERIED, error="no UniProt accession to query")
        try:
            payload = self.client.get_json(
                "afdb_prediction",
                AFDB_PREDICTION.format(accession=accession),
                allow_statuses=(400, 404),
            )
        except HttpError as exc:
            self.failures.append({"accession": accession, "error": str(exc)})
            return AlphaFoldModel(accession=accession, status=STATUS_FAILED, error=str(exc))
        return parse_prediction(accession, payload)

    def lookup_many(self, accessions: Iterable[str], *, workers: int = 4) -> dict[str, AlphaFoldModel]:
        unique = sorted({a for a in accessions if a})
        if not unique:
            return {}
        if workers <= 1:
            return {accession: self.lookup(accession) for accession in unique}
        from concurrent.futures import ThreadPoolExecutor

        with ThreadPoolExecutor(max_workers=workers) as pool:
            models = list(pool.map(self.lookup, unique))
        return dict(zip(unique, models))

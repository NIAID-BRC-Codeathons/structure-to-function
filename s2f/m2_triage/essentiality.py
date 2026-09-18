"""Essentiality by orthology to public relatives (issue #29).

BV-BRC precomputes `Essential Gene` rows from flux-balance analysis, but **only for public
genomes**: the public G37 record has 148, a fresh CGA run has none. A blinded genome therefore
gets no essentiality at all, and M2's score has a component for it.

So we transfer it. Build a reference set from FBA-essential proteins in public relatives, align
our proteome against it with DIAMOND, and call a protein essential when it has a strong ortholog.
Every call names **which relative and what identity** justified it, because that is the only way a
reviewer can judge it.

What this is and is not:

- FBA essentiality is **computational**, not experimental — a metabolic model's prediction that
  removing the gene stops growth in silico. Transferring it by homology adds a second inference.
  Both are recorded, and `evidence` says `FBA` so nothing downstream can mistake it for a knockout
  experiment.
- The transfer threshold (40% identity, 70% coverage) is the conventional bar for functional
  transfer between bacteria. It is a judgement, recorded in the run manifest and adjustable.
"""

from __future__ import annotations

import csv
import subprocess
import tempfile
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Iterable, Sequence

from .human_homology import OUTFMT, diamond_version, write_fasta

BVBRC_API = "https://www.bv-brc.org/api"

#: Fetch the whole relative reference set. Capping it costs recall: 5,000 of the 12,742
#: Mycoplasma FBA proteins recovered 113 of G37's 148 known essential genes; all 12,742
#: recovered 148 of 148.
DEFAULT_REFERENCE_LIMIT = 15000

DEFAULT_MIN_IDENTITY = 40.0
DEFAULT_MIN_COVERAGE = 0.7
DEFAULT_EVALUE = 1e-10
DEFAULT_SENSITIVITY = "--very-sensitive"
BATCH = 100

STATUS_ESSENTIAL = "essential-ortholog"
STATUS_NO_CALL = "no-call"
STATUS_NOT_RUN = "not-run"


@dataclass
class EssentialityCall:
    """One protein's essentiality evidence, or an explicit absence of it."""

    feature_id: str
    status: str = STATUS_NOT_RUN
    reference_patric_id: str = ""
    reference_genome: str = ""
    reference_gene: str = ""
    reference_product: str = ""
    identity: float | None = None
    coverage: float | None = None
    evalue: float | None = None
    evidence: str = "FBA"
    note: str = ""

    @property
    def essential(self) -> bool:
        return self.status == STATUS_ESSENTIAL

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)

    def as_annotation(self, *, retrieved_at: str) -> dict[str, Any]:
        return {
            "source": "essentiality_ortholog",
            "hit": self.reference_patric_id or None,
            "description": (
                f"ortholog of {self.reference_gene or self.reference_product or 'an essential gene'}"
                f" in {self.reference_genome}"
                if self.essential
                else "no essential-gene ortholog"
            ),
            "identity": round(self.identity / 100.0, 4) if self.identity is not None else None,
            "coverage": round(self.coverage, 4) if self.coverage is not None else None,
            "evalue": self.evalue,
            "retrieved_at": retrieved_at,
            "evidence": self.evidence,
            "note": self.note or (
                "FBA essentiality in a relative, transferred by homology: two inferences, not an "
                "experimental knockout"
                if self.essential
                else None
            ),
        }


@dataclass
class ReferenceProtein:
    patric_id: str
    genome_name: str
    gene: str
    product: str
    sequence: str


@dataclass
class EssentialitySearch:
    calls: dict[str, EssentialityCall] = field(default_factory=dict)
    reference_proteins: int = 0
    reference_genomes: int = 0
    reference_query: str = ""
    min_identity: float = DEFAULT_MIN_IDENTITY
    min_coverage: float = DEFAULT_MIN_COVERAGE
    elapsed_seconds: float = 0.0
    diamond_version: str = ""

    def summary(self) -> dict[str, Any]:
        return {
            "reference": "BV-BRC sp_gene, property 'Essential Gene', evidence FBA",
            "reference_query": self.reference_query,
            "reference_proteins": self.reference_proteins,
            "reference_genomes": self.reference_genomes,
            "min_identity": self.min_identity,
            "min_coverage": self.min_coverage,
            "diamond_version": self.diamond_version,
            "essential_calls": sum(1 for c in self.calls.values() if c.essential),
            "no_call": sum(1 for c in self.calls.values() if c.status == STATUS_NO_CALL),
            # Distinguishes "we searched and found no ortholog" from "the search never ran".
            # Without this an empty reference set reports all-zero and reads as a completed
            # transfer that found nothing (observed on lambda0, 2026-09-18, with an empty
            # --essentiality-keyword sending `keyword()` to BV-BRC).
            "not_run": sum(1 for c in self.calls.values() if c.status == STATUS_NOT_RUN),
            "elapsed_seconds": round(self.elapsed_seconds, 2),
            "caveat": (
                "FBA essentiality is a metabolic-model prediction, transferred here by homology. "
                "Not an experimental knockout."
            ),
        }


def _api_get(client, core: str, rql: str) -> list[dict[str, Any]]:
    """BV-BRC uses RQL in the query string. Values containing spaces need %22 quotes."""
    url = f"{BVBRC_API}/{core}/?{rql}&http_accept=application/json"
    payload = client.get_json(f"bvbrc_{core}", url)
    return payload if isinstance(payload, list) else []


def fetch_reference_set(
    client,
    *,
    keyword: str,
    limit: int = DEFAULT_REFERENCE_LIMIT,
) -> list[ReferenceProtein]:
    """FBA-essential proteins from public relatives, with their sequences.

    ``keyword`` selects the relatives — a genus or species name, matched by BV-BRC's keyword
    search. Three calls: the essential-gene rows, the features' sequence checksums, then the
    sequences themselves, batched.
    """
    rows = _api_get(
        client,
        "sp_gene",
        f"and(eq(evidence,FBA),keyword({keyword}))"
        f"&select(patric_id,genome_id,genome_name,gene,product)&limit({limit})",
    )
    by_id = {r["patric_id"]: r for r in rows if r.get("patric_id")}
    if not by_id:
        return []

    checksums: dict[str, str] = {}
    identifiers = sorted(by_id)
    for start in range(0, len(identifiers), BATCH):
        batch = identifiers[start : start + BATCH]
        features = _api_get(
            client,
            "genome_feature",
            f"in(patric_id,({','.join(batch)}))&select(patric_id,aa_sequence_md5)&limit({BATCH})",
        )
        for feature in features:
            if feature.get("aa_sequence_md5"):
                checksums[feature["patric_id"]] = feature["aa_sequence_md5"]

    sequences: dict[str, str] = {}
    unique_md5 = sorted(set(checksums.values()))
    for start in range(0, len(unique_md5), BATCH):
        batch = unique_md5[start : start + BATCH]
        for record in _api_get(
            client, "feature_sequence", f"in(md5,({','.join(batch)}))&select(md5,sequence)&limit({BATCH})"
        ):
            if record.get("md5") and record.get("sequence"):
                sequences[record["md5"]] = record["sequence"]

    reference: list[ReferenceProtein] = []
    for patric_id, row in by_id.items():
        md5 = checksums.get(patric_id)
        sequence = sequences.get(md5 or "")
        if not sequence:
            continue
        reference.append(
            ReferenceProtein(
                patric_id=patric_id,
                genome_name=row.get("genome_name") or "",
                gene=row.get("gene") or "",
                product=row.get("product") or "",
                sequence=sequence,
            )
        )
    return reference


def resolve_reference_set(
    client,
    *,
    keyword: str = "",
    taxonomy: Any = None,
    limit: int = DEFAULT_REFERENCE_LIMIT,
) -> tuple[str, list[ReferenceProtein]]:
    """Find a keyword that actually returns relatives, and the set it returns.

    An explicit ``keyword`` is used as given. Otherwise walk the query genome's lineage outwards
    from the nearest rank, because the caller does not know what the organism is — that is the
    premise of the pipeline — and because BV-BRC's keyword index does not track NCBI's current
    names. Measured 2026-09-18: this genome resolves to genus *Mycoplasmoides*, which returns 0
    rows, while its former genus *Mycoplasma* returns 12,742 rows across 96 genomes.
    """
    if keyword:
        return keyword, fetch_reference_set(client, keyword=_rql_value(keyword), limit=limit)
    candidates = lineage_candidates(taxonomy)
    for candidate in candidates:
        try:
            reference = fetch_reference_set(client, keyword=_rql_value(candidate), limit=limit)
        except (RuntimeError, OSError, ValueError):
            # One rank failing must not abort the walk: the whole point is that we do not know
            # which name this index carries, so try the next rank out.
            continue
        if reference:
            return candidate, reference
    return (candidates[0] if candidates else ""), []


def _rql_value(value: str) -> str:
    """RQL needs %22 quotes around any value containing a space (see ``_api_get``)."""
    return f"%22{value}%22" if " " in value else value


def lineage_candidates(taxonomy: Any) -> list[str]:
    """Names to try as a keyword, nearest rank first.

    M1 writes each lineage entry as ``[name, taxon_id, rank]``; a bare string is accepted too so
    a hand-built taxonomy still works.
    """
    names: list[str] = []
    for entry in getattr(taxonomy, "lineage_names", None) or []:
        if isinstance(entry, (list, tuple)):
            entry = entry[0] if entry else ""
        if isinstance(entry, str) and entry.strip():
            names.append(entry.strip())
    candidates = list(reversed(names))[:4]
    scientific = (getattr(taxonomy, "scientific_name", "") or "").strip()
    if scientific and scientific not in candidates:
        candidates.insert(0, scientific)
    return candidates


def parse_calls(
    rows: Iterable[Sequence[str]],
    reference: dict[str, ReferenceProtein],
    *,
    min_identity: float = DEFAULT_MIN_IDENTITY,
    min_coverage: float = DEFAULT_MIN_COVERAGE,
) -> dict[str, EssentialityCall]:
    """Best ortholog per protein, gated on identity and coverage."""
    best: dict[str, tuple[float, EssentialityCall]] = {}
    for row in rows:
        if len(row) < len(OUTFMT):
            continue
        record = dict(zip(OUTFMT, row))
        try:
            bitscore = float(record["bitscore"])
            identity = float(record["pident"])
            qlen = int(record["qlen"])
            qstart, qend = int(record["qstart"]), int(record["qend"])
        except (TypeError, ValueError):
            continue
        if qlen <= 0:
            continue
        coverage = max(0.0, min(1.0, (qend - qstart + 1) / qlen))
        feature_id = record["qseqid"]
        if feature_id in best and bitscore <= best[feature_id][0]:
            continue

        source = reference.get(record["sseqid"])
        passes = identity >= min_identity and coverage >= min_coverage
        call = EssentialityCall(
            feature_id=feature_id,
            status=STATUS_ESSENTIAL if passes else STATUS_NO_CALL,
            reference_patric_id=record["sseqid"] if passes else "",
            reference_genome=(source.genome_name if source and passes else ""),
            reference_gene=(source.gene if source and passes else ""),
            reference_product=(source.product if source and passes else ""),
            identity=identity,
            coverage=coverage,
            evalue=float(record["evalue"]),
            note=(
                ""
                if passes
                else f"best ortholog {identity:.0f}% identity over {coverage:.0%} of the query, "
                f"below the {min_identity:.0f}%/{min_coverage:.0%} transfer threshold"
            ),
        )
        best[feature_id] = (bitscore, call)
    return {feature_id: call for feature_id, (_score, call) in best.items()}


def search(
    proteins: Iterable[tuple[str, str]],
    reference: Sequence[ReferenceProtein],
    *,
    reference_query: str = "",
    min_identity: float = DEFAULT_MIN_IDENTITY,
    min_coverage: float = DEFAULT_MIN_COVERAGE,
    evalue: float = DEFAULT_EVALUE,
    sensitivity: str = DEFAULT_SENSITIVITY,
    threads: int = 4,
    workdir: Path | None = None,
) -> EssentialitySearch:
    """Align our proteome against the reference set and transfer essentiality."""
    proteins = list(proteins)
    result = EssentialitySearch(
        min_identity=min_identity, min_coverage=min_coverage, reference_query=reference_query
    )
    result.diamond_version = diamond_version()
    result.reference_proteins = len(reference)
    result.reference_genomes = len({r.genome_name for r in reference})
    if not reference:
        result.calls = {
            feature_id: EssentialityCall(
                feature_id=feature_id, status=STATUS_NOT_RUN, note="no reference set available"
            )
            for feature_id, _ in proteins
        }
        return result

    by_id = {r.patric_id: r for r in reference}
    started = time.monotonic()
    with tempfile.TemporaryDirectory(dir=workdir) as tmp:
        tmpdir = Path(tmp)
        query_path = tmpdir / "query.faa"
        reference_path = tmpdir / "reference.faa"
        database = tmpdir / "essential"
        hits_path = tmpdir / "hits.tsv"

        write_fasta(query_path, proteins)
        write_fasta(reference_path, ((r.patric_id, r.sequence) for r in reference))
        subprocess.run(
            ["diamond", "makedb", "--in", str(reference_path), "-d", str(database), "--quiet"],
            check=True, capture_output=True, text=True, timeout=1800,
        )
        process = subprocess.run(
            [
                "diamond", "blastp", "-q", str(query_path), "-d", str(database),
                "-o", str(hits_path), "--quiet", "--outfmt", "6", *OUTFMT,
                "--evalue", str(evalue), "--max-target-seqs", "5",
                "--threads", str(threads), sensitivity,
            ],
            capture_output=True, text=True, timeout=7200,
        )
        if process.returncode != 0:
            raise RuntimeError(f"diamond blastp failed: {(process.stderr or process.stdout)[:500]}")
        with hits_path.open(encoding="utf-8") as handle:
            result.calls = parse_calls(
                csv.reader(handle, delimiter="\t"), by_id,
                min_identity=min_identity, min_coverage=min_coverage,
            )

    for feature_id, _ in proteins:
        if feature_id not in result.calls:
            result.calls[feature_id] = EssentialityCall(
                feature_id=feature_id,
                status=STATUS_NO_CALL,
                note=f"no ortholog in the reference set at e-value <= {evalue:g}",
            )
    result.elapsed_seconds = time.monotonic() - started
    return result

"""Human-homology search for every protein (issue #29).

BV-BRC precomputes `Human Homolog` rows, but **only for public genomes**: a fresh CGA run on the
G37 test genome returns 17 specialty rows where the public record has 174. So on a blinded
genome the human-homolog penalty in M2's score silently never fires — and that penalty is the
only thing standing between us and recommending a target whose close human counterpart makes it
a toxicity risk (pitfall #3/#4). This module computes the evidence instead of borrowing it.

Method: DIAMOND blastp against the reviewed human proteome (UniProt UP000005640, Swiss-Prot,
~20,400 sequences), **`--very-sensitive`**.

Two settings decided from measurement, not preference:

- **Sensitivity.** DIAMOND's default mode misses 40% of the proteins that have human homologs
  (107 vs 179 on G37) and misses `rpsL` outright — a genuine 50.6%-identity, 83-residue, gapless
  alignment, confirmed independently with a pairwise aligner and matching BV-BRC's own 50%.
  `--very-sensitive` recovers all five of the genome's known human homologs at identities within
  ~1% of BV-BRC's, and costs 1.4 s for 542 proteins.
- **Coverage gate.** A hit counts as a homolog only when it covers enough of the query
  (default 0.5). BV-BRC's own rows record an identity with no coverage at all, which would let a
  50% match over 20 residues penalise a target as heavily as a full-length one.

Every protein gets a record: a hit with its numbers, or an explicit "no hit". Never a silent zero.
"""

from __future__ import annotations

import csv
import shutil
import subprocess
import tempfile
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Iterable

HUMAN_PROTEOME_QUERY = "proteome:UP000005640 AND reviewed:true"
HUMAN_PROTEOME_URL = "https://rest.uniprot.org/uniprotkb/stream"
HUMAN_PROTEOME_ID = "UP000005640"

DEFAULT_SENSITIVITY = "--very-sensitive"
DEFAULT_EVALUE = 1e-5
DEFAULT_MIN_COVERAGE = 0.5
DEFAULT_MAX_TARGETS = 5

STATUS_HIT = "hit"
STATUS_NO_HIT = "no-hit"
STATUS_NOT_RUN = "not-run"

OUTFMT = [
    "qseqid", "sseqid", "pident", "length", "qlen", "slen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
]


class DiamondMissing(RuntimeError):
    """DIAMOND is not installed. See the README's external-tools table."""


@dataclass
class HumanHit:
    """One protein's best human homolog, or an explicit absence."""

    feature_id: str
    status: str = STATUS_NOT_RUN
    accession: str = ""
    entry_name: str = ""
    identity: float | None = None          # percent, as BV-BRC reports it
    coverage: float | None = None          # fraction of the query aligned
    evalue: float | None = None
    bitscore: float | None = None
    alignment_length: int | None = None
    query_start: int | None = None
    query_end: int | None = None
    counted: bool = False                  # passed the coverage gate
    note: str = ""

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)

    def as_annotation(self, *, retrieved_at: str, source_version: str = "") -> dict[str, Any]:
        """The `annotations[]` row for report.json (schema requires `source`)."""
        return {
            "source": "human_homology_diamond",
            "hit": self.accession or None,
            "description": self.entry_name or ("no human homolog" if self.status == STATUS_NO_HIT else None),
            "identity": round(self.identity / 100.0, 4) if self.identity is not None else None,
            "coverage": round(self.coverage, 4) if self.coverage is not None else None,
            "evalue": self.evalue,
            "retrieved_at": retrieved_at,
            "reference": f"UniProt {HUMAN_PROTEOME_ID} (reviewed)",
            "reference_version": source_version or None,
            "counted_as_homolog": self.counted,
            "note": self.note or None,
        }


def _parse_subject(sseqid: str) -> tuple[str, str]:
    """`sp|O15235|RT12_HUMAN` -> ('O15235', 'RT12_HUMAN')."""
    parts = sseqid.split("|")
    if len(parts) >= 3:
        return parts[1], parts[2]
    return sseqid, ""


def parse_diamond_hits(
    rows: Iterable[list[str]],
    *,
    min_coverage: float = DEFAULT_MIN_COVERAGE,
) -> dict[str, HumanHit]:
    """Best hit per query by bitscore, with the coverage gate applied."""
    best: dict[str, HumanHit] = {}
    for row in rows:
        if len(row) < len(OUTFMT):
            continue
        record = dict(zip(OUTFMT, row))
        try:
            bitscore = float(record["bitscore"])
            qlen = int(record["qlen"])
            qstart, qend = int(record["qstart"]), int(record["qend"])
        except (TypeError, ValueError):
            continue
        if qlen <= 0:
            continue
        coverage = max(0.0, min(1.0, (qend - qstart + 1) / qlen))
        feature_id = record["qseqid"]
        existing = best.get(feature_id)
        if existing is not None and existing.bitscore is not None and bitscore <= existing.bitscore:
            continue
        accession, entry_name = _parse_subject(record["sseqid"])
        counted = coverage >= min_coverage
        best[feature_id] = HumanHit(
            feature_id=feature_id,
            status=STATUS_HIT,
            accession=accession,
            entry_name=entry_name,
            identity=float(record["pident"]),
            coverage=coverage,
            evalue=float(record["evalue"]),
            bitscore=bitscore,
            alignment_length=int(record["length"]),
            query_start=qstart,
            query_end=qend,
            counted=counted,
            note=(
                ""
                if counted
                else f"alignment covers {coverage:.0%} of the query, below the {min_coverage:.0%} gate"
            ),
        )
    return best


@dataclass
class HumanHomologySearch:
    """Result of one search over a proteome."""

    hits: dict[str, HumanHit] = field(default_factory=dict)
    reference_sequences: int = 0
    sensitivity: str = DEFAULT_SENSITIVITY
    evalue: float = DEFAULT_EVALUE
    min_coverage: float = DEFAULT_MIN_COVERAGE
    elapsed_seconds: float = 0.0
    diamond_version: str = ""
    error: str = ""

    def summary(self) -> dict[str, Any]:
        counted = [h for h in self.hits.values() if h.counted]
        below_gate = [h for h in self.hits.values() if h.status == STATUS_HIT and not h.counted]
        return {
            "reference": f"UniProt {HUMAN_PROTEOME_ID} (reviewed)",
            "reference_sequences": self.reference_sequences,
            "diamond_version": self.diamond_version,
            "sensitivity": self.sensitivity,
            "evalue": self.evalue,
            "min_coverage": self.min_coverage,
            "proteins_with_hit": sum(1 for h in self.hits.values() if h.status == STATUS_HIT),
            "counted_as_homolog": len(counted),
            "below_coverage_gate": len(below_gate),
            "elapsed_seconds": round(self.elapsed_seconds, 2),
        }


def diamond_version() -> str:
    """Version string, or raise if DIAMOND is not installed."""
    binary = shutil.which("diamond")
    if binary is None:
        raise DiamondMissing(
            "diamond not found. Install it (`brew install diamond` or "
            "`conda install -c bioconda diamond`) — see the README's external-tools table."
        )
    result = subprocess.run([binary, "--version"], capture_output=True, text=True, timeout=60)
    return result.stdout.strip() or result.stderr.strip()


def write_fasta(path: Path, records: Iterable[tuple[str, str]]) -> int:
    """Write (id, sequence) pairs. Returns how many were written."""
    count = 0
    with path.open("w", encoding="utf-8") as handle:
        for identifier, sequence in records:
            if not identifier or not sequence:
                continue
            handle.write(f">{identifier}\n")
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start : start + 60] + "\n")
            count += 1
    return count


def search(
    proteins: Iterable[tuple[str, str]],
    reference_fasta: Path,
    *,
    sensitivity: str = DEFAULT_SENSITIVITY,
    evalue: float = DEFAULT_EVALUE,
    min_coverage: float = DEFAULT_MIN_COVERAGE,
    max_targets: int = DEFAULT_MAX_TARGETS,
    threads: int = 4,
    workdir: Path | None = None,
) -> HumanHomologySearch:
    """Run DIAMOND blastp of `proteins` against the human proteome.

    Every protein that produced no hit gets an explicit ``no-hit`` record, so "we looked and found
    nothing" is distinguishable from "we never looked".
    """
    proteins = list(proteins)
    result = HumanHomologySearch(
        sensitivity=sensitivity, evalue=evalue, min_coverage=min_coverage
    )
    result.diamond_version = diamond_version()

    reference_fasta = Path(reference_fasta)
    if not reference_fasta.exists():
        raise FileNotFoundError(f"human proteome FASTA not found: {reference_fasta}")

    started = time.monotonic()
    with tempfile.TemporaryDirectory(dir=workdir) as tmp:
        tmpdir = Path(tmp)
        query_path = tmpdir / "query.faa"
        write_fasta(query_path, proteins)
        database = tmpdir / "human_sp"
        hits_path = tmpdir / "hits.tsv"

        subprocess.run(
            ["diamond", "makedb", "--in", str(reference_fasta), "-d", str(database), "--quiet"],
            check=True, capture_output=True, text=True, timeout=1800,
        )
        # Output goes to a real file: DIAMOND writes temp files beside it, so /dev/stdout fails.
        command = [
            "diamond", "blastp",
            "-q", str(query_path), "-d", str(database), "-o", str(hits_path),
            "--quiet", "--outfmt", "6", *OUTFMT,
            "--evalue", str(evalue), "--max-target-seqs", str(max_targets),
            "--threads", str(threads), sensitivity,
        ]
        process = subprocess.run(command, capture_output=True, text=True, timeout=7200)
        if process.returncode != 0:
            result.error = (process.stderr or process.stdout).strip()[:500]
            raise RuntimeError(f"diamond blastp failed: {result.error}")

        with hits_path.open(encoding="utf-8") as handle:
            result.hits = parse_diamond_hits(
                csv.reader(handle, delimiter="\t"), min_coverage=min_coverage
            )

    for feature_id, _ in proteins:
        if feature_id not in result.hits:
            result.hits[feature_id] = HumanHit(
                feature_id=feature_id,
                status=STATUS_NO_HIT,
                note=f"no human hit at e-value <= {evalue:g} ({sensitivity})",
            )

    result.elapsed_seconds = time.monotonic() - started
    with reference_fasta.open(encoding="utf-8") as handle:
        result.reference_sequences = sum(1 for line in handle if line.startswith(">"))
    return result


def fetch_human_proteome(client, destination: Path) -> Path:
    """Download the reviewed human proteome once and keep it on disk.

    Not routed through the JSON cache: this is a 14 MB FASTA, not an API record.
    """
    destination = Path(destination)
    if destination.exists() and destination.stat().st_size > 0:
        return destination
    destination.parent.mkdir(parents=True, exist_ok=True)
    response = client.session.get(
        HUMAN_PROTEOME_URL,
        params={"query": HUMAN_PROTEOME_QUERY, "format": "fasta"},
        timeout=900,
        stream=True,
    )
    response.raise_for_status()
    with destination.open("wb") as handle:
        for chunk in response.iter_content(chunk_size=1 << 20):
            handle.write(chunk)
    return destination

"""Load M1-shaped BV-BRC output: protein FASTA plus feature and specialty tables.

Canonical protein key is the BV-BRC ``patric_id`` (``feature_id`` in
``docs/00-architecture.md``). Column names vary between the BV-BRC web UI and the ``p3-`` CLI,
so headers are matched through aliases rather than assumed.
"""

from __future__ import annotations

import csv
import hashlib
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

ID_ALIASES = ("patric_id", "feature_id", "brc_id", "genome_feature_id")
PRODUCT_ALIASES = ("product", "function")
GENE_ALIASES = ("gene", "gene_name", "gene_symbol")
LOCUS_ALIASES = ("refseq_locus_tag", "locus_tag", "alt_locus_tag")
PGFAM_ALIASES = ("pgfam_id", "pgfam", "pgfams")
PLFAM_ALIASES = ("plfam_id", "plfam", "plfams")
PROPERTY_ALIASES = ("property", "specialty_gene_property", "type")
SOURCE_ALIASES = ("source", "database", "specialty_source")
IDENTITY_ALIASES = ("identity", "pct_identity", "percent_identity")
COVERAGE_ALIASES = ("query_coverage", "coverage", "query_cov")

VIRULENCE_PROPERTIES = {"virulence factor", "virulance factor"}
AMR_PROPERTIES = {"antibiotic resistance"}
ESSENTIAL_PROPERTIES = {"essential gene"}
DRUG_TARGET_PROPERTIES = {"drug target"}
HUMAN_HOMOLOG_PROPERTIES = {"human homolog"}
TRANSPORTER_PROPERTIES = {"transporter"}
METAL_PROPERTIES = {"metal resistance"}

UNCHARACTERIZED = re.compile(
    r"hypothetical|uncharacteri[sz]ed|putative|unknown function|unnamed|"
    r"\bDUF\d+|protein of unknown|predicted protein",
    re.IGNORECASE,
)


@dataclass
class SpecialtyHit:
    """One row of the specialty-gene table."""

    property_name: str
    source: str
    product: str
    identity: float | None
    query_coverage: float | None


@dataclass
class Protein:
    """One protein with the M1 metadata this module scores against."""

    feature_id: str
    sequence: str
    product: str = ""
    gene: str = ""
    locus_tag: str = ""
    pgfam: str = ""
    plfam: str = ""
    specialty: list[SpecialtyHit] = field(default_factory=list)

    @property
    def sequence_hash(self) -> str:
        return hashlib.sha256(self.sequence.encode("ascii")).hexdigest()

    def _properties(self) -> set[str]:
        return {hit.property_name.strip().lower() for hit in self.specialty}

    @property
    def has_virulence_or_amr(self) -> bool:
        return bool(self._properties() & (VIRULENCE_PROPERTIES | AMR_PROPERTIES))

    @property
    def is_essential_ortholog(self) -> bool:
        return bool(self._properties() & ESSENTIAL_PROPERTIES)

    @property
    def is_drug_target(self) -> bool:
        return bool(self._properties() & DRUG_TARGET_PROPERTIES)

    @property
    def is_transporter(self) -> bool:
        return bool(self._properties() & TRANSPORTER_PROPERTIES)

    @property
    def has_metal_resistance(self) -> bool:
        return bool(self._properties() & METAL_PROPERTIES)

    @property
    def human_homolog_identity(self) -> float | None:
        """Highest reported identity to a human protein, as a percentage."""
        identities = [
            hit.identity
            for hit in self.specialty
            if hit.property_name.strip().lower() in HUMAN_HOMOLOG_PROPERTIES
            and hit.identity is not None
        ]
        return max(identities) if identities else None

    @property
    def is_uncharacterized(self) -> bool:
        return bool(UNCHARACTERIZED.search(self.product or ""))


@dataclass
class InputBundle:
    """Loaded proteins plus what could not be joined, for the run manifest."""

    proteins: list[Protein]
    specialty_without_protein: list[str] = field(default_factory=list)
    table_rows_without_sequence: list[str] = field(default_factory=list)
    sequences_without_table_row: list[str] = field(default_factory=list)

    def unique_sequences(self) -> dict[str, list[Protein]]:
        """Group proteins by sequence so each distinct sequence is searched once."""
        groups: dict[str, list[Protein]] = defaultdict(list)
        for protein in self.proteins:
            groups[protein.sequence_hash].append(protein)
        return dict(groups)


def _pick(row: dict[str, str], aliases: tuple[str, ...]) -> str:
    for alias in aliases:
        if alias in row and row[alias] is not None:
            value = str(row[alias]).strip()
            if value:
                return value
    return ""


def _normalize_header(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", (name or "").strip().lower()).strip("_")


def _read_table(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        sample = handle.read(8192)
        handle.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
        except csv.Error:
            dialect = csv.excel
        reader = csv.DictReader(handle, dialect=dialect)
        rows = []
        for raw in reader:
            rows.append({_normalize_header(k): v for k, v in raw.items() if k is not None})
        return rows


def _read_fasta(path: Path) -> list[tuple[str, str]]:
    """Return (header, sequence) pairs. Kept dependency-free; headers stay raw."""
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header = line[1:]
                chunks = []
            else:
                chunks.append(line.replace(" ", "").upper())
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def _feature_id_from_header(header: str) -> tuple[str, str]:
    """Split ``fig|1125630.4.peg.5196  DNA replication protein`` into id and product."""
    parts = header.split(None, 1)
    return parts[0].strip(), (parts[1].strip() if len(parts) > 1 else "")


def _float_or_none(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def load_input(m1_dir: Path) -> InputBundle:
    """Load ``proteins.faa`` plus the feature and (optional) specialty tables."""
    m1_dir = Path(m1_dir)
    fasta_path = m1_dir / "proteins.faa"
    if not fasta_path.exists():
        raise FileNotFoundError(f"missing protein FASTA: {fasta_path}")

    table_rows: list[dict[str, str]] = []
    for name in ("genes_proteins.csv", "features.csv", "proteins.tsv", "features.tsv"):
        candidate = m1_dir / name
        if candidate.exists():
            table_rows = _read_table(candidate)
            break

    by_id: dict[str, dict[str, str]] = {}
    for row in table_rows:
        feature_id = _pick(row, ID_ALIASES)
        if feature_id:
            by_id[feature_id] = row

    proteins: list[Protein] = []
    sequences_without_row: list[str] = []
    seen_ids: set[str] = set()
    for header, sequence in _read_fasta(fasta_path):
        feature_id, header_product = _feature_id_from_header(header)
        if not feature_id or not sequence:
            continue
        if feature_id in seen_ids:
            continue
        seen_ids.add(feature_id)
        row = by_id.get(feature_id)
        if row is None:
            sequences_without_row.append(feature_id)
            row = {}
        proteins.append(
            Protein(
                feature_id=feature_id,
                sequence=sequence,
                product=_pick(row, PRODUCT_ALIASES) or header_product,
                gene=_pick(row, GENE_ALIASES),
                locus_tag=_pick(row, LOCUS_ALIASES),
                pgfam=_pick(row, PGFAM_ALIASES),
                plfam=_pick(row, PLFAM_ALIASES),
            )
        )

    protein_index = {protein.feature_id: protein for protein in proteins}

    specialty_orphans: list[str] = []
    specialty_path = m1_dir / "specialty_genes_all.csv"
    if specialty_path.exists():
        for row in _read_table(specialty_path):
            feature_id = _pick(row, ID_ALIASES)
            target = protein_index.get(feature_id)
            if target is None:
                if feature_id:
                    specialty_orphans.append(feature_id)
                continue
            target.specialty.append(
                SpecialtyHit(
                    property_name=_pick(row, PROPERTY_ALIASES),
                    source=_pick(row, SOURCE_ALIASES),
                    product=_pick(row, PRODUCT_ALIASES),
                    identity=_float_or_none(_pick(row, IDENTITY_ALIASES)),
                    query_coverage=_float_or_none(_pick(row, COVERAGE_ALIASES)),
                )
            )

    return InputBundle(
        proteins=proteins,
        specialty_without_protein=specialty_orphans,
        table_rows_without_sequence=sorted(set(by_id) - seen_ids),
        sequences_without_table_row=sequences_without_row,
    )

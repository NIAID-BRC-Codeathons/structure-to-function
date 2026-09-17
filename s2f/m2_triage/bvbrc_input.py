"""Load M1-shaped BV-BRC output: protein FASTA plus feature and specialty tables.

Canonical protein key is the BV-BRC ``patric_id`` (``feature_id`` in
``docs/00-architecture.md``). Column names vary between the BV-BRC web UI and the ``p3-`` CLI,
so headers are matched through aliases rather than assumed.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

ID_ALIASES = ("patric_id", "feature_id", "brc_id", "genome_feature_id")
PRODUCT_ALIASES = ("product", "function")
GENE_ALIASES = ("gene", "gene_name", "gene_symbol")
LOCUS_ALIASES = ("refseq_locus_tag", "locus_tag", "alt_locus_tag")
PGFAM_ALIASES = ("pgfam_id", "pgfam", "pgfams")
PLFAM_ALIASES = ("plfam_id", "plfam", "plfams")
PROPERTY_ALIASES = ("property", "specialty_gene_property", "type")
CLASSIFICATION_ALIASES = ("classification", "amr_classification", "assertion")
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


#: An "Antibiotic Resistance" row does not always mean resistance. BV-BRC classifies most of
#: them as "antibiotic target in susceptible species" — gyrA and rpoB are the *targets* of
#: fluoroquinolones and rifamycins, not resistance genes (issue #30). Scoring those as
#: resistance evidence is the specific mistake cmmann21 flagged, and on G37 it would apply to
#: 10 of 14 AMR rows.
TARGET_IN_SUSCEPTIBLE = "antibiotic target in susceptible species"
RESISTANCE_MARKERS = (
    "conferring antibiotic resistance", "antibiotic efflux", "antibiotic inactivation",
    "efflux pump", "antibiotic target alteration", "antibiotic target replacement",
    "antibiotic target protection", "resistance via absence", "reduced permeability",
)


@dataclass
class SpecialtyHit:
    """One row of the specialty-gene table."""

    property_name: str
    source: str
    product: str
    identity: float | None
    query_coverage: float | None
    classification: str = ""

    @property
    def is_resistance_mechanism(self) -> bool:
        text = (self.classification or "").lower()
        return any(marker in text for marker in RESISTANCE_MARKERS)

    @property
    def is_target_in_susceptible_species(self) -> bool:
        return TARGET_IN_SUSCEPTIBLE in (self.classification or "").lower()


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
    def is_virulence_factor(self) -> bool:
        return bool(self._properties() & VIRULENCE_PROPERTIES)

    def _amr_rows(self) -> list[SpecialtyHit]:
        return [h for h in self.specialty if h.property_name.strip().lower() in AMR_PROPERTIES]

    @property
    def amr_evidence(self) -> tuple[float, str]:
        """How much an AMR row supports *resistance*, and on what basis.

        1.0 a classified resistance mechanism; 0.5 an AMR row whose classification is missing,
        so the export cannot say which it is; 0.0 rows that are only "antibiotic target in
        susceptible species" — those are drug targets and are scored as such instead.
        """
        rows = self._amr_rows()
        if not rows:
            return 0.0, ""
        if any(hit.is_resistance_mechanism for hit in rows):
            mechanisms = sorted({h.classification for h in rows if h.is_resistance_mechanism})
            return 1.0, "; ".join(mechanisms)[:120]
        if all(hit.is_target_in_susceptible_species for hit in rows):
            return 0.0, TARGET_IN_SUSCEPTIBLE
        return 0.5, "AMR row with no classification recorded"

    @property
    def is_antibiotic_target(self) -> bool:
        """Target of an antibiotic in a susceptible species — a drug target, not resistance."""
        return any(hit.is_target_in_susceptible_species for hit in self._amr_rows())

    @property
    def is_essential_ortholog(self) -> bool:
        return bool(self._properties() & ESSENTIAL_PROPERTIES)

    @property
    def is_drug_target(self) -> bool:
        return bool(self._properties() & DRUG_TARGET_PROPERTIES) or self.is_antibiotic_target

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
    #: The query genome's organism, read from the FASTA headers. Used to tell a structural hit
    #: against our own organism from a genuine cross-organism transfer.
    organism: str = ""
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


def _read_records(path: Path) -> list[dict[str, str]]:
    """Read a feature or specialty table as CSV/TSV or as BV-BRC JSON."""
    if path.suffix.lower() == ".json":
        payload = json.loads(path.read_text(encoding="utf-8"))
        if isinstance(payload, dict):
            for key in ("features", "items", "records", "response", "docs"):
                if isinstance(payload.get(key), list):
                    payload = payload[key]
                    break
            else:
                payload = [payload]
        return [
            {_normalize_header(k): ("" if v is None else str(v)) for k, v in row.items()}
            for row in payload
            if isinstance(row, dict)
        ]
    return _read_table(path)


def _find_input(m1_dir: Path, names: tuple[str, ...], patterns: tuple[str, ...]) -> Path | None:
    """Exact filenames first, then glob patterns, so both export layouts work."""
    for name in names:
        candidate = m1_dir / name
        if candidate.exists():
            return candidate
    for pattern in patterns:
        matches = sorted(m1_dir.glob(pattern))
        if matches:
            return matches[0]
    return None


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


FIG_ID = re.compile(r"^(fig\|[^|\s]+)")
ORGANISM_SUFFIX = re.compile(r"\s*\[[^\]]*\]\s*$")


def organism_key(name: str) -> tuple[str, str]:
    """(genus, species) from an organism name, lowercased and stripped of strain detail.

    Genus names get renamed — BV-BRC writes "Mycoplasma genitalium" where the PDB and UniProt
    now say "Mycoplasmoides genitalium" — so comparison falls back to the species epithet when
    the genus disagrees.
    """
    tokens = [t for t in re.split(r"[\s_]+", (name or "").strip()) if t]
    if not tokens:
        return "", ""
    genus = tokens[0].lower().rstrip(".,")
    species = tokens[1].lower().rstrip(".,") if len(tokens) > 1 else ""
    if species in {"sp", "sp.", "subsp", "subsp."} and len(tokens) > 2:
        species = tokens[2].lower().rstrip(".,")
    return genus, species


def same_species(query: str, hit: str) -> bool:
    """Both names denote the same species. Tolerates the genus renaming above."""
    q_genus, q_species = organism_key(query)
    h_genus, h_species = organism_key(hit)
    if not q_species or not h_species:
        return False
    if q_species != h_species:
        return False
    # Same species epithet: accept when the genus matches, or when one genus name is a
    # lengthened form of the other (Mycoplasma -> Mycoplasmoides).
    return q_genus == h_genus or q_genus.startswith(h_genus[:6]) or h_genus.startswith(q_genus[:6])


def same_genus(query: str, hit: str) -> bool:
    q_genus, _ = organism_key(query)
    h_genus, _ = organism_key(hit)
    if not q_genus or not h_genus:
        return False
    return q_genus == h_genus or q_genus.startswith(h_genus[:6]) or h_genus.startswith(q_genus[:6])


def _organism_from_header(header: str) -> str:
    """`... product [Mycoplasma genitalium G37 | 243273.25]` -> the organism name."""
    match = re.search(r"\[([^\]]+)\]\s*$", header or "")
    if not match:
        return ""
    return match.group(1).split("|")[0].strip()


def _feature_id_from_header(header: str) -> tuple[str, str]:
    """Split a BV-BRC FASTA header into feature ID and product.

    Two layouts are in circulation:
    ``fig|1125630.4.peg.5196  DNA replication protein`` (p3-CLI/web export) and
    ``fig|243273.25.peg.308|MG_267|VBIMycGen98045_0308| hypothetical protein [Organism | id]``
    (the fixture on main). Only the ``fig|`` part is the canonical ID.
    """
    header = header.strip()
    first, _, rest = header.partition(" ")
    match = FIG_ID.match(header)
    feature_id = match.group(1) if match else first.strip()
    product = rest.strip() if match is None else header[len(first):].strip()
    product = product.lstrip("| ").strip()
    return feature_id, ORGANISM_SUFFIX.sub("", product).strip()


def _float_or_none(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def load_input(m1_dir: Path) -> InputBundle:
    """Load ``proteins.faa`` plus the feature and (optional) specialty tables."""
    m1_dir = Path(m1_dir)
    # *.faa only: a genome directory also holds nucleotide *.fna, which is not our input.
    fasta_path = _find_input(m1_dir, ("proteins.faa",), ("*.faa",))
    if fasta_path is None:
        raise FileNotFoundError(f"missing protein FASTA (proteins.faa or *.faa) in {m1_dir}")

    table_path = _find_input(
        m1_dir,
        ("genes_proteins.csv", "features.csv", "proteins.tsv", "features.tsv", "features.json"),
        ("*.features.json", "*features*.csv", "*features*.tsv"),
    )
    table_rows = _read_records(table_path) if table_path is not None else []

    by_id: dict[str, dict[str, str]] = {}
    for row in table_rows:
        feature_id = _pick(row, ID_ALIASES)
        if feature_id:
            by_id[feature_id] = row

    proteins: list[Protein] = []
    sequences_without_row: list[str] = []
    seen_ids: set[str] = set()
    organisms: list[str] = []
    for header, sequence in _read_fasta(fasta_path):
        feature_id, header_product = _feature_id_from_header(header)
        organism = _organism_from_header(header)
        if organism:
            organisms.append(organism)
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
    specialty_path = _find_input(
        m1_dir,
        ("specialty_genes_all.csv", "specialty_genes.csv", "sp_gene.json"),
        ("*.sp_gene.json", "*specialty*.csv", "*specialty*.tsv"),
    )
    if specialty_path is not None:
        for row in _read_records(specialty_path):
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
                    classification=_pick(row, CLASSIFICATION_ALIASES),
                )
            )

    return InputBundle(
        proteins=proteins,
        organism=(Counter(organisms).most_common(1)[0][0] if organisms else ""),
        specialty_without_protein=specialty_orphans,
        table_rows_without_sequence=sorted(set(by_id) - seen_ids),
        sequences_without_table_row=sequences_without_row,
    )

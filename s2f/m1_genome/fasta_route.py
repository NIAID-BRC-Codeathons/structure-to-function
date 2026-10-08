"""M1's third route: protein sequences in, no genome and no BV-BRC (issue #89).

    python -m s2f.m1_genome --run runs/q1 --from-fasta query.faa [--organism NAME]

The CGA and API routes describe a genome. This one describes a set of proteins someone handed
us, so everything a genome would have supplied — contigs, coordinates, specialty genes,
subsystems, the taxonomy tree — is recorded as **not available**, never as empty. An empty
`specialty` list on a genome route means "BV-BRC looked and found nothing"; here nobody looked,
and the report must be able to tell the two apart (#90).

FASTA headers are untrusted free text. The first token is the protein id and the rest is kept
as the product, but no organism is read out of them (#55): the organism comes from
`--organism` / `--taxon-id` or not at all.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

#: What a genome route supplies and this route cannot. Recorded in `genome.not_available` and
#: per protein, so a reader of a sequence-only report sees what was never looked at.
GENOME_ONLY = {
    "contigs": "no assembly: the input is protein sequences",
    "coordinates": "no assembly: proteins have no contig, start, end or strand",
    "specialty": "BV-BRC specialty genes need a BV-BRC genome",
    "subsystems": "BV-BRC subsystems need a BV-BRC genome",
    "closest_genomes": "no genome to compare",
    "tree": "no genome to place on a tree",
}

#: UniProt and NCBI headers append the organism and other fields to the product. Dropping them
#: keeps the product readable; the full header is kept in `proteins[].header`.
_HEADER_TAIL = re.compile(r"\s+(?:OS=.*|\[[^\]]*\]\s*)$")
_ID_UNSAFE = re.compile(r"[^A-Za-z0-9_.|:-]")
_AMINO = re.compile(r"^[A-Za-z*.-]+$")
_NUCLEOTIDE = re.compile(r"^[ACGTUN]+$")


def read_protein_fasta(path: str | Path) -> list[tuple[str, str]]:
    """`(header, sequence)` pairs, sequences uppercased with gaps and stop marks removed."""
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks)))
            header, chunks = line[1:].strip(), []
        elif header is None:
            raise ValueError(f"{path}: sequence data before the first '>' header")
        elif not _AMINO.match(line):
            raise ValueError(f"{path}: not a protein sequence line under {header[:60]!r}")
        else:
            chunks.append(re.sub(r"[*.-]", "", line).upper())
    if header is not None:
        records.append((header, "".join(chunks)))
    # A, C, G and T are amino-acid letters too, so a nucleotide FASTA parses cleanly and
    # would be searched as nonsense proteins. Refuse it rather than score it.
    for header, sequence in records:
        if len(sequence) >= 30 and _NUCLEOTIDE.match(sequence):
            raise ValueError(f"{path}: {header[:60]!r} looks like nucleotide sequence; "
                             "--from-fasta takes protein sequences")
    return records


def query_proteins(records: list[tuple[str, str]]) -> list[dict[str, Any]]:
    """One record per sequence: a unique id, the product text, the sequence and its header."""
    proteins: list[dict[str, Any]] = []
    seen: set[str] = set()
    for index, (header, sequence) in enumerate(records, start=1):
        if not sequence:
            continue
        first, _, rest = header.partition(" ")
        base = _ID_UNSAFE.sub("_", first) or f"query_{index}"
        feature_id, n = base, 2
        while feature_id in seen:  # duplicate ids in user FASTA are common; never merge them
            feature_id, n = f"{base}_{n}", n + 1
        seen.add(feature_id)
        proteins.append({
            "feature_id": feature_id,
            "product": _HEADER_TAIL.sub("", rest).strip() or None,
            "header": header,
            "sequence": sequence,
        })
    return proteins


def write_fasta_m1_dir(proteins: list[dict[str, Any]], m1_dir: str | Path) -> int:
    """Write `proteins.faa`, the one file M2 needs. Headers carry no organism bracket."""
    m1_dir = Path(m1_dir)
    m1_dir.mkdir(parents=True, exist_ok=True)
    with (m1_dir / "proteins.faa").open("w", encoding="utf-8") as handle:
        for protein in proteins:
            handle.write(f">{protein['feature_id']} {protein['product'] or ''}".rstrip() + "\n")
            sequence = protein["sequence"]
            for i in range(0, len(sequence), 60):
                handle.write(sequence[i:i + 60] + "\n")
    return len(proteins)


def proteins_section(proteins: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """`proteins[]` for report.json. Sequences stay in the FASTA, as on the other routes."""
    return [{
        "feature_id": p["feature_id"],
        "product": p["product"],
        "header": p["header"],
        "aa_length": len(p["sequence"]),
        "not_available": ["specialty", "subsystems", "coordinates"],
    } for p in proteins]


def genome_section(*, source: str, n_proteins: int, organism: str | None,
                   taxon_id: int | None) -> dict[str, Any]:
    """`genome` for report.json: what the user told us, and what no genome could supply."""
    return {
        "annotation_route": "fasta",
        "input_fasta": source,
        "protein_count": n_proteins,
        "taxonomy": {"scientific_name": organism, "source": "user"} if organism else None,
        "taxon_id": taxon_id,
        "closest_genomes": [],
        "not_available": dict(GENOME_ONLY),
    }

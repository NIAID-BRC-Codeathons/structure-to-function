"""Parse a BV-BRC Comprehensive Genome Analysis output directory (issue #5).

Pure parsing: no network, no `s2f.common` imports, so it can be tested against a
downloaded CGA directory on its own. Layout and field meanings are documented in
`docs/01a-cga-coverage.md`.

Produces three things:

* the `proteins[]` and `genome` sections of `report.json`
* `run` manifest fields (job id, runtime, tool versions)
* the `<run>/m1/` files M2 reads: `proteins.faa`, `genes_proteins.csv`,
  `specialty_genes_all.csv`
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

FEATURE_COLUMNS = [
    "patric_id", "refseq_locus_tag", "gene", "product", "aa_length",
    "start", "end", "strand", "plfam_id", "pgfam_id",
]
SPECIALTY_COLUMNS = [
    "property", "gene", "product", "source", "classification",
    "antibiotics_class", "identity", "query_coverage", "patric_id",
]

#: Genome-quality values that must not be handed downstream without a decision.
POOR_QUALITY = {"poor"}


class CgaLayoutError(RuntimeError):
    """The directory does not look like CGA output."""


@dataclass
class CgaRun:
    """One parsed CGA output directory."""

    root: Path
    job: dict[str, Any]
    genome: dict[str, Any]
    annotated: dict[str, Any]
    features: list[dict[str, Any]]
    specialty: list[dict[str, Any]]
    subsystem_rows: list[dict[str, Any]]
    sequences: dict[str, str]
    tree_newick: str | None = None
    ingroup: list[str] = field(default_factory=list)

    # -- quality ---------------------------------------------------------
    @property
    def genome_id(self) -> str:
        return str(self.genome.get("genome_id") or self.annotated.get("id") or "")

    @property
    def quality(self) -> str:
        return str(self.genome.get("genome_quality") or "")

    @property
    def quality_flags(self) -> list[str]:
        return list((self.annotated.get("quality") or {}).get("genome_quality_flags") or [])

    @property
    def is_poor(self) -> bool:
        return self.quality.strip().lower() in POOR_QUALITY or bool(self.quality_flags)

    def quality_reason(self) -> str:
        bits = []
        if self.quality.strip().lower() in POOR_QUALITY:
            bits.append(f"genome_quality={self.quality}")
        if self.quality_flags:
            bits.append("flags=" + ",".join(self.quality_flags))
        return "; ".join(bits)


def _as_int(value: Any) -> int | None:
    """BV-BRC writes coordinates as strings; the schema wants integers or null."""
    if value is None or value == "":
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _read_text(path: Path) -> str:
    """Read a CGA file whose encoding is not guaranteed.

    BV-BRC writes bytes that are not valid UTF-8 into some tables - a 0xa0
    (Latin-1 non-breaking space) inside a real `sp_gene.json` killed a strict
    read. Fall back rather than fail: the field values matter, the exact
    codepoint of a stray space does not.
    """
    data = path.read_bytes()
    for encoding in ("utf-8", "utf-8-sig", "cp1252", "latin-1"):
        try:
            return data.decode(encoding)
        except UnicodeDecodeError:
            continue
    return data.decode("utf-8", errors="replace")


def _as_float(value: Any) -> float | None:
    """Numeric fields arrive as strings from some sources and numbers from others.

    In one USA300 run, DIAMOND rows carried `identity: '99'` while AMRFinderPlus and
    RGI carried `identity: 99.77`. Downstream comparisons need one type. None stays
    None: a k-mer row has no identity, which is not the same as zero.
    """
    if value is None or value == "":
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _load_json(path: Path) -> Any:
    return json.loads(_read_text(path))


def _read_fasta(path: Path) -> dict[str, str]:
    """feature_id -> sequence, keyed on the `fig|...` part of the header."""
    out: dict[str, str] = {}
    key: str | None = None
    chunks: list[str] = []
    for line in _read_text(path).splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if key:
                out[key] = "".join(chunks)
            # Two header shapes are in circulation:
            #   >fig|243273.147.peg.4 DNA gyrase subunit A [Organism ...]
            #   >fig|243273.25.peg.308|MG_267|VBIMycGen98045_0308| hypothetical protein [...]
            # Only the fig| id is canonical, and it ends at the first space or pipe.
            token = line[1:].split()[0]
            key = "|".join(token.split("|")[:2]) if token.startswith("fig|") else token
            chunks = []
        else:
            chunks.append(line)
    if key:
        out[key] = "".join(chunks)
    return out


def load_cga(cga_dir: str | Path) -> CgaRun:
    """Read a CGA output directory. Accepts the `.name_cga` dir or its parent."""
    root = Path(cga_dir)
    if not (root / "annotation").exists():
        nested = [p for p in root.iterdir() if p.is_dir() and (p / "annotation").exists()]
        if len(nested) != 1:
            raise CgaLayoutError(
                f"no CGA 'annotation' job record under {root} "
                f"(found {len(nested)} candidate subdirectories)"
            )
        root = nested[0]

    ann = root / ".annotation"
    load = ann / "load_files"
    for required in (root / "annotation", ann / "annotation.genome", load / "genome.json"):
        if not required.exists():
            raise CgaLayoutError(f"missing expected CGA file: {required}")

    genome = _load_json(load / "genome.json")
    if isinstance(genome, list):
        genome = genome[0] if genome else {}

    tree = root / "codontree_tree.nwk"
    ingroup_file = root / "tree_ingroup.txt"
    return CgaRun(
        root=root,
        job=_load_json(root / "annotation"),
        genome=genome,
        annotated=_load_json(ann / "annotation.genome"),
        features=_load_json(load / "genome_feature.json"),
        specialty=_load_json(load / "sp_gene.json"),
        subsystem_rows=_load_json(load / "subsystem.json"),
        sequences=_read_fasta(ann / "annotation.feature_protein.fasta"),
        tree_newick=_read_text(tree).strip() if tree.exists() else None,
        ingroup=ingroup_file.read_text().split() if ingroup_file.exists() else [],
    )


def _subsystems_by_feature(run: CgaRun) -> dict[str, list[dict[str, Any]]]:
    out: dict[str, list[dict[str, Any]]] = {}
    for sub in run.annotated.get("subsystems") or []:
        name = sub.get("name")
        classification = sub.get("classification") or []
        for binding in sub.get("role_bindings") or []:
            for fid in binding.get("features") or []:
                out.setdefault(fid, []).append(
                    {"name": name, "classification": classification,
                     "role": binding.get("role_id")}
                )
    return out


def _specialty_by_feature(run: CgaRun) -> dict[str, list[dict[str, Any]]]:
    """Specialty rows keyed by feature, keeping `classification` and the evidence type.

    `K-mer Search` rows carry no identity or coverage; their confidence is the
    per-feature k-mer count in `annotation.genome` (docs/01a-cga-coverage.md).
    Nulls stay null rather than becoming zero.
    """
    out: dict[str, list[dict[str, Any]]] = {}
    for row in run.specialty:
        fid = row.get("patric_id")
        if not fid:
            continue
        antibiotics = row.get("antibiotics") or []
        out.setdefault(fid, []).append(
            {
                "type": row.get("property"),
                "classification": row.get("classification"),
                "database": row.get("source"),
                "evidence": row.get("evidence"),
                "gene": row.get("gene"),
                "hit": row.get("source_id") or row.get("product"),
                "identity": _as_float(row.get("identity")),
                "coverage": _as_float(row.get("query_coverage")),
                "subject_coverage": _as_float(row.get("subject_coverage")),
                "e_value": _as_float(row.get("e_value")),
                "same_species": row.get("same_species"),
                "same_genus": row.get("same_genus"),
                "antibiotics": antibiotics if isinstance(antibiotics, list) else [antibiotics],
                "pmid": row.get("pmid") or [],
            }
        )
    return out


def _kmer_confidence(run: CgaRun) -> dict[str, dict[str, Any]]:
    out: dict[str, dict[str, Any]] = {}
    for feature in run.annotated.get("features") or []:
        quality = feature.get("quality") or {}
        if quality:
            out[feature.get("id")] = {
                "hit_count": quality.get("hit_count"),
                "weighted_hit_count": quality.get("weighted_hit_count"),
                "priority": quality.get("priority"),
            }
    return out


def proteins_section(run: CgaRun) -> list[dict[str, Any]]:
    """`proteins[]` for report.json. Sequences are NOT included (pitfall: keep them in the FASTA)."""
    subs = _subsystems_by_feature(run)
    spec = _specialty_by_feature(run)
    kmer = _kmer_confidence(run)
    proteins = []
    for row in run.features:
        if row.get("feature_type") != "CDS":
            continue
        fid = row.get("patric_id")
        if not fid:
            continue
        proteins.append(
            {
                "feature_id": fid,
                "locus_tag": row.get("refseq_locus_tag"),
                "gene": row.get("gene"),
                "product": row.get("product"),
                "contig": row.get("sequence_id") or row.get("accession"),
                "start": _as_int(row.get("start")),
                "end": _as_int(row.get("end")),
                "strand": row.get("strand"),
                "aa_length": _as_int(row.get("aa_length")),
                "subsystems": subs.get(fid, []),
                "specialty": spec.get(fid, []),
                "families": {"pgfam": row.get("pgfam_id"), "plfam": row.get("plfam_id")},
                "kmer_confidence": kmer.get(fid),
            }
        )
    return proteins


def genome_section(run: CgaRun, taxon_call: dict[str, Any] | None = None) -> dict[str, Any]:
    """`genome` for report.json.

    `closest_genomes` comes from the Minhash call when one was made; CGA's own
    `close_genomes` is empty in practice (docs/01a-cga-coverage.md). The codon-tree
    ingroup is reported separately because it is chosen from sequence, not taxonomy.

    `annotation_route` is `"cga"`: these features describe the *submitted assembly*.
    The API route writes `"api"` for the same key, where they describe a reference
    genome for the same organism instead, and gene presence or absence is not a
    property of the isolate.
    """
    quality = run.annotated.get("quality") or {}
    closest = []
    if taxon_call:
        for hit in taxon_call.get("hits") or []:
            closest.append(
                {
                    "genome_id": hit.get("genome_id"),
                    "name": hit.get("genome_name"),
                    "mash_distance": hit.get("distance"),
                    "pvalue": hit.get("pvalue"),
                    "shared_kmers": hit.get("counts"),
                    "ani": None,
                    "snp_distance": None,
                }
            )
    return {
        "genome_id": run.genome_id,
        "taxon_id": int(run.genome["taxon_id"]) if run.genome.get("taxon_id") else None,
        "taxonomy": {
            "scientific_name": run.annotated.get("scientific_name"),
            "lineage_names": run.annotated.get("ncbi_lineage"),
            "genetic_code": run.annotated.get("genetic_code"),
            "called_by": (taxon_call or {}).get("called_by", "supplied"),
            "called_rank": (taxon_call or {}).get("called_rank"),
            "top_hit_distance": ((taxon_call or {}).get("top_hit") or {}).get("distance"),
        },
        "closest_genomes": closest,
        "tree_ingroup": run.ingroup,
        "tree_newick": run.tree_newick,
        "cga_job_id": str(run.job.get("id")) if run.job.get("id") else None,
        # Stated, not inferred. The API route sets this to "api"; a consumer asking
        # `annotation_route == "cga"` used to match nothing at all, because only the
        # other route declared itself. The renderers' `or ("cga" if cga_job_id ...)`
        # fallback stays for runs written before this, but it is no longer the only way
        # to identify a CGA run.
        "annotation_route": "cga",
        "quality": {
            "genome_quality": run.quality,
            "genome_quality_flags": run.quality_flags,
            "coarse_consistency": run.genome.get("coarse_consistency"),
            "fine_consistency": run.genome.get("fine_consistency"),
            "cds_ratio": run.genome.get("cds_ratio"),
            "hypothetical_cds_ratio": run.genome.get("hypothetical_cds_ratio"),
            "contigs": run.genome.get("contigs"),
            "genome_length": run.genome.get("genome_length"),
            "feature_summary": quality.get("feature_summary"),
            "protein_summary": quality.get("protein_summary"),
        },
    }


def run_fields(run: CgaRun) -> dict[str, Any]:
    """Job identity, wall-clock and tool versions for the `run` manifest."""
    tools = sorted({
        str(event.get("tool_name"))
        for event in run.annotated.get("analysis_events") or []
        if event.get("tool_name")
    })
    return {
        "cga": {
            "job_id": str(run.job.get("id")),
            "elapsed_seconds": run.job.get("elapsed_time"),
            "hostname": run.job.get("hostname"),
            "parameters": run.job.get("parameters"),
            "genetic_code": run.annotated.get("genetic_code"),
        },
        "tool_versions": {"cga_analysis_events": tools},
    }


def write_m1_dir(run: CgaRun, m1_dir: str | Path) -> dict[str, int]:
    """Write the files M2 reads: proteins.faa + the two tables."""
    m1_dir = Path(m1_dir)
    m1_dir.mkdir(parents=True, exist_ok=True)

    written = 0
    with (m1_dir / "proteins.faa").open("w", encoding="utf-8") as fh:
        for row in run.features:
            if row.get("feature_type") != "CDS":
                continue
            fid = row.get("patric_id")
            seq = run.sequences.get(fid)
            if not fid or not seq:
                continue
            fh.write(f">{fid} {row.get('product') or ''}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")
            written += 1

    with (m1_dir / "genes_proteins.csv").open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FEATURE_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        rows = 0
        for row in run.features:
            if row.get("feature_type") != "CDS" or not row.get("patric_id"):
                continue
            writer.writerow({k: row.get(k, "") for k in FEATURE_COLUMNS})
            rows += 1

    with (m1_dir / "specialty_genes_all.csv").open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=SPECIALTY_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        spec_rows = 0
        for row in run.specialty:
            antibiotics = row.get("antibiotics_class") or row.get("antibiotics") or ""
            writer.writerow({
                "property": row.get("property", ""),
                "gene": row.get("gene", ""),
                "product": row.get("product", ""),
                "source": row.get("source", ""),
                "classification": row.get("classification", ""),
                "antibiotics_class": ";".join(antibiotics) if isinstance(antibiotics, list) else antibiotics,
                "identity": "" if row.get("identity") is None else row.get("identity"),
                "query_coverage": "" if row.get("query_coverage") is None else row.get("query_coverage"),
                "patric_id": row.get("patric_id", ""),
            })
            spec_rows += 1

    return {"proteins_faa": written, "feature_rows": rows, "specialty_rows": spec_rows}

#!/usr/bin/env python3
"""Report which CGA analyses produced data, for one or more retrieved CGA directories.

    python scripts/cga_coverage_report.py runs/usa300/cga [runs/mgen_G37/cga ...]

Answers issue #33: several CGA analyses run but return nothing for a
poorly-characterized organism, and we need to know which of them populate for a
well-characterized pathogen before M2 and M5 depend on them. Output is markdown,
to paste into the issue.

Reads only files already on disk. No network.
"""

from __future__ import annotations

import collections
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from s2f.m1_genome.parse import _read_text, load_cga  # noqa: E402


def _json(path: Path, default=None):
    try:
        return json.loads(_read_text(path))
    except (OSError, json.JSONDecodeError):
        return default


def _data_rows(path: Path) -> int:
    """Rows in a TSV beyond the header; -1 when the file is missing."""
    if not path.exists():
        return -1
    lines = [l for l in _read_text(path).splitlines() if l.strip()]
    return max(0, len(lines) - 1)


def survey(cga_dir: str) -> dict:
    run = load_cga(cga_dir)
    ann = run.root / ".annotation"
    load = ann / "load_files"

    specialty = collections.Counter(
        (r.get("property"), r.get("source"), r.get("evidence")) for r in run.specialty
    )
    by_property = collections.Counter(r.get("property") for r in run.specialty)
    # BV-BRC carries both spellings; treat them as one thing.
    virulence = sum(n for p, n in by_property.items() if p and p.lower().startswith("virul"))

    typing = _json(load / "genome_typing.json", []) or []
    cgmlst = None
    for row in typing:
        if str(row.get("method", "")).lower() == "cgmlst":
            try:
                cgmlst = float(row.get("pct_called"))
            except (TypeError, ValueError):
                cgmlst = None
    mlst_events = [t for t in (run.annotated.get("typing") or []) if t.get("typing_method") == "MLST"]
    mlst = None
    if mlst_events:
        # The analysis ran; sequence_type stays null when it called nothing.
        mlst = mlst_events[0].get("sequence_type") or 0

    # The tree table lists every ingroup genome, ours among them, not necessarily first.
    genes_per_genome = run.root / "codontree.genesPerGenome.txt"
    single_copy = None
    if genes_per_genome.exists():
        rows = [r.split("\t") for r in genes_per_genome.read_text().splitlines() if r.strip()]
        if len(rows) > 1:
            header = rows[0]
            col = header.index("Filtered_SingleCopy") if "Filtered_SingleCopy" in header else -1
            ours = [r for r in rows[1:] if r and r[0] == run.genome_id]
            row = ours[0] if ours else max(rows[1:], key=lambda r: int(r[col]) if r[col].isdigit() else 0)
            single_copy = int(row[col]) if row[col].lstrip("-").isdigit() else row[col]

    return {
        "dir": str(run.root),
        "genome_id": run.genome_id,
        "name": run.annotated.get("scientific_name"),
        "taxon_id": run.genome.get("taxon_id"),
        "genetic_code": run.annotated.get("genetic_code"),
        "length": run.genome.get("genome_length"),
        "contigs": run.genome.get("contigs"),
        "plasmids": (run.annotated.get("quality") or {}).get("plasmids"),
        "cds": run.genome.get("patric_cds"),
        "quality": run.quality,
        "flags": run.quality_flags,
        "runtime_s": round(run.job.get("elapsed_time") or 0, 1),
        "specialty_total": len(run.specialty),
        "specialty_by_property": dict(by_property),
        "specialty_detail": {" | ".join(str(x) for x in k): v for k, v in sorted(specialty.items())},
        "virulence_rows": virulence,
        "genome_amr_rows": len(_json(load / "genome_amr.json", []) or []),
        "taxonomy_rows": len(_json(load / "taxonomy.json", []) or []),
        "amrfinder_rows": _data_rows(ann / "specialty-amrfinder.txt"),
        "rgi_rows": _data_rows(ann / "specialty-rgi.txt"),
        "specialty_blast_rows": _data_rows(ann / "specialty-blast.txt"),
        "mlst": mlst,
        "cgmlst_pct_called": cgmlst,
        "subsystem_rows": len(run.subsystem_rows),
        "pathway_rows": len(_json(load / "pathway.json", []) or []),
        "tree_ingroup": len(run.ingroup),
        "tree_single_copy_genes": single_copy,
    }


def verdict(value) -> str:
    """Distinguish "did not run" from "ran and found nothing" - the whole point here."""
    if value is None:
        return "not run"
    if value == -1:
        return "file missing"
    if isinstance(value, (int, float)) and not value:
        return "**ran, empty**"
    return f"populated ({value})"


ROWS = [
    ("Virulence (VFDB, Victors)", "virulence_rows"),
    ("Specialty genes, all sources", "specialty_total"),
    ("AMR phenotype (genome_amr: MIC/SIR)", "genome_amr_rows"),
    ("AMRFinderPlus table", "amrfinder_rows"),
    ("RGI / CARD table", "rgi_rows"),
    ("DIAMOND specialty hits (identity + coverage)", "specialty_blast_rows"),
    ("MLST", "mlst"),
    ("cgMLST, percent of loci called", "cgmlst_pct_called"),
    ("taxonomy.json", "taxonomy_rows"),
    ("Subsystem rows", "subsystem_rows"),
    ("Pathway rows", "pathway_rows"),
    ("Codon-tree ingroup genomes", "tree_ingroup"),
    ("Codon-tree single-copy genes used", "tree_single_copy_genes"),
]


def main(argv: list[str]) -> int:
    if not argv:
        print(__doc__)
        return 2
    surveys = [survey(d) for d in argv]

    print("| | " + " | ".join(f"{s['name']} ({s['genome_id']})" for s in surveys) + " |")
    print("| --- |" + " --- |" * len(surveys))
    for label, key in [("Genome length", "length"), ("Contigs", "contigs"), ("CDS", "cds"),
                       ("Genetic code", "genetic_code"), ("CGA compute (s)", "runtime_s"),
                       ("Genome quality", "quality")]:
        print(f"| {label} | " + " | ".join(str(s[key]) for s in surveys) + " |")
    print()
    print("| Analysis | " + " | ".join(s["genome_id"] for s in surveys) + " |")
    print("| --- |" + " --- |" * len(surveys))
    for label, key in ROWS:
        print(f"| {label} | " + " | ".join(verdict(s[key]) for s in surveys) + " |")

    for s in surveys:
        print(f"\n**{s['genome_id']} specialty-gene detail** "
              f"({s['specialty_total']} rows, quality {s['quality']}, flags {s['flags'] or 'none'}):\n")
        if not s["specialty_detail"]:
            print("- none")
        for key, count in sorted(s["specialty_detail"].items(), key=lambda kv: -kv[1]):
            print(f"- {key}: {count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))

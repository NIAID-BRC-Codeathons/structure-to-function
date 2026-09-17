"""Plain-markdown summary of an M1 run, rendered from `report.json`.

The HTML report is the one people read; this is the one they can diff, grep, paste into a
pull request or an issue, and read over SSH. It is deliberately a summary rather than a
second full report: counts, the top candidates, and the caveats. Anything exhaustive lives
in `report.json` and the CSVs.

    python -m s2f.m1_genome.report_md --run runs/mgen_G37

Like `report_html`, it reads the contract rather than a collector's in-memory result, so it
renders either M1 route and runs made before it existed.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any

from ..common.io import read_report
from .pathogens import resolve_disease_profile
from .priority import CAT_LABEL
from .report_html import _mapping, _pairs, _rows, collect_priorities

TOP_N = 15


def _clean(text: Any) -> str:
    """Flatten a value onto one line and neutralise the pipe that would split a cell."""
    return str("" if text is None else text).replace("|", "\\|").replace("\n", " ").strip()


def build_markdown(report: dict[str, Any], *, top: int = TOP_N) -> str:
    """Render the summary."""
    genome = _mapping(report.get("genome"))
    run = _mapping(report.get("run"))
    quality = _mapping(genome.get("quality"))
    nutrition = _mapping(genome.get("nutrition"))
    isolation = _mapping(genome.get("isolation"))
    route = genome.get("annotation_route") or ("cga" if genome.get("cga_job_id") else "unknown")
    name = _mapping(genome.get("taxonomy")).get("scientific_name") or genome.get("genome_id")

    priorities = collect_priorities(report)
    selected = [p for p in priorities if p["selected"]]
    profile = resolve_disease_profile({
        "species": genome.get("species"), "genus": genome.get("genus"),
        "disease": genome.get("disease"),
    })

    out: list[str] = []
    w = out.append

    w(f"# M1 genome report — {_clean(name)}")
    w("")
    w(f"- **Run:** `{_clean(run.get('run_id'))}`  ·  "
      f"**generated:** {_clean(run.get('created_at'))}")
    w(f"- **BV-BRC genome:** `{_clean(genome.get('genome_id'))}`  ·  "
      f"**taxon:** `{_clean(genome.get('taxon_id'))}`")
    w(f"- **Acquisition route:** `{route}`  ·  "
      f"**annotation:** {_clean(genome.get('annotation_source') or 'CGA')}")
    facts = []
    if quality.get("genome_length"):
        facts.append(f"{_clean(quality['genome_length'])} bp")
    if genome.get("gc_content") not in (None, ""):
        facts.append(f"{_clean(genome['gc_content'])}% GC")
    facts.append(f"{len(_rows(report.get('proteins')))} proteins")
    facts.append(f"quality {_clean(quality.get('genome_quality') or 'unknown')}")
    w(f"- **Genome:** {', '.join(facts)}")
    if genome.get("bvbrc_url"):
        w(f"- **Live BV-BRC report:** {genome['bvbrc_url']}")
    w("")

    if route == "api":
        w("> **This report describes a reference genome, not the submitted assembly.** The")
        w("> API route characterises a BV-BRC reference or representative genome for the")
        w("> same organism, so an isolate-specific gene is invisible here and a")
        w("> reference-only gene is a false positive. See `docs/01b-m1-api-mode.md`.")
        if genome.get("resolution"):
            w(">")
            w(f"> Resolved by: {_clean(genome['resolution'])}")
        w("")

    w("## Proteins and priority ranking")
    w("")
    w(f"- Proteins ranked: **{len(priorities)}**")
    w(f"- Selected for Module 2: **{len(selected)}**")
    w(f"- With a host-interaction mechanism: "
      f"**{sum(1 for p in priorities if p['categories'])}**")
    w(f"- Virulence-flagged: "
      f"**{sum(1 for p in priorities if 'virulence' in p['specialty_types'])}**")
    w("")
    if priorities:
        w(f"Top {min(top, len(priorities))} by M1 pathogenesis priority:")
        w("")
        w("| Rank | Score | Gene | Product | aa | Mechanism |")
        w("| ---: | ----: | --- | --- | ---: | --- |")
        for record in priorities[:top]:
            mechanisms = ", ".join(CAT_LABEL[c] for c in record["categories"]
                                   if c in CAT_LABEL) or "—"
            w(f"| {_clean(record['rank'])} | {_clean(record['score'])} "
              f"| {_clean(record.get('gene') or '—')} "
              f"| {_clean(record.get('product'))} "
              f"| {_clean(record.get('aa_length'))} | {_clean(mechanisms)} |")
        w("")
        w("Full ranking with per-protein score breakdowns: `report.json` "
          "(`proteins[].m1_priority`) and `m1/proteins_ranked.csv`.")
        w("")

    counts = {str(k): v for k, v in _mapping(genome.get("specialty_gene_counts")).items()}
    if counts:
        w("## Specialty genes")
        w("")
        for key, value in sorted(counts.items(), key=lambda kv: str(kv[1]), reverse=True):
            w(f"- {_clean(key)}: **{_clean(value)}**")
        w("")

    phenotypes = _rows(genome.get("amr_phenotypes"))
    if phenotypes:
        resistant = sorted({p["antibiotic"] for p in phenotypes
                            if p.get("resistant_phenotype") == "Resistant"
                            and p.get("antibiotic")})
        w("## Antimicrobial resistance")
        w("")
        w(f"- Laboratory/computational phenotype records: **{len(phenotypes)}**")
        w(f"- Recorded resistant to: {_clean(', '.join(resistant)) or '_none recorded_'}")
        w("- Absence of a phenotype record is not evidence of susceptibility.")
        w("")

    growth = _rows(genome.get("growth"))
    if growth:
        w("## Growth conditions")
        w("")
        for row in growth:
            w(f"- **{_clean(row.get('property'))}:** {_clean(row.get('value'))} "
              f"_({_clean(row.get('source'))})_")
        w("")

    own = _rows(isolation.get("genome"))
    distribution = _mapping(isolation.get("species_distribution"))
    if own or distribution:
        w("## Isolation and origin")
        w("")
        for row in own:
            w(f"- **{_clean(row.get('property'))}:** {_clean(row.get('value'))}")
        for label, items in distribution.items():
            rendered = ", ".join(f"{_clean(v)} ({_clean(c)})" for v, c in _pairs(items))
            if rendered:
                w(f"- **{_clean(label)}** (species-wide): {rendered}")
        w("")

    if nutrition:
        w("## Nutrition requirement (inferred)")
        w("")
        w(f"- **Oxygen requirement:** "
          f"{_clean(nutrition.get('oxygen_requirement') or 'unknown')}")
        w(f"- Amino-acid biosynthesis subsystems encoded: "
          f"**{_clean(nutrition.get('amino_acid_biosynthesis_count'))}**")
        w(f"- Cofactor/vitamin biosynthesis subsystems encoded: "
          f"**{_clean(nutrition.get('cofactor_vitamin_biosynthesis_count'))}**")
        w(f"- _{_clean(nutrition.get('basis') or 'inferred from encoded pathways')}._")
        w("")

    close = _mapping(genome.get("close_human_pathogens"))
    if close.get("species"):
        w("## Related human-associated species")
        w("")
        for row in _rows(close["species"]):
            w(f"- *{_clean(row.get('species'))}* — "
              f"{_clean(row.get('human_associated_genomes'))} human-host genomes")
        w(f"- _Basis: {_clean(close.get('basis'))}. A species here means BV-BRC holds "
          f"human isolates of it, not that it is a validated human pathogen._")
        w("")

    if profile.get("diseases") or profile.get("symptoms"):
        w("## Disease context")
        w("")
        w(f"_Source: {_clean(profile.get('source'))}. {_clean(profile.get('scope'))}._")
        w("")
        for title, key in [("Diseases", "diseases"), ("Symptoms", "symptoms"),
                           ("Transmission", "spread"), ("Invasion strategy", "invasion")]:
            values = profile.get(key) or []
            if values:
                w(f"**{title}**")
                w("")
                for value in values:
                    w(f"- {_clean(value)}")
                w("")

    w("## Hand-off and limitations")
    w("")
    w("Module 2 reads `report.json` plus `m1/proteins.faa`, `m1/genes_proteins.csv` and")
    w("`m1/specialty_genes_all.csv`. The M1 priority score ranks on annotation evidence;")
    w("M2's `triage` score ranks on structural tractability. They are different questions")
    w("and are expected to disagree.")
    w("")
    w("Nutrition and growth statements are inferences from encoded pathways, not assays.")
    w("Disease text describes the species from curated literature and BV-BRC metadata, not")
    w("this isolate. Mechanism hypotheses are annotation-driven and need experimental")
    w("validation. Research use only; no clinical or treatment implication is intended.")
    w("")
    return "\n".join(out)


def render_run(run_dir: str | Path, *, top: int = TOP_N) -> str:
    return build_markdown(read_report(Path(run_dir)), top=top)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m s2f.m1_genome.report_md", description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run", required=True, help="run directory holding report.json")
    parser.add_argument("--out", help="output path (default: <run>/m1/report.md)")
    parser.add_argument("--top", type=int, default=TOP_N,
                        help=f"proteins in the ranking table (default {TOP_N})")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    run_dir = Path(args.run)
    if not (run_dir / "report.json").exists():
        print(f"no report.json in {run_dir}; run M1 first", file=sys.stderr)
        return 2
    text = render_run(run_dir, top=args.top)
    out = Path(args.out) if args.out else run_dir / "m1" / "report.md"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text, encoding="utf-8")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Run the pipeline end to end.

    python -m s2f.run --run runs/mgen_G37 --contigs assembly.fna \\
        --ws-dir /USER@bvbrc/home/s2f

Stages run in order and each writes its own section of ``runs/<id>/report.json``.
A stage that is not implemented yet is skipped and recorded as such, so this is
runnable today with M1 and M2 and grows as the other modules land.

Offline, from a CGA directory already retrieved:

    python -m s2f.run --run runs/mgen_G37 --from-cga-dir data/cga --skip m2
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

#: stage -> module, in execution order. M3-M6 appear here as they are implemented.
STAGES: list[tuple[str, str]] = [
    ("m1", "s2f.m1_genome"),
    ("m2", "s2f.m2_triage"),
    ("m3", "s2f.m3_fold"),
    ("m6", "s2f.m6_report"),
    # Last, and deliberately after m6: it reads the whole run directory, including what
    # the other stages just wrote into their own manifests.
    ("full", "s2f.full_report"),
]


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="python -m s2f.run", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--run", required=True, help="run directory, e.g. runs/mgen_G37")
    p.add_argument("--contigs", help="assembly FASTA for M1")
    p.add_argument("--ws-dir", help="BV-BRC workspace dir, e.g. /USER@bvbrc/home/s2f")
    p.add_argument("--from-cga-dir", help="M1 parses this directory instead of submitting a job")
    p.add_argument("--from-fasta", metavar="FAA",
                   help="M1 starts from protein sequences: no genome, no BV-BRC (#89)")
    p.add_argument("--organism", help="M1 --from-fasta route: organism of the sequences")
    p.add_argument("--from-bvbrc-api", action="store_true",
                   help="M1 uses the BV-BRC Data API route instead of the CGA service")
    p.add_argument("--species", help="M1 API route: species name to resolve")
    p.add_argument("--genome-id", help="M1 API route: explicit BV-BRC genome id")
    p.add_argument("--html", action="store_true",
                   help="M1 also writes the interactive HTML report")
    p.add_argument("--figures", action="store_true",
                   help="M1 also writes publication figures (needs matplotlib)")
    p.add_argument("--annotation", choices=["PATRIC", "RefSeq"],
                   help="M1 API route: annotation source for CDS features")
    p.add_argument("--max-features", type=int, help="M1 API route: cap on CDS features")
    p.add_argument("--seq-cap", type=int, help="M1 API route: cap on protein sequences")
    p.add_argument("--tree-leaves", type=int,
                   help="M1 API route: max genomes on the gene-content phylogeny")
    p.add_argument("--no-tree", action="store_true",
                   help="M1 API route: skip the gene-content phylogeny")
    p.add_argument("--taxon-id", type=int, help="M1 CGA route: skip the Minhash call")
    p.add_argument("--genetic-code", type=int,
                   help="M1 CGA route: translation table; required with --taxon-id")
    p.add_argument("--skip", action="append", default=[], metavar="STAGE",
                   help="skip a stage; repeatable")
    p.add_argument("--only", action="append", default=[], metavar="STAGE",
                   help="run only these stages; repeatable")
    # One flag per stage (#80). --limit used to reach M3 as well, so capping structure
    # downloads silently capped how many proteins M2 scored.
    p.add_argument("--limit", type=int, help="M2: score only the first N proteins")
    p.add_argument("--structure-limit", type=int,
                   help="M3: collect structures for only the first N selected proteins")
    p.add_argument("--map-ids", action="store_true", help="M2: resolve UniProt/ChEMBL xrefs")
    p.add_argument("--kg", action="store_true", help="M2: assemble the knowledge subgraph")
    p.add_argument("--annotate", action="store_true",
                   help="M2: function, localization and membrane flags (docs/02c)")
    p.add_argument("--interproscan", metavar="PATH_OR_AUTO",
                   help="M2: InterProScan output to ingest, or 'auto' to find a local install")
    for tool in ("eggnog", "deeptmhmm", "signalp", "psortb"):
        p.add_argument(f"--{tool}", metavar="PATH", help=f"M2: {tool} output file")
    # The two selectivity/essentiality components carry 0.40 of the weight mass (#81), so
    # a pipeline run must be able to reach them, and must say so when they are missing.
    p.add_argument("--full-evidence", action="store_true",
                   help="M2: turn on every component that needs external data "
                        "(--human-homology and --essentiality)")
    p.add_argument("--human-homology", action="store_true",
                   help="M2: compute human homologs with DIAMOND (human_homolog_penalty)")
    p.add_argument("--human-proteome", metavar="PATH", help="M2: reviewed human proteome FASTA")
    p.add_argument("--human-min-coverage", type=float,
                   help="M2: minimum query coverage for a human homolog")
    p.add_argument("--essentiality", action="store_true",
                   help="M2: transfer FBA essentiality from public relatives (essential)")
    p.add_argument("--essentiality-keyword",
                   help="M2: genus or species selecting the public relatives")
    p.add_argument("--essentiality-limit", type=int, help="M2: max reference genomes")
    p.add_argument("--essentiality-min-identity", type=float,
                   help="M2: minimum percent identity for an ortholog")
    p.add_argument("--essentiality-min-coverage", type=float,
                   help="M2: minimum coverage for an ortholog")
    p.add_argument("--offline", action="store_true",
                   help="replay from cache with no network (M2, and M1's API route)")
    p.add_argument("--allow-poor", action="store_true", help="passed to M1")
    p.add_argument("--dry-run", action="store_true",
                   help="print the command for each stage without running it")
    return p


def _taxon_id(run_dir: str) -> int | None:
    """M1 wrote the taxon; M2 needs it for id mapping and STRING. Read it rather than ask twice."""
    path = Path(run_dir) / "report.json"
    if not path.exists():
        return None
    try:
        return (json.loads(path.read_text()).get("genome") or {}).get("taxon_id")
    except (json.JSONDecodeError, OSError):
        return None


def _organism(run_dir: str) -> str | None:
    """M1 determined the organism; M2 needs it for the same-species flags. Read it, do not
    make M2 infer it from FASTA header brackets (issue #55)."""
    path = Path(run_dir) / "report.json"
    if not path.exists():
        return None
    try:
        genome = json.loads(path.read_text()).get("genome") or {}
    except (json.JSONDecodeError, OSError):
        return None
    return ((genome.get("taxonomy") or {}).get("scientific_name") or "").strip() or None


def stage_args(stage: str, args: argparse.Namespace) -> list[str]:
    common = ["--run", args.run]
    if stage == "m1":
        out = list(common)
        for flag, value in (("--contigs", args.contigs), ("--ws-dir", args.ws_dir),
                            ("--from-cga-dir", args.from_cga_dir),
                            ("--from-fasta", args.from_fasta),
                            ("--organism", args.organism),
                            ("--species", args.species), ("--genome-id", args.genome_id),
                            ("--annotation", args.annotation),
                            ("--max-features", args.max_features),
                            ("--seq-cap", args.seq_cap),
                            ("--tree-leaves", args.tree_leaves),
                            ("--taxon-id", args.taxon_id),
                            ("--genetic-code", args.genetic_code)):
            if value is not None:
                out += [flag, str(value)]
        for flag, on in (("--allow-poor", args.allow_poor), ("--html", args.html),
                         ("--figures", args.figures), ("--no-tree", args.no_tree)):
            if on:
                out.append(flag)
        if args.from_bvbrc_api:
            out.append("--from-bvbrc-api")
            # --offline only means anything on the API route; the CGA route's network is
            # the p3 CLI, which has no cache to replay from.
            if args.offline:
                out.append("--offline")
        return out
    if stage == "m2":
        # --report is not optional here: the point of the runner is one report.json
        # per run, written through the shared contract (00-architecture.md).
        out = common + ["--report"]
        if args.limit:
            out += ["--limit", str(args.limit)]
        for flag, on in (("--map-ids", args.map_ids), ("--kg", args.kg), ("--offline", args.offline),
                         ("--annotate", args.annotate),
                         ("--human-homology", args.human_homology or args.full_evidence),
                         ("--essentiality", args.essentiality or args.full_evidence)):
            if on:
                out.append(flag)
        for flag, value in (("--human-proteome", args.human_proteome),
                            ("--human-min-coverage", args.human_min_coverage),
                            ("--essentiality-keyword", args.essentiality_keyword),
                            ("--essentiality-limit", args.essentiality_limit),
                            ("--essentiality-min-identity", args.essentiality_min_identity),
                            ("--essentiality-min-coverage", args.essentiality_min_coverage)):
            if value is not None:
                out += [flag, str(value)]
        # A provider path implies --annotate: passing one and getting nothing would be a trap.
        for tool in ("interproscan", "eggnog", "deeptmhmm", "signalp", "psortb"):
            value = getattr(args, tool, None)
            if value:
                if "--annotate" not in out:
                    out.append("--annotate")
                out += [f"--{tool}", value]
        taxon = _taxon_id(args.run)
        if taxon:
            out += ["--taxon", str(taxon)]
        organism = _organism(args.run)
        if organism:
            out += ["--organism", organism]
        return out
    if stage == "m3":
        out = list(common)
        if args.structure_limit:
            out += ["--limit", str(args.structure_limit)]
        if args.offline:
            out.append("--offline")
        return out
    return common


#: M2 components that need external data, and the run.json section recording whether they ran.
EVIDENCE_COMPONENTS = (("human_homolog_penalty", "human_homology", "--human-homology"),
                       ("essential", "essentiality", "--essentiality"))


def missing_evidence(run_dir: str) -> list[str]:
    """One line per weighted M2 component that had no evidence behind it in this run.

    A ranking without them looks exactly as finished as one with them, so the runner says
    so instead of leaving it to `m2_pdb/run.json`.
    """
    path = Path(run_dir) / "m2_pdb" / "run.json"
    try:
        manifest = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return []
    weights = manifest.get("weights") or {}
    lines = []
    for component, section, flag in EVIDENCE_COMPONENTS:
        summary = manifest.get(section) or {}
        if summary.get("enabled") or not weights.get(component):
            continue
        reason = f"failed: {summary['error']}" if summary.get("error") else f"from BV-BRC's precomputed rows only (pass {flag})"
        lines.append(f"{component} (weight {weights[component]:+.2f}) {reason}")
    for component, entry in (manifest.get("components_available") or {}).items():
        if weights.get(component) and not entry.get("measured_for"):
            lines.append(f"{component} (weight {weights[component]:+.2f}) measured for 0 proteins")
    return lines


def stage_caps(args: argparse.Namespace, stages: list[str]) -> dict[str, dict]:
    """The caps the runner applied, per stage, for `run.caps` in report.json.

    A truncated run validates exactly like a complete one, so the caps have to be in the
    report, not only on the console. `None` means the stage ran uncapped or with its own
    default.
    """
    caps = {"m1": {"max_features": args.max_features, "seq_cap": args.seq_cap},
            "m2": {"limit": args.limit or None},
            "m3": {"limit": args.structure_limit or None}}
    return {stage: caps[stage] for stage in stages if stage in caps}


def _record_caps(args: argparse.Namespace, stages: list[str]) -> None:
    from .common.io import read_report, update_section

    run = dict(read_report(args.run).get("run") or {})
    if not run:
        return  # M1 creates the run section; nothing to annotate before it exists
    run["caps"] = {**(run.get("caps") or {}), **stage_caps(args, stages)}
    update_section(args.run, "run", run)


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    Path(args.run).mkdir(parents=True, exist_ok=True)

    selected = [(name, mod) for name, mod in STAGES
                if name not in args.skip and (not args.only or name in args.only)]
    if not selected:
        print("nothing to run", file=sys.stderr)
        return 2

    print(f"pipeline: {' -> '.join(name for name, _ in selected)}  (run {args.run})\n")
    for name, module in selected:
        cmd = [sys.executable, "-m", module, *stage_args(name, args)]
        print(f"=== {name}: {' '.join(cmd)}", flush=True)
        if args.dry_run:
            continue
        started = time.time()
        code = subprocess.call(cmd)
        elapsed = time.time() - started
        if code != 0:
            print(f"=== {name} failed (exit {code}) after {elapsed:.1f}s", file=sys.stderr)
            return code
        _record_caps(args, [name for name, _ in selected])
        print(f"=== {name} ok in {elapsed:.1f}s\n", flush=True)

    if not args.dry_run and any(name == "m2" for name, _ in selected):
        missing = missing_evidence(args.run)
        if missing:
            print("WARNING: the M2 ranking is missing weighted components:", file=sys.stderr)
            for line in missing:
                print(f"  {line}", file=sys.stderr)
    print(f"done. report: {Path(args.run) / 'report.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

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
]


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="python -m s2f.run", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--run", required=True, help="run directory, e.g. runs/mgen_G37")
    p.add_argument("--contigs", help="assembly FASTA for M1")
    p.add_argument("--ws-dir", help="BV-BRC workspace dir, e.g. /USER@bvbrc/home/s2f")
    p.add_argument("--from-cga-dir", help="M1 parses this directory instead of submitting a job")
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
    p.add_argument("--limit", type=int, help="passed to M2 as --limit")
    p.add_argument("--map-ids", action="store_true", help="M2: resolve UniProt/ChEMBL xrefs")
    p.add_argument("--kg", action="store_true", help="M2: assemble the knowledge subgraph")
    p.add_argument("--annotate", action="store_true",
                   help="M2: function, localization and membrane flags (docs/02c)")
    p.add_argument("--interproscan", metavar="PATH_OR_AUTO",
                   help="M2: InterProScan output to ingest, or 'auto' to find a local install")
    for tool in ("eggnog", "deeptmhmm", "signalp", "psortb"):
        p.add_argument(f"--{tool}", metavar="PATH", help=f"M2: {tool} output file")
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
                         ("--annotate", args.annotate)):
            if on:
                out.append(flag)
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
        if args.limit:
            out += ["--limit", str(args.limit)]
        if args.offline:
            out.append("--offline")
        return out
    return common


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
        print(f"=== {name} ok in {elapsed:.1f}s\n", flush=True)

    print(f"done. report: {Path(args.run) / 'report.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

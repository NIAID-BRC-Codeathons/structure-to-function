"""M1: assembly in, annotated genome out (issue #5).

    python -m s2f.m1_genome --run runs/<id> --contigs assembly.fna \\
        --ws-dir /USER@bvbrc/home/s2f

Order matters and is the point of this module: Minhash predicts the taxon *before*
CGA runs, because CGA's annotation quality depends on it and CGA cannot infer it.
Annotating this project's test genome with no taxon doubles the CDS count and halves
mean protein length (`docs/01a-cga-coverage.md`).

Offline re-parse of a job already retrieved:

    python -m s2f.m1_genome --run runs/<id> --from-cga-dir data/cga
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from ..common.io import init_report, update_section
from . import cga
from .parse import load_cga, genome_section, proteins_section, run_fields, write_m1_dir
from .taxon import TaxonCallError, call_taxon


def _now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="python -m s2f.m1_genome", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--run", required=True, help="run directory, e.g. runs/mgen_G37")
    p.add_argument("--contigs", help="assembly FASTA (required unless --from-cga-dir)")
    p.add_argument("--ws-dir", help="BV-BRC workspace dir, e.g. /USER@bvbrc/home/s2f")
    p.add_argument("--out-name", default=None, help="CGA output name (default: the run id)")
    p.add_argument("--from-cga-dir", help="parse an already-retrieved CGA directory; no network")
    p.add_argument("--job-id", help="resume: skip submission and poll this job id")
    p.add_argument("--taxon-id", type=int, help="skip the Minhash call and use this taxon")
    p.add_argument("--genetic-code", type=int, help="translation table; required with --taxon-id")
    p.add_argument("--scientific-name", help="overrides the name from the taxon call")
    p.add_argument("--max-distance", type=float, default=0.2, help="Minhash cut-off (default 0.2)")
    p.add_argument("--max-hits", type=int, default=50)
    p.add_argument("--reference-scope", action="store_true",
                   help="search only reference/representative sketches (usually too narrow)")
    p.add_argument("--no-blind", action="store_true", help="keep original FASTA headers")
    p.add_argument("--poll-interval", type=int, default=30)
    p.add_argument("--poll-timeout", type=int, default=7200)
    p.add_argument("--allow-poor", action="store_true",
                   help="continue even if CGA reports Poor quality or sets a quality flag")
    p.add_argument("--dry-run", action="store_true", help="print the CGA submission and stop")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    run_dir = Path(args.run)
    run_id = run_dir.name
    m1_dir = run_dir / "m1"
    m1_dir.mkdir(parents=True, exist_ok=True)
    started = time.time()

    init_report(run_dir, run_id, input_file=args.contigs, modules={})
    taxon_call: dict | None = None
    cga_dir = args.from_cga_dir

    if cga_dir is None:
        if not (args.contigs and args.ws_dir):
            print("--contigs and --ws-dir are required unless --from-cga-dir is given",
                  file=sys.stderr)
            return 2

        contigs = Path(args.contigs)
        if not args.no_blind:
            # Named for the run: the workspace upload keeps this basename, and p3-cp
            # will not overwrite, so a shared name would hand Minhash a stale genome.
            contigs = cga.blind_contigs(contigs, m1_dir / f"{run_id}.contigs.fna")
            print(f"blinded contigs -> {contigs}")

        if args.taxon_id:
            if not args.genetic_code:
                print("--genetic-code is required with --taxon-id", file=sys.stderr)
                return 2
            taxon_call = {"called_by": "supplied", "taxon_id": args.taxon_id,
                          "genetic_code": args.genetic_code,
                          "taxon_name": args.scientific_name or str(args.taxon_id),
                          "called_rank": None, "hits": []}
        else:
            ws_fasta = cga.upload(contigs, args.ws_dir)
            print(f"uploaded -> {ws_fasta}\ncalling taxon via Minhash ...")
            try:
                taxon_call = call_taxon(ws_fasta, max_distance=args.max_distance,
                                        max_hits=args.max_hits,
                                        all_public=not args.reference_scope)
            except TaxonCallError as exc:
                print(f"taxon call failed: {exc}", file=sys.stderr)
                return 3
            print(f"  {taxon_call['n_hits']} hits; top {taxon_call['top_hit']['genome_id']} "
                  f"at d={taxon_call['top_hit']['distance']:.5f}")
            print(f"  called {taxon_call['called_rank']}: {taxon_call['taxon_name']} "
                  f"(taxon {taxon_call['taxon_id']}), genetic code {taxon_call['genetic_code']}")
        (m1_dir / "taxon_call.json").write_text(json.dumps(taxon_call, indent=1))

        out_name = args.out_name or run_id
        job_id = args.job_id
        if job_id is None:
            job_id = cga.submit(
                contigs, ws_dir=args.ws_dir, out_name=out_name,
                scientific_name=args.scientific_name or taxon_call["taxon_name"],
                taxon_id=taxon_call["taxon_id"], genetic_code=taxon_call["genetic_code"],
                label=run_id, dry_run=args.dry_run,
            )
            if args.dry_run:
                return 0
            print(f"submitted CGA job {job_id}")
        print("polling ...")
        print(f"  final status: {cga.wait(job_id, interval=args.poll_interval, timeout=args.poll_timeout)}")
        cga_dir = str(run_dir / "cga")
        cga.fetch(args.ws_dir, out_name, cga_dir)
        print(f"retrieved -> {cga_dir}")

    parsed = load_cga(cga_dir)
    proteins = proteins_section(parsed)
    counts = write_m1_dir(parsed, m1_dir)
    update_section(run_dir, "genome", genome_section(parsed, taxon_call))
    update_section(run_dir, "proteins", proteins)

    report = init_report(run_dir, run_id)
    run_section = dict(report.get("run") or {})
    run_section.update(run_fields(parsed))
    if taxon_call:
        run_section["taxon_call"] = {k: taxon_call.get(k) for k in
                                     ("called_by", "called_rank", "taxon_id", "taxon_name",
                                      "genetic_code", "n_hits")}
    modules = dict(run_section.get("modules") or {})
    modules["m1"] = {"status": "ok", "finished_at": _now(),
                     "elapsed_seconds": round(time.time() - started, 1), **counts}
    run_section["modules"] = modules
    update_section(run_dir, "run", run_section)

    print(f"\nwrote {len(proteins)} proteins, {counts['specialty_rows']} specialty rows -> {m1_dir}")
    print(f"genome {parsed.genome_id}, quality {parsed.quality or 'unknown'}")

    if parsed.is_poor:
        print(f"\nQUALITY GATE: {parsed.quality_reason()}", file=sys.stderr)
        print("A Poor annotation usually means the taxon or genetic code was wrong; "
              "the CDS count inflates and proteins fragment (docs/01a-cga-coverage.md).",
              file=sys.stderr)
        if not args.allow_poor:
            print("Refusing to hand this downstream. Re-run with --allow-poor to override.",
                  file=sys.stderr)
            return 4
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

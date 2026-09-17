"""Collect existing structures for proteins selected by M2.

    python -m s2f.m3_fold --run runs/<run_id> [--offline] [--dry-run]

This first issue-13 slice collects PDB/AlphaFold DB structures. Prediction and confidence
gating are intentionally not part of this command yet.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from ..common.http import CachedJsonClient
from ..common.io import read_report, report_path, update_section
from .collect import StructureCandidate, collect_existing

REPO_ROOT = Path(__file__).resolve().parents[2]
REPORT_FIXTURE = REPO_ROOT / "fixtures" / "report.fixture.json"
STRUCTURE_FIXTURES = REPO_ROOT / "fixtures" / "m3"


def _seed_dry_run(run_dir: Path) -> None:
    if report_path(run_dir).exists():
        return
    fixture = json.loads(REPORT_FIXTURE.read_text(encoding="utf-8"))
    for key, value in fixture.items():
        if key != "schema_version":
            update_section(run_dir, key, value)


def _fixture_fetch(candidate: StructureCandidate) -> bytes:
    name = candidate.accession.split("_", 1)[0] if candidate.source == "pdb" else candidate.accession
    path = STRUCTURE_FIXTURES / f"{name}.cif"
    if not path.exists():
        raise OSError(f"no offline structure fixture for {candidate.accession}")
    return path.read_bytes()


def run(args: argparse.Namespace) -> int:
    run_dir = Path(args.run)
    if args.dry_run:
        _seed_dry_run(run_dir)

    report = read_report(run_dir)
    if "proteins" not in report:
        raise ValueError(f"{report_path(run_dir)} has no proteins section")

    if args.dry_run:
        fetch = _fixture_fetch
    else:
        client = CachedJsonClient(offline=args.offline)
        cache_dir = Path(args.cache) if args.cache else run_dir / "cache" / "files"

        def fetch(candidate: StructureCandidate) -> bytes:
            return client.get_bytes("m3_structure", candidate.url, cache_dir=cache_dir)

    summary = collect_existing(report, run_dir, fetch, limit=args.limit)
    update_section(run_dir, "structures", summary.structures)

    print(
        f"{summary.selected} selected proteins: {summary.collected} existing structures collected, "
        f"{summary.prediction_required} require prediction, {summary.failed} failed."
    )
    return 0 if summary.failed == 0 else 1


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="python -m s2f.m3_fold", description=__doc__)
    parser.add_argument("--run", default="runs/dev", help="run directory containing report.json")
    parser.add_argument("--limit", type=int, default=0, help="process only the first N selected proteins")
    parser.add_argument("--cache", default="", help="override the downloaded-file cache directory")
    parser.add_argument("--offline", action="store_true", help="use cached files; never call the network")
    parser.add_argument("--dry-run", action="store_true", help="use committed report and structure fixtures")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.dry_run and args.run == "runs/dev":
        args.run = "runs/dry-run"
    try:
        return run(args)
    except (OSError, ValueError) as exc:
        print(f"Error: {exc}")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

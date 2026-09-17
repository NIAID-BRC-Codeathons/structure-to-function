"""Collect existing structures for proteins selected by M2.

    python -m s2f.m3_fold --run runs/<run_id> [--offline] [--dry-run]

This issue-13 slice collects PDB/AlphaFold DB structures and applies a provisional,
configurable quality gate. Prediction is intentionally not part of this command yet.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import sys
import time
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path

from ..common.http import CachedJsonClient, DownloadedFile
from ..common.io import read_report, report_path, update_section
from .collect import StructureCandidate, collect_existing
from .quality import POLICY_VERSION, QualityThresholds

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


def _fixture_fetch(candidate: StructureCandidate) -> DownloadedFile:
    name = candidate.accession.split("_", 1)[0] if candidate.source == "pdb" else candidate.accession
    path = STRUCTURE_FIXTURES / f"{name}.cif"
    if not path.exists():
        raise OSError(f"no offline structure fixture for {candidate.accession}")
    retrieved_at = datetime.fromtimestamp(path.stat().st_mtime, UTC)
    return DownloadedFile(path.read_bytes(), retrieved_at, from_cache=False)


def run(args: argparse.Namespace) -> int:
    started = datetime.now(UTC)
    clock = time.monotonic()
    run_dir = Path(args.run)
    if args.dry_run:
        _seed_dry_run(run_dir)

    input_report = report_path(run_dir)
    report = read_report(run_dir)
    if "proteins" not in report:
        raise ValueError(f"{input_report} has no proteins section")
    input_report_sha256 = hashlib.sha256(input_report.read_bytes()).hexdigest()

    cache_dir: Path | None = None
    if args.dry_run:
        fetch = _fixture_fetch
    else:
        client = CachedJsonClient(offline=args.offline)
        cache_dir = Path(args.cache) if args.cache else run_dir / "cache" / "files"

        def fetch(candidate: StructureCandidate) -> DownloadedFile:
            assert cache_dir is not None
            return client.get_file("m3_structure", candidate.url, cache_dir=cache_dir)

    quality_thresholds = QualityThresholds(
        min_pdb_identity=args.min_pdb_identity,
        min_pdb_coverage=args.min_pdb_coverage,
        max_pdb_resolution=args.max_pdb_resolution,
        min_mean_plddt=args.min_mean_plddt,
        min_afdb_coverage=args.min_afdb_coverage,
        min_local_plddt=args.min_local_plddt,
    )
    summary = collect_existing(
        report, run_dir, fetch, limit=args.limit, quality_thresholds=quality_thresholds
    )
    update_section(run_dir, "structures", summary.structures)
    _write_manifest(
        run_dir,
        args,
        summary,
        started=started,
        elapsed_seconds=time.monotonic() - clock,
        input_report_sha256=input_report_sha256,
        cache_dir=cache_dir,
        quality_thresholds=quality_thresholds,
    )

    print(
        f"{summary.selected} selected proteins: {summary.collected} existing structures collected, "
        f"{summary.prediction_required} require prediction, {summary.failed} failed; "
        f"{sum(record.get('usable_for_docking') is True for record in summary.structures)} "
        "passed the provisional quality gate."
    )
    return 0 if summary.failed == 0 else 1


def _write_manifest(
    run_dir: Path,
    args: argparse.Namespace,
    summary,
    *,
    started: datetime,
    elapsed_seconds: float,
    input_report_sha256: str,
    cache_dir: Path | None,
    quality_thresholds: QualityThresholds,
) -> None:
    finished = datetime.now(UTC)
    source_counts = Counter(record["source"] for record in summary.structures)
    manifest = {
        "module": "m3_fold",
        "docs": "docs/03-m3-fold.md",
        "issue": 13,
        "started_utc": started.isoformat(),
        "finished_utc": finished.isoformat(),
        "elapsed_seconds": round(elapsed_seconds, 6),
        "command": getattr(args, "command", None),
        "input_report": str(report_path(run_dir)),
        "input_report_sha256": input_report_sha256,
        "output_report_sha256": hashlib.sha256(report_path(run_dir).read_bytes()).hexdigest(),
        "output_dir": str(run_dir / "structures"),
        "cache_dir": str(cache_dir) if cache_dir is not None else None,
        "offline": args.offline,
        "dry_run": args.dry_run,
        "parameters": {
            "limit": args.limit,
            "selection_policy": "experimental_homolog_then_afdb",
            "confidence_gate_applied": True,
            "quality_gate": {
                "policy_version": POLICY_VERSION,
                "thresholds": quality_thresholds.as_dict(),
            },
        },
        "environment": {
            "python": platform.python_version(),
            "python_executable": sys.executable,
            "platform": platform.platform(),
        },
        "git": _git_state(),
        "counts": {
            "selected": summary.selected,
            "collected": summary.collected,
            "prediction_required": summary.prediction_required,
            "failed": summary.failed,
            "cache_hits": sum(record.get("cache_hit") is True for record in summary.structures),
            "sources": dict(sorted(source_counts.items())),
            "usable_for_docking": sum(
                record.get("usable_for_docking") is True for record in summary.structures
            ),
            "not_usable_for_docking": sum(
                record.get("usable_for_docking") is False for record in summary.structures
            ),
            "quality_pending": sum(
                record.get("usable_for_docking") is None for record in summary.structures
            ),
        },
        "failures": [
            {"feature_id": record["feature_id"], "reason": record.get("reason")}
            for record in summary.structures
            if record.get("collection_status") == "failed"
        ],
    }
    out = run_dir / "m3_fold" / "run.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def _git_state() -> dict[str, object]:
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    status = subprocess.run(
        ["git", "status", "--porcelain", "--untracked-files=no"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    return {
        "commit": commit.stdout.strip() or None,
        "dirty": bool(status.stdout.strip()),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="python -m s2f.m3_fold", description=__doc__)
    parser.add_argument("--run", default="runs/dev", help="run directory containing report.json")
    parser.add_argument("--limit", type=int, default=0, help="process only the first N selected proteins")
    parser.add_argument("--cache", default="", help="override the downloaded-file cache directory")
    parser.add_argument("--offline", action="store_true", help="use cached files; never call the network")
    parser.add_argument("--dry-run", action="store_true", help="use committed report and structure fixtures")
    gate = parser.add_argument_group("provisional quality gate")
    gate.add_argument(
        "--min-pdb-identity",
        type=float,
        default=0.25,
        help="minimum PDB-homolog sequence identity (default: 0.25)",
    )
    gate.add_argument(
        "--min-pdb-coverage",
        type=float,
        default=0.50,
        help="minimum PDB-homolog query coverage (default: 0.50)",
    )
    gate.add_argument(
        "--max-pdb-resolution",
        type=float,
        default=3.50,
        help="maximum X-ray/EM resolution in angstroms (default: 3.50)",
    )
    gate.add_argument(
        "--min-mean-plddt",
        type=float,
        default=70.0,
        help="minimum AFDB mean pLDDT (default: 70)",
    )
    gate.add_argument(
        "--min-afdb-coverage",
        type=float,
        default=0.80,
        help="minimum AFDB sequence coverage (default: 0.80)",
    )
    gate.add_argument(
        "--min-local-plddt",
        type=float,
        default=70.0,
        help="mask AFDB residues below this pLDDT (default: 70)",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    raw_argv = list(sys.argv[1:] if argv is None else argv)
    args = build_parser().parse_args(raw_argv)
    args.command = [sys.executable, "-m", "s2f.m3_fold", *raw_argv]
    if args.dry_run and args.run == "runs/dev":
        args.run = "runs/dry-run"
    try:
        return run(args)
    except (OSError, ValueError) as exc:
        print(f"Error: {exc}")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

"""The full report — every section of a run, as one self-contained HTML page.

    python -m s2f.full_report --run runs/<id>

Runs last and reads the whole run directory. Owns `full_report` in the contract and
never writes another module's key; `s2f.m6_report` keeps `report` and its own compact
rendering. A section that is missing produces a tab saying so rather than a tab that
quietly disappears.
"""

from __future__ import annotations

import argparse
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from ..common.io import read_report, report_path, update_section
from .render import TABS, build_page, render_run

REPO_ROOT = Path(__file__).resolve().parents[2]
REPORT_FIXTURE = REPO_ROOT / "fixtures" / "report.fixture.json"


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m s2f.full_report",
                                     description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run", required=True, help="run directory holding report.json")
    parser.add_argument("--out", default="",
                        help="output path (default: <run>/full_report.html)")
    parser.add_argument("--dry-run", action="store_true",
                        help="render the committed fixture instead of a run; no network, "
                             "no run directory needed")
    args = parser.parse_args(argv)

    started = time.time()
    run_dir = Path(args.run)

    if args.dry_run:
        # The contract asks every module for a --dry-run that needs no network. M6 needs
        # none at any time, so this only proves the renderer survives the fixture.
        import json as _json
        fixture = _json.loads(REPORT_FIXTURE.read_text(encoding="utf-8"))
        page = build_page(fixture, {})
        target = Path(args.out) if args.out else run_dir / "full_report.html"
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(page, encoding="utf-8")
        print(f"dry run: fixture rendered -> {target} ({len(page) / 1024:.0f} KB)")
        return 0
    if not report_path(run_dir).exists():
        print(f"no report.json in {run_dir}", file=sys.stderr)
        return 2

    page = render_run(run_dir)
    target = Path(args.out) if args.out else run_dir / "full_report.html"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(page, encoding="utf-8")

    report = read_report(run_dir)
    sections = [key for key in ("run", "genome", "proteins", "kg", "structures",
                                "ligands", "docking", "disease") if report.get(key)]
    counts = {
        "proteins": len(report.get("proteins") or []),
        "selected": sum(1 for protein in (report.get("proteins") or [])
                        if (protein.get("triage") or {}).get("selected")),
        "structures": len(report.get("structures") or []),
        "ligands": len(report.get("ligands") or []),
    }
    update_section(run_dir, "full_report", {
        "rendered_path": str(target),
        "rendered_at": _now(),
        "bytes": len(page.encode("utf-8")),
        "sections_present": sections,
        "sections_missing": [key for key in ("structures", "ligands", "docking",
                                             "disease", "kg") if not report.get(key)],
        "tabs": [label for _key, label, _empty in TABS],
        "counts": counts,
        "elapsed_seconds": round(time.time() - started, 2),
    })
    print(f"full report -> {target} ({len(page) / 1024:.0f} KB)")
    print(f"  sections present: {', '.join(sections) or 'none'}")
    print(f"  {counts['proteins']} proteins, {counts['selected']} selected, "
          f"{counts['structures']} structures, {counts['ligands']} ligands")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

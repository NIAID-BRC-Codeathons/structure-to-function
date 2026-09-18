"""Render `report.json` into one self-contained HTML file (issue #66).

    python -m s2f.m6_report --run runs/<run_id>
    python -m s2f.m6_report --dry-run          # from fixtures/report.fixture.json

Writes `<run>/report.html` and records the render in the report's own `report` section, which
is the one key M6 owns (`docs/00-architecture.md`).

The file has no external references — no CDN, no fonts, no images — because it has to open
from disk with the network off, which is also how the demo runs. 3Dmol.js views and RDKit
ligand drawings from the original scope need external assets or a build step, so they are
deferred rather than bolted on; the sections they would fill say "not run" instead.
"""

from __future__ import annotations

import argparse
import json
from datetime import UTC, datetime
from pathlib import Path

from jinja2 import Environment, FileSystemLoader, StrictUndefined

from ..common.io import read_report, update_section
from . import view as view_mod

TEMPLATE_DIR = Path(__file__).resolve().parent / "templates"
FIXTURE = Path(__file__).resolve().parents[2] / "fixtures" / "report.fixture.json"


def render(report: dict, *, top: int = 25) -> str:
    """Report dict -> HTML string. StrictUndefined so a missing field fails loudly here."""
    environment = Environment(
        loader=FileSystemLoader(TEMPLATE_DIR),
        autoescape=True,
        undefined=StrictUndefined,
        trim_blocks=True,
        lstrip_blocks=True,
    )
    template = environment.get_template("report.html.j2")
    return template.render(view=view_mod.build(report, top=top))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m s2f.m6_report", description=__doc__)
    parser.add_argument("--run", default="runs/dev", help="run directory holding report.json")
    parser.add_argument("--out", default="", help="output path (default <run>/report.html)")
    parser.add_argument("--top", type=int, default=25, help="how many selected proteins to show")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="render the committed fixture instead of a run, proving it works before any module has",
    )
    args = parser.parse_args(argv)

    run_dir = Path(args.run)
    if args.dry_run:
        report = json.loads(FIXTURE.read_text(encoding="utf-8"))
    else:
        report = read_report(run_dir)
        if not any(key in report for key in view_mod.SECTIONS):
            print(f"Error: {run_dir / 'report.json'} has no sections to render.")
            return 2

    html = render(report, top=args.top)
    out_path = Path(args.out) if args.out else run_dir / "report.html"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(html, encoding="utf-8")

    built = view_mod.build(report, top=args.top)
    status = {section.key: section.status for section in built.sections}
    if not args.dry_run:
        # The one key M6 owns. Written through the shared writer like every other section.
        update_section(
            run_dir,
            "report",
            {
                "rendered_path": str(out_path),
                "rendered_at": datetime.now(UTC).isoformat(timespec="seconds"),
                "counts": built.counts,
                "sections": status,
            },
        )

    ran = [key for key, value in status.items() if value == "ran"]
    missing = [key for key, value in status.items() if value != "ran"]
    print(f"Wrote {out_path} ({out_path.stat().st_size / 1024:.0f} KB)")
    print(f"  sections rendered: {', '.join(ran) or 'none'}")
    if missing:
        print(f"  shown as not run:  {', '.join(missing)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

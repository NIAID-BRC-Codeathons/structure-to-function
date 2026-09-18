"""Full report tests.

The acceptance checks are the same ones `docs/06-m6-report.md` sets for any report: it
renders before the modules that feed it exist, it renders from a partial file, it opens
with no network, and a section that did not run says so instead of looking empty.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

from s2f.full_report.__main__ import REPORT_FIXTURE, main
from s2f.full_report.render import TABS, build_page

FIXTURE_REPORT = json.loads(REPORT_FIXTURE.read_text(encoding="utf-8"))


def test_renders_from_the_committed_fixture() -> None:
    """Day-one check: the report exists before any other module produces data."""
    html = build_page(FIXTURE_REPORT, {})

    assert html.startswith("<!doctype html>")
    assert "</html>" in html.strip()[-16:]


def test_renders_from_a_partial_report() -> None:
    html = build_page({"run": {"run_id": "partial"}, "genome": {"genome_id": "1.1"}}, {})

    assert "partial" in html
    assert "has not run" in html


def test_an_empty_report_still_renders() -> None:
    html = build_page({}, {})

    assert html.startswith("<!doctype html>")
    for _key, label, _empty in TABS:
        assert label in html


def test_every_tab_is_present_even_when_its_section_is_absent() -> None:
    """A missing section and an empty one are different facts."""
    html = build_page({}, {})

    for key, _label, empty_text in TABS:
        assert f"id='tab-{key}'" in html
        if empty_text != "Nothing recorded yet.":
            assert empty_text in html


def test_tabs_are_named_for_content_not_for_modules() -> None:
    labels = [label for _key, label, _empty in TABS]

    assert not [label for label in labels if re.fullmatch(r"M\d.*", label)]
    assert "Genome Overview" in labels
    assert "Candidate targets" in labels


def test_the_file_has_no_external_references() -> None:
    """It has to open offline and by itself: one file, nothing fetched to render it."""
    html = build_page(FIXTURE_REPORT, {})

    assert not re.search(r"<script[^>]+src=", html)
    assert not re.search(r"<link[^>]*>", html)
    assert not re.search(r"<img[^>]*>", html)
    assert "url(" not in html
    assert "@import" not in html
    for token in ("fetch(", "XMLHttpRequest", "localStorage"):
        assert token not in html


def test_borrowed_metadata_names_the_genome_it_came_from() -> None:
    """'This isolate came from a human' and 'a relative did' are different claims."""
    report = {
        "run": {"run_id": "borrowed"},
        "genome": {
            "genome_id": "2097.118",
            "metadata_provenance": {
                "basis": "relative", "genome_id": "243273.25",
                "genome_name": "Mycoplasma genitalium G37",
                "mash_distance": 0.0, "note": "borrowed from the closest public match",
            },
        },
    }
    html = build_page(report, {})

    assert "243273.25" in html
    assert "Mycoplasma genitalium G37" in html
    assert "not this assembly" in html


def test_the_limitations_are_always_present() -> None:
    html = build_page({}, {})

    assert "Research use only" in html
    assert "not comparable across proteins" in html


def test_cli_writes_the_file_and_records_its_own_section(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    (run_dir / "report.json").write_text(json.dumps(FIXTURE_REPORT), encoding="utf-8")

    assert main(["--run", str(run_dir)]) == 0

    target = run_dir / "full_report.html"
    assert target.is_file()
    written = json.loads((run_dir / "report.json").read_text(encoding="utf-8"))
    assert written["full_report"]["rendered_path"] == str(target)
    assert written["full_report"]["bytes"] > 0
    # It owns `full_report` and never another module's key.
    assert "report" not in written or written.get("report") != written["full_report"]


def test_cli_dry_run_needs_no_run_directory(tmp_path: Path) -> None:
    target = tmp_path / "out.html"

    assert main(["--run", str(tmp_path / "absent"), "--dry-run", "--out", str(target)]) == 0
    assert target.is_file()


def test_cli_refuses_a_run_directory_with_no_report(tmp_path: Path) -> None:
    assert main(["--run", str(tmp_path)]) == 2

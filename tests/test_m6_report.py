"""M6 report renderer tests (issue #66).

The acceptance checks in `docs/06-m6-report.md` are the tests: renders from the fixture before
any module exists, renders from a partial file, nothing computed in the template, opens offline.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from s2f.common.io import read_report
from s2f.m6_report import view as view_mod
from s2f.m6_report.__main__ import FIXTURE, main, render

FIXTURE_REPORT = json.loads(FIXTURE.read_text(encoding="utf-8"))


def test_renders_from_the_committed_fixture() -> None:
    """Day-one check: the report exists before any other module produces data."""
    html = render(FIXTURE_REPORT)

    assert html.startswith("<!DOCTYPE html>")
    assert "Structure to Function" in html
    assert "Mycoplasmoides genitalium" in html


def test_renders_from_a_partial_report() -> None:
    """M3-M5 missing must render as 'not run', not as an error and not as zero."""
    partial = {k: FIXTURE_REPORT[k] for k in ("schema_version", "run", "genome", "proteins")}

    html = render(partial)

    assert "not run" in html
    # The absent modules say so in words, so a reader cannot mistake it for "none found".
    assert "M4 did not contribute" in html
    assert "M5 did not contribute" in html


def test_an_empty_report_still_renders() -> None:
    html = render({"schema_version": "0.1.0"})

    assert "<html" in html
    assert "not determined" in html  # organism, rather than a crash or a blank


def test_the_file_has_no_external_references() -> None:
    """It has to open from disk with the network off — same rule as the slide deck."""
    html = render(FIXTURE_REPORT)

    assert "http://" not in html
    assert "https://" not in html
    assert "<script" not in html.lower()
    assert re.search(r"<(img|link)\b", html, re.IGNORECASE) is None


def test_the_limitations_block_is_always_present() -> None:
    """Fixed text, not conditional on which modules ran (pitfalls #1, #5, #6, #9)."""
    for report in (FIXTURE_REPORT, {"schema_version": "0.1.0"}):
        html = render(report)
        assert "not binding affinities" in html
        assert "no cofactors" in html
        assert "Nothing here is experimentally validated" in html
        assert "No clinical or treatment implication" in html


def test_ligands_are_grouped_by_protein_and_never_ranked_across_them() -> None:
    """Pitfall #1: the single most likely way this produces a confidently wrong headline."""
    html = render(FIXTURE_REPORT)

    assert "Grouped by protein" in html
    assert "never comparable across proteins" in html
    built = view_mod.build(FIXTURE_REPORT)
    # Every candidate belongs to exactly one protein group; there is no global ligand table.
    assert built.candidates_by_protein
    assert all(group["feature_id"] for group in built.candidates_by_protein)


def test_nothing_is_computed_in_the_template() -> None:
    """Formatting and arithmetic live in view.py, so the template prints prepared strings."""
    source = (Path(view_mod.__file__).parent / "templates" / "report.html.j2").read_text()
    body = "\n".join(line for line in source.splitlines() if not line.strip().startswith("{#"))

    for forbidden in ("|sum", "|round", "|length >", "* 100", "/ ", " + "):
        assert forbidden not in body, f"template does arithmetic: {forbidden}"


def test_a_taxonomy_string_renders_like_a_taxonomy_object() -> None:
    """The schema allows both; M1 writes a dict, the fixture carries a string."""
    as_dict = view_mod.build(
        {"genome": {"taxonomy": {"scientific_name": "Escherichia coli"}}}
    )
    as_string = view_mod.build({"genome": {"taxonomy": "Escherichia coli"}})

    assert as_dict.organism == as_string.organism == "Escherichia coli"


def test_missing_sections_are_reported_as_not_run() -> None:
    built = view_mod.build({"schema_version": "0.1.0", "proteins": []})
    status = {section.key: section.status for section in built.sections}

    assert set(status) == set(view_mod.SECTIONS)
    assert all(value == "not run" for value in status.values())
    assert any("not the same as empty" in caveat for caveat in built.caveats)


def test_a_same_species_run_says_the_annotation_was_recovered_not_discovered() -> None:
    report = {
        "proteins": [
            {"feature_id": "fig|1.1.peg.1", "flags": {"same_species_hit": True},
             "triage": {"selected": True, "rank": 1, "score": 0.9, "components": {}, "reason": "x"}}
        ]
    }
    built = view_mod.build(report)

    assert built.counts["same_species_hits"] == 1
    assert any("recovered, not discovered" in caveat for caveat in built.caveats)


def test_cli_writes_the_file_and_records_the_report_section(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    (run_dir / "report.json").write_text(json.dumps(FIXTURE_REPORT), encoding="utf-8")

    assert main(["--run", str(run_dir)]) == 0

    out = run_dir / "report.html"
    assert out.exists() and out.stat().st_size > 2000
    # M6 owns exactly one key, written through the shared section writer.
    section = read_report(run_dir)["report"]
    assert section["rendered_path"].endswith("report.html")
    assert section["rendered_at"]
    assert section["sections"]["proteins"] == "ran"


def test_cli_dry_run_needs_no_run_directory(tmp_path: Path) -> None:
    out = tmp_path / "fixture.html"

    assert main(["--dry-run", "--out", str(out)]) == 0
    assert out.exists()
    assert not (tmp_path / "report.json").exists()  # a dry run writes no report section


def test_cli_refuses_an_empty_run_directory(tmp_path: Path) -> None:
    assert main(["--run", str(tmp_path / "nothing")]) == 2

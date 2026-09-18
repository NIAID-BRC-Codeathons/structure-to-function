"""The evaluation CLI, end to end, from committed fixtures (issue #49).

`00-architecture.md` requires `--dry-run` to work against `fixtures/` with no network and
no credentials. Nothing here touches either.
"""

from __future__ import annotations

import json

import pytest

from s2f.common.io import read_report
from s2f.common.schema import validate_report
from s2f.evaluation.__main__ import REPORT_FIXTURE, build_parser, main

OUTPUTS = ("metrics.json", "truth_matches.tsv", "summary.md", "run.json")


@pytest.fixture
def run_dir(tmp_path):
    return tmp_path / "eval-run"


def run_cli(run_dir, *extra: str) -> int:
    return main(["--dry-run", "--run", str(run_dir), *extra])


def test_dry_run_succeeds_and_writes_every_output(run_dir, capsys):
    assert run_cli(run_dir) == 0
    for name in OUTPUTS:
        assert (run_dir / "eval" / name).is_file(), name
    out = capsys.readouterr().out
    assert "precision" in out and "recall" in out


def test_the_fixture_is_committed_and_is_a_valid_report():
    """The dry-run input has to satisfy the same schema every other module writes."""
    assert REPORT_FIXTURE.is_file()
    validate_report(json.loads(REPORT_FIXTURE.read_text(encoding="utf-8")))


def test_it_never_writes_a_report_section(run_dir):
    """The contract gives every top-level key an owner and evaluation is not one of them.

    `common.schema.validate` rejects an unknown section, so inventing an `evaluation` key
    would be the improvised schema change `00-architecture.md` forbids. The module reads
    the report and writes its own directory; this is the test that keeps it that way.
    """
    assert run_cli(run_dir) == 0
    before = json.loads(REPORT_FIXTURE.read_text(encoding="utf-8"))
    after = read_report(run_dir)
    assert set(after) == set(before)
    validate_report(after)


def test_metrics_json_carries_the_provenance_a_rerun_needs(run_dir):
    assert run_cli(run_dir) == 0
    metrics = json.loads((run_dir / "eval" / "metrics.json").read_text(encoding="utf-8"))
    assert metrics["truthset"] == "kpneumoniae_hs11286"
    assert len(metrics["weights_fingerprint"]) == 12
    assert metrics["weights"] and metrics["measured"]["k"] > 0
    assert metrics["positives_in_genome"] > 0
    assert set(metrics["baselines"]) == {"keyword", "length", "random"}


def test_manifest_records_the_truthset_hash_and_the_scorer(run_dir):
    assert run_cli(run_dir) == 0
    manifest = json.loads((run_dir / "eval" / "run.json").read_text(encoding="utf-8"))
    assert manifest["issue"] == 49
    assert manifest["dry_run"] is True
    assert len(manifest["truthset"]["sha256"]) == 64
    assert manifest["scorer"]["module"] == "s2f.m1_genome.priority"
    assert manifest["parameters"]["top_k"] > 0


def test_matches_tsv_has_one_row_per_matched_protein(run_dir):
    assert run_cli(run_dir) == 0
    lines = (run_dir / "eval" / "truth_matches.tsv").read_text(encoding="utf-8").splitlines()
    metrics = json.loads((run_dir / "eval" / "metrics.json").read_text(encoding="utf-8"))
    assert len(lines) - 1 == metrics["matched"]
    assert lines[0].split("\t")[:4] == ["rank", "score", "feature_id", "symbol"]


def test_the_fixture_exercises_both_failure_modes_the_issue_reports(run_dir):
    """The fixture is only useful if the known-bad proteins actually score."""
    assert run_cli(run_dir) == 0
    metrics = json.loads((run_dir / "eval" / "metrics.json").read_text(encoding="utf-8"))
    by_symbol = {row["symbol"]: row for row in metrics["rows"]}
    # envelope vocabulary scoring as secretion, and a plasmid addiction toxin as a toxin
    assert by_symbol["skp"]["label"] == "negative" and by_symbol["skp"]["score"] > 0
    assert by_symbol["pemK"]["label"] == "negative" and by_symbol["pemK"]["score"] > 0
    # and the ceiling tie the issue reports
    assert metrics["ties"]["at_max"] > 1


def test_tonb_dependent_receptor_is_not_matched_to_tonb(run_dir):
    """End-to-end guard for the #49 regex trap, with the weak tier switched on."""
    assert run_cli(run_dir, "--include-product-matches") == 0
    metrics = json.loads((run_dir / "eval" / "metrics.json").read_text(encoding="utf-8"))
    tonb = [row for row in metrics["rows"] if row["symbol"] == "tonB"]
    assert [row["gene"] for row in tonb] == ["tonB"]
    assert all("dependent" not in (row["product"] or "") for row in tonb)


def test_summary_states_the_precision_denominator(run_dir):
    """A precision figure without its denominator is not a measurement."""
    assert run_cli(run_dir) == 0
    summary = (run_dir / "eval" / "summary.md").read_text(encoding="utf-8")
    assert "labelled" in summary
    assert "pitfall #12" in summary
    assert "Leakage audit" in summary


def test_top_k_is_honoured(run_dir):
    assert run_cli(run_dir, "--top-k", "5") == 0
    metrics = json.loads((run_dir / "eval" / "metrics.json").read_text(encoding="utf-8"))
    assert metrics["measured"]["k"] == 5


def test_unknown_truthset_exits_two_and_lists_the_bundled_ones(run_dir, capsys):
    assert main(["--dry-run", "--run", str(run_dir), "--truthset", "nope"]) == 2
    assert "kpneumoniae_hs11286" in capsys.readouterr().out


def test_a_run_without_a_report_exits_two(tmp_path, capsys):
    assert main(["--run", str(tmp_path / "empty")]) == 2
    assert "does not exist" in capsys.readouterr().out


def test_defaults_match_the_issue(run_dir):
    args = build_parser().parse_args([])
    assert args.top_k == 50, "issue #49 asks for the top 50"
    assert args.include_product_matches is False, "the weak tier is opt-in"
    assert args.truthset == "kpneumoniae_hs11286"


def test_rerunning_is_idempotent(run_dir):
    assert run_cli(run_dir) == 0
    first = (run_dir / "eval" / "truth_matches.tsv").read_text(encoding="utf-8")
    assert run_cli(run_dir) == 0
    assert (run_dir / "eval" / "truth_matches.tsv").read_text(encoding="utf-8") == first

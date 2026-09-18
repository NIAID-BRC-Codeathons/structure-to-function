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
from s2f.m2_triage.score import WEIGHTS as TRIAGE_WEIGHTS

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


# ---------------------------------------------------------------------------------------
# The three analyses that only exist once M2 and M3 have run. The committed fixture carries
# triage components and a structures section so the dry run exercises all of them.

def test_triage_can_be_measured_as_well_as_m1_priority(run_dir):
    assert run_cli(run_dir, "--score", "triage") == 0
    metrics = json.loads((run_dir / "eval" / "metrics.json").read_text(encoding="utf-8"))
    assert metrics["score"] == "triage"
    assert metrics["weights"]["pdb_evidence"] == 0.4


def test_the_default_score_is_m1_priority(run_dir):
    assert run_cli(run_dir) == 0
    metrics = json.loads((run_dir / "eval" / "metrics.json").read_text(encoding="utf-8"))
    assert metrics["score"] == "m1_priority"


def test_a_score_the_report_does_not_carry_exits_two(tmp_path, capsys):
    """Asking for triage on an M1-only run names what is actually there."""
    run_dir = tmp_path / "m1only"
    assert run_cli(run_dir) == 0
    report = json.loads((run_dir / "report.json").read_text(encoding="utf-8"))
    for protein in report["proteins"]:
        protein.pop("triage", None)
    (run_dir / "report.json").write_text(json.dumps(report), encoding="utf-8")
    assert main(["--run", str(run_dir), "--score", "triage"]) == 2
    assert "m1_priority" in capsys.readouterr().out


def test_ablation_is_written_and_sweeps_every_component(run_dir):
    assert run_cli(run_dir, "--score", "triage", "--top-k", "20") == 0
    ablation = json.loads((run_dir / "eval" / "ablation.json").read_text(encoding="utf-8"))
    assert len(ablation["effects"]) == 8
    assert ablation["k"] == 20
    # the fixture never fires membrane_penalty, which must read as a data gap
    assert "membrane_penalty" in ablation["never_fired"]


def test_ablation_changes_no_production_weight(run_dir):
    """An ablation reports counterfactuals; it must not rewrite what M2 scored."""
    assert run_cli(run_dir, "--score", "triage") == 0
    ablation = json.loads((run_dir / "eval" / "ablation.json").read_text(encoding="utf-8"))
    assert ablation["weights"] == dict(TRIAGE_WEIGHTS)
    report = read_report(run_dir)
    fixture = json.loads(REPORT_FIXTURE.read_text(encoding="utf-8"))
    assert [p["triage"]["score"] for p in report["proteins"]] == \
           [p["triage"]["score"] for p in fixture["proteins"]]


def test_outcome_is_scored_against_m3_and_excludes_technical_failures(run_dir):
    assert run_cli(run_dir, "--score", "triage", "--top-k", "20") == 0
    outcome = json.loads((run_dir / "eval" / "outcome.json").read_text(encoding="utf-8"))
    assert outcome["counts"]["dockable"] > 0
    assert outcome["counts"]["failed"] > 0
    assert outcome["policy_version"] == "provisional-2026-09-17"
    assert outcome["measured"]["precision"] is not None


def test_agreement_is_written_when_both_scores_are_present(run_dir):
    assert run_cli(run_dir, "--top-k", "20") == 0
    pair = json.loads((run_dir / "eval" / "agreement.json").read_text(encoding="utf-8"))
    assert pair["left"] == "m1_priority" and pair["right"] == "triage"
    assert -1.0 <= pair["spearman"] <= 1.0
    assert 0 <= pair["overlap_at_k"] <= pair["k"]
    assert pair["biggest_disagreements"]


def test_the_summary_carries_all_three_analyses(run_dir):
    assert run_cli(run_dir, "--score", "triage", "--top-k", "20") == 0
    summary = (run_dir / "eval" / "summary.md").read_text(encoding="utf-8")
    assert "Against M3's own outcome" in summary
    assert "Which triage components decide the top" in summary
    assert "m1_priority vs triage" in summary
    assert "No weight is changed" in summary


def test_the_extra_analyses_are_skipped_when_their_inputs_are_missing(tmp_path):
    """An M1-only report still produces the calibration, and no empty extras."""
    run_dir = tmp_path / "m1only"
    report = json.loads(REPORT_FIXTURE.read_text(encoding="utf-8"))
    for protein in report["proteins"]:
        protein.pop("triage", None)
    report["structures"] = []
    run_dir.mkdir(parents=True)
    (run_dir / "report.json").write_text(json.dumps(report), encoding="utf-8")
    assert main(["--run", str(run_dir)]) == 0
    assert (run_dir / "eval" / "metrics.json").is_file()
    assert not (run_dir / "eval" / "ablation.json").is_file()
    assert not (run_dir / "eval" / "agreement.json").is_file()

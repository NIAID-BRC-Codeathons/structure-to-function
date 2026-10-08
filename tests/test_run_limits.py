"""A cap meant for one stage never reaches another, and the caps are in the report (#80)."""

from __future__ import annotations

from s2f.common.io import init_report, read_report
from s2f.m2_triage.__main__ import build_parser as m2_parser
from s2f.run import _record_caps, build_parser, stage_args


def _args(tmp_path, *argv):
    return build_parser().parse_args(["--run", str(tmp_path), *argv])


def test_a_structure_limit_leaves_protein_scoring_untouched(tmp_path) -> None:
    args = _args(tmp_path, "--structure-limit", "10")
    assert "--limit" not in stage_args("m2", args)
    assert m2_parser().parse_args(stage_args("m2", args)).limit == 0
    m3 = stage_args("m3", args)
    assert m3[m3.index("--limit") + 1] == "10"


def test_the_m2_limit_does_not_reach_m3(tmp_path) -> None:
    args = _args(tmp_path, "--limit", "10")
    assert "--limit" not in stage_args("m3", args)
    m2 = stage_args("m2", args)
    assert m2[m2.index("--limit") + 1] == "10"


def test_the_caps_in_force_are_recorded_in_the_run_section(tmp_path) -> None:
    init_report(tmp_path, "r")
    _record_caps(_args(tmp_path, "--structure-limit", "10", "--seq-cap", "500"),
                 ["m1", "m2", "m3"])
    caps = read_report(tmp_path)["run"]["caps"]
    assert caps == {"m1": {"max_features": None, "seq_cap": 500},
                    "m2": {"limit": None}, "m3": {"limit": 10}}


def test_caps_are_not_invented_before_m1_has_written_the_run_section(tmp_path) -> None:
    _record_caps(_args(tmp_path, "--limit", "5"), ["m2"])
    assert "run" not in read_report(tmp_path)

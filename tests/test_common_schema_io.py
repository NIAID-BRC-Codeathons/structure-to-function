"""Data-contract tests (issue #2): schema validation and the section writer."""

from __future__ import annotations

import copy
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest

from s2f.common.io import (
    init_report,
    read_report,
    report_path,
    update_proteins,
    update_section,
)
from s2f.common.schema import (
    SCHEMA_VERSION,
    SECTION_SCHEMAS,
    SchemaError,
    is_valid,
    validate,
    validate_report,
)

FIXTURE = Path(__file__).resolve().parents[1] / "fixtures" / "report.fixture.json"


def load_fixture() -> dict:
    return json.loads(FIXTURE.read_text(encoding="utf-8"))


# --- the definition of done: accepts the fixture, rejects a mangled copy -------------


def test_fixture_validates() -> None:
    validate_report(load_fixture())


def test_fixture_has_the_shape_the_issue_asks_for() -> None:
    report = load_fixture()
    assert len(report["proteins"]) == 5
    assert len(report["structures"]) == 1
    assert len(report["ligands"]) == 3
    assert len(report["docking"]) == 2
    assert len(report["disease"]["claims"]) == 1


@pytest.mark.parametrize(
    ("section", "mangle", "reason"),
    [
        ("proteins", lambda r: r["proteins"][0].pop("feature_id"), "canonical key removed"),
        ("proteins", lambda r: r["proteins"][0]["annotations"][0].pop("source"), "annotation without a source"),
        ("proteins", lambda r: r["proteins"][0]["triage"].pop("components"), "score not reproducible"),
        ("proteins", lambda r: r["proteins"][0]["annotations"][0].update(identity=7.0), "identity outside 0-1"),
        ("structures", lambda r: r["structures"][0].update(source="vibes"), "source not an allowed value"),
        ("kg", lambda r: r["kg"]["edges"][0].pop("provenance"), "edge with no provenance"),
        ("kg", lambda r: r["kg"]["edges"][0]["provenance"].pop("retrieved_at"), "provenance without a date"),
        ("ligands", lambda r: r["ligands"][0].pop("provenance"), "ligand with no provenance"),
        ("docking", lambda r: r["docking"][0].pop("feature_id"), "pose not tied to a protein"),
        ("disease", lambda r: r["disease"]["claims"][0].pop("claim"), "empty claim"),
        ("run", lambda r: r["run"].pop("run_id"), "run with no id"),
    ],
)
def test_mangled_sections_are_rejected(section, mangle, reason) -> None:
    report = load_fixture()
    mangle(report)

    assert not is_valid(section, report[section]), f"should have been rejected: {reason}"
    with pytest.raises(SchemaError):
        validate(section, report[section])


def test_unknown_section_is_an_error() -> None:
    with pytest.raises(SchemaError, match="unknown section"):
        validate("nonsense", {})

    report = load_fixture()
    report["surprise"] = {}
    with pytest.raises(SchemaError, match="unknown top-level key"):
        validate_report(report)


def test_error_message_names_every_problem_not_just_the_first() -> None:
    report = load_fixture()
    report["proteins"][0].pop("feature_id")
    report["proteins"][1].pop("feature_id")

    with pytest.raises(SchemaError) as excinfo:
        validate("proteins", report["proteins"])
    assert "2 problem(s)" in str(excinfo.value)


def test_every_section_has_a_schema() -> None:
    expected = {
        "run", "genome", "proteins", "structures", "kg",
        "ligands", "docking", "analyses", "disease", "report",
    }
    assert set(SECTION_SCHEMAS) == expected


def test_schema_is_permissive_about_unknown_fields() -> None:
    """Modules must be able to ship a field before the shape is settled."""
    report = load_fixture()
    report["proteins"][0]["something_new"] = {"from": "a future module"}
    validate("proteins", report["proteins"])


# --- the section writer ---------------------------------------------------------------


def test_update_section_writes_and_reads_back(tmp_path: Path) -> None:
    report = update_section(tmp_path, "genome", {"taxon_id": 243273})

    assert report["genome"]["taxon_id"] == 243273
    assert report["schema_version"] == SCHEMA_VERSION
    assert read_report(tmp_path)["genome"]["taxon_id"] == 243273


def test_update_section_validates_before_writing_anything(tmp_path: Path) -> None:
    with pytest.raises(SchemaError):
        update_section(tmp_path, "structures", [{"source": "pdb"}])  # no feature_id

    assert not report_path(tmp_path).exists()  # nothing was written


def test_a_module_never_clobbers_another_ones_section(tmp_path: Path) -> None:
    """The definition of done: two modules writing different sections concurrently."""
    fixture = load_fixture()
    update_section(tmp_path, "run", fixture["run"])

    sections = {
        "genome": fixture["genome"],
        "proteins": fixture["proteins"],
        "structures": fixture["structures"],
        "kg": fixture["kg"],
        "ligands": fixture["ligands"],
        "docking": fixture["docking"],
        "disease": fixture["disease"],
    }
    with ThreadPoolExecutor(max_workers=len(sections)) as pool:
        list(pool.map(lambda item: update_section(tmp_path, item[0], item[1]), sections.items()))

    final = read_report(tmp_path)
    for key in sections:
        assert key in final, f"{key} was clobbered by a concurrent writer"
    assert final["run"]["run_id"] == fixture["run"]["run_id"]
    validate_report(final)


def test_repeated_writes_of_one_section_keep_the_last(tmp_path: Path) -> None:
    for taxon in (1, 2, 3):
        update_section(tmp_path, "genome", {"taxon_id": taxon})
    assert read_report(tmp_path)["genome"]["taxon_id"] == 3


def test_a_failed_write_leaves_the_previous_version_intact(tmp_path: Path) -> None:
    update_section(tmp_path, "genome", {"taxon_id": 243273})
    before = report_path(tmp_path).read_text()

    with pytest.raises(SchemaError):
        update_section(tmp_path, "proteins", [{"locus_tag": "no feature_id here"}])

    assert report_path(tmp_path).read_text() == before
    assert not list(tmp_path.glob(".report-*.tmp"))  # no temp files left behind


def test_init_report_is_idempotent(tmp_path: Path) -> None:
    first = init_report(tmp_path, "run-1", input_file="x.fna")
    created = first["run"]["created_at"]

    second = init_report(tmp_path, "run-2")

    assert second["run"]["run_id"] == "run-1"  # does not overwrite an existing run
    assert second["run"]["created_at"] == created


def test_read_report_of_a_new_run_is_an_empty_skeleton(tmp_path: Path) -> None:
    assert read_report(tmp_path) == {"schema_version": SCHEMA_VERSION}


def test_read_report_rejects_a_corrupt_file(tmp_path: Path) -> None:
    report_path(tmp_path).parent.mkdir(parents=True, exist_ok=True)
    report_path(tmp_path).write_text("{not json")

    with pytest.raises(ValueError, match="not valid JSON"):
        read_report(tmp_path)


# --- per-protein enrichment -----------------------------------------------------------


def test_update_proteins_merges_enrichment_into_m1s_list(tmp_path: Path) -> None:
    fixture = load_fixture()
    bare = [{"feature_id": p["feature_id"], "product": p["product"]} for p in fixture["proteins"]]
    update_section(tmp_path, "proteins", bare)

    target = bare[0]["feature_id"]
    update_proteins(
        tmp_path,
        {target: {"xrefs": {"uniprot": "Q49434"}, "annotations": [{"source": "rcsb_sequence"}]}},
    )

    record = next(p for p in read_report(tmp_path)["proteins"] if p["feature_id"] == target)
    assert record["product"] == bare[0]["product"]  # M1's field survives
    assert record["xrefs"]["uniprot"] == "Q49434"
    assert record["annotations"] == [{"source": "rcsb_sequence"}]


def test_annotations_accumulate_across_modules(tmp_path: Path) -> None:
    """M2, Foldseek (#9) and eggNOG (#10) all write here; `source` keeps them apart."""
    update_section(tmp_path, "proteins", [{"feature_id": "fig|1.1.peg.1"}])

    update_proteins(tmp_path, {"fig|1.1.peg.1": {"annotations": [{"source": "rcsb_sequence", "identity": 0.9}]}})
    update_proteins(tmp_path, {"fig|1.1.peg.1": {"annotations": [{"source": "foldseek", "tm_score": 0.7}]}})

    record = read_report(tmp_path)["proteins"][0]
    assert [a["source"] for a in record["annotations"]] == ["rcsb_sequence", "foldseek"]


def test_enriching_an_unknown_feature_id_is_refused(tmp_path: Path) -> None:
    update_section(tmp_path, "proteins", [{"feature_id": "fig|1.1.peg.1"}])

    with pytest.raises(KeyError, match="not present in proteins"):
        update_proteins(tmp_path, {"fig|9.9.peg.9": {"xrefs": {"uniprot": "P00001"}}})


def test_concurrent_protein_enrichment_does_not_lose_updates(tmp_path: Path) -> None:
    ids = [f"fig|1.1.peg.{i}" for i in range(1, 21)]
    update_section(tmp_path, "proteins", [{"feature_id": fid} for fid in ids])

    def enrich(fid: str) -> None:
        update_proteins(tmp_path, {fid: {"xrefs": {"uniprot": f"P{fid.split('.')[-1]:0>5}"}}})

    with ThreadPoolExecutor(max_workers=8) as pool:
        list(pool.map(enrich, ids))

    records = read_report(tmp_path)["proteins"]
    assert all(r.get("xrefs", {}).get("uniprot") for r in records), "a concurrent update was lost"


# --- the guarantee across processes, which is what modules actually are -----------------

WRITER = """
import sys, time
sys.path.insert(0, {repo!r})
from s2f.common.io import update_section
run_dir, key, start_at = sys.argv[1], sys.argv[2], float(sys.argv[3])
payload = [{{"method": key, "parameters": {{"pad": "y" * 400}}}} for _ in range(3000)]
while time.time() < start_at:
    pass
update_section(run_dir, key, payload, validate_section=False)
"""


def test_separate_processes_writing_at_the_same_instant_keep_every_section(tmp_path: Path) -> None:
    """Modules are processes, not threads, so the lock has to hold across them.

    Verified to be a real detector: with the flock removed, only one section of six survives.
    The payload is large and the writers are synchronised on a wall-clock instant, because a
    small or staggered write does not reproduce the race.
    """
    import subprocess
    import sys
    import time

    script = tmp_path / "writer.py"
    script.write_text(WRITER.format(repo=str(Path(__file__).resolve().parents[1])))
    run_dir = tmp_path / "run"
    run_dir.mkdir()

    keys = ["analyses", "a1", "a2", "a3", "a4", "a5"]
    start_at = time.time() + 1.0
    procs = [
        subprocess.Popen([sys.executable, str(script), str(run_dir), key, str(start_at)])
        for key in keys
    ]
    for proc in procs:
        assert proc.wait(timeout=60) == 0

    written = {k for k in read_report(run_dir) if k != "schema_version"}
    assert written == set(keys), f"lost sections: {sorted(set(keys) - written)}"


# --- the contract file must be readable by more than its owner ----------------------
def test_report_json_is_not_written_owner_only(tmp_path: Path) -> None:
    """`NamedTemporaryFile` creates its file 0600 and `os.replace` preserves the mode.

    Without an explicit chmod, report.json lands owner-only. That is invisible on a laptop
    and breaks collaboration on a shared filesystem: no other account in the run
    directory's group can read the one file every module is supposed to join on. Observed
    on a shared compute node, 2026-09-17.
    """
    from s2f.common import io as io_module

    init_report(tmp_path, "run_perms")
    mode = report_path(tmp_path).stat().st_mode & 0o777
    expected = 0o666 & ~io_module._UMASK

    assert mode == expected, (
        f"report.json is mode {mode:04o}, expected {expected:04o} under umask "
        f"{io_module._UMASK:04o}. The chmod in _atomic_write has been removed, so the "
        f"file keeps NamedTemporaryFile's 0600."
    )
    if expected != 0o600:
        assert mode & 0o044, (
            f"report.json is mode {mode:04o}: unreadable outside its owner even though "
            f"umask {io_module._UMASK:04o} permits it."
        )

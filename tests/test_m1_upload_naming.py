"""The workspace-upload path: the caller's side of `cga.upload`, and the `name=` trap.

`TestUploadSafety` in `test_m1_parse.py` covers `upload` itself — that it verifies the
workspace copy by size and refuses a stale or missing one. These cover the two things it
cannot:

1. **The caller.** `upload` can only *detect* a collision. What *prevents* one is M1 giving
   each run's blinded contigs a distinct filename, and nothing tested that. A constant name
   has been in the tree at least once.
2. **The path/argv correspondence.** The tests in `test_m1_parse.py` monkeypatch
   `remote_size`, so they would still pass if the path `upload` returns were decoupled from
   the file `p3-cp` actually wrote. A `name=` parameter that renamed only the return value
   existed at one point and did exactly that.

Both failures are silent: the run exits 0, writes a plausible `report.json`, and CGA has
annotated with a taxon from the wrong organism (`docs/pitfalls.md` #23).
"""

from __future__ import annotations

import inspect
from pathlib import Path

from s2f.m1_genome import cga


def _fake_p3(monkeypatch, *, ls_size: int | None):
    """Stand in for the `p3-` CLI, recording the argv of every call."""
    calls: list[list[str]] = []

    def fake_run(cmd, timeout=900):
        calls.append(list(cmd))
        if cmd[0] == "p3-ls":
            if ls_size is None:
                raise cga.CgaError("p3-ls: No such file or directory")
            return f"-rw-- someone@bvbrc {ls_size} Sep 17 09:32 {Path(cmd[-1]).name}\n"
        return ""

    monkeypatch.setattr(cga, "_run", fake_run)
    return calls


def test_the_returned_path_is_where_p3_cp_actually_wrote_the_file(tmp_path,
                                                                 monkeypatch) -> None:
    """`p3-cp SRC ws:DIR` writes `DIR/<local basename>`; the return value must match it.

    Asserted against the recorded argv rather than against a literal, so any future
    change that decouples the two fails here instead of surfacing as Minhash reading a
    workspace path with no file at it.
    """
    local = tmp_path / "run_zeta.contigs.fna"
    local.write_text(">contig_1\nACGT\n")
    calls = _fake_p3(monkeypatch, ls_size=local.stat().st_size)

    remote = cga.upload(local, "/someone@bvbrc/home/s2f")

    copy_call = next(call for call in calls if call[0] == "p3-cp")
    source, destination = Path(copy_call[1]), copy_call[2]
    # Sliced rather than str.removeprefix, which is 3.9+: the repo targets 3.11, but a
    # test has no reason to be the thing that fails first on an older interpreter.
    target = destination[len("ws:"):] if destination.startswith("ws:") else destination
    expected = f"{target}/{source.name}"
    assert remote == expected, (
        f"upload returned {remote} but p3-cp wrote {expected}; nothing may decouple the "
        f"returned path from the local basename (docs/pitfalls.md #23)"
    )


def test_upload_offers_no_rename_parameter() -> None:
    """A `name` parameter cannot work here, so it must not come back.

    `p3-cp` writes under the local basename, so a `name` that differed from it would
    change only the string `upload` returns. Kept as a named check so the failure
    explains itself rather than surfacing months later as a dangling workspace path.
    """
    assert "name" not in inspect.signature(cga.upload).parameters, (
        "cga.upload has grown a `name` parameter. p3-cp writes the file under its local "
        "basename, so a differing `name` returns a workspace path with no file at it. "
        "Make the local filename unique instead (docs/pitfalls.md #23)."
    )


def test_m1_gives_each_run_a_distinct_workspace_path(tmp_path, monkeypatch) -> None:
    """Two runs must not upload to the same workspace path.

    This is the guard `upload`'s size check cannot provide: it detects a collision, it
    does not prevent one. M1 names the blinded contigs for the run id; with a constant
    name the second run's `p3-cp` silently no-ops and Minhash calls the first run's
    taxon, which `docs/01a-cga-coverage.md` measures as roughly doubling the CDS count
    and halving mean protein length.
    """
    from s2f.m1_genome.__main__ import main

    assembly = tmp_path / "assembly.fna"
    assembly.write_text(">c1 Klebsiella pneumoniae\nACGTACGT\n")

    uploaded: list[str] = []

    def fake_upload(local, ws_dir, **kwargs):
        # Mirrors p3-cp: the workspace path comes from the LOCAL basename.
        path = f"{ws_dir}/uploads/{Path(local).name}"
        uploaded.append(path)
        return path

    monkeypatch.setattr("s2f.m1_genome.cga.upload", fake_upload)
    monkeypatch.setattr("s2f.m1_genome.cga.submit", lambda *a, **kw: "")
    monkeypatch.setattr(
        "s2f.m1_genome.__main__.call_taxon",
        lambda ws_fasta, **kw: {"called_by": "minhash", "taxon_id": 573,
                                "genetic_code": 11, "taxon_name": "Klebsiella pneumoniae",
                                "called_rank": "species", "n_hits": 3,
                                "top_hit": {"genome_id": "573.1", "distance": 0.001},
                                "hits": []})

    for run_id in ("run_one", "run_two"):
        assert main(["--run", str(tmp_path / run_id), "--contigs", str(assembly),
                     "--ws-dir", "/u@bvbrc/home/s2f", "--dry-run"]) == 0

    assert len(uploaded) == 2
    assert len(set(uploaded)) == 2, (
        f"both runs uploaded to {uploaded[0]}; the second p3-cp will no-op and Minhash "
        f"will call the first run's taxon. M1 must name the blinded contigs for the run "
        f"(docs/pitfalls.md #23)."
    )

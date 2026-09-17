"""Reading and writing `runs/<run_id>/report.json` (issue #2).

One state file per run; each module owns one top-level key and never edits another's
(``00-architecture.md``). Every write goes through :func:`update_section`.

Why not just ``json.dump``: two modules finishing at the same moment both read the file, each
replaces its own key, and the second writer saves a copy of what it read *before* the first
writer's change — silently deleting a section. So a write here is:

1. take an exclusive lock on ``report.json.lock`` (across processes, not just threads),
2. read the current file,
3. replace exactly one key,
4. write a temp file in the same directory and ``os.replace`` it, which is atomic on POSIX.

A reader therefore never sees a half-written file, and a crash mid-write leaves the previous
version intact.
"""

from __future__ import annotations

import fcntl
import json
import os
import tempfile
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterator

from .schema import SCHEMA_VERSION, SchemaError, validate

REPORT_NAME = "report.json"
LOCK_NAME = "report.json.lock"

#: Read at import, which is single-threaded: `os.umask` can only be *read* by setting it,
#: so doing this lazily would race. `NamedTemporaryFile` creates its file 0600 and
#: `os.replace` preserves the mode, so without the chmod below every report.json landed
#: owner-only -- unreadable to collaborators on a shared filesystem, and to any other
#: account in the run directory's group. Observed on a shared node, 2026-09-17.
_UMASK = os.umask(0)
os.umask(_UMASK)


def report_path(run_dir: str | Path) -> Path:
    return Path(run_dir) / REPORT_NAME


@contextmanager
def _locked(run_dir: Path) -> Iterator[None]:
    """Exclusive inter-process lock for the whole read-modify-write."""
    run_dir.mkdir(parents=True, exist_ok=True)
    lock_path = run_dir / LOCK_NAME
    with open(lock_path, "w", encoding="utf-8") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def read_report(run_dir: str | Path) -> dict[str, Any]:
    """Read the report, or an empty skeleton if the run has not written one yet."""
    path = report_path(run_dir)
    if not path.exists():
        return {"schema_version": SCHEMA_VERSION}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"{path} is not valid JSON: {exc}") from exc


def _atomic_write(path: Path, payload: dict[str, Any]) -> None:
    """Write via a temp file in the same directory, then rename. Never a partial file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", dir=path.parent, prefix=".report-", suffix=".tmp", delete=False
    )
    try:
        with handle:
            json.dump(payload, handle, indent=2, sort_keys=False)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        # Before the rename, never after: a reader that opens the new inode must not be
        # able to observe a 0600 window.
        os.chmod(handle.name, 0o666 & ~_UMASK)
        os.replace(handle.name, path)
    except BaseException:
        Path(handle.name).unlink(missing_ok=True)
        raise


def init_report(run_dir: str | Path, run_id: str, **run_fields: Any) -> dict[str, Any]:
    """Create the report with its `run` section if it does not exist yet."""
    run_dir = Path(run_dir)
    report = read_report(run_dir)
    if "run" in report:
        return report
    run_section = {
        "run_id": run_id,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "schema_version": SCHEMA_VERSION,
        **run_fields,
    }
    return update_section(run_dir, "run", run_section)


def update_section(
    run_dir: str | Path, key: str, obj: Any, *, validate_section: bool = True
) -> dict[str, Any]:
    """Replace exactly one top-level key and return the whole report.

    Locked and atomic, so concurrent writers of *different* sections cannot clobber each other.
    Validation is on by default: a module that writes an invalid section is broken, not the
    schema. Pass ``validate_section=False`` only to stage something deliberately partial.
    """
    run_dir = Path(run_dir)
    if validate_section:
        validate(key, obj)  # raises SchemaError before anything is written
    with _locked(run_dir):
        report = read_report(run_dir)
        report.setdefault("schema_version", SCHEMA_VERSION)
        report[key] = obj
        _atomic_write(report_path(run_dir), report)
        return report


def update_proteins(
    run_dir: str | Path, updates: dict[str, dict[str, Any]], *, validate_section: bool = True
) -> dict[str, Any]:
    """Merge per-protein enrichment into `proteins[]`, keyed by feature_id.

    M1 owns the list; M2 adds `xrefs`, `annotations`, `flags` and `triage` to existing records.
    Merging inside the lock keeps that from being a read-then-write race. A feature_id that M1
    never wrote is reported rather than invented, since the canonical key is the one thing every
    module joins on (pitfall #11).
    """
    run_dir = Path(run_dir)
    with _locked(run_dir):
        report = read_report(run_dir)
        proteins = report.get("proteins") or []
        index = {p.get("feature_id"): p for p in proteins if isinstance(p, dict)}

        unknown = sorted(set(updates) - set(index))
        if unknown:
            raise KeyError(
                f"{len(unknown)} feature_id(s) not present in proteins[]: "
                + ", ".join(unknown[:5])
                + (" ..." if len(unknown) > 5 else "")
            )

        for feature_id, fields in updates.items():
            record = index[feature_id]
            for field, value in fields.items():
                if field == "annotations":
                    record.setdefault("annotations", [])
                    record["annotations"].extend(value)
                else:
                    record[field] = value

        if validate_section:
            validate("proteins", proteins)
        report["proteins"] = proteins
        report.setdefault("schema_version", SCHEMA_VERSION)
        _atomic_write(report_path(run_dir), report)
        return report


__all__ = [
    "REPORT_NAME",
    "SchemaError",
    "init_report",
    "read_report",
    "report_path",
    "update_proteins",
    "update_section",
]

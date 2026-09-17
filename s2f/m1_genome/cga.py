"""Submit, poll and retrieve a BV-BRC Comprehensive Genome Analysis job (issue #5).

Thin wrappers over the `p3-` CLI, which must be installed and logged in
(`p3-login`). Everything here shells out rather than calling the AppService API
directly, so the credentials and retry behaviour stay BV-BRC's problem.

Measured runtimes on this project's smoke-test genome (580 kb, 1 contig):
165 s with an exact taxon, 220 s at species level, 156 s with no useful taxon.
"""

from __future__ import annotations

import re
import subprocess
import time
from pathlib import Path

JOB_ID = re.compile(r"\b(\d{6,})\b")
TERMINAL_OK = {"completed", "complete", "succeeded", "success", "ok"}
TERMINAL_BAD = {"failed", "error", "deleted", "cancelled", "canceled", "killed"}


class CgaError(RuntimeError):
    pass


def _run(cmd: list[str], timeout: int = 900) -> str:
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    if proc.returncode != 0:
        raise CgaError(f"{cmd[0]} failed: {(proc.stderr or proc.stdout).strip()[:500]}")
    return proc.stdout


def blind_contigs(source: str | Path, destination: str | Path) -> Path:
    """Strip FASTA headers to contig_1, contig_2, ... (the blinding the charter asks for).

    Blinding headers does not blind the run: the taxon and genetic code passed to CGA
    are what actually tell the annotator what the organism is.
    """
    source, destination = Path(source), Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with source.open(encoding="utf-8") as src, destination.open("w", encoding="utf-8") as dst:
        for line in src:
            if line.startswith(">"):
                n += 1
                dst.write(f">contig_{n}\n")
            else:
                dst.write(line)
    if n == 0:
        raise CgaError(f"no FASTA records in {source}")
    return destination


def upload(local: str | Path, ws_dir: str, *, name: str | None = None) -> str:
    """Copy a local file into the workspace and return its workspace path.

    Minhash reads the FASTA from the workspace, so it has to be there before the
    taxon can be called - `p3-submit-CGA` uploads its own copy later.
    """
    local = Path(local)
    target = f"{ws_dir}/uploads"
    _run(["p3-cp", str(local), f"ws:{target}"], timeout=1800)
    return f"{target}/{name or local.name}"


def submit(contigs: str | Path, *, ws_dir: str, out_name: str, scientific_name: str,
           taxon_id: int, genetic_code: int, label: str = "s2f", domain: str = "Bacteria",
           overwrite: bool = True, dry_run: bool = False) -> str:
    """Submit a CGA job; returns the job id printed by `p3-submit-CGA`."""
    cmd = [
        "p3-submit-CGA",
        "--contigs", str(contigs),
        "--workspace-upload-path", f"{ws_dir}/uploads",
        "--scientific-name", scientific_name,
        "--taxonomy-id", str(taxon_id),
        "--code", str(genetic_code),
        "--domain", domain,
        "--label", label,
    ]
    if overwrite:
        cmd.append("--overwrite")
    if dry_run:
        cmd.insert(1, "--dry-run")
    cmd += [ws_dir, out_name]

    out = _run(cmd)
    if dry_run:
        print(out)
        return ""
    match = JOB_ID.search(out)
    if not match:
        raise CgaError(
            "could not find a job id in p3-submit-CGA output; pass --job-id to resume:\n" + out
        )
    return match.group(1)


def job_status(job_id: str) -> str:
    return _run(["p3-job-status", job_id], timeout=120).strip()


def wait(job_id: str, *, interval: int = 30, timeout: int = 7200, verbose: bool = True) -> str:
    """Poll until the status is terminal. Returns the final status text."""
    started = time.time()
    while True:
        status = job_status(job_id)
        low = status.lower()
        if any(word in low for word in TERMINAL_BAD):
            raise CgaError(f"CGA job {job_id} did not succeed: {status}")
        if any(word in low for word in TERMINAL_OK):
            return status
        if time.time() - started > timeout:
            raise CgaError(f"CGA job {job_id} still running after {timeout}s (status: {status})")
        if verbose:
            print(f"  [{int(time.time() - started):>5}s] {status}", flush=True)
        time.sleep(interval)


def fetch(ws_dir: str, out_name: str, destination: str | Path) -> Path:
    """Copy the job output directory out of the workspace. Needs the `ws:` prefix."""
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    _run(["p3-cp", "-r", f"ws:{ws_dir}/.{out_name}", str(destination)], timeout=1800)
    return destination

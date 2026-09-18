"""ANI between the analysed assembly and its closest genomes (issue #7).

`genome.closest_genomes[]` arrives from Similar Genome Finder carrying a Mash distance
and a shared-k-mer count. Neither is an identity: Mash distance approximates ``1 - ANI``
over a 1000-hash sketch, which is enough to call a taxon and not enough to put a number
in a report. This module computes the identity itself.

**skani, not fastANI.** skani estimates ANI from sparse chaining rather than from
fragment-level BLAST-like mapping, which makes it roughly two orders of magnitude faster
and — the part that matters here — accurate on fragmented drafts, where fastANI degrades.
It is a single static binary with no database, so it runs anywhere the pipeline runs.

**What "no ANI" means.** skani drops a pair whose aligned fraction falls below
``--min-af`` (15% by default) rather than reporting a low identity, because below that
the estimate is not meaningful. A closest genome can therefore come back with no row at
all, and that is an answer: `ani` stays ``None`` and the reason is recorded in
`ani.json`, never guessed at from the Mash distance.

**The self-match.** A genome we submit blind is usually already in BV-BRC, so the
nearest hit is the assembly matching itself at ~100% ANI. Identity cannot be used to
detect this — CGA mints a fresh genome id for the submission (`2097.118`), which matches
no public record — so the row is flagged from the number: ANI at or above
``self_ani`` over most of the genome. It is kept, not dropped: "this assembly is already
public" is a finding. `ani.json` counts the non-self rows separately so a run that found
*only* itself is visible rather than silent (issue #7's definition of done).

**SNP distance is not here.** Issue #7 also asked for a SNP distance per closest genome.
That was deferred to its own issue: it needs a second aligner and a defensible answer to
"how many SNPs" for pairs below the species boundary, where the question stops meaning
anything. The `closest_genomes[]` rows keep their `snp_distance: null` and `ani.json`
records the deferral rather than leaving a silent gap.
"""

from __future__ import annotations

import re
import shutil
import subprocess
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence

from ..common.http import HttpError, OfflineCacheMiss

#: Reference genomes come from the BV-BRC Data API, not from `ftp.bvbrc.org`. That host
#: resolves (140.221.78.70) but refuses both 80 and 443, so it serves FTP only; the Data
#: API answers over HTTPS, needs no `p3-` CLI and no token, and is the same host the rest
#: of M1 already talks to. `http_accept` is a query parameter rather than a request
#: header so the recorded provenance URL is pasteable and reproduces the exact bytes.
#:
#: Verified 2026-09-18 from lambda0: `243273.25` returns 589,848 bytes in one record,
#: byte-identical to `fixtures/genomes/mgen_G37/mgen_G37.fna`.
BASE = "https://www.bv-brc.org/api"
GENOME_FASTA_URL = (BASE + "/genome_sequence/?eq(genome_id,{gid})"
                    "&limit(10000)&http_accept=application/dna+fasta")
#: One row of the same query, used only for its `Content-Range` total: the authoritative
#: contig count, which the downloaded FASTA is checked against.
GENOME_COUNT_URL = BASE + "/genome_sequence/?eq(genome_id,{gid})&select(sequence_id)"

NAMESPACE = "bvbrc_genome_fna"

DEFAULT_MAX_GENOMES = 10
DEFAULT_THREADS = 3
DEFAULT_MIN_AF = 15.0
#: ANI at or above this, over `SELF_MIN_AF` of the query, is the assembly matching itself.
#: 99.95% is well above the ~99.5% seen between the closest distinct strains of a species
#: and below the float noise of a true self-comparison, which reports exactly 100.00.
SELF_ANI = 99.95
SELF_MIN_AF = 90.0

#: skani presets, as `--fast` / `--medium` / `--slow`. "default" passes none.
PRESETS = ("default", "fast", "medium", "slow")

SNP_DEFERRED = ("deferred to a follow-up issue; issue #7 landed ANI only "
                "(docs/01-m1-genome.md)")


class SkaniError(RuntimeError):
    """skani is missing, or exited non-zero."""


@dataclass(frozen=True)
class SkaniRow:
    """One line of `skani dist` output."""

    ref_file: str
    query_file: str
    ani: float
    align_fraction_ref: float
    align_fraction_query: float
    ref_name: str
    query_name: str


# --- skani ------------------------------------------------------------------

#: The columns `skani dist` writes without `--detailed`, in order.
_SKANI_COLUMNS = ("Ref_file", "Query_file", "ANI", "Align_fraction_ref",
                  "Align_fraction_query", "Ref_name", "Query_name")


def parse_skani(text: str) -> list[SkaniRow]:
    """Parse a `skani dist` table.

    Keyed off the header line rather than column position, so a skani release that adds a
    column does not silently shift the numbers. Rows are matched back to genomes by
    ``Ref_file`` — a path this module wrote — and never by ``Ref_name``, which is the
    reference's FASTA defline and belongs to BV-BRC.
    """
    lines = [line for line in text.splitlines() if line.strip()]
    if not lines:
        return []
    header = lines[0].split("\t")
    missing = [c for c in _SKANI_COLUMNS if c not in header]
    if missing:
        raise SkaniError(f"skani output is missing columns {missing}; got {header}")
    index = {name: header.index(name) for name in _SKANI_COLUMNS}
    rows: list[SkaniRow] = []
    for line in lines[1:]:
        parts = line.split("\t")
        if len(parts) < len(header):
            continue
        try:
            rows.append(SkaniRow(
                ref_file=parts[index["Ref_file"]],
                query_file=parts[index["Query_file"]],
                ani=float(parts[index["ANI"]]),
                align_fraction_ref=float(parts[index["Align_fraction_ref"]]),
                align_fraction_query=float(parts[index["Align_fraction_query"]]),
                ref_name=parts[index["Ref_name"]],
                query_name=parts[index["Query_name"]],
            ))
        except ValueError:
            continue
    return rows


def skani_version(exe: str = "skani") -> str:
    """`skani 0.3.2`, or raise if the binary is not usable."""
    path = shutil.which(exe) or (exe if Path(exe).exists() else None)
    if path is None:
        raise SkaniError(
            f"skani not found on PATH as {exe!r}. Install it with one of:\n"
            "  curl -sSL -o ~/bin/skani "
            "https://github.com/bluenote-1577/skani/releases/download/latest/skani "
            "&& chmod +x ~/bin/skani\n"
            "  conda install -c bioconda skani\n"
            "  cargo install skani\n"
            "then pass --skani PATH if it is not on PATH."
        )
    proc = subprocess.run([path, "-V"], capture_output=True, text=True)
    if proc.returncode != 0:
        raise SkaniError(f"{path} -V exited {proc.returncode}: {proc.stderr.strip()[:200]}")
    return (proc.stdout or proc.stderr).strip()


def run_skani(query: str | Path, references: Sequence[str | Path], *,
              exe: str = "skani", threads: int = DEFAULT_THREADS,
              preset: str = "default", min_af: float = DEFAULT_MIN_AF,
              ) -> tuple[list[SkaniRow], list[str]]:
    """Run `skani dist` for one query against many references.

    Returns the parsed rows and the argv that produced them, so the run manifest can
    record the exact command rather than a prose description of it.
    """
    if preset not in PRESETS:
        raise SkaniError(f"unknown skani preset {preset!r}; pick one of {PRESETS}")
    if not references:
        return [], []
    path = shutil.which(exe) or exe
    argv = [str(path), "dist", "-t", str(threads), "--min-af", str(min_af),
            "-q", str(query), "-r", *[str(r) for r in references]]
    if preset != "default":
        argv.insert(2, f"--{preset}")
    proc = subprocess.run(argv, capture_output=True, text=True)
    if proc.returncode != 0:
        raise SkaniError(f"skani dist exited {proc.returncode}: {proc.stderr.strip()[:400]}")
    return parse_skani(proc.stdout), argv


# --- reference genomes ------------------------------------------------------

_DEFLINE = re.compile(rb"^>", re.MULTILINE)


def count_fasta_records(payload: bytes) -> int:
    """How many `>` records a FASTA body contains."""
    return len(_DEFLINE.findall(payload))


def _looks_like_fasta(payload: bytes) -> bool:
    return payload.lstrip()[:1] == b">"


def _short(exc: BaseException, limit: int = 160) -> str:
    """One readable line from an exception whose message may be an HTML error page."""
    text = " ".join(str(exc).split())
    return text if len(text) <= limit else text[:limit].rstrip() + " …"


def expected_contig_count(client, gid: str, *, cache_dir: Path | None = None) -> int | None:
    """The genome's contig count from `Content-Range`, or None if it cannot be read.

    One row is requested and discarded; only the result-set total is wanted. This exists
    because the FASTA route gives no way to tell a complete download from a truncated one:
    the Data API pages at 25 rows by default, `limit()` did not visibly raise it in
    testing, and a reference missing contigs would deflate every ANI computed against it
    without failing anything.

    None means "could not be checked", which is not the same as a mismatch and does not
    reject the download — `items 0-24/*` is a legal answer meaning the total is unknown.
    """
    url = GENOME_COUNT_URL.format(gid=gid)
    try:
        envelope = client.get_envelope(NAMESPACE + "_count", url,
                                       headers={"Range": "items=0-0"})
    except (HttpError, OfflineCacheMiss, OSError, ValueError):
        return None
    raw = str((envelope.get("headers") or {}).get("content-range", ""))
    match = re.search(r"/\s*(\d+)\s*$", raw)
    return int(match.group(1)) if match else None


def fetch_reference_fasta(client, gid: str, *, genomes_dir: Path, cache_dir: Path,
                          ) -> tuple[Path | None, dict[str, Any]]:
    """Get one BV-BRC genome's contigs as a local FASTA, cached.

    The named copy under ``genomes_dir`` is what skani is handed: the shared HTTP cache
    stores bodies under a hash of the URL, and a tool that reads a file called ``a3f9…``
    gives an output table nobody can read. The named copy is hard-linked to the cached
    body where the filesystem allows it, so the bytes are stored once.

    Two failure modes the Data API has and an FTP server does not, both handled here:

    * **An unknown genome id answers 200 with an empty body**, not 404. So "did the
      download work" cannot be read off the status code and is decided from the body.
    * **A truncated result is indistinguishable from a complete one.** The record count is
      checked against `Content-Range`, and a short FASTA is refused rather than silently
      used — half a reference genome produces a plausible, wrong ANI.
    """
    url = GENOME_FASTA_URL.format(gid=gid)
    provenance: dict[str, Any] = {"genome_id": gid, "url": url}
    target = genomes_dir / f"{gid}.fna"
    if target.exists() and target.stat().st_size > 0:
        provenance["from_cache"] = True
        provenance["contigs"] = count_fasta_records(target.read_bytes())
        return target, provenance

    genomes_dir.mkdir(parents=True, exist_ok=True)
    expected = expected_contig_count(client, gid, cache_dir=cache_dir)
    provenance["contigs_expected"] = expected

    try:
        downloaded = client.get_file(NAMESPACE, url, cache_dir=cache_dir)
    except (HttpError, OfflineCacheMiss, OSError) as exc:
        provenance["error"] = f"{type(exc).__name__}: {_short(exc)}"
        return None, provenance

    payload = downloaded.content
    if not _looks_like_fasta(payload):
        provenance["error"] = (
            f"no sequence returned for genome id {gid} "
            f"({len(payload)} bytes; the Data API answers 200 with an empty body for a "
            f"genome id it does not hold)")
        return None, provenance

    contigs = count_fasta_records(payload)
    provenance["contigs"] = contigs
    if expected is not None and contigs != expected:
        provenance["error"] = (
            f"truncated reference: {contigs} contigs downloaded, {expected} in BV-BRC. "
            f"Refusing it — a reference missing contigs deflates ANI without failing.")
        return None, provenance

    cached_path = _cached_body_path(cache_dir, NAMESPACE, client, url)
    _link_or_write(cached_path, target, payload)
    provenance["retrieved_at"] = downloaded.retrieved_at.isoformat()
    provenance["from_cache"] = downloaded.from_cache
    provenance["bytes"] = len(payload)
    return target, provenance


def _cached_body_path(cache_dir: Path, namespace: str, client, url: str) -> Path:
    """Where `CachedJsonClient.get_file` put the body for this URL."""
    return Path(cache_dir) / namespace / client.cache_key({"url": url})


def _link_or_write(source: Path, target: Path, payload: bytes) -> None:
    """Hard-link `target` to the cached body, or write a copy if that is not possible."""
    try:
        if source.exists():
            target.hardlink_to(source)
            return
    except (OSError, AttributeError):
        pass
    target.write_bytes(payload)


# --- the step ---------------------------------------------------------------

def _as_float(value: Any) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def add_ani(closest: Iterable[dict[str, Any]], query_fasta: str | Path, *,
            client, run_dir: Path, exe: str = "skani",
            threads: int = DEFAULT_THREADS, preset: str = "default",
            min_af: float = DEFAULT_MIN_AF, max_genomes: int = DEFAULT_MAX_GENOMES,
            self_ani: float = SELF_ANI,
            ) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Compute ANI for the nearest `max_genomes` rows and return them plus a summary.

    Rows are returned in their original order, whether or not they got a number, and
    every row past ``max_genomes`` is returned untouched — a cap is a budget, not a
    filter, and dropping the tail would change what `closest_genomes` means.
    """
    rows = [dict(row or {}) for row in closest]
    query = Path(query_fasta)
    if not query.exists():
        raise SkaniError(f"query assembly not found: {query}")

    version = skani_version(exe)
    genomes_dir = run_dir / "m1" / "genomes"
    cache_dir = run_dir / "cache"

    considered = [row for row in rows if str(row.get("genome_id") or "")][:max_genomes]
    fetched: dict[str, Path] = {}
    fetch_provenance: list[dict[str, Any]] = []
    for row in considered:
        gid = str(row["genome_id"])
        path, provenance = fetch_reference_fasta(
            client, gid, genomes_dir=genomes_dir, cache_dir=cache_dir)
        fetch_provenance.append(provenance)
        if path is not None:
            fetched[gid] = path

    skani_rows, argv = run_skani(query, list(fetched.values()), exe=exe, threads=threads,
                                 preset=preset, min_af=min_af)
    by_path = {Path(r.ref_file).resolve(): r for r in skani_rows}

    reasons: dict[str, str] = {}
    self_match: dict[str, Any] | None = None
    scored = 0
    for row in considered:
        gid = str(row["genome_id"])
        provenance = next((p for p in fetch_provenance if p["genome_id"] == gid), {})
        # Every considered row carries the keys whether or not it got a number, so a
        # reader can tell "compared, no answer" from "never compared" — which is the
        # difference between a row past the cap and a row skani declined to score.
        row["ani"] = None
        row["ani_align_fraction_query"] = None
        row["ani_align_fraction_ref"] = None
        row["ani_source"] = None
        row["self_match"] = False
        row["ani_reference_url"] = provenance.get("url")
        row["ani_retrieved_at"] = provenance.get("retrieved_at")
        if gid not in fetched:
            reasons[gid] = provenance.get("error") or "reference FASTA not retrieved"
            continue
        hit = by_path.get(fetched[gid].resolve())
        if hit is None:
            reasons[gid] = (f"no skani row: aligned fraction below --min-af {min_af}, "
                            "so ANI is not estimable for this pair")
            continue
        row["ani"] = hit.ani
        row["ani_align_fraction_query"] = hit.align_fraction_query
        row["ani_align_fraction_ref"] = hit.align_fraction_ref
        row["ani_source"] = "skani"
        row["self_match"] = bool(hit.ani >= self_ani
                                 and hit.align_fraction_query >= SELF_MIN_AF)
        scored += 1
        if row["self_match"] and self_match is None:
            self_match = {"genome_id": gid, "name": row.get("name"), "ani": hit.ani}

    non_self = [row for row in considered
                if row.get("ani") is not None and not row.get("self_match")]
    meta = {
        "tool": "skani",
        "version": version,
        "command": argv,
        "params": {"threads": threads, "preset": preset, "min_af": min_af,
                   "max_genomes": max_genomes, "self_ani": self_ani},
        "query_fasta": str(query),
        "computed_at": datetime.now(timezone.utc).isoformat(),
        "genomes_in_section": len(rows),
        "genomes_considered": len(considered),
        "genomes_fetched": len(fetched),
        "genomes_with_ani": scored,
        "non_self_with_ani": len(non_self),
        "self_match": self_match,
        "no_ani_reason": reasons,
        "references": fetch_provenance,
        "snp_distance": SNP_DEFERRED,
    }
    if self_match is not None and not non_self:
        meta["warning"] = ("the only genome with an ANI is the assembly matching itself; "
                           "no distinct relative was measured (issue #7 asks for both)")
    return rows, meta

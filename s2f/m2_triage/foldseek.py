"""Foldseek structure search for proteins whose sequence matches nothing (issue #9).

Our PDB route is a *sequence* search, and its floor is real: below roughly 20–25% identity,
alignment cannot separate homology from chance. Structure survives far longer than sequence, so
a protein can share a fold, an active site and a function at 15% identity and be invisible to
every sequence method.

Measured on ten G37 proteins that our sequence search found nothing for: **seven gained an
informative fold match**, every one at 10–22% identity — all below the 25% cutoff the sequence
search uses. Examples: a "hypothetical protein" matching ArgX/LysX amino-acid ligases at
e-value 5e-19, another matching HU nucleoid-associated protein across three independent
structures, and two "putative esterases" confirmed against bacterial esterases at 1e-12.

**Same fold is not same function.** TIM barrels and Rossmann folds are scaffolds shared by
enzymes doing unrelated chemistry, so a fold match is a lead, not an assignment. Every hit
therefore carries its e-value, alignment coverage, identity and the target's own description —
and a hit whose target is itself an uncharacterized protein is recorded as exactly that, because
it advances nothing.

Transport: the Foldseek web server's ticket API (submit, poll, fetch), so no multi-gigabyte
database download is needed and teammates can run this as-is. Results are cached by structure
checksum, so a rerun costs nothing and ``--offline`` replays with no network. A local
``foldseek`` binary with a downloaded database would be faster for a final full run and is the
only way to get a true TM-score — the web API does not return one (see ``TM_SCORE_NOTE``).
"""

from __future__ import annotations

import hashlib
import json
import re
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Iterable, Sequence

import requests

from ..common.http import JsonCache, OfflineCacheMiss

FOLDSEEK_API = "https://search.foldseek.com/api"
DEFAULT_DATABASES = ("pdb100",)
DEFAULT_MODE = "3diaa"

#: From the pilot: 1e-4 separates consistent, reproducible fold matches from noise. A hit at
#: e-value 1.89 looked like a named protein and was meaningless.
DEFAULT_MAX_EVALUE = 1e-4
DEFAULT_MIN_COVERAGE = 0.4
POLL_INTERVAL_SECONDS = 5.0
MAX_POLLS = 60

#: This is a shared public server with an opaque rate limit. Measured behaviour: batches of
#: 10-25 searches succeed completely; a 163-search run is refused with HTTP 429 for most of it
#: (102 of 163 even with 2 s pacing and three attempts). The 429 carries **no Retry-After and no
#: rate-limit headers**, so the correct wait cannot be computed — only backed off blindly.
#:
#: Practical consequence, and the boundary between the two routes: the web API is for
#: development, spot checks and teammates without a database. A genome-scale run needs the local
#: `foldseek` binary with a downloaded database, which is also the only way to get a TM-score.
MIN_SUBMIT_INTERVAL_SECONDS = 2.0
MAX_ATTEMPTS = 3
#: Blind backoff for 429, in seconds. Long, because the server will not say how long to wait.
RATE_LIMIT_BACKOFF_SECONDS = (30.0, 90.0)
#: Stop hammering a shared server once it is clearly refusing.
CONSECUTIVE_RATE_LIMITS_BEFORE_GIVING_UP = 5

TM_SCORE_NOTE = (
    "The Foldseek web API returns e-value, probability, identity and alignment spans, but not a "
    "TM-score; a true TM-score needs the local binary with a downloaded database."
)

STATUS_FOUND = "found"
STATUS_NO_HIT = "no-hit"
STATUS_FAILED = "query-failed"
STATUS_NO_STRUCTURE = "not-queried-no-structure"

#: Targets whose own name tells us nothing. A fold match to one of these is a recorded
#: non-result, not a functional lead.
UNINFORMATIVE_TARGET = re.compile(
    r"hypothetical|uncharacteri[sz]ed|unknown function|\bDUF\d+|putative protein|"
    r"structural genomics|protein of unknown",
    re.IGNORECASE,
)

# Per-hit fields worth keeping. The API also returns tCa (target C-alpha coordinates) and full
# alignment strings, which are large and useless downstream — dropped before caching.
KEEP_FIELDS = (
    "target", "eval", "prob", "seqId", "score", "alnLength",
    "qStartPos", "qEndPos", "qLen", "dbStartPos", "dbEndPos", "dbLen", "taxId", "taxName",
)


@dataclass
class FoldseekHit:
    """One structural neighbour."""

    target: str
    database: str = ""
    evalue: float | None = None
    probability: float | None = None
    identity: float | None = None       # percent, as the API reports it
    score: float | None = None
    query_coverage: float | None = None
    target_coverage: float | None = None
    alignment_length: int | None = None
    taxonomy: str = ""

    @property
    def pdb_id(self) -> str:
        """`3ia2-assembly2.cif.gz_D ...` -> `3ia2`."""
        return self.target.split("-")[0].split(".")[0][:4].lower()

    @property
    def description(self) -> str:
        parts = self.target.split(" ", 1)
        return parts[1].strip() if len(parts) > 1 else ""

    @property
    def informative(self) -> bool:
        """Does the target's own name say anything about function?"""
        return bool(self.description) and not UNINFORMATIVE_TARGET.search(self.description)

    def as_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["pdb_id"] = self.pdb_id
        payload["informative"] = self.informative
        return payload


@dataclass
class FoldseekResult:
    """Search outcome for one protein."""

    feature_id: str
    status: str = STATUS_NO_STRUCTURE
    structure_source: str = ""
    structure_id: str = ""
    hits: list[FoldseekHit] = field(default_factory=list)
    databases: tuple[str, ...] = DEFAULT_DATABASES
    error: str = ""
    note: str = ""

    @property
    def best(self) -> FoldseekHit | None:
        return self.hits[0] if self.hits else None

    @property
    def best_informative(self) -> FoldseekHit | None:
        for hit in self.hits:
            if hit.informative:
                return hit
        return None

    def as_annotations(self, *, retrieved_at: str, limit: int = 3) -> list[dict[str, Any]]:
        """`annotations[]` rows. `source` is what keeps these apart from sequence hits."""
        if self.status != STATUS_FOUND:
            return [
                {
                    "source": "foldseek",
                    "hit": None,
                    "description": {
                        STATUS_NO_HIT: "no structural neighbour above the significance gate",
                        STATUS_FAILED: "structure search failed",
                        STATUS_NO_STRUCTURE: "no structure available to search with",
                    }.get(self.status, self.status),
                    "retrieved_at": retrieved_at,
                    "databases": list(self.databases),
                    "note": self.note or self.error or None,
                }
            ]
        return [
            {
                "source": "foldseek",
                "hit": hit.pdb_id or hit.target,
                "description": hit.description or None,
                "identity": round(hit.identity / 100.0, 4) if hit.identity is not None else None,
                "coverage": round(hit.query_coverage, 4) if hit.query_coverage is not None else None,
                "evalue": hit.evalue,
                "probability": hit.probability,
                "tm_score": None,  # see TM_SCORE_NOTE
                "retrieved_at": retrieved_at,
                "query_structure": self.structure_id,
                "query_structure_source": self.structure_source,
                "organism": hit.taxonomy or None,
                "informative_target": hit.informative,
                "note": TM_SCORE_NOTE,
            }
            for hit in self.hits[:limit]
        ]


def structure_key(structure: bytes, databases: Sequence[str], mode: str) -> str:
    digest = hashlib.sha256(structure).hexdigest()
    return f"{mode}:{','.join(sorted(databases))}:{digest}"


def parse_alignments(payload: dict[str, Any]) -> list[tuple[str, dict[str, Any]]]:
    """Flatten the API response into (database, trimmed alignment) pairs."""
    out: list[tuple[str, dict[str, Any]]] = []
    for result_set in payload.get("results") or []:
        database = str(result_set.get("db") or "")
        alignments = result_set.get("alignments") or []
        flat = alignments[0] if alignments and isinstance(alignments[0], list) else alignments
        for alignment in flat or []:
            if isinstance(alignment, dict):
                out.append((database, {k: alignment.get(k) for k in KEEP_FIELDS}))
    return out


def _to_float(value: Any) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def build_hits(
    pairs: Iterable[tuple[str, dict[str, Any]]],
    *,
    max_evalue: float = DEFAULT_MAX_EVALUE,
    min_coverage: float = DEFAULT_MIN_COVERAGE,
) -> list[FoldseekHit]:
    """Significant hits, best first. The gate is the pilot's 1e-4, plus a coverage floor."""
    hits: list[FoldseekHit] = []
    for database, alignment in pairs:
        evalue = _to_float(alignment.get("eval"))
        if evalue is None or evalue > max_evalue:
            continue
        qlen = _to_float(alignment.get("qLen")) or 0.0
        qstart, qend = _to_float(alignment.get("qStartPos")), _to_float(alignment.get("qEndPos"))
        coverage = None
        if qlen > 0 and qstart is not None and qend is not None:
            coverage = max(0.0, min(1.0, (qend - qstart + 1) / qlen))
            if coverage < min_coverage:
                continue
        dblen = _to_float(alignment.get("dbLen")) or 0.0
        dbstart, dbend = _to_float(alignment.get("dbStartPos")), _to_float(alignment.get("dbEndPos"))
        target_coverage = None
        if dblen > 0 and dbstart is not None and dbend is not None:
            target_coverage = max(0.0, min(1.0, (dbend - dbstart + 1) / dblen))
        hits.append(
            FoldseekHit(
                target=str(alignment.get("target") or ""),
                database=database,
                evalue=evalue,
                probability=_to_float(alignment.get("prob")),
                identity=_to_float(alignment.get("seqId")),
                score=_to_float(alignment.get("score")),
                query_coverage=coverage,
                target_coverage=target_coverage,
                alignment_length=int(_to_float(alignment.get("alnLength")) or 0) or None,
                taxonomy=str(alignment.get("taxName") or ""),
            )
        )
    hits.sort(key=lambda h: (h.evalue if h.evalue is not None else 1.0))
    return hits


class FoldseekClient:
    """Ticket-based Foldseek web search, cached by structure checksum."""

    def __init__(
        self,
        cache: JsonCache | None = None,
        *,
        session: requests.Session | None = None,
        offline: bool = False,
        databases: Sequence[str] = DEFAULT_DATABASES,
        mode: str = DEFAULT_MODE,
        poll_interval: float = POLL_INTERVAL_SECONDS,
        max_polls: int = MAX_POLLS,
        sleep=time.sleep,
        api: str = FOLDSEEK_API,
        min_submit_interval: float = MIN_SUBMIT_INTERVAL_SECONDS,
        max_attempts: int = MAX_ATTEMPTS,
    ) -> None:
        self.cache = cache
        self.session = session or requests.Session()
        self.offline = offline
        self.databases = tuple(databases)
        self.mode = mode
        self.poll_interval = poll_interval
        self.max_polls = max_polls
        self._sleep = sleep
        self.api = api
        self.min_submit_interval = min_submit_interval
        self.max_attempts = max_attempts
        self._last_submit = 0.0
        self._consecutive_rate_limits = 0
        self.rate_limited = False
        self.failures: list[dict[str, str]] = []

    def _throttle(self) -> None:
        if self.min_submit_interval <= 0:
            return
        wait = self.min_submit_interval - (time.monotonic() - self._last_submit)
        if wait > 0:
            self._sleep(wait)
        self._last_submit = time.monotonic()

    def search_structure(self, structure: bytes, filename: str = "query.cif") -> dict[str, Any]:
        """Submit, poll, fetch. Returns the trimmed payload; caches it by structure checksum."""
        key = structure_key(structure, self.databases, self.mode)
        if self.cache is not None:
            cached = self.cache.get("foldseek_result", key)
            if cached is not None:
                return cached
        if self.offline:
            raise OfflineCacheMiss(f"no cached Foldseek result for {key[-16:]}")

        payload = self._submit_with_retry(structure, filename)
        # Keep only what we use: the raw response carries target coordinates per hit.
        trimmed = {
            "results": [
                {"db": db, "alignments": [alignment]} for db, alignment in parse_alignments(payload)
            ]
        }
        if self.cache is not None:
            self.cache.set("foldseek_result", key, trimmed)
        return trimmed

    def _submit_with_retry(self, structure: bytes, filename: str) -> dict[str, Any]:
        """Submit, poll and fetch, retrying the whole cycle with backoff."""
        data = [("mode", self.mode)] + [("database[]", db) for db in self.databases]
        last_error = ""
        for attempt in range(self.max_attempts):
            try:
                self._throttle()
                response = self.session.post(
                    f"{self.api}/ticket", files={"q": (filename, structure)}, data=data, timeout=300
                )
                response.raise_for_status()
                ticket = response.json().get("id")
                if not ticket:
                    raise RuntimeError(f"no ticket returned: {response.text[:200]}")

                status = ""
                for _ in range(self.max_polls):
                    poll = self.session.get(f"{self.api}/ticket/{ticket}", timeout=120)
                    poll.raise_for_status()
                    status = str(poll.json().get("status") or "")
                    if status in {"COMPLETE", "ERROR"}:
                        break
                    self._sleep(self.poll_interval)
                if status != "COMPLETE":
                    raise RuntimeError(f"ticket {ticket} ended as {status or 'TIMEOUT'}")

                result = self.session.get(f"{self.api}/result/{ticket}/0", timeout=300)
                result.raise_for_status()
                self._consecutive_rate_limits = 0
                return result.json()
            except (requests.RequestException, RuntimeError, ValueError) as exc:
                last_error = f"{type(exc).__name__}: {exc}"
                throttled = "429" in last_error or "Too Many Requests" in last_error
                if throttled:
                    self._consecutive_rate_limits += 1
                    if self._consecutive_rate_limits >= CONSECUTIVE_RATE_LIMITS_BEFORE_GIVING_UP:
                        self.rate_limited = True
                        raise RuntimeError(
                            "Foldseek public server is rate limiting this run "
                            f"({self._consecutive_rate_limits} refusals in a row). It sends no "
                            "Retry-After, so the wait cannot be computed. Use the local foldseek "
                            "binary with a downloaded database for a genome-scale run, or lower "
                            "--foldseek-limit."
                        ) from exc
                if attempt < self.max_attempts - 1:
                    delay = (
                        RATE_LIMIT_BACKOFF_SECONDS[min(attempt, len(RATE_LIMIT_BACKOFF_SECONDS) - 1)]
                        if throttled
                        else min(2.0**attempt * self.min_submit_interval, 30.0)
                    )
                    self._sleep(delay)
        raise RuntimeError(f"Foldseek failed after {self.max_attempts} attempts — {last_error}")

    def search_protein(
        self,
        feature_id: str,
        structure: bytes | None,
        *,
        structure_source: str = "",
        structure_id: str = "",
        max_evalue: float = DEFAULT_MAX_EVALUE,
        min_coverage: float = DEFAULT_MIN_COVERAGE,
    ) -> FoldseekResult:
        """Search one protein. A protein with no structure is `not-queried`, never `no-hit`."""
        result = FoldseekResult(
            feature_id=feature_id,
            structure_source=structure_source,
            structure_id=structure_id,
            databases=self.databases,
        )
        if not structure:
            result.status = STATUS_NO_STRUCTURE
            result.note = "no usable structure; needs a prediction before a structure search"
            return result
        try:
            payload = self.search_structure(structure, filename=f"{structure_id or feature_id}.cif")
        except (requests.RequestException, RuntimeError, OfflineCacheMiss) as exc:
            result.status = STATUS_FAILED
            result.error = str(exc)[:300]
            self.failures.append({"feature_id": feature_id, "error": result.error})
            return result

        result.hits = build_hits(
            parse_alignments(payload), max_evalue=max_evalue, min_coverage=min_coverage
        )
        result.status = STATUS_FOUND if result.hits else STATUS_NO_HIT
        if not result.hits:
            result.note = f"no neighbour at e-value <= {max_evalue:g} and coverage >= {min_coverage:.0%}"
        return result


def fetch_structure(session: requests.Session, url: str, destination: Path) -> bytes:
    """Download a model (AlphaFold mmCIF) once and keep it under the run directory."""
    destination = Path(destination)
    if destination.exists() and destination.stat().st_size > 0:
        return destination.read_bytes()
    destination.parent.mkdir(parents=True, exist_ok=True)
    response = session.get(url, timeout=300)
    response.raise_for_status()
    destination.write_bytes(response.content)
    return response.content


def summarize(results: Iterable[FoldseekResult]) -> dict[str, Any]:
    results = list(results)
    statuses: dict[str, int] = {}
    for result in results:
        statuses[result.status] = statuses.get(result.status, 0) + 1
    found = [r for r in results if r.status == STATUS_FOUND]
    informative = [r for r in found if r.best_informative is not None]
    identities = [
        r.best.identity for r in found if r.best is not None and r.best.identity is not None
    ]
    return {
        "proteins": len(results),
        "status": dict(sorted(statuses.items())),
        "with_hits": len(found),
        "with_informative_hit": len(informative),
        "uninformative_only": len(found) - len(informative),
        "median_identity_of_best_hit": (
            round(sorted(identities)[len(identities) // 2], 1) if identities else None
        ),
        "databases": list(DEFAULT_DATABASES),
        "tm_score": TM_SCORE_NOTE,
        "caveat": (
            "A shared fold is a lead, not a function: common scaffolds carry unrelated chemistry. "
            "Hits whose target is itself uncharacterized are marked informative_target=false."
        ),
    }

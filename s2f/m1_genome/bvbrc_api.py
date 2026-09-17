"""BV-BRC Data API client and genome resolution — the no-account M1 route (issue #5).

`cga.py` submits an assembly to the Comprehensive Genome Analysis *service*, which needs a
BV-BRC account, the `p3-*` CLI and a job that runs for hours. This module is the other
route: query the public Data REST API for a genome BV-BRC has **already** annotated. No
account, no job, seconds instead of hours.

The two are not interchangeable and the difference is recorded in `genome.annotation_route`
so nothing downstream has to guess:

* `cga` annotates *this* assembly. The feature calls belong to the sample in hand.
* `api` describes a *reference or representative* genome for the same organism. Gene content
  is the reference's, not the isolate's, so an isolate-specific gene is invisible and a
  reference-only gene is a false positive. Fine for orienting a run and for the report;
  never a substitute for annotating the sample (`docs/01b-m1-api-mode.md`).

Every request goes through :class:`s2f.common.http.CachedJsonClient`, so a run is cached,
retried and replayable offline like the rest of the pipeline (issue #4).
"""

from __future__ import annotations

import gzip
import re
from typing import Any
from urllib.parse import quote

from ..common.http import CachedJsonClient, HttpError

BASE = "https://www.bv-brc.org/api"
GENOME_URL = "https://www.bv-brc.org/view/Genome/{gid}"
FEATURE_URL = "https://www.bv-brc.org/view/Feature/{fid}"

#: Values BV-BRC uses for "we did not record this". Treated as absent, not as data.
EMPTYISH = {"", "-", "unknown", "na", "none", "not collected", "not provided", "missing"}


class GenomeResolutionError(RuntimeError):
    """No BV-BRC genome could be identified for the given input."""


def rql_value(value: Any) -> str:
    """URL-encode one RQL value, quoting it when it contains whitespace.

    An unquoted multi-word value matches loosely in RQL, so `eq(species,Klebsiella
    pneumoniae)` silently becomes a prefix-ish match on the first token. Quoting first,
    then percent-encoding, is what makes it an exact match.
    """
    text = str(value)
    if re.search(r"\s", text):
        text = f'"{text}"'
    return quote(text, safe="")


class BvbrcApi:
    """RQL-over-HTTP reads from the BV-BRC Data API, cached through the common client."""

    def __init__(self, client: CachedJsonClient, *, base: str = BASE) -> None:
        self.client = client
        self.base = base.rstrip("/")
        #: Non-fatal failures, so one empty metadata facet cannot abort a whole run.
        self.failures: list[dict[str, str]] = []

    # --- raw access ---------------------------------------------------------
    def _url(self, core: str, rql: str) -> str:
        """RQL goes in the raw query string, never through `params`.

        RQL is built from parentheses and commas that are already percent-encoded where
        they need to be. Handing it to requests as a parameter mapping would re-encode
        the structural characters and turn `eq(genome_id,x)` into a literal field name.
        """
        return f"{self.base}/{core}/?{rql}"

    def query(self, core: str, rql: str) -> list[dict[str, Any]]:
        """One page of results. `rql` is already-encoded RQL, passed through verbatim."""
        envelope = self.client.get_envelope(
            f"bvbrc_{core}", self._url(core, rql),
            headers={"Accept": "application/json"},
        )
        body = envelope.get("body")
        return body if isinstance(body, list) else []

    def query_all(self, core: str, rql: str, *, page: int = 25000,
                  cap: int = 50000) -> list[dict[str, Any]]:
        """Page through a result set with the `Range` header until the total or `cap`.

        BV-BRC reports the result-set size only in `Content-Range`, which makes two
        failure modes possible and silent, so both are checked rather than assumed:

        * **A missing or unparseable total** used to end the loop after one page and
          report success — an intermediary that strips `Content-Range` would quietly halve
          a proteome. A full page with no total is now a hard error, because a truncated
          proteome that looks complete is worse than a failed run. A *short* page with no
          total is the genuine end of the data and is accepted.
        * **A server that ignores `Range`** returns page 1 every time, appending
          duplicates until `cap` and inflating every count in the report. Progress is
          checked by feature identity, and a page that adds nothing new ends the loop.
        """
        out: list[dict[str, Any]] = []
        seen: set[str] = set()
        start = 0
        while start < cap:
            fetch = min(page, cap - start)
            envelope = self.client.get_envelope(
                f"bvbrc_{core}", self._url(core, rql),
                headers={"Accept": "application/json",
                         "Range": f"items={start}-{start + fetch - 1}"},
            )
            chunk = envelope.get("body")
            if not isinstance(chunk, list) or not chunk:
                break

            added = 0
            for row in chunk:
                key = _row_key(row)
                if key is not None:
                    if key in seen:
                        continue
                    seen.add(key)
                out.append(row)
                added += 1
            if added == 0:
                break  # the server is replaying a page it already gave us

            total = _content_range_total(envelope)
            if total is None:
                if len(chunk) >= fetch:
                    raise HttpError(
                        f"{core}: no usable Content-Range total after {len(out)} rows, and "
                        f"the page was full — refusing to report a possibly truncated "
                        f"result set as complete"
                    )
                break  # short page, no total: the end of the data
            if len(out) >= total or len(out) >= cap:
                break
            start += len(chunk)
        return out[:cap]

    def facet(self, core: str, base_rql: str, field: str, *, mincount: int = 1,
              limit: int = 25) -> list[tuple[str, int]]:
        """Solr facet counts for one field as `[(value, count), ...]`, commonest first.

        **Returns a list, not a dict, and sorts explicitly.** Solr happens to emit facets
        in descending count order and every consumer here wants the commonest value, but
        relying on that order through a dict is a trap: `JsonCache` serialises with
        `sort_keys=True`, so a cached facet comes back *alphabetically*. A live run and an
        `--offline` replay of the same run would then report different top values —
        "the top host for this species" could change from a 5000-genome *Homo sapiens* to
        a 4-genome *Ailuropoda*. Sorting by count here makes the answer identical either
        way; the value is the tie-break so equal counts are still deterministic.

        Facets are context, not findings — a missing one degrades the report rather than
        failing the run, so a transport error is recorded and an empty list returned.
        """
        rql = (f"{base_rql}&limit(1)&facet((field,{field}),(mincount,{mincount}),"
               f"(limit,{limit}))&json(nl,map)")
        try:
            envelope = self.client.get_envelope(
                f"bvbrc_{core}_facet", self._url(core, rql),
                headers={"Accept": "application/solr+json"},
            )
        except HttpError as exc:
            self.failures.append({"facet": f"{core}.{field}", "error": str(exc)})
            return []
        body = envelope.get("body") or {}
        if not isinstance(body, dict):
            return []
        counts = body.get("facet_counts", {}).get("facet_fields", {}).get(field, {})
        if not isinstance(counts, dict):
            return []
        pairs = [(str(value), int(count)) for value, count in counts.items()
                 if isinstance(count, (int, float))]
        pairs.sort(key=lambda item: (-item[1], item[0]))
        return pairs

    # --- genome resolution --------------------------------------------------
    def resolve_genome(self, *, genome_id: str | None = None, species: str | None = None,
                       fasta: str | None = None) -> tuple[dict[str, Any], str]:
        """Return `(genome_record, how_it_was_resolved)`.

        The resolution string is kept and written into `genome.resolution` because the
        answer to "which genome is this report actually about" is not recoverable later:
        a species name can resolve to a different reference as BV-BRC is recurated.
        """
        if genome_id:
            recs = self.query("genome", f"eq(genome_id,{rql_value(genome_id)})")
            if not recs:
                raise GenomeResolutionError(f"genome_id {genome_id} not found in BV-BRC")
            return recs[0], f"explicit genome_id={genome_id}"

        name = species
        note = f"--species {species!r}" if species else ""
        if not name and fasta:
            defline, _, _ = read_fasta_defline(fasta)
            name = organism_from_defline(defline)
            note = f"organism parsed from FASTA header: {name!r}"
            if not name:
                raise GenomeResolutionError(
                    "could not detect an organism from the FASTA header; "
                    "pass --species 'Genus species' or --genome-id, or annotate the "
                    "assembly itself with the CGA route"
                )
        if not name:
            raise GenomeResolutionError("no genome_id, species or FASTA to resolve from")

        # Preference order is deliberate: a Reference/Representative genome is curated and
        # stable, a Complete+Good one is at least trustworthy, and only then do we fall
        # back to "the biggest thing with this name" — which is recorded, because it is a
        # much weaker claim about the organism.
        fields = ("select(genome_id,genome_name,reference_genome,genome_quality,"
                  "genome_status,patric_cds)")
        attempts = [
            (f"and(eq(species,{rql_value(name)}),or(eq(reference_genome,Reference),"
             f"eq(reference_genome,Representative)))", "reference/representative genome"),
            (f"and(eq(species,{rql_value(name)}),eq(genome_status,Complete),"
             f"eq(genome_quality,Good))", "complete genome of Good quality"),
            (f"eq(species,{rql_value(name)})", "largest genome for the species"),
            (f"eq(genus,{rql_value(name.split()[0])})", "largest genome in the genus only"),
        ]
        for rql, how in attempts:
            hits = self.query("genome", f"{rql}&{fields}&sort(-patric_cds)&limit(1)")
            if not hits:
                continue
            full = self.query("genome", f"eq(genome_id,{rql_value(hits[0]['genome_id'])})")
            if not full:
                continue
            record = full[0]
            detail = (f"{note}; matched the {how} -> BV-BRC genome "
                      f"{record['genome_id']} ({record.get('genome_name')})")
            return record, detail.lstrip("; ")
        raise GenomeResolutionError(
            f"no BV-BRC genome found for {name!r}; try --genome-id, or annotate the "
            f"assembly itself with the CGA route"
        )


#: Fields that identify a BV-BRC row, best first. Used only to detect a server replaying
#: a page, so any stable per-row id will do.
_ROW_ID_FIELDS = ("patric_id", "genome_id", "feature_id", "md5", "taxon_id", "id")


def _row_key(row: Any) -> str | None:
    """A stable identity for one row, or None when it has no usable id."""
    if not isinstance(row, dict):
        return None
    for field in _ROW_ID_FIELDS:
        value = row.get(field)
        if value:
            return f"{field}={value}"
    return None


def _content_range_total(envelope: dict[str, Any]) -> int | None:
    """The result-set total from `Content-Range`, or None if it does not carry one.

    `items 0-24/*` is legal and means "total unknown", so it must come back as None
    rather than matching a digit run elsewhere in the header.
    """
    raw = str((envelope.get("headers") or {}).get("content-range", ""))
    match = re.search(r"/\s*(\d+)\s*$", raw)
    return int(match.group(1)) if match else None


def read_fasta_defline(path: str) -> tuple[str, int, int]:
    """Return `(first_defline, n_sequences, total_bases)` from a FASTA, gzipped or not."""
    opener = gzip.open if str(path).endswith(".gz") else open
    first: str | None = None
    n_seq = 0
    n_bases = 0
    with opener(path, "rt", errors="ignore") as handle:  # type: ignore[operator]
        for line in handle:
            if line.startswith(">"):
                n_seq += 1
                if first is None:
                    first = line[1:].strip()
            else:
                n_bases += len(line.strip())
    if first is None:
        raise ValueError(f"no FASTA header (>) found in {path}")
    return first, n_seq, n_bases


def organism_from_defline(defline: str) -> str | None:
    """Best-effort `Genus species` from a FASTA header.

    Deliberately conservative: it wants a capitalised genus followed by a lower-case
    epithet, and returns None rather than a guess. A wrong organism resolves to a wrong
    reference genome and every downstream claim inherits the error, so failing loudly and
    asking for --species is the cheaper outcome.
    """
    cleaned = re.sub(r"^\S+\s+", "", defline)  # drop the leading accession token
    match = re.search(r"\b([A-Z][a-z]{2,})\s+([a-z]{3,})\b", cleaned)
    return f"{match.group(1)} {match.group(2)}" if match else None


def nonempty_facets(pairs: list[tuple[str, int]]) -> list[tuple[str, int]]:
    """Drop BV-BRC's placeholders for missing data, keeping the commonest-first order."""
    return [(value, count) for value, count in pairs or []
            if value and str(value).strip().lower() not in EMPTYISH]


def top_nonempty(pairs: list[tuple[str, int]]) -> tuple[str | None, int]:
    """The commonest facet value that is not a placeholder for missing data."""
    real = nonempty_facets(pairs)
    return real[0] if real else (None, 0)

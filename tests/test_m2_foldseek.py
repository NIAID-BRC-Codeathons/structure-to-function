"""Foldseek structure-search tests (issue #9). Payloads are trimmed real API responses."""

from __future__ import annotations

import json

import pytest

from s2f.common.http import JsonCache
from s2f.m2_triage.foldseek import (
    STATUS_FAILED,
    STATUS_FOUND,
    STATUS_NO_HIT,
    STATUS_NO_STRUCTURE,
    FoldseekClient,
    FoldseekHit,
    build_hits,
    parse_alignments,
    structure_key,
    summarize,
)

# Shape of a real response, including the coordinate blob we must not keep.
API_PAYLOAD = {
    "results": [
        {
            "db": "pdb100",
            "alignments": [
                [
                    {
                        "target": "3ia2-assembly2.cif.gz_D Pseudomonas fluorescens esterase complexed to a transition state analog",
                        "eval": 1.637e-12, "prob": 1, "seqId": 15.9, "score": 365, "alnLength": 200,
                        "qStartPos": 5, "qEndPos": 205, "qLen": 220,
                        "dbStartPos": 1, "dbEndPos": 198, "dbLen": 210,
                        "taxId": 294, "taxName": "Pseudomonas fluorescens",
                        "tCa": "x" * 5000, "tSeq": "y" * 500, "qAln": "z" * 200, "dbAln": "w" * 200,
                    },
                    {
                        "target": "1ybx-assembly1.cif.gz_A Conserved hypothetical protein Cth-383",
                        "eval": 3.4e-05, "prob": 0.9, "seqId": 19.7, "score": 120, "alnLength": 150,
                        "qStartPos": 10, "qEndPos": 160, "qLen": 220,
                        "dbStartPos": 1, "dbEndPos": 150, "dbLen": 160,
                        "taxId": 1515, "taxName": "Acetivibrio thermocellus",
                        "tCa": "x" * 5000,
                    },
                    {   # noise: looks like a named protein, e-value says otherwise
                        "target": "2co5-assembly1.cif.gz_A F93 winged-helix DNA-binding protein",
                        "eval": 1.894, "prob": 0.2, "seqId": 23.8, "score": 30, "alnLength": 40,
                        "qStartPos": 100, "qEndPos": 140, "qLen": 220,
                        "dbStartPos": 1, "dbEndPos": 40, "dbLen": 95,
                        "taxName": "Sulfolobus virus",
                    },
                ]
            ],
        }
    ]
}


def test_parse_alignments_flattens_and_drops_coordinate_blobs() -> None:
    pairs = parse_alignments(API_PAYLOAD)

    assert len(pairs) == 3
    database, alignment = pairs[0]
    assert database == "pdb100"
    # tCa/tSeq/qAln/dbAln are megabytes across a genome and useless downstream.
    for dropped in ("tCa", "tSeq", "qAln", "dbAln"):
        assert dropped not in alignment
    assert alignment["eval"] == 1.637e-12


def test_the_evalue_gate_rejects_plausible_looking_noise() -> None:
    """A hit at e-value 1.89 named a real protein and meant nothing."""
    hits = build_hits(parse_alignments(API_PAYLOAD))

    targets = [h.target for h in hits]
    assert len(hits) == 2
    assert not any("winged-helix" in t for t in targets)
    assert hits[0].evalue == 1.637e-12  # sorted best first


def test_coverage_gate_and_thresholds_are_adjustable() -> None:
    assert len(build_hits(parse_alignments(API_PAYLOAD), max_evalue=10, min_coverage=0.0)) == 3
    # The esterase hit spans 5-205 of 220 = 91%; a 95% floor drops everything.
    assert build_hits(parse_alignments(API_PAYLOAD), min_coverage=0.95) == []


def test_hit_fields_are_derived_from_the_target_string() -> None:
    hit = build_hits(parse_alignments(API_PAYLOAD))[0]

    assert hit.pdb_id == "3ia2"
    assert hit.description.startswith("Pseudomonas fluorescens esterase")
    assert hit.informative is True
    assert hit.query_coverage == pytest.approx(201 / 220, abs=1e-3)
    assert hit.target_coverage == pytest.approx(198 / 210, abs=1e-3)
    assert hit.taxonomy == "Pseudomonas fluorescens"


def test_a_hit_to_another_unknown_protein_is_marked_uninformative() -> None:
    """Matching a 'conserved hypothetical protein' advances nothing, and must say so."""
    hits = build_hits(parse_alignments(API_PAYLOAD))
    unknown = next(h for h in hits if "Cth-383" in h.target)

    assert unknown.evalue == 3.4e-05  # significant
    assert unknown.informative is False  # but useless


class FakeResponse:
    def __init__(self, payload, status_code=200):
        self.status_code = status_code
        self._payload = payload
        self.content = json.dumps(payload).encode()
        self.text = self.content.decode()

    def raise_for_status(self):
        if self.status_code >= 400:
            raise RuntimeError(f"HTTP {self.status_code}")

    def json(self):
        return self._payload


class TicketSession:
    """Mimics submit -> poll -> fetch, counting calls."""

    def __init__(self, *, statuses=("COMPLETE",), payload=None):
        self.statuses = list(statuses)
        self.payload = payload if payload is not None else API_PAYLOAD
        self.posts = 0
        self.gets: list[str] = []
        self.headers: dict[str, str] = {}

    def post(self, url, files=None, data=None, timeout=None):
        self.posts += 1
        return FakeResponse({"id": "TICKET1", "status": "PENDING"})

    def get(self, url, timeout=None, **kwargs):
        self.gets.append(url)
        if "/ticket/" in url:
            return FakeResponse({"status": self.statuses.pop(0) if self.statuses else "COMPLETE"})
        return FakeResponse(self.payload)


def make_client(tmp_path, session, **kwargs) -> FoldseekClient:
    return FoldseekClient(
        cache=JsonCache(tmp_path / "cache.sqlite"), session=session,
        poll_interval=0, sleep=lambda _s: None, **kwargs,
    )


def test_search_submits_polls_and_parses(tmp_path) -> None:
    session = TicketSession(statuses=("RUNNING", "COMPLETE"))
    client = make_client(tmp_path, session)

    result = client.search_protein("fig|1.1.peg.1", b"data_CIF", structure_source="afdb", structure_id="AF-X-F1")

    assert result.status == STATUS_FOUND
    assert session.posts == 1
    assert sum(1 for u in session.gets if "/ticket/" in u) == 2  # polled until COMPLETE
    assert result.best.pdb_id == "3ia2"
    assert result.best_informative.informative is True


def test_results_are_cached_by_structure(tmp_path) -> None:
    session = TicketSession()
    client = make_client(tmp_path, session)

    client.search_protein("fig|1.1.peg.1", b"data_CIF")
    client.search_protein("fig|1.1.peg.2", b"data_CIF")  # same structure bytes

    assert session.posts == 1  # second search served from cache


def test_offline_records_a_failure_rather_than_pretending_there_is_no_hit(tmp_path) -> None:
    client = make_client(tmp_path, TicketSession(), offline=True)

    result = client.search_protein("fig|1.1.peg.1", b"data_CIF")

    assert result.status == STATUS_FAILED
    assert result.status != STATUS_NO_HIT
    assert client.failures and client.failures[0]["feature_id"] == "fig|1.1.peg.1"


def test_a_protein_with_no_structure_is_not_queried(tmp_path) -> None:
    client = make_client(tmp_path, TicketSession())

    result = client.search_protein("fig|1.1.peg.9", None)

    assert result.status == STATUS_NO_STRUCTURE
    assert "needs a prediction" in result.note
    rows = result.as_annotations(retrieved_at="2026-09-17T12:00:00Z")
    assert rows[0]["source"] == "foldseek"
    assert rows[0]["hit"] is None


def test_empty_results_are_no_hit_not_found(tmp_path) -> None:
    client = make_client(tmp_path, TicketSession(payload={"results": []}))

    result = client.search_protein("fig|1.1.peg.1", b"data_CIF")

    assert result.status == STATUS_NO_HIT
    assert "no neighbour" in result.note


def test_annotation_rows_carry_source_and_the_tm_score_caveat(tmp_path) -> None:
    client = make_client(tmp_path, TicketSession())
    result = client.search_protein("fig|1.1.peg.1", b"data_CIF", structure_id="AF-X-F1")

    rows = result.as_annotations(retrieved_at="2026-09-17T12:00:00Z", limit=2)

    assert [r["source"] for r in rows] == ["foldseek", "foldseek"]
    assert rows[0]["identity"] == pytest.approx(0.159, abs=1e-3)
    assert rows[0]["tm_score"] is None  # the web API does not provide one
    assert "TM-score" in rows[0]["note"]
    assert rows[0]["query_structure"] == "AF-X-F1"
    assert rows[1]["informative_target"] is False


def test_structure_key_distinguishes_databases_and_structures() -> None:
    a = structure_key(b"x", ["pdb100"], "3diaa")
    assert a == structure_key(b"x", ["pdb100"], "3diaa")
    assert a != structure_key(b"y", ["pdb100"], "3diaa")
    assert a != structure_key(b"x", ["pdb100", "afdb50"], "3diaa")


def test_summary_separates_informative_from_merely_significant() -> None:
    from s2f.m2_triage.foldseek import FoldseekResult

    informative = FoldseekResult("a", status=STATUS_FOUND, hits=[FoldseekHit("1abc_A Real enzyme", evalue=1e-9, identity=18.0)])
    vague = FoldseekResult("b", status=STATUS_FOUND, hits=[FoldseekHit("2xyz_A hypothetical protein", evalue=1e-9, identity=20.0)])
    nothing = FoldseekResult("c", status=STATUS_NO_HIT)

    summary = summarize([informative, vague, nothing])

    assert summary["with_hits"] == 2
    assert summary["with_informative_hit"] == 1
    assert summary["uninformative_only"] == 1
    assert summary["status"] == {STATUS_FOUND: 2, STATUS_NO_HIT: 1}


# --- the public server's rate limit ---------------------------------------------------


class RateLimitedSession(TicketSession):
    """Refuses every submission the way the real server does: 429, no Retry-After."""

    def post(self, url, files=None, data=None, timeout=None):
        self.posts += 1
        return FakeResponse({"error": "Too Many Requests"}, status_code=429)


def test_rate_limiting_stops_the_run_with_an_actionable_message(tmp_path) -> None:
    """The 429 carries no Retry-After, so blind backoff is all there is — then stop."""
    session = RateLimitedSession()
    client = make_client(tmp_path, session)

    results = [client.search_protein(f"fig|1.1.peg.{i}", b"cif" + str(i).encode()) for i in range(6)]

    assert all(r.status == STATUS_FAILED for r in results)
    assert client.rate_limited is True
    # It gives up rather than hammering a shared server indefinitely.
    last = results[-1]
    assert "rate limiting" in last.error
    assert "local foldseek binary" in last.error


def test_a_success_resets_the_rate_limit_counter(tmp_path) -> None:
    class Flaky(TicketSession):
        def __init__(self):
            super().__init__()
            self.calls = 0

        def post(self, url, files=None, data=None, timeout=None):
            self.calls += 1
            self.posts += 1
            if self.calls == 1:
                return FakeResponse({"error": "Too Many Requests"}, status_code=429)
            return FakeResponse({"id": "TICKET1", "status": "PENDING"})

    client = make_client(tmp_path, Flaky())
    result = client.search_protein("fig|1.1.peg.1", b"cif")

    assert result.status == STATUS_FOUND  # retried past the 429
    assert client.rate_limited is False

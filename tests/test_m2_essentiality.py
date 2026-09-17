"""Essentiality-by-orthology tests (issue #29)."""

from __future__ import annotations

import json

import pytest

from s2f.common.http import CachedJsonClient, JsonCache
from s2f.m2_triage.essentiality import (
    STATUS_ESSENTIAL,
    STATUS_NO_CALL,
    EssentialityCall,
    ReferenceProtein,
    fetch_reference_set,
    parse_calls,
)

REFERENCE = {
    "fig|1.1.peg.1": ReferenceProtein("fig|1.1.peg.1", "Mycoplasma pirum ATCC 25960", "atpD",
                                      "ATP synthase beta chain", "MKALIV"),
    "fig|2.2.peg.2": ReferenceProtein("fig|2.2.peg.2", "Mycoplasma alvi ATCC 29626", "",
                                      "Pyruvate kinase", "MSTKNK"),
}

# qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore
ROWS = [
    # clean ortholog: 78% identity over 98% of the query
    ["fig|9.9.peg.1", "fig|1.1.peg.1", "78.0", "480", "490", "500", "5", "484", "1", "480", "1e-200", "700"],
    # identity fine, coverage too low (0.66 < 0.70)
    ["fig|9.9.peg.2", "fig|1.1.peg.1", "58.0", "100", "150", "500", "10", "109", "1", "100", "1e-30", "120"],
    # coverage fine, identity too low (38.4 < 40)
    ["fig|9.9.peg.3", "fig|2.2.peg.2", "38.4", "180", "190", "200", "5", "183", "1", "179", "1e-25", "110"],
    # two hits for one protein: the higher bitscore wins
    ["fig|9.9.peg.4", "fig|2.2.peg.2", "45.0", "190", "200", "200", "1", "195", "1", "190", "1e-40", "150"],
    ["fig|9.9.peg.4", "fig|1.1.peg.1", "80.0", "190", "200", "500", "1", "195", "1", "190", "1e-90", "300"],
]


def test_a_strong_ortholog_becomes_an_essentiality_call() -> None:
    calls = parse_calls(ROWS, REFERENCE)

    call = calls["fig|9.9.peg.1"]
    assert call.status == STATUS_ESSENTIAL
    assert call.essential is True
    assert call.identity == 78.0
    assert call.coverage == pytest.approx(480 / 490, abs=1e-3)
    # The DoD: a call must name its source and the identity behind it.
    assert call.reference_genome == "Mycoplasma pirum ATCC 25960"
    assert call.reference_gene == "atpD"
    assert call.reference_patric_id == "fig|1.1.peg.1"
    assert call.evidence == "FBA"


def test_partial_alignments_do_not_become_calls() -> None:
    calls = parse_calls(ROWS, REFERENCE)

    low_coverage = calls["fig|9.9.peg.2"]
    assert low_coverage.status == STATUS_NO_CALL
    assert low_coverage.essential is False
    assert low_coverage.identity == 58.0          # evidence is kept
    assert "below the 40%/70% transfer threshold" in low_coverage.note
    assert low_coverage.reference_genome == ""    # but no source is claimed

    low_identity = calls["fig|9.9.peg.3"]
    assert low_identity.status == STATUS_NO_CALL
    assert low_identity.coverage == pytest.approx(179 / 190, abs=1e-3)


def test_best_hit_wins_by_bitscore() -> None:
    call = parse_calls(ROWS, REFERENCE)["fig|9.9.peg.4"]

    assert call.identity == 80.0
    assert call.reference_patric_id == "fig|1.1.peg.1"


def test_thresholds_are_adjustable() -> None:
    relaxed = parse_calls(ROWS, REFERENCE, min_identity=30.0, min_coverage=0.6)
    assert relaxed["fig|9.9.peg.2"].essential is True
    assert relaxed["fig|9.9.peg.3"].essential is True

    strict = parse_calls(ROWS, REFERENCE, min_identity=90.0, min_coverage=0.9)
    assert strict["fig|9.9.peg.1"].essential is False


def test_annotation_row_states_the_inference_chain() -> None:
    """A homology transfer of a metabolic-model prediction must not read as an experiment."""
    call = parse_calls(ROWS, REFERENCE)["fig|9.9.peg.1"]

    row = call.as_annotation(retrieved_at="2026-09-17T12:00:00Z")

    assert row["source"] == "essentiality_ortholog"
    assert row["evidence"] == "FBA"
    assert "Mycoplasma pirum" in row["description"]
    assert "not an" in row["note"] and "experimental knockout" in row["note"]
    assert row["identity"] == pytest.approx(0.78, abs=1e-3)


def test_a_protein_with_no_ortholog_is_an_explicit_no_call() -> None:
    call = EssentialityCall(feature_id="fig|9.9.peg.8", status=STATUS_NO_CALL, note="no ortholog")
    row = call.as_annotation(retrieved_at="2026-09-17T12:00:00Z")

    assert row["hit"] is None
    assert row["description"] == "no essential-gene ortholog"
    assert call.identity is None  # absence, not zero


class RecordingSession:
    """Replays canned BV-BRC responses and records the RQL that was sent."""

    def __init__(self) -> None:
        self.urls: list[str] = []
        self.headers: dict[str, str] = {}

    def get(self, url, params=None, timeout=None):
        self.urls.append(url)

        class R:
            status_code = 200
            headers: dict[str, str] = {}

            def __init__(self, payload):
                self.content = json.dumps(payload).encode()
                self.text = self.content.decode()
                self._payload = payload

            @property
            def ok(self):
                return True

            def json(self):
                return self._payload

        if "/sp_gene/" in url:
            return R([{"patric_id": "fig|1.1.peg.1", "genome_id": "1.1", "genome_name": "Mycoplasma pirum",
                       "gene": "atpD", "product": "ATP synthase beta chain"}])
        if "/genome_feature/" in url:
            return R([{"patric_id": "fig|1.1.peg.1", "aa_sequence_md5": "abc123"}])
        if "/feature_sequence/" in url:
            return R([{"md5": "abc123", "sequence": "MKALIVMKALIV"}])
        raise AssertionError(f"unexpected url {url}")

    def post(self, url, json=None, timeout=None):  # pragma: no cover
        raise AssertionError("essentiality should not POST")


def test_reference_set_is_fetched_and_joined_to_sequences(tmp_path) -> None:
    session = RecordingSession()
    client = CachedJsonClient(
        cache=JsonCache(tmp_path / "cache.sqlite"), session=session,
        min_interval_seconds=0, sleep=lambda _s: None,
    )

    reference = fetch_reference_set(client, keyword="Mycoplasma", limit=10)

    assert len(reference) == 1
    assert reference[0].sequence == "MKALIVMKALIV"
    assert reference[0].genome_name == "Mycoplasma pirum"
    # BV-BRC's RQL: FBA evidence, keyword-selected relatives, in three joined calls.
    assert any("eq(evidence,FBA)" in u and "keyword(Mycoplasma)" in u for u in session.urls)
    assert any("/genome_feature/" in u and "in(patric_id," in u for u in session.urls)
    assert any("/feature_sequence/" in u and "in(md5," in u for u in session.urls)

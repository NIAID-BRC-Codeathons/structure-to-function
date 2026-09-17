"""AlphaFold DB lookup tests (issue #8). Payloads are trimmed real API responses."""

from __future__ import annotations

import json

from s2f.common.http import CachedJsonClient, JsonCache
from s2f.m2_triage.afdb import (
    STATUS_FAILED,
    STATUS_FOUND,
    STATUS_INVALID,
    STATUS_NO_MODEL,
    STATUS_NOT_QUERIED,
    AlphaFoldClient,
    AlphaFoldModel,
    confidence_band,
    parse_prediction,
)

PREDICTION = [
    {
        "entryId": "AF-P47641-F1",
        "latestVersion": 6,
        "globalMetricValue": 92.5,
        "uniprotStart": 1,
        "uniprotEnd": 518,
        "uniprotSequence": "M" * 518,
        "cifUrl": "https://alphafold.ebi.ac.uk/files/AF-P47641-F1-model_v6.cif",
        "pdbUrl": "https://alphafold.ebi.ac.uk/files/AF-P47641-F1-model_v6.pdb",
        "paeDocUrl": "https://alphafold.ebi.ac.uk/files/AF-P47641-F1-predicted_aligned_error_v6.json",
        "modelCreatedDate": "2025-08-01T00:00:00Z",
    }
]

# Real shape of a long protein: the model covers a fragment, not the whole accession.
FRAGMENT_PREDICTION = [
    {
        "entryId": "AF-0000000365840311",
        "latestVersion": 1,
        "globalMetricValue": 96.08,
        "uniprotStart": 1368,
        "uniprotEnd": 1493,
        "uniprotSequence": "M" * 7096,
        "cifUrl": "https://example.invalid/frag.cif",
    }
]


class FakeResponse:
    def __init__(self, payload, status_code: int = 200) -> None:
        self.status_code = status_code
        self._payload = payload
        self.content = b"" if payload is None else json.dumps(payload).encode()
        self.text = self.content.decode()
        self.headers: dict[str, str] = {}

    @property
    def ok(self) -> bool:
        return 200 <= self.status_code < 300

    def json(self):
        return self._payload


class AfdbSession:
    def __init__(self, responses: dict[str, FakeResponse]) -> None:
        self.responses = responses
        self.calls: list[str] = []
        self.headers: dict[str, str] = {}

    def get(self, url, params=None, timeout=None):
        self.calls.append(url)
        accession = url.rstrip("/").split("/")[-1]
        return self.responses.get(accession, FakeResponse({}, status_code=404))

    def post(self, url, json=None, timeout=None):  # pragma: no cover
        raise AssertionError("afdb client should not POST")


def make_client(tmp_path, session) -> AlphaFoldClient:
    return AlphaFoldClient(
        CachedJsonClient(
            cache=JsonCache(tmp_path / "cache.sqlite"),
            session=session,
            min_interval_seconds=0,
            sleep=lambda _s: None,
        )
    )


def test_confidence_bands_follow_alphafolds_own_thresholds() -> None:
    assert confidence_band(95.0) == "very high"
    assert confidence_band(90.0) == "very high"
    assert confidence_band(70.0) == "confident"
    assert confidence_band(69.9) == "low"
    assert confidence_band(49.9) == "very low"
    assert confidence_band(None) == ""


def test_parse_prediction_reads_a_full_length_model() -> None:
    model = parse_prediction("P47641", PREDICTION)

    assert model.status == STATUS_FOUND
    assert model.found is True
    assert model.entry_id == "AF-P47641-F1"
    assert model.mean_plddt == 92.5
    assert model.confidence_band == "very high"
    assert model.covered_residues == 518
    assert model.coverage == 1.0
    assert model.cif_url.endswith(".cif")


def test_a_fragment_model_reports_its_real_coverage() -> None:
    """A model of 2% of a protein is not a model of that protein."""
    model = parse_prediction("P0DTD1", FRAGMENT_PREDICTION)

    assert model.found is True
    assert model.mean_plddt == 96.08  # high confidence...
    assert model.covered_residues == 126
    assert model.coverage == round(126 / 7096, 4)
    usable, reason = model.usable_for_docking()
    assert usable is False  # ...but it covers almost nothing
    assert "1368-1493" in reason


def test_http_statuses_become_distinct_explicit_answers() -> None:
    assert parse_prediction("Q8WZ42", {"_http_status": 404}).status == STATUS_NO_MODEL
    assert parse_prediction("NOTANACC", {"_http_status": 400}).status == STATUS_INVALID
    assert parse_prediction("P00000", {"_http_status": 503}).status == STATUS_FAILED
    assert parse_prediction("P00000", []).status == STATUS_NO_MODEL


def test_the_highest_confidence_fragment_wins() -> None:
    entries = [
        dict(PREDICTION[0], entryId="AF-X-F1", globalMetricValue=55.0),
        dict(PREDICTION[0], entryId="AF-X-F2", globalMetricValue=88.0),
    ]
    model = parse_prediction("X", entries)

    assert model.entry_id == "AF-X-F2"
    assert model.fragments == 2


def test_usable_for_docking_gates_on_confidence_and_coverage() -> None:
    good = parse_prediction("P47641", PREDICTION)
    assert good.usable_for_docking()[0] is True

    weak = parse_prediction("X", [dict(PREDICTION[0], globalMetricValue=56.47)])
    usable, reason = weak.usable_for_docking()
    assert usable is False
    assert "low" in reason

    missing = AlphaFoldModel(accession="X", status=STATUS_NO_MODEL)
    usable, reason = missing.usable_for_docking()
    assert usable is False
    assert "no AlphaFold model" in reason


def test_lookup_treats_404_as_an_answer_not_a_failure(tmp_path) -> None:
    session = AfdbSession({"Q8WZ42": FakeResponse({}, status_code=404)})
    client = make_client(tmp_path, session)

    model = client.lookup("Q8WZ42")

    assert model.status == STATUS_NO_MODEL
    assert client.failures == []  # a missing model is not a broken query


def test_lookup_records_a_real_failure(tmp_path) -> None:
    session = AfdbSession({"P00001": FakeResponse({"error": "boom"}, status_code=500)})
    client = make_client(tmp_path, session)

    model = client.lookup("P00001")

    assert model.status == STATUS_FAILED
    assert client.failures and client.failures[0]["accession"] == "P00001"


def test_lookup_without_an_accession_is_not_queried(tmp_path) -> None:
    client = make_client(tmp_path, AfdbSession({}))

    model = client.lookup("")

    assert model.status == STATUS_NOT_QUERIED
    assert "no UniProt accession" in model.error


def test_lookup_many_deduplicates_and_caches(tmp_path) -> None:
    session = AfdbSession({"P47641": FakeResponse(PREDICTION)})
    client = make_client(tmp_path, session)

    models = client.lookup_many(["P47641", "P47641", ""], workers=1)

    assert set(models) == {"P47641"}
    assert models["P47641"].found is True
    assert len(session.calls) == 1

    client.lookup("P47641")
    assert len(session.calls) == 1  # served from cache

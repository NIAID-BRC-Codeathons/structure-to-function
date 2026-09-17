"""RCSB client tests. Response shapes are trimmed copies of real API responses."""

from __future__ import annotations

import json

import pytest
import requests

from s2f.common.http import CachedJsonClient, HttpError, JsonCache, OfflineCacheMiss
from s2f.m2_triage.pdb_evidence import (
    PdbEvidenceClient,
    SequenceHit,
    apply_metadata,
    classify_components,
    parse_entity_metadata,
    parse_search_response,
    search_payload,
)

SEARCH_RESPONSE = {
    "query_id": "test",
    "result_type": "polymer_entity",
    "total_count": 2,
    "result_set": [
        {
            "identifier": "6M0J_2",
            "score": 1.0,
            "services": [
                {
                    "service_type": "sequence",
                    "nodes": [
                        {
                            "match_context": [
                                {
                                    "sequence_identity": 1.0,
                                    "evalue": 2.869e-155,
                                    "bitscore": 484,
                                    "alignment_length": 223,
                                    "query_beg": 1,
                                    "query_end": 223,
                                    "query_length": 223,
                                    "subject_length": 229,
                                }
                            ]
                        }
                    ],
                }
            ],
        },
        {
            "identifier": "7A94_1",
            "score": 0.5,
            "services": [
                {
                    "service_type": "sequence",
                    "nodes": [
                        {
                            "match_context": [
                                {
                                    "sequence_identity": 0.4,
                                    "evalue": 1e-20,
                                    "bitscore": 90,
                                    "query_beg": 10,
                                    "query_end": 120,
                                    "query_length": 223,
                                    "subject_length": 200,
                                },
                                {
                                    "sequence_identity": 0.9,
                                    "evalue": 1e-30,
                                    "bitscore": 180,
                                    "query_beg": 1,
                                    "query_end": 200,
                                    "query_length": 223,
                                    "subject_length": 200,
                                },
                            ]
                        }
                    ],
                }
            ],
        },
    ],
}

GRAPHQL_RESPONSE = {
    "data": {
        "polymer_entities": [
            {
                "rcsb_id": "6M0J_2",
                "rcsb_polymer_entity": {"pdbx_description": "Spike protein S1", "rcsb_ec_lineage": None},
                "rcsb_polymer_entity_container_identifiers": {"uniprot_ids": ["P0DTC2"]},
                "rcsb_entity_source_organism": [
                    {"ncbi_scientific_name": "Severe acute respiratory syndrome coronavirus 2", "ncbi_taxonomy_id": 2697049}
                ],
                "rcsb_polymer_entity_annotation": [{"type": "GO"}, {"type": "Pfam"}],
                "entry": {
                    "rcsb_entry_info": {"resolution_combined": [2.45], "experimental_method": "X-ray"},
                    "nonpolymer_entities": [
                        {"pdbx_entity_nonpoly": {"comp_id": "NAG"}},
                        {"pdbx_entity_nonpoly": {"comp_id": "ZN"}},
                        {"pdbx_entity_nonpoly": {"comp_id": "CL"}},
                    ],
                    "polymer_entities": [{"rcsb_id": "6M0J_1"}, {"rcsb_id": "6M0J_2"}],
                },
            },
            None,
        ]
    }
}


class FakeResponse:
    def __init__(self, status_code: int, payload=None, headers=None) -> None:
        self.status_code = status_code
        self._payload = payload
        self.headers = headers or {}
        self.content = b"" if payload is None else json.dumps(payload).encode()
        self.text = self.content.decode()

    @property
    def ok(self) -> bool:
        return 200 <= self.status_code < 300

    def json(self):
        return self._payload


class FakeSession:
    """Stands in for requests.Session, returning queued responses."""

    def __init__(self, responses) -> None:
        self.responses = list(responses)
        self.calls = []
        self.headers: dict[str, str] = {}

    def post(self, url, json=None, timeout=None):  # noqa: A002 - mirrors requests' signature
        self.calls.append((url, json))
        response = self.responses.pop(0)
        if isinstance(response, Exception):
            raise response
        return response


def make_client(tmp_path, responses, **kwargs) -> tuple[CachedJsonClient, FakeSession]:
    session = FakeSession(responses)
    client = CachedJsonClient(
        cache=JsonCache(tmp_path / "cache.sqlite"),
        session=session,
        min_interval_seconds=0,
        sleep=lambda _seconds: None,
        **kwargs,
    )
    return client, session


def test_search_payload_carries_documented_parameters() -> None:
    payload = search_payload("MKALIV")
    parameters = payload["query"]["parameters"]

    assert payload["return_type"] == "polymer_entity"
    assert parameters["evalue_cutoff"] == 1e-3
    assert parameters["identity_cutoff"] == 0.25
    assert payload["request_options"]["results_verbosity"] == "verbose"


def test_parse_search_response_keeps_the_best_match_context() -> None:
    hits = parse_search_response(SEARCH_RESPONSE)

    assert [hit.entity_id for hit in hits] == ["6M0J_2", "7A94_1"]
    assert hits[0].query_coverage == 1.0
    # The second result has two contexts; the higher-bitscore one wins.
    assert hits[1].identity == 0.9
    assert hits[1].bitscore == 180


def test_classify_components_splits_ligands_metals_and_additives() -> None:
    ligands, metals = classify_components(["MK1", "ZN", "HOH", "GOL", "NAG", "", "mg"])

    assert ligands == ["MK1"]
    assert metals == ["ZN", "MG"]


def test_parse_entity_metadata_maps_by_rcsb_id() -> None:
    metadata = parse_entity_metadata(GRAPHQL_RESPONSE)

    assert set(metadata) == {"6M0J_2"}
    entity = metadata["6M0J_2"]
    assert entity.resolution == 2.45
    assert entity.experimental_method == "X-ray"
    assert entity.organism.startswith("Severe acute")
    assert entity.uniprot_ids == ["P0DTC2"]
    assert entity.ligands == []  # NAG is a glycan, ZN a metal, CL an additive
    assert entity.metals == ["ZN"]
    assert entity.partner_entities == 1
    assert entity.has_go is True
    assert entity.has_ec is False


def test_apply_metadata_leaves_unmatched_hits_alone() -> None:
    hits = [
        SequenceHit("6M0J_2", 1.0, 1e-100, 400, 1, 223, 223, 229),
        SequenceHit("9ZZZ_1", 0.5, 1e-10, 80, 1, 100, 223, 120),
    ]
    apply_metadata(hits, parse_entity_metadata(GRAPHQL_RESPONSE))

    assert hits[0].resolution == 2.45
    assert hits[1].resolution is None
    assert hits[1].experimental_method == ""


def test_search_reports_no_hit_for_204_and_found_for_results(tmp_path) -> None:
    client, _ = make_client(tmp_path, [FakeResponse(204), FakeResponse(200, SEARCH_RESPONSE)])
    pdb = PdbEvidenceClient(client)

    assert pdb.search_sequence("MKALIV").status == "no-hit"
    found = pdb.search_sequence("MKALIW")
    assert found.status == "found"
    assert len(found.hits) == 2


def test_search_failure_is_query_failed_not_no_hit(tmp_path) -> None:
    client, _ = make_client(tmp_path, [FakeResponse(400, {"error": "bad request"})])
    pdb = PdbEvidenceClient(client)

    result = pdb.search_sequence("MKALIV")
    assert result.status == "query-failed"
    assert "400" in result.error
    assert result.hits == []


def test_retries_429_then_succeeds(tmp_path) -> None:
    client, session = make_client(
        tmp_path,
        [
            FakeResponse(429, {"error": "slow down"}, headers={"Retry-After": "0"}),
            FakeResponse(500, {"error": "boom"}),
            FakeResponse(200, SEARCH_RESPONSE),
        ],
    )
    pdb = PdbEvidenceClient(client)

    assert pdb.search_sequence("MKALIV").status == "found"
    assert len(session.calls) == 3


def test_gives_up_after_max_attempts(tmp_path) -> None:
    client, session = make_client(tmp_path, [FakeResponse(503, {}) for _ in range(4)], max_attempts=4)

    with pytest.raises(HttpError):
        client.post_json("rcsb_sequence", "https://example.invalid", {"q": 1})
    assert len(session.calls) == 4


def test_network_errors_are_retried(tmp_path) -> None:
    client, session = make_client(
        tmp_path,
        [requests.ConnectionError("reset"), FakeResponse(200, SEARCH_RESPONSE)],
    )

    payload = client.post_json("rcsb_sequence", "https://example.invalid", {"q": 1})
    assert payload["total_count"] == 2
    assert len(session.calls) == 2


def test_cache_prevents_a_second_request(tmp_path) -> None:
    client, session = make_client(tmp_path, [FakeResponse(200, SEARCH_RESPONSE)])

    first = client.post_json("rcsb_sequence", "https://example.invalid", {"q": 1})
    second = client.post_json("rcsb_sequence", "https://example.invalid", {"q": 1})

    assert first == second
    assert len(session.calls) == 1


def test_offline_mode_never_calls_the_network(tmp_path) -> None:
    warm, session = make_client(tmp_path, [FakeResponse(200, SEARCH_RESPONSE)])
    warm.post_json("rcsb_sequence", "https://example.invalid", {"q": 1})

    offline = CachedJsonClient(
        cache=JsonCache(tmp_path / "cache.sqlite"),
        session=FakeSession([]),
        offline=True,
        min_interval_seconds=0,
    )
    assert offline.post_json("rcsb_sequence", "https://example.invalid", {"q": 1})["total_count"] == 2
    with pytest.raises(OfflineCacheMiss):
        offline.post_json("rcsb_sequence", "https://example.invalid", {"q": 2})


def test_metadata_batch_failures_are_recorded(tmp_path) -> None:
    client, _ = make_client(tmp_path, [FakeResponse(500, {}) for _ in range(4)])
    pdb = PdbEvidenceClient(client)

    metadata = pdb.fetch_entity_metadata(["6M0J_2"])
    assert metadata == {}
    assert pdb.metadata_failures and "6M0J_2" in pdb.metadata_failures[0]["entity_ids"]


def test_graphql_errors_are_recorded_even_when_data_returns(tmp_path) -> None:
    payload = dict(GRAPHQL_RESPONSE, errors=[{"message": "FieldUndefined"}])
    client, _ = make_client(tmp_path, [FakeResponse(200, payload)])
    pdb = PdbEvidenceClient(client)

    metadata = pdb.fetch_entity_metadata(["6M0J_2"])
    assert "6M0J_2" in metadata
    assert pdb.metadata_failures

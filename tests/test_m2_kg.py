"""Knowledge-graph assembly tests (issue #11). Payloads are trimmed real API responses."""

from __future__ import annotations

import json

from s2f.common.http import CachedJsonClient, JsonCache
from s2f.m2_triage.kg import (
    EDGE_CHEMBL_TARGET,
    EDGE_HOMOLOG,
    EDGE_INTERACTS,
    EDGE_TARGET_OF,
    KnowledgeGraphBuilder,
    protein_node_id,
)

FEATURE_ID = "fig|243273.25.peg.405"

STRING_PARTNERS = [
    {"stringId_B": "243273.MG_218", "preferredName_B": "MG_218", "score": 0.782, "escore": 0, "dscore": 0.4, "tscore": 0.6},
    {"stringId_B": "243273.MG_219", "preferredName_B": "MG_219", "score": 0.9, "escore": 0.5, "dscore": 0, "tscore": 0.3},
    {"stringId_B": "243273.MG_220", "preferredName_B": "MG_220", "score": 0.71, "escore": 0, "dscore": 0.2, "tscore": 0.5},
]

CHEMBL_TARGETS = {
    "page_meta": {"total_count": 1},
    "targets": [
        {
            "target_chembl_id": "CHEMBL1982",
            "pref_name": "Isoleucine--tRNA ligase",
            "organism": "Staphylococcus aureus",
            "target_type": "SINGLE PROTEIN",
            "target_components": [{"accession": "P41972"}],
        }
    ],
}

CHEMBL_ACTIVITIES = {
    "page_meta": {"total_count": 2},
    "activities": [
        {
            "molecule_chembl_id": "CHEMBL538163",
            "standard_type": "IC50",
            "standard_value": "300.0",
            "standard_units": "nM",
            "pchembl_value": "6.52",
            "assay_chembl_id": "CHEMBL880123",
        },
        {
            "molecule_chembl_id": "CHEMBL538164",
            "standard_type": "MIC",
            "standard_value": "1000.0",
            "standard_units": "nM",
            "pchembl_value": "6.0",
            "assay_chembl_id": "CHEMBL880124",
        },
    ],
}


class FakeResponse:
    def __init__(self, payload, status_code: int = 200) -> None:
        self.status_code = status_code
        self._payload = payload
        self.content = json.dumps(payload).encode()
        self.text = self.content.decode()
        self.headers: dict[str, str] = {}

    @property
    def ok(self) -> bool:
        return 200 <= self.status_code < 300

    def json(self):
        return self._payload


class KgSession:
    def __init__(self, *, partners=None, targets=None, activities=None, fail=()) -> None:
        self.partners = STRING_PARTNERS if partners is None else partners
        self.targets = CHEMBL_TARGETS if targets is None else targets
        self.activities = CHEMBL_ACTIVITIES if activities is None else activities
        self.fail = set(fail)
        self.calls: list[tuple[str, dict]] = []
        self.headers: dict[str, str] = {}

    def get(self, url, params=None, timeout=None):
        params = params or {}
        self.calls.append((url, params))
        if "version" in url:
            return FakeResponse([{"string_version": "12.0"}])
        if "interaction_partners" in url:
            if "string" in self.fail:
                return FakeResponse({"error": "boom"}, status_code=500)
            return FakeResponse(self.partners)
        if "target.json" in url:
            if "chembl_target" in self.fail:
                return FakeResponse({"error": "boom"}, status_code=500)
            return FakeResponse(self.targets)
        if "activity.json" in url:
            return FakeResponse(self.activities)
        raise AssertionError(f"unexpected URL {url}")

    def post(self, url, json=None, timeout=None):  # pragma: no cover
        raise AssertionError("kg builder should not POST")


def make_builder(tmp_path, session, **kwargs) -> KnowledgeGraphBuilder:
    client = CachedJsonClient(
        cache=JsonCache(tmp_path / "cache.sqlite"),
        session=session,
        min_interval_seconds=0,
        sleep=lambda _s: None,
    )
    return KnowledgeGraphBuilder(client, taxon_id=243273, **kwargs)


def test_string_partners_become_capped_interacts_edges(tmp_path) -> None:
    builder = make_builder(tmp_path, KgSession(), max_partners=2)

    added = builder.add_string_partners(FEATURE_ID, "MG_405")

    assert added == 2  # cap applied
    edges = [e for e in builder.graph.edges if e.type == EDGE_INTERACTS]
    assert len(edges) == 2
    assert edges[0].evidence["combined_score"] == 0.782
    assert "not a measured interaction" in edges[0].evidence["note"]
    assert edges[0].provenance["source"] == "STRING"
    assert edges[0].provenance["source_version"] == "12.0"
    assert edges[0].provenance["retrieved_at"]


def test_homolog_edge_carries_the_identity_that_justifies_it(tmp_path) -> None:
    builder = make_builder(tmp_path, KgSession())

    builder.add_homolog(
        FEATURE_ID, "P41972", identity=37.6, coverage=0.99,
        organism="Staphylococcus aureus", via="PDB sequence hit 1FFY_2",
    )

    edge = next(e for e in builder.graph.edges if e.type == EDGE_HOMOLOG)
    assert edge.source_id == protein_node_id(FEATURE_ID)
    assert edge.target_id == "uniprot:P41972"
    assert edge.evidence["percent_identity"] == 37.6
    assert edge.evidence["organism"] == "Staphylococcus aureus"
    assert "not this protein's function" in edge.evidence["note"]


def test_human_homolog_edge_records_identity_without_inventing_an_accession(tmp_path) -> None:
    builder = make_builder(tmp_path, KgSession())

    node_id = builder.add_named_human_homolog(FEATURE_ID, "SSU ribosomal protein S12p", 50.0)

    node = builder.graph.nodes[node_id]
    assert node.attributes["accession"] is None
    assert node.attributes["organism"] == "Homo sapiens"
    edge = next(e for e in builder.graph.edges if e.type == EDGE_HOMOLOG)
    assert edge.evidence["percent_identity"] == 50.0
    assert "never a bonus" in edge.evidence["note"]


def test_chembl_bridges_compounds_to_the_homolog_not_to_us(tmp_path) -> None:
    builder = make_builder(tmp_path, KgSession())

    added = builder.add_chembl_for_accession("P41972")

    assert added == 2
    represented = next(e for e in builder.graph.edges if e.type == EDGE_CHEMBL_TARGET)
    assert represented.source_id == "uniprot:P41972"
    assert represented.target_id == "chembl:CHEMBL1982"

    target_edges = [e for e in builder.graph.edges if e.type == EDGE_TARGET_OF]
    assert {e.source_id for e in target_edges} == {"compound:CHEMBL538163", "compound:CHEMBL538164"}
    # Every compound edge points at the ChEMBL target, never at our protein directly.
    assert all(e.target_id == "chembl:CHEMBL1982" for e in target_edges)
    assert all(e.source_id != protein_node_id(FEATURE_ID) for e in target_edges)
    assert target_edges[0].evidence["pchembl_value"] == "6.52"
    assert "homolog of our protein" in target_edges[0].evidence["note"]


def test_chembl_target_is_rejected_when_the_filter_did_not_apply(tmp_path) -> None:
    """ChEMBL silently ignores some filter names; a target that does not list the accession is dropped."""
    wrong = {
        "page_meta": {"total_count": 1},
        "targets": [
            {
                "target_chembl_id": "CHEMBL999",
                "pref_name": "Maltase-glucoamylase",
                "organism": "Homo sapiens",
                "target_components": [{"accession": "O43451"}],
            }
        ],
    }
    builder = make_builder(tmp_path, KgSession(targets=wrong))

    added = builder.add_chembl_for_accession("P41972")

    assert added == 0
    assert not builder.graph.edges
    assert builder.graph.failures
    assert "filter did not apply" in builder.graph.failures[0]["error"]


def test_source_failures_are_recorded_not_silently_dropped(tmp_path) -> None:
    builder = make_builder(tmp_path, KgSession(fail={"string", "chembl_target"}))

    assert builder.add_string_partners(FEATURE_ID, "MG_405") == 0
    assert builder.add_chembl_for_accession("P41972") == 0

    steps = {f["step"] for f in builder.graph.failures}
    assert steps == {"interaction_partners", "target"}


def test_build_assembles_the_whole_subgraph(tmp_path) -> None:
    builder = make_builder(tmp_path, KgSession(), max_partners=2)

    graph = builder.build(
        [
            {
                "feature_id": FEATURE_ID,
                "locus_tag": "MG_405",
                "product": "Isoleucyl-tRNA synthetase",
                "homologs": [
                    {
                        "accession": "P41972",
                        "identity": 37.6,
                        "coverage": 0.99,
                        "organism": "Staphylococcus aureus",
                        "via": "PDB sequence hit 1FFY_2",
                    }
                ],
                "human_homolog_identity": 50.0,
            }
        ]
    )

    summary = graph.summary()
    assert summary["edge_types"] == {
        EDGE_HOMOLOG: 2,  # the S. aureus homolog and the named human one
        EDGE_INTERACTS: 2,
        EDGE_CHEMBL_TARGET: 1,
        EDGE_TARGET_OF: 2,
    }
    # Every edge is attributable.
    assert all(e.provenance.get("source") and e.provenance.get("retrieved_at") for e in graph.edges)

    payload = graph.to_dict(caps=builder.caps())
    assert payload["caps"]["max_string_partners_per_protein"] == 2
    assert payload["source_versions"]["STRING"] == "12.0"
    assert "no knowledge graph is ingested or hosted" in payload["note"]


def test_graph_converts_to_networkx(tmp_path) -> None:
    builder = make_builder(tmp_path, KgSession(), max_partners=1)
    graph = builder.build(
        [{"feature_id": FEATURE_ID, "locus_tag": "MG_405", "product": "x", "homologs": []}]
    )

    nx_graph = graph.to_networkx()

    assert nx_graph.number_of_nodes() == len(graph.nodes)
    assert nx_graph.number_of_edges() == len(graph.edges)
    assert nx_graph.nodes[protein_node_id(FEATURE_ID)]["type"] == "protein"

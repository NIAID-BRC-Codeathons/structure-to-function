"""Identifier mapping tests (issue #3). Response shapes are trimmed real UniProt responses."""

from __future__ import annotations

import json

import pytest

from s2f.common.http import CachedJsonClient, JsonCache
from s2f.common.ids import (
    ROUTE_LOCUS,
    ROUTE_PROTEOME_LOCUS,
    ROUTE_PROTEOME_SEQUENCE,
    ROUTE_PATRIC,
    ROUTE_PDB_HIT,
    ROUTE_SEQUENCE,
    ROUTE_UNMAPPED,
    IdMapper,
    Xrefs,
    build_taxon_index,
    _select_uniprot_crossref,
    crc64,
)

FEATURE_ID = "fig|243273.25.peg.468"
# Real P47641 sequence, so the checksum test is a genuine round-trip against UniProt.
ATPA = "MADKLNEYVALIKTEIKKYSKKIFNSEIGQVISVADGIAKVSGLENALLNELIQFENNIQGIVLNLEQNTVGIALFGDYSSLREGSTAKRTHSVMKTPVGDVMLGRIVNALGEAIDGRGDIKATEYDQIEKIAPGVMKRKSVNQPLETGILTIDALFPIGKGQRELIVGDRQTGKTAIAIDTIINQKDKDVYCVYVAIGQKNSSVAQIVHQLEVNDSMKYTTVVCATASDSDSMVYLSPFTGITIAEYWLKKGKDVLIVFDDLSKHAVAYRTLSLLLKRPPGREAFPGDVFYLHSRLLERACKLNDENGGGSITALPIIETQAGDISAYIPTNVISITDGQLFMVSSLFNAGQRPAIQIGLSVSRVGSAAQTKAIKQQTGSLKLELAQYSELDSFSQFGSDLDENTKKVLEHGKRVMEMIKQPNGKPYSQVHEALFLFAINKAFIKFIPVDEIAKFKQRITEEFNGSHPLFKELSNKKEFTEDLESKTKTAFKMLVKRFISTLTDYDITKFGSIEELN"

UNIPARC_ENTRY = {
    "uniParcId": "UPI00001263F8",
    "uniParcCrossReferences": [
        {
            "database": "UniProtKB/Swiss-Prot",
            "id": "P47641",
            "active": True,
            "geneName": "atpA",
            "organism": {"scientificName": "Mycoplasmoides genitalium G37", "taxonId": 243273},
        },
        {
            "database": "UniProtKB/TrEMBL",
            "id": "J3TUZ1",
            "active": False,
            "organism": {"taxonId": 662946},
        },
        {"database": "PATRIC", "id": FEATURE_ID, "active": True, "organism": {"taxonId": 243273}},
        {"database": "RefSeq", "id": "WP_009885622", "active": True, "organism": {"taxonId": 2097}},
    ],
}

UNIPROTKB_ENTRY = {
    "results": [
        {
            "primaryAccession": "P47641",
            "uniProtkbId": "ATPA_MYCGE",
            "proteinDescription": {"recommendedName": {"fullName": {"value": "ATP synthase subunit alpha"}}},
            "genes": [{"geneName": {"value": "atpA"}, "orderedLocusNames": [{"value": "MG401"}]}],
            "organism": {"taxonId": 243273},
            "uniProtKBCrossReferences": [
                {"database": "RefSeq", "id": "WP_009885622.1"},
                {"database": "PDB", "id": "1ABC"},
                {"database": "ChEMBL", "id": "CHEMBL1234"},
            ],
        }
    ]
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


class RoutingSession:
    """Answers by URL and query, the way the real endpoints split up."""

    def __init__(self, *, uniparc_hits=None, uniprotkb=None, stream=None) -> None:
        self.uniparc_hits = uniparc_hits or {}
        self.uniprotkb = uniprotkb if uniprotkb is not None else UNIPROTKB_ENTRY
        self.stream = stream if stream is not None else {"results": []}
        self.calls: list[tuple[str, dict]] = []
        self.headers: dict[str, str] = {}

    def get(self, url, params=None, timeout=None):
        params = params or {}
        self.calls.append((url, params))
        if "uniprotkb/stream" in url:
            return FakeResponse(self.stream)
        if "uniparc/search" in url:
            query = params.get("query", "")
            for needle, upi in self.uniparc_hits.items():
                if needle in query:
                    return FakeResponse({"results": [{"uniParcId": upi}]})
            return FakeResponse({"results": []})
        if "uniparc/" in url:
            return FakeResponse(UNIPARC_ENTRY)
        if "uniprotkb/search" in url:
            return FakeResponse(self.uniprotkb)
        raise AssertionError(f"unexpected URL {url}")

    def post(self, url, json=None, timeout=None):  # pragma: no cover - mapper only GETs
        raise AssertionError("IdMapper should not POST")


def make_mapper(tmp_path, session, **kwargs) -> IdMapper:
    client = CachedJsonClient(
        cache=JsonCache(tmp_path / "cache.sqlite"),
        session=session,
        min_interval_seconds=0,
        sleep=lambda _s: None,
    )
    return IdMapper(client, **kwargs)


def test_crc64_matches_uniprots_own_checksum() -> None:
    # UniProt publishes sequence.crc64; this is the value from P47641.
    assert crc64(ATPA) == "2CBBDAF094BCDEEB"
    assert crc64("") == "0000000000000000"


def test_sequence_route_maps_and_enriches(tmp_path) -> None:
    session = RoutingSession(uniparc_hits={"checksum:2CBBDAF094BCDEEB": "UPI00001263F8"})
    mapper = make_mapper(tmp_path, session, taxon_id=243273)

    xrefs = mapper.map_protein(FEATURE_ID, ATPA, locus_tag="MG_401")

    assert xrefs.uniprot == "P47641"
    assert xrefs.route == ROUTE_SEQUENCE
    assert xrefs.mapped is True
    assert xrefs.uniparc == "UPI00001263F8"
    assert xrefs.gene_name == "atpA"
    assert xrefs.uniprot_entry_name == "ATPA_MYCGE"
    assert xrefs.pdb == ["1ABC"]
    assert xrefs.chembl == ["CHEMBL1234"]
    assert xrefs.refseq == ["WP_009885622.1"]
    assert "lists this feature ID" in xrefs.note


def test_falls_back_to_feature_id_then_locus_tag(tmp_path) -> None:
    # No checksum hit, but the feature ID is indexed in UniParc.
    session = RoutingSession(uniparc_hits={'dbid:"' + FEATURE_ID + '"': "UPI00001263F8"})
    mapper = make_mapper(tmp_path, session, taxon_id=243273)
    assert mapper.map_protein(FEATURE_ID, "MKALIV").route == ROUTE_PATRIC

    # Neither UniParc route hits: the locus tag is verified against the entry's OLN.
    session = RoutingSession(uniparc_hits={})
    mapper = make_mapper(tmp_path / "b", session, taxon_id=243273)
    xrefs = mapper.map_protein(FEATURE_ID, "MKALIV", locus_tag="MG_401")
    assert xrefs.route == ROUTE_LOCUS
    assert xrefs.uniprot == "P47641"


def test_locus_tag_comparison_ignores_punctuation(tmp_path) -> None:
    # BV-BRC writes MG_401, UniProt writes MG401.
    session = RoutingSession(uniparc_hits={})
    mapper = make_mapper(tmp_path, session, taxon_id=243273)
    assert mapper.map_protein(FEATURE_ID, "", locus_tag="MG401").uniprot == "P47641"
    assert mapper.map_protein(FEATURE_ID, "", locus_tag="MG_401").uniprot == "P47641"


def test_locus_tag_route_rejects_an_unverified_match(tmp_path) -> None:
    # Free-text search can return an entry that does not actually carry the locus.
    other = {"results": [dict(UNIPROTKB_ENTRY["results"][0], genes=[{"geneName": {"value": "other"}}])]}
    session = RoutingSession(uniparc_hits={}, uniprotkb=other)
    mapper = make_mapper(tmp_path, session, taxon_id=243273)

    xrefs = mapper.map_protein(FEATURE_ID, "", locus_tag="MG_401")
    assert xrefs.uniprot is None
    assert xrefs.route == ROUTE_UNMAPPED


def test_unmapped_is_explicit_never_a_guess(tmp_path) -> None:
    session = RoutingSession(uniparc_hits={}, uniprotkb={"results": []})
    mapper = make_mapper(tmp_path, session)

    xrefs = mapper.map_protein("fig|9.9.peg.1", "MKALIV", locus_tag="ZZ_1")

    assert xrefs.uniprot is None
    assert xrefs.mapped is False
    assert xrefs.route == ROUTE_UNMAPPED
    assert xrefs.note


def test_pdb_homolog_accessions_are_related_not_our_identity(tmp_path) -> None:
    session = RoutingSession(uniparc_hits={}, uniprotkb={"results": []})
    mapper = make_mapper(tmp_path, session)

    xrefs = mapper.map_protein(
        "fig|9.9.peg.1", "MKALIV", pdb_hit_uniprot_ids=["P00001", "P00002", "P00001"]
    )

    assert xrefs.uniprot is None
    assert xrefs.mapped is False
    assert xrefs.related_uniprot == ["P00001", "P00002"]
    assert xrefs.route == ROUTE_PDB_HIT
    assert "not ours" in xrefs.note


def test_crossref_selection_prefers_feature_match_taxon_then_swissprot() -> None:
    accession, note = _select_uniprot_crossref(
        UNIPARC_ENTRY["uniParcCrossReferences"], feature_id=FEATURE_ID, taxon_id=243273
    )
    assert accession == "P47641"
    assert "lists this feature ID" in note

    inactive_only = [{"database": "UniProtKB/TrEMBL", "id": "X", "active": False}]
    accession, note = _select_uniprot_crossref(inactive_only)
    assert accession is None
    assert "no active" in note

    two_active = [
        {"database": "UniProtKB/TrEMBL", "id": "T1", "active": True, "organism": {"taxonId": 999}},
        {"database": "UniProtKB/Swiss-Prot", "id": "S1", "active": True, "organism": {"taxonId": 243273}},
    ]
    accession, note = _select_uniprot_crossref(two_active, taxon_id=243273)
    assert accession == "S1"
    assert "2 active accessions" in note


def test_second_run_makes_no_network_calls(tmp_path) -> None:
    """Definition of done: cached, so a rerun is free."""
    session = RoutingSession(uniparc_hits={"checksum:2CBBDAF094BCDEEB": "UPI00001263F8"})
    mapper = make_mapper(tmp_path, session, taxon_id=243273)
    first = mapper.map_protein(FEATURE_ID, ATPA, locus_tag="MG_401")
    calls_after_first = len(session.calls)
    assert calls_after_first > 0

    again = mapper.map_protein(FEATURE_ID, ATPA, locus_tag="MG_401")
    assert len(session.calls) == calls_after_first
    assert again.uniprot == first.uniprot

    offline_session = RoutingSession(uniparc_hits={})
    offline = CachedJsonClient(
        cache=JsonCache(tmp_path / "cache.sqlite"),
        session=offline_session,
        offline=True,
        min_interval_seconds=0,
    )
    replayed = IdMapper(offline, taxon_id=243273).map_protein(FEATURE_ID, ATPA, locus_tag="MG_401")
    assert replayed.uniprot == "P47641"
    assert offline_session.calls == []


def test_map_many_keys_by_feature_id(tmp_path) -> None:
    session = RoutingSession(uniparc_hits={"checksum:2CBBDAF094BCDEEB": "UPI00001263F8"})
    mapper = make_mapper(tmp_path, session, taxon_id=243273)

    results = mapper.map_many(
        [
            {"feature_id": FEATURE_ID, "sequence": ATPA, "locus_tag": "MG_401"},
            {"feature_id": "fig|9.9.peg.1", "sequence": "MKALIV"},
        ],
        workers=1,
    )

    assert set(results) == {FEATURE_ID, "fig|9.9.peg.1"}
    assert results[FEATURE_ID].mapped is True


def test_xrefs_serializes_with_an_explicit_mapped_flag() -> None:
    payload = Xrefs(feature_id="fig|1.1.peg.1").as_dict()
    assert payload["mapped"] is False
    assert payload["uniprot"] is None
    assert payload["route"] == ROUTE_UNMAPPED


# --- taxon-index fast path ------------------------------------------------

STREAM_ENTRY = dict(
    UNIPROTKB_ENTRY["results"][0],
    sequence={"crc64": "2CBBDAF094BCDEEB", "length": 518},
)
STREAM_PAYLOAD = {"results": [STREAM_ENTRY]}


def test_build_taxon_index_keys_by_checksum_and_locus() -> None:
    index = build_taxon_index(243273, STREAM_PAYLOAD["results"])

    assert index.entries == 1
    assert index.by_checksum["2CBBDAF094BCDEEB"]["primaryAccession"] == "P47641"
    # BV-BRC writes MG_401, UniProt MG401: the index key is punctuation-free.
    assert index.by_locus["MG401"]["primaryAccession"] == "P47641"


def test_proteome_checksum_maps_without_a_per_protein_request(tmp_path) -> None:
    session = RoutingSession(stream=STREAM_PAYLOAD)
    mapper = make_mapper(tmp_path, session, taxon_id=243273)

    xrefs = mapper.map_protein(FEATURE_ID, ATPA, locus_tag="MG_401")

    assert xrefs.uniprot == "P47641"
    assert xrefs.route == ROUTE_PROTEOME_SEQUENCE
    assert xrefs.gene_name == "atpA"
    assert xrefs.pdb == ["1ABC"]
    # One bulk fetch, and nothing per protein.
    assert len(session.calls) == 1
    assert "uniprotkb/stream" in session.calls[0][0]

    mapper.map_protein("fig|243273.25.peg.999", ATPA)
    assert len(session.calls) == 1  # index reused in-process


def test_proteome_locus_route_when_the_sequence_differs(tmp_path) -> None:
    session = RoutingSession(stream=STREAM_PAYLOAD)
    mapper = make_mapper(tmp_path, session, taxon_id=243273)

    xrefs = mapper.map_protein(FEATURE_ID, "MKALIVDIFFERENT", locus_tag="MG_401")

    assert xrefs.route == ROUTE_PROTEOME_LOCUS
    assert xrefs.uniprot == "P47641"


def test_falls_back_to_uniparc_when_not_in_the_index(tmp_path) -> None:
    session = RoutingSession(
        stream={"results": []},
        uniparc_hits={"checksum:2CBBDAF094BCDEEB": "UPI00001263F8"},
    )
    mapper = make_mapper(tmp_path, session, taxon_id=243273)

    xrefs = mapper.map_protein(FEATURE_ID, ATPA)

    assert xrefs.route == ROUTE_SEQUENCE
    assert xrefs.uniprot == "P47641"


def test_a_failed_index_fetch_does_not_block_the_other_routes(tmp_path) -> None:
    class FailingStream(RoutingSession):
        def get(self, url, params=None, timeout=None):
            if "uniprotkb/stream" in url:
                self.calls.append((url, params or {}))
                return FakeResponse({"error": "boom"}, status_code=500)
            return super().get(url, params, timeout)

    session = FailingStream(uniparc_hits={"checksum:2CBBDAF094BCDEEB": "UPI00001263F8"})
    mapper = make_mapper(tmp_path, session, taxon_id=243273)

    xrefs = mapper.map_protein(FEATURE_ID, ATPA)

    assert xrefs.uniprot == "P47641"
    assert xrefs.route == ROUTE_SEQUENCE
    assert mapper.failures and mapper.failures[0]["feature_id"] == "taxon:243273"

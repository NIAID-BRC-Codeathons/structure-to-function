"""Shared test doubles for the BV-BRC Data API route.

These live here rather than in one of the `test_m1_*` modules because a test module is
not a reliable import target: whether `from tests.test_m1_bvbrc_api import ...` resolves
depends on the pytest version, its import mode, and whether `tests/` is a package. It
worked locally and broke under a different pytest. `conftest.py` is loaded by pytest
itself on every version, so fixtures defined here are always available by name.

Nothing here touches the network. Responses come from
`fixtures/m1/bvbrc_api/klebsiella_hs11286.json`.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from s2f.common.http import CachedJsonClient, JsonCache
from s2f.m1_genome.bvbrc_api import BvbrcApi

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
KP_FIXTURE = FIXTURES / "m1/bvbrc_api/klebsiella_hs11286.json"
GENOME_ID = "1125630.4"


class FakeResponse:
    """Enough of `requests.Response` for `CachedJsonClient`."""

    def __init__(self, body, *, status: int = 200, headers: dict | None = None) -> None:
        self._body = body
        self.status_code = status
        self.headers = headers or {}
        self.content = b"x" if body is not None else b""
        self.text = ""

    @property
    def ok(self) -> bool:
        return self.status_code < 400

    def json(self):
        return self._body


class FakeSession:
    """Routes a BV-BRC URL to the right slice of the fixture.

    Matching on the core plus the distinguishing RQL fragment, rather than on an exact
    URL, is deliberate: the tests should keep passing when a `select(...)` list gains a
    field, and should fail when the *query* changes meaning.
    """

    def __init__(self, payload: dict, *, fail_facets: bool = False) -> None:
        self.payload = payload
        self.fail_facets = fail_facets
        self.calls: list[str] = []
        self.headers: dict[str, str] = {}

    def get(self, url, params=None, timeout=None, headers=None):
        self.calls.append(url)
        headers = headers or {}
        body, total = self._route(url)
        if body is None:
            return FakeResponse(None, status=404)
        if isinstance(body, list):
            start, end = self._range(headers.get("Range"))
            page = body[start:end + 1] if end is not None else body
            return FakeResponse(
                page,
                headers={"Content-Range": f"items {start}-{start + len(page) - 1}/{total}"},
            )
        return FakeResponse(body)

    @staticmethod
    def _range(header: str | None) -> tuple[int, int | None]:
        if not header:
            return 0, None
        match = re.search(r"items=(\d+)-(\d+)", header)
        return (int(match.group(1)), int(match.group(2))) if match else (0, None)

    def _route(self, url: str):
        data = self.payload
        if "facet(" in url:
            if self.fail_facets:
                return None, 0
            field = re.search(r"\(field,([a-z_]+)\)", url)
            name = field.group(1) if field else ""
            if name == "species":
                return {"facet_counts": {"facet_fields":
                                         {"species": data["close_pathogens_facet"]}}}, 0
            if name == "disease":
                return {"facet_counts": {"facet_fields":
                                         {"disease": {"Nosocomial infection": 12}}}}, 0
            if name == "oxygen_requirement":
                return {"facet_counts": {"facet_fields":
                                         {"oxygen_requirement": {"Facultative": 90}}}}, 0
            if name in {"superclass", "class"}:
                return {"facet_counts": {"facet_fields": {name: {"Metabolism": 40}}}}, 0
            # growth / isolation facets: an empty-ish first value must be skipped
            return {"facet_counts": {"facet_fields": {name: {"": 5, "Negative": 57}}}}, 0
        if "/taxonomy/" in url:
            return data["taxonomy"], len(data["taxonomy"])
        if "/genome_feature/" in url:
            if "select(pgfam_id)" in url:
                rows = [{"pgfam_id": f["pgfam_id"]} for f in data["features"]
                        if f.get("pgfam_id")]
                return rows, len(rows)
            return data["features"], len(data["features"])
        if "/sp_gene/" in url:
            return data["specialty"], len(data["specialty"])
        if "/genome_amr/" in url:
            return data["amr_phenotypes"], len(data["amr_phenotypes"])
        if "/subsystem/" in url:
            return data["subsystems"], len(data["subsystems"])
        if "/feature_sequence/" in url:
            rows = [{"md5": k, "sequence": v} for k, v in data["sequences"].items()]
            return rows, len(rows)
        if "/genome/" in url:
            if "eq(genus," in url:
                return data["neighbors"], len(data["neighbors"])
            if f"eq(genome_id,{GENOME_ID})" in url:
                return [data["genome"]], 1
            if "eq(species," in url:
                return [data["genome"]], 1
            return [], 0
        return None, 0


class HeaderlessSession(FakeSession):
    """A server (or intermediary) that returns no usable `Content-Range`."""

    def __init__(self, payload, *, total_marker: str = "") -> None:
        super().__init__(payload)
        self.total_marker = total_marker

    def get(self, url, params=None, timeout=None, headers=None):
        response = super().get(url, params, timeout, headers)
        if isinstance(response._body, list):
            response.headers = ({"Content-Range": self.total_marker}
                                if self.total_marker else {})
        return response


class RangeIgnoringSession(FakeSession):
    """A server that ignores `Range` and replays page one forever."""

    def get(self, url, params=None, timeout=None, headers=None):
        return super().get(url, params, timeout, headers={"Range": "items=0-3"})


# --- fixtures ---------------------------------------------------------------
@pytest.fixture(scope="session")
def kp_payload() -> dict:
    """The canned Klebsiella HS11286 API responses, read once per session."""
    return json.loads(KP_FIXTURE.read_text())


@pytest.fixture
def fake_session(kp_payload):
    """Factory: `fake_session()` or `fake_session(cls=..., **kwargs)`."""
    def build(cls=FakeSession, **kwargs):
        return cls(kp_payload, **kwargs)
    return build


@pytest.fixture
def make_api(kp_payload, tmp_path):
    """Factory returning `(BvbrcApi, session)` backed by the fixture and a real cache."""
    def build(session=None, *, cache_name: str = "cache.sqlite", offline: bool = False,
              **session_kwargs):
        session = session if session is not None else FakeSession(kp_payload,
                                                                  **session_kwargs)
        client = CachedJsonClient(cache=JsonCache(tmp_path / cache_name),
                                  session=session, offline=offline,
                                  min_interval_seconds=0, sleep=lambda _: None)
        return BvbrcApi(client), session
    return build


@pytest.fixture
def api_bundle(make_api, kp_payload):
    """The collected section bundle for the fixture genome."""
    from s2f.m1_genome.collect import collect_all

    api, _ = make_api()
    return collect_all(api, kp_payload["genome"], cap=500)

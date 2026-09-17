"""Binary-file caching on the shared HTTP client."""

from __future__ import annotations

from pathlib import Path

import pytest

from s2f.common.http import CachedJsonClient, OfflineCacheMiss


class FileResponse:
    status_code = 200
    content = b"structure bytes"
    text = "structure bytes"
    headers: dict[str, str] = {}
    ok = True


class FileSession:
    def __init__(self) -> None:
        self.calls: list[str] = []
        self.headers: dict[str, str] = {}

    def get(self, url, params=None, timeout=None):
        self.calls.append(url)
        return FileResponse()

    def post(self, url, json=None, timeout=None):  # pragma: no cover
        raise AssertionError("file download should not POST")


def test_file_download_is_cached_and_replays_offline(tmp_path: Path) -> None:
    url = "https://example.invalid/1ABC.cif"
    session = FileSession()
    client = CachedJsonClient(session=session, min_interval_seconds=0)

    first = client.get_bytes("structures", url, cache_dir=tmp_path)
    second = client.get_bytes("structures", url, cache_dir=tmp_path)

    assert first == second == b"structure bytes"
    assert session.calls == [url]

    offline_session = FileSession()
    offline = CachedJsonClient(session=offline_session, offline=True, min_interval_seconds=0)
    assert offline.get_bytes("structures", url, cache_dir=tmp_path) == first
    assert offline_session.calls == []


def test_offline_file_cache_miss_is_explicit(tmp_path: Path) -> None:
    client = CachedJsonClient(session=FileSession(), offline=True, min_interval_seconds=0)

    with pytest.raises(OfflineCacheMiss, match="no cached file"):
        client.get_bytes("structures", "https://example.invalid/missing.cif", cache_dir=tmp_path)

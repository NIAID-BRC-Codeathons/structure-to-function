"""Cached HTTP with retry and polite rate limiting (issue #4).

Every outbound call in the pipeline should go through this client so responses are cached on
disk, failures are retried, and a run can be replayed offline. Cache and retry behaviour are
adapted from ``AutoPDB/src/autopdb/cache.py`` and ``rcsb.py`` (a separate MIT-licensed
repository by the same author), kept here rather than imported so this repo has no cross-repo
dependency.

Provenance stamping (the other half of #4) is not here yet; records currently carry their own
source fields.
"""

from __future__ import annotations

import hashlib
import json
import os
import sqlite3
import tempfile
import threading
import time
from datetime import UTC, datetime, timedelta
from pathlib import Path
from typing import Any

import requests

USER_AGENT = "s2f/0.1 (NIAID-BRC Codeathon Project 9; +https://github.com/NIAID-BRC-Codeathons/structure-to-function)"
RETRY_STATUS = {429, 500, 502, 503, 504}


class HttpError(RuntimeError):
    """Raised when a request fails after its retries."""


class OfflineCacheMiss(HttpError):
    """Raised when offline mode cannot satisfy a request from the cache."""


class JsonCache:
    """SQLite store of JSON responses, keyed by namespace and request hash."""

    def __init__(self, path: str | Path) -> None:
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._lock = threading.Lock()
        self._connection = sqlite3.connect(self.path, check_same_thread=False)
        self._connection.execute(
            """
            CREATE TABLE IF NOT EXISTS api_cache (
                namespace TEXT NOT NULL,
                cache_key TEXT NOT NULL,
                payload_json TEXT NOT NULL,
                fetched_at TEXT NOT NULL,
                PRIMARY KEY (namespace, cache_key)
            )
            """
        )
        self._connection.commit()

    def get(
        self, namespace: str, cache_key: str, *, max_age: timedelta | None = None
    ) -> Any | None:
        with self._lock:
            row = self._connection.execute(
                "SELECT payload_json, fetched_at FROM api_cache WHERE namespace = ? AND cache_key = ?",
                (namespace, cache_key),
            ).fetchone()
        if row is None:
            return None
        payload_json, fetched_at = row
        if max_age is not None:
            fetched = datetime.fromisoformat(fetched_at)
            if datetime.now(UTC) - fetched > max_age:
                return None
        return json.loads(payload_json)

    def set(self, namespace: str, cache_key: str, payload: Any) -> None:
        with self._lock:
            self._connection.execute(
                """
                INSERT INTO api_cache (namespace, cache_key, payload_json, fetched_at)
                VALUES (?, ?, ?, ?)
                ON CONFLICT(namespace, cache_key) DO UPDATE SET
                    payload_json = excluded.payload_json,
                    fetched_at = excluded.fetched_at
                """,
                (namespace, cache_key, json.dumps(payload, sort_keys=True), datetime.now(UTC).isoformat()),
            )
            self._connection.commit()

    def close(self) -> None:
        with self._lock:
            self._connection.close()

    def __enter__(self) -> "JsonCache":
        return self

    def __exit__(self, *_: object) -> None:
        self.close()


class CachedJsonClient:
    """POST JSON with cache, retry/backoff and a shared minimum request interval."""

    def __init__(
        self,
        *,
        cache: JsonCache | None = None,
        session: requests.Session | None = None,
        offline: bool = False,
        timeout_seconds: float = 60.0,
        max_age: timedelta | None = timedelta(days=30),
        max_attempts: int = 4,
        min_interval_seconds: float = 0.05,
        sleep = time.sleep,
    ) -> None:
        self.cache = cache
        self.session = session or requests.Session()
        self.offline = offline
        self.timeout_seconds = timeout_seconds
        self.max_age = max_age
        self.max_attempts = max_attempts
        self.min_interval_seconds = min_interval_seconds
        self._sleep = sleep
        self._rate_lock = threading.Lock()
        self._last_request = 0.0
        self.session.headers.setdefault("User-Agent", USER_AGENT)

    @staticmethod
    def cache_key(payload: Any) -> str:
        return hashlib.sha256(json.dumps(payload, sort_keys=True).encode("utf-8")).hexdigest()

    def post_json(self, namespace: str, url: str, payload: dict[str, Any]) -> dict[str, Any]:
        key = self.cache_key({"url": url, "payload": payload})
        if self.cache is not None:
            cached = self.cache.get(namespace, key, max_age=None if self.offline else self.max_age)
            if cached is not None:
                return cached
        if self.offline:
            raise OfflineCacheMiss(f"no cached response for {namespace}:{key[:12]} ({url})")

        response = self._post_with_retry(url, payload)
        body = response.json() if response.content else {}
        if self.cache is not None:
            self.cache.set(namespace, key, body)
        return body

    def get_json(
        self,
        namespace: str,
        url: str,
        params: dict[str, Any] | None = None,
        *,
        allow_statuses: tuple[int, ...] = (),
    ) -> dict[str, Any]:
        """GET with the same cache, retry and offline semantics as post_json.

        ``allow_statuses`` names statuses that are a legitimate answer rather than a failure —
        an AlphaFold 404 means "no model for this accession", not "the request broke". Those
        come back as ``{"_http_status": <code>}`` and are cached, since they are stable answers.
        """
        params = dict(params or {})
        key = self.cache_key({"url": url, "params": params})
        if self.cache is not None:
            cached = self.cache.get(namespace, key, max_age=None if self.offline else self.max_age)
            if cached is not None:
                return cached
        if self.offline:
            raise OfflineCacheMiss(f"no cached response for {namespace}:{key[:12]} ({url})")

        response = self._request_with_retry("GET", url, params=params, allow_statuses=allow_statuses)
        if not response.ok:
            body: dict[str, Any] = {"_http_status": response.status_code}
        else:
            body = response.json() if response.content else {}
        if self.cache is not None:
            self.cache.set(namespace, key, body)
        return body

    def get_bytes(self, namespace: str, url: str, *, cache_dir: str | Path) -> bytes:
        """Download a file through the shared retry/rate-limit path and cache it on disk.

        Structure files are not JSON, so they do not belong in :class:`JsonCache`. The URL is
        hashed into a small file cache instead. Keeping this operation on the shared client
        preserves the repository rule that all outbound HTTP uses one retry, rate-limit, and
        user-agent implementation.
        """
        key = self.cache_key({"url": url})
        path = Path(cache_dir) / namespace / key
        if path.exists():
            return path.read_bytes()
        if self.offline:
            raise OfflineCacheMiss(f"no cached file for {namespace}:{key[:12]} ({url})")

        response = self._request_with_retry("GET", url)
        payload = response.content
        path.parent.mkdir(parents=True, exist_ok=True)
        handle = tempfile.NamedTemporaryFile("wb", dir=path.parent, prefix=".download-", delete=False)
        try:
            with handle:
                handle.write(payload)
                handle.flush()
                os.fsync(handle.fileno())
            os.replace(handle.name, path)
        except BaseException:
            Path(handle.name).unlink(missing_ok=True)
            raise
        return payload

    def _throttle(self) -> None:
        if self.min_interval_seconds <= 0:
            return
        with self._rate_lock:
            elapsed = time.monotonic() - self._last_request
            wait = self.min_interval_seconds - elapsed
            if wait > 0:
                self._sleep(wait)
            self._last_request = time.monotonic()

    def _post_with_retry(self, url: str, payload: dict[str, Any]) -> requests.Response:
        return self._request_with_retry("POST", url, json=payload)

    def _request_with_retry(
        self,
        method: str,
        url: str,
        *,
        json: dict[str, Any] | None = None,
        params: dict[str, Any] | None = None,
        allow_statuses: tuple[int, ...] = (),
    ) -> requests.Response:
        last_error: str = ""
        for attempt in range(self.max_attempts):
            self._throttle()
            try:
                if method == "POST":
                    response = self.session.post(url, json=json, timeout=self.timeout_seconds)
                else:
                    response = self.session.get(url, params=params, timeout=self.timeout_seconds)
            except requests.RequestException as exc:
                last_error = f"{type(exc).__name__}: {exc}"
            else:
                if response.status_code == 204:
                    return response
                if response.ok:
                    return response
                if response.status_code in allow_statuses:
                    return response
                last_error = f"HTTP {response.status_code}: {response.text[:200]}"
                if response.status_code not in RETRY_STATUS:
                    raise HttpError(f"{method} {url} failed — {last_error}")
                retry_after = _retry_after_seconds(response)
                if retry_after is not None:
                    self._sleep(retry_after)
                    continue
            if attempt < self.max_attempts - 1:
                self._sleep(min(2.0**attempt, 8.0))
        raise HttpError(f"{method} {url} failed after {self.max_attempts} attempts — {last_error}")


def _retry_after_seconds(response: requests.Response) -> float | None:
    raw = response.headers.get("Retry-After")
    if not raw:
        return None
    try:
        return max(0.0, min(float(raw), 60.0))
    except ValueError:
        return None

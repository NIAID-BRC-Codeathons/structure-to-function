"""Predict a taxon for a blinded assembly with BV-BRC's Minhash service (issue #5).

Similar Genome Finder is not an AppService app and has no `p3-` command: it is a
synchronous JSON-RPC service (`Minhash`), which the web UI calls at
``https://p3.theseed.org/services/minhash_service``. Methods and parameter order come
from the service spec, `BV-BRC/p3_minhash`.

The taxon is called at the rank the Mash distance justifies, not at the rank of the
top hit. Mash distance approximates 1 - ANI, so 0.05 is about the 95% ANI species
boundary. A novel organism lands at genus or at nothing, and nothing is a deferred
case (see `docs/issues.md`).

The genetic code is read from the *called* taxon's own record. That matters: the
class Mollicutes reports code 11 while its genera report 4, so calling at too coarse
a rank silently yields the wrong translation table.
"""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path
from typing import Any

MINHASH_URL = "https://p3.theseed.org/services/minhash_service"
DATA_API = "https://www.bv-brc.org/api"

#: (max Mash distance, rank to call). First match wins; past the last, no call.
RANK_BY_DISTANCE: tuple[tuple[float, str], ...] = ((0.05, "species"), (0.20, "genus"))

TOKEN_PATH = "~/.patric_token"


class TaxonCallError(RuntimeError):
    pass


def _token() -> str:
    path = Path(os.path.expanduser(TOKEN_PATH))
    if not path.exists():
        raise TaxonCallError(f"no BV-BRC token at {path}; run p3-login first")
    return path.read_text().strip()


def _curl(url: str, payload: dict[str, Any] | None = None, auth: bool = False,
          timeout: int = 180) -> Any:
    """Fetch JSON with curl.

    Not requests: the BV-BRC data API rejects some default user agents with a 403,
    and curl is present wherever the `p3-` CLI is.
    """
    cmd = ["curl", "-sS", "--fail", "--max-time", str(timeout), url, "-H", "Accept: application/json"]
    if auth:
        cmd += ["-H", "Authorization: " + _token()]
    if payload is not None:
        cmd += ["-H", "Content-Type: application/json", "-d", json.dumps(payload)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise TaxonCallError(f"request failed ({url}): {proc.stderr.strip()[:300]}")
    try:
        return json.loads(proc.stdout)
    except json.JSONDecodeError as exc:
        raise TaxonCallError(f"non-JSON response from {url}: {proc.stdout[:200]}") from exc


def minhash_hits(ws_fasta_path: str, *, max_pvalue: float = 1e-5, max_distance: float = 0.2,
                 max_hits: int = 50, all_public: bool = True,
                 bacterial: bool = True, viral: bool = False) -> list[dict[str, Any]]:
    """Closest BV-BRC genomes for a FASTA already in the user's workspace.

    ``all_public`` false restricts the search to reference and representative
    sketches, which is what the web form's Scope control does; for an unknown
    organism that is usually too narrow.
    """
    scope = 0 if all_public else 1
    payload = {
        "id": 1,
        "version": "1.1",
        "method": "Minhash.compute_genome_distance_for_fasta2",
        "params": [ws_fasta_path, max_pvalue, max_distance, max_hits,
                   scope, scope, int(bacterial), int(viral)],
    }
    response = _curl(MINHASH_URL, payload, auth=True)
    if "error" in response:
        raise TaxonCallError(f"Minhash error: {json.dumps(response['error'])[:300]}")
    return [
        {"genome_id": g, "distance": float(d), "pvalue": float(p), "counts": c}
        for g, d, p, c in response["result"][0]
    ]


def genome_metadata(genome_ids: list[str]) -> dict[str, dict[str, Any]]:
    if not genome_ids:
        return {}
    query = (
        f"{DATA_API}/genome/?in(genome_id,({','.join(genome_ids)}))"
        "&select(genome_id,genome_name,taxon_id,genome_status,genome_quality)"
        f"&limit({len(genome_ids) + 5})&http_accept=application/json"
    )
    return {row["genome_id"]: row for row in _curl(query)}


def taxonomy(taxon_id: int | str) -> dict[str, Any]:
    query = (
        f"{DATA_API}/taxonomy/?eq(taxon_id,{taxon_id})"
        "&select(taxon_id,taxon_name,taxon_rank,genetic_code,lineage_ids,lineage_names,lineage_ranks)"
        "&http_accept=application/json"
    )
    rows = _curl(query)
    if not rows:
        raise TaxonCallError(f"no taxonomy record for taxon {taxon_id}")
    return rows[0]


def call_taxon(ws_fasta_path: str, **minhash_kwargs: Any) -> dict[str, Any]:
    """Minhash -> closest genomes -> taxon and genetic code at a justified rank."""
    hits = minhash_hits(ws_fasta_path, **minhash_kwargs)
    if not hits:
        raise TaxonCallError(
            "Minhash returned no hits; behaviour for this case is deferred "
            "(see 'define behaviour when Similar Genome Finder returns no usable hit')"
        )
    meta = genome_metadata([h["genome_id"] for h in hits[:25]])
    for hit in hits[:25]:
        hit.update({k: meta.get(hit["genome_id"], {}).get(k)
                    for k in ("genome_name", "taxon_id", "genome_status", "genome_quality")})

    top = hits[0]
    rank = next((r for cut, r in RANK_BY_DISTANCE if top["distance"] <= cut), None)
    if rank is None:
        raise TaxonCallError(
            f"top hit at Mash distance {top['distance']:.4f} exceeds every cut-off "
            f"{RANK_BY_DISTANCE}; no taxon call made"
        )
    if not top.get("taxon_id"):
        raise TaxonCallError(f"top hit {top['genome_id']} has no taxon_id")

    lineage = taxonomy(top["taxon_id"])
    ranks = lineage.get("lineage_ranks") or []
    if rank not in ranks:
        raise TaxonCallError(f"lineage of taxon {top['taxon_id']} has no {rank} rank: {ranks}")
    called = taxonomy(lineage["lineage_ids"][ranks.index(rank)])

    return {
        "called_by": "minhash",
        "service": MINHASH_URL,
        "ws_fasta": ws_fasta_path,
        "n_hits": len(hits),
        "hits": hits,
        "top_hit": top,
        "rank_rule": [list(pair) for pair in RANK_BY_DISTANCE],
        "called_rank": rank,
        "taxon_id": int(called["taxon_id"]),
        "taxon_name": called["taxon_name"],
        "genetic_code": int(called["genetic_code"]),
        "lineage": list(zip(lineage.get("lineage_ranks") or [],
                            lineage.get("lineage_names") or [],
                            lineage.get("lineage_ids") or [])),
    }

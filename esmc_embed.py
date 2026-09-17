"""ESMC embeddings (and optional SAE features) for M2 proteins with no PDB hit (issue #9, partial).

Reads ``<run>/m2_pdb/no_pdb_hit.tsv`` (the output of ``s2f.m2_triage`` / issue #8) and computes
an ESMC embedding for every protein whose retrieval was a genuine ``no-hit`` — never for
``query-failed`` rows, which are recorded but skipped, and never for ``found`` rows. See
"Eligibility" below for why this distinction matters.

Relationship to the rest of the pipeline
-----------------------------------------
``docs/02-m2-triage.md`` defines two different output shapes:

- ``proteins[].annotations[]``: ``source``, ``hit``, ``identity/coverage or TM-score``,
  ``description`` — this is a *homology hit* shape (PDB sequence search, Foldseek structure
  search). A raw embedding vector has no "hit," so it does not belong here.
- ``proteins[].xrefs.esm_atlas[]``: this is where an ESMC embedding (or an ESM Atlas search
  result derived from one) belongs.

This module therefore writes a new, additive artifact (``esmc_embeddings.jsonl`` /
``esmc_embeddings.tsv``) rather than annotation rows. If a later step searches the ESM Atlas
with these embeddings and gets back a homolog, *that* result is an annotations-shaped row with
``source: "esm_atlas"``, parallel to ``source: "foldseek"`` (issue #9's Foldseek half) and
``source: "rcsb_sequence"`` (issue #8) — kept distinguishable by ``source`` alone, per the
project's convention. Nothing in this module writes such a row yet.

Eligibility
-----------
``no_pdb_hit.tsv`` is every protein with ``proteins.tsv``'s ``no_pdb_hit`` flag set, and that
flag is true for **both** ``retrieval_status == "no-hit"`` and ``retrieval_status ==
"query-failed"`` (see ``score.py``: ``"no_pdb_hit": not qualifying``, which is also true when a
search failed and produced no hits to qualify). Treating a failed search as "confirmed no
homolog, go compute an embedding" would silently launder a network failure into a structural
claim. So only ``no-hit`` rows are embedded here; ``query-failed`` and ``found`` rows are kept in
the output with a distinct ``esmc_status`` so they stay visible rather than silently dropped.

Caching
-------
``no_pdb_hit.tsv`` carries no sequence column, so sequences are re-joined from the M1 input
(``proteins.faa``) by ``feature_id``, the same input directory ``s2f.m2_triage`` reads.

The ESMC embedding call goes through Biohub's ``esm`` SDK client, which manages its own HTTP
transport — it cannot be handed ``s2f.common.http.CachedJsonClient`` directly. What *is* reused
from ``s2f.common.http`` is ``JsonCache`` (the on-disk SQLite cache) and the ``HttpError`` /
``OfflineCacheMiss`` exception vocabulary, so ``--offline`` gives the same "replay from cache,
never touch the network" guarantee the rest of the pipeline has. The retry/backoff *pattern*
mirrors ``common/http.py`` (bounded attempts, exponential backoff) but is a separate
implementation around the SDK call, not the same code path.

SAE feature extraction is optional (``--sae``) and, per the caveats in
``fetch_esmc_features.py``, uses local Hugging-Face-weights inference rather than a confirmed
hosted API call — this is slower and needs real compute, so it defaults to off.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import time
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Callable, Iterable

from ..common.http import HttpError, JsonCache, OfflineCacheMiss
from .bvbrc_input import load_input

FIXTURE_DIR = Path(__file__).resolve().parents[2] / "fixtures" / "m2"

DEFAULT_MODEL = "esmc-600m-2024-12"
DEFAULT_SAE_LAYER = 60
MAX_ATTEMPTS = 4

# esmc_status values. Distinct from the upstream PDB retrieval_status (never overwritten).
NOT_APPLICABLE_FOUND = "not-applicable-found"
SKIPPED_QUERY_FAILED = "skipped-upstream-query-failed"
SKIPPED_NO_SEQUENCE = "skipped-no-sequence"
EMBEDDED = "embedded"
QUERY_FAILED = "query-failed"

# sae_status values.
SAE_NOT_REQUESTED = "not-requested"
SAE_EXTRACTED = "extracted"
SAE_FAILED = "query-failed"

SPECIALTY_COLUMNS = ("virulence_amr", "essential", "drug_target", "annotation_gap", "human_homolog_penalty")

METADATA_COLUMNS = [
    "feature_id", "product", *SPECIALTY_COLUMNS, "retrieval_status", "esmc_status",
    "embedding_model", "embedding_dim", "sae_status", "sae_model", "sae_layer", "sae_dim", "error",
]


@dataclass
class EmbeddingResult:
    """Outcome of one embedding call, as returned by an injected ``embed_fn``."""

    status: str  # "embedded" | "query-failed"
    vector: list[float] = field(default_factory=list)
    model: str = ""
    error: str = ""


@dataclass
class SaeResult:
    """Outcome of one SAE feature extraction, as returned by an injected ``sae_fn``."""

    status: str = SAE_NOT_REQUESTED
    vector: list[float] = field(default_factory=list)
    model: str = ""
    layer: int | None = None
    error: str = ""


EmbedFn = Callable[[str], EmbeddingResult]
SaeFn = Callable[[str], SaeResult]


def sequence_cache_key(sequence: str, model: str) -> str:
    digest = hashlib.sha256(sequence.strip().upper().encode("ascii", errors="ignore")).hexdigest()
    return f"{model}:{digest}"


def cached_embed(
    cache: JsonCache | None, offline: bool, model: str, sequence: str, embed_fn: EmbedFn
) -> EmbeddingResult:
    """Cache wrapper around one embedding call. Raises OfflineCacheMiss, never calls the
    network, when ``offline`` and nothing is cached — the same contract as CachedJsonClient."""
    key = sequence_cache_key(sequence, model)
    if cache is not None:
        cached = cache.get("esmc_embedding", key)
        if cached is not None:
            return EmbeddingResult(**cached)
    if offline:
        raise OfflineCacheMiss(f"no cached ESMC embedding for {key[:16]}")

    result = embed_fn(sequence)
    if cache is not None and result.status == EMBEDDED:
        cache.set(
            "esmc_embedding",
            key,
            {"status": result.status, "vector": result.vector, "model": result.model, "error": result.error},
        )
    return result


def cached_sae(
    cache: JsonCache | None,
    offline: bool,
    model: str,
    layer: int,
    sequence: str,
    sae_fn: SaeFn,
) -> SaeResult:
    key = sequence_cache_key(sequence, f"{model}:layer{layer}")
    if cache is not None:
        cached = cache.get("esmc_sae", key)
        if cached is not None:
            return SaeResult(**cached)
    if offline:
        raise OfflineCacheMiss(f"no cached SAE features for {key[:16]}")

    result = sae_fn(sequence)
    if cache is not None and result.status == SAE_EXTRACTED:
        cache.set(
            "esmc_sae",
            key,
            {
                "status": result.status,
                "vector": result.vector,
                "model": result.model,
                "layer": result.layer,
                "error": result.error,
            },
        )
    return result


@dataclass
class EsmcEmbedRow:
    """One protein's outcome: metadata carried through, plus the embedding/SAE vectors."""

    feature_id: str
    product: str
    specialty: dict[str, str]
    retrieval_status: str
    esmc_status: str
    embedding_model: str = ""
    embedding_dim: int = 0
    embedding: list[float] = field(default_factory=list)
    sae_status: str = SAE_NOT_REQUESTED
    sae_model: str = ""
    sae_layer: int | None = None
    sae_dim: int = 0
    sae: list[float] = field(default_factory=list)
    error: str = ""

    def metadata_row(self) -> dict[str, Any]:
        row: dict[str, Any] = {
            "feature_id": self.feature_id,
            "product": self.product,
            "retrieval_status": self.retrieval_status,
            "esmc_status": self.esmc_status,
            "embedding_model": self.embedding_model,
            "embedding_dim": self.embedding_dim,
            "sae_status": self.sae_status,
            "sae_model": self.sae_model,
            "sae_layer": self.sae_layer if self.sae_layer is not None else "",
            "sae_dim": self.sae_dim,
            "error": self.error,
        }
        row.update(self.specialty)
        return row

    def vector_row(self) -> dict[str, Any]:
        payload: dict[str, Any] = {"feature_id": self.feature_id, "embedding": self.embedding}
        if self.sae:
            payload["sae"] = self.sae
        return payload


def load_no_pdb_hit_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def embed_eligible_proteins(
    rows: Iterable[dict[str, str]],
    sequence_by_id: dict[str, str],
    *,
    embed_fn: EmbedFn,
    model: str = DEFAULT_MODEL,
    cache: JsonCache | None = None,
    offline: bool = False,
    sae_fn: SaeFn | None = None,
    sae_model: str = "",
    sae_layer: int = DEFAULT_SAE_LAYER,
) -> list[EsmcEmbedRow]:
    """Core logic. Network-free at this call site: ``embed_fn``/``sae_fn`` are injected, so
    tests exercise real control flow (eligibility, caching, offline behaviour) without a real
    ESMC client. See the module docstring's "Eligibility" section for the status rules.
    """
    out: list[EsmcEmbedRow] = []
    for row in rows:
        feature_id = row.get("feature_id", "")
        retrieval_status = row.get("retrieval_status", "")
        specialty = {name: row.get(name, "") for name in SPECIALTY_COLUMNS}
        product = row.get("product", "")

        if retrieval_status == "found":
            out.append(EsmcEmbedRow(feature_id, product, specialty, retrieval_status, NOT_APPLICABLE_FOUND))
            continue
        if retrieval_status == "query-failed":
            out.append(EsmcEmbedRow(feature_id, product, specialty, retrieval_status, SKIPPED_QUERY_FAILED))
            continue
        if retrieval_status != "no-hit":
            out.append(
                EsmcEmbedRow(
                    feature_id, product, specialty, retrieval_status,
                    f"skipped-unrecognized-status:{retrieval_status}",
                )
            )
            continue

        sequence = sequence_by_id.get(feature_id, "")
        if not sequence:
            out.append(
                EsmcEmbedRow(
                    feature_id, product, specialty, retrieval_status, SKIPPED_NO_SEQUENCE,
                    error="no sequence found in M1 input (proteins.faa) for this feature_id",
                )
            )
            continue

        try:
            embedding = cached_embed(cache, offline, model, sequence, embed_fn)
        except (OfflineCacheMiss, HttpError) as exc:
            out.append(
                EsmcEmbedRow(
                    feature_id, product, specialty, retrieval_status, QUERY_FAILED,
                    embedding_model=model, error=str(exc),
                )
            )
            continue

        record = EsmcEmbedRow(
            feature_id=feature_id,
            product=product,
            specialty=specialty,
            retrieval_status=retrieval_status,
            esmc_status=embedding.status,
            embedding_model=embedding.model or model,
            embedding_dim=len(embedding.vector),
            embedding=embedding.vector,
            error=embedding.error,
        )

        if embedding.status == EMBEDDED and sae_fn is not None:
            try:
                sae = cached_sae(cache, offline, sae_model or model, sae_layer, sequence, sae_fn)
            except (OfflineCacheMiss, HttpError) as exc:
                sae = SaeResult(status=SAE_FAILED, error=str(exc))
            record.sae_status = sae.status
            record.sae_model = sae.model or sae_model or model
            record.sae_layer = sae.layer if sae.layer is not None else sae_layer
            record.sae_dim = len(sae.vector)
            record.sae = sae.vector

        out.append(record)
    return out


# --------------------------------------------------------------------------
# Real ESMC client (lazy-imports `esm`; not exercised by unit tests)
# --------------------------------------------------------------------------

def make_esmc_embed_fn(model: str, url: str, token: str) -> EmbedFn:
    """Builds a real embed_fn against Biohub's ESMC platform API.

    Mirrors common/http.py's retry pattern (bounded attempts, exponential backoff) around the
    SDK call, since the SDK's own client cannot be handed our CachedJsonClient/session.
    """
    from esm.sdk.forge import ESMCForgeInferenceClient
    from esm.sdk.api import ESMProtein, ESMProteinError, LogitsConfig

    client = ESMCForgeInferenceClient(model=model, url=url, token=token)

    def embed(sequence: str) -> EmbeddingResult:
        last_error = ""
        for attempt in range(MAX_ATTEMPTS):
            try:
                protein = ESMProtein(sequence=sequence)
                protein_tensor = client.encode(protein)
                if isinstance(protein_tensor, ESMProteinError):
                    raise RuntimeError(str(protein_tensor))
                output = client.logits(protein_tensor, LogitsConfig(sequence=True, return_embeddings=True))
                emb = output.embeddings
                pooled = emb.mean(dim=0) if hasattr(emb, "mean") else emb
                vector = pooled.tolist() if hasattr(pooled, "tolist") else list(pooled)
                return EmbeddingResult(status=EMBEDDED, vector=vector, model=model)
            except Exception as exc:  # noqa: BLE001 - SDK exception types are not all documented
                last_error = f"{type(exc).__name__}: {exc}"
                if attempt < MAX_ATTEMPTS - 1:
                    time.sleep(min(2.0**attempt, 8.0))
        return EmbeddingResult(status=QUERY_FAILED, model=model, error=last_error)

    return embed


def make_local_sae_fn(model: str, layer: int, device: str = "cpu") -> SaeFn:
    """Local Hugging-Face-weights SAE extraction. See the module docstring's SAE caveat."""
    import torch
    from esm.models.esmc import EsmcForMaskedLM, EsmcSaeModel, EsmcTokenizer

    hf_model = EsmcForMaskedLM.from_pretrained(model, device=device).eval()
    tokenizer = EsmcTokenizer()
    sae = EsmcSaeModel.from_pretrained(
        f"{model}-sae-k64-codebook16384",
        allow_patterns=["config.json", f"layer_{layer}.safetensors"],
        device=hf_model.device,
    )
    sae.initialize_layers([layer])
    hf_model.add_sae_models([sae.layers[str(layer)]])

    def extract(sequence: str) -> SaeResult:
        try:
            inputs = tokenizer(sequence, return_tensors="pt", padding=True)
            inputs = {k: v.to(hf_model.device) for k, v in inputs.items()}
            with torch.inference_mode():
                output = hf_model(**inputs)
            sparse = output.sae_outputs[f"layer{layer}"]
            dense = sparse.to_dense() if hasattr(sparse, "to_dense") else sparse
            pooled = dense.mean(dim=0) if dense.ndim > 1 else dense
            return SaeResult(status=SAE_EXTRACTED, vector=pooled.tolist(), model=model, layer=layer)
        except Exception as exc:  # noqa: BLE001
            return SaeResult(status=SAE_FAILED, model=model, layer=layer, error=f"{type(exc).__name__}: {exc}")

    return extract


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def _write_metadata_tsv(path: Path, rows: list[EsmcEmbedRow]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=METADATA_COLUMNS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row.metadata_row())


def _write_vectors_jsonl(path: Path, rows: list[EsmcEmbedRow]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for row in rows:
            if row.esmc_status == EMBEDDED:
                handle.write(json.dumps(row.vector_row()) + "\n")


def run(args: argparse.Namespace) -> int:
    run_dir = Path(args.run)
    m1_dir = FIXTURE_DIR if args.dry_run else (Path(args.m1) if args.m1 else run_dir / "m1")
    out_dir = run_dir / "m2_pdb"
    no_hit_path = Path(args.no_pdb_hit) if args.no_pdb_hit else out_dir / "no_pdb_hit.tsv"

    started = datetime.now(UTC)
    clock = time.monotonic()

    rows = load_no_pdb_hit_rows(no_hit_path) if no_hit_path.exists() else []
    bundle = load_input(m1_dir)
    if args.limit:
        rows = rows[: args.limit]
    sequence_by_id = {protein.feature_id: protein.sequence for protein in bundle.proteins}

    offline = args.offline or args.dry_run
    cache_path = Path(args.cache) if args.cache else (FIXTURE_DIR / "cache.sqlite" if args.dry_run else run_dir / "cache" / "m2_esmc.sqlite")

    embed_fn = args.embed_fn or make_esmc_embed_fn(args.model, args.url, args.token or "")
    sae_fn = None
    if args.sae:
        sae_fn = args.sae_fn or make_local_sae_fn(args.model, args.sae_layer, device=args.device)

    with JsonCache(cache_path) as cache:
        results = embed_eligible_proteins(
            rows,
            sequence_by_id,
            embed_fn=embed_fn,
            model=args.model,
            cache=cache,
            offline=offline,
            sae_fn=sae_fn,
            sae_layer=args.sae_layer,
        )

        out_dir.mkdir(parents=True, exist_ok=True)
        _write_metadata_tsv(out_dir / "esmc_embeddings.tsv", results)
        _write_vectors_jsonl(out_dir / "esmc_embeddings.jsonl", results)

        status_counts: dict[str, int] = {}
        for row in results:
            status_counts[row.esmc_status] = status_counts.get(row.esmc_status, 0) + 1

        manifest = {
            "module": "m2_esmc_embed",
            "issues": [9],
            "started_utc": started.isoformat(),
            "finished_utc": datetime.now(UTC).isoformat(),
            "elapsed_seconds": round(time.monotonic() - clock, 2),
            "input_no_pdb_hit": str(no_hit_path),
            "input_m1": str(m1_dir),
            "output_dir": str(out_dir),
            "cache_path": str(cache_path),
            "offline": offline,
            "dry_run": args.dry_run,
            "model": args.model,
            "sae_requested": args.sae,
            "sae_model": args.model if args.sae else "",
            "sae_layer": args.sae_layer if args.sae else None,
            "counts": {
                "input_rows": len(rows),
                "embedded": status_counts.get(EMBEDDED, 0),
                "query_failed": status_counts.get(QUERY_FAILED, 0),
                "skipped_upstream_query_failed": status_counts.get(SKIPPED_QUERY_FAILED, 0),
                "not_applicable_found": status_counts.get(NOT_APPLICABLE_FOUND, 0),
                "skipped_no_sequence": status_counts.get(SKIPPED_NO_SEQUENCE, 0),
                "esmc_status": status_counts,
            },
        }
        (out_dir / "esmc_run.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print(
        f"{manifest['counts']['embedded']} embedded, "
        f"{manifest['counts']['query_failed']} failed, "
        f"{manifest['counts']['skipped_upstream_query_failed']} skipped (upstream query-failed), "
        f"{manifest['counts']['not_applicable_found']} not applicable (found)."
    )
    print(f"Wrote {out_dir / 'esmc_embeddings.tsv'} and {out_dir / 'esmc_embeddings.jsonl'}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="python -m s2f.m2_triage.esmc_embed", description=__doc__)
    parser.add_argument("--run", default="runs/dev", help="run directory (reads <run>/m2_pdb/no_pdb_hit.tsv, <run>/m1)")
    parser.add_argument("--m1", default="", help="override the M1 input directory (for sequence join)")
    parser.add_argument("--no-pdb-hit", default="", help="override the no_pdb_hit.tsv path")
    parser.add_argument("--limit", type=int, default=0, help="process only the first N rows")
    parser.add_argument("--model", default=DEFAULT_MODEL, help="ESMC model name")
    parser.add_argument("--url", default="https://biohub.ai", help="Biohub platform base URL")
    parser.add_argument("--token", default="", help="Biohub API token (else reads BIOHUB_API_TOKEN)")
    parser.add_argument("--cache", default="", help="override the SQLite cache path")
    parser.add_argument("--offline", action="store_true", help="replay from cache; never call the network")
    parser.add_argument("--dry-run", action="store_true", help="run from fixtures/m2 with no network")
    parser.add_argument("--sae", action="store_true", help="also extract SAE features (local HF inference)")
    parser.add_argument("--sae-layer", type=int, default=DEFAULT_SAE_LAYER, help="SAE layer to extract")
    parser.add_argument("--device", default="cpu", help="device for local SAE inference (cpu or cuda)")
    # Test-only injection points; not exposed as real CLI flags in --help-driven usage.
    parser.set_defaults(embed_fn=None, sae_fn=None)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.dry_run and args.run == "runs/dev":
        args.run = "runs/dry-run"
    import os

    if not args.token:
        args.token = os.environ.get("BIOHUB_API_TOKEN", "")
    try:
        return run(args)
    except (OSError, ValueError) as exc:
        print(f"Error: {exc}")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

"""Tests for s2f.m2_triage.esmc_embed (issue #9, ESMC embedding half).

Real network / real esm-SDK / real model weights are never touched here: embed_fn and sae_fn
are always fakes. What's under test is the control flow this module owns: eligibility
filtering (no-hit vs query-failed vs found), sequence join against M1 input, caching, and
--offline behaviour.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from s2f.common.http import JsonCache, OfflineCacheMiss
from s2f.m2_triage.esmc_embed import (
    EMBEDDED,
    NOT_APPLICABLE_FOUND,
    QUERY_FAILED,
    SAE_EXTRACTED,
    SAE_FAILED,
    SAE_NOT_REQUESTED,
    SKIPPED_QUERY_FAILED,
    EmbeddingResult,
    SaeResult,
    embed_eligible_proteins,
    load_no_pdb_hit_rows,
)

GENOME_FIXTURE = Path(__file__).resolve().parents[1] / "fixtures" / "genomes" / "mgen_G37"


def make_row(feature_id: str, retrieval_status: str, product: str = "hypothetical protein", **extra) -> dict[str, str]:
    row = {
        "feature_id": feature_id,
        "product": product,
        "retrieval_status": retrieval_status,
        "virulence_amr": "0.0",
        "essential": "0.0",
        "drug_target": "0.0",
        "annotation_gap": "1.0",
        "human_homolog_penalty": "0.0",
    }
    row.update(extra)
    return row


def fake_embed_ok(sequence: str) -> EmbeddingResult:
    # A tiny deterministic "embedding" so tests can assert on it without a real model.
    return EmbeddingResult(status=EMBEDDED, vector=[float(len(sequence)), 1.0, 2.0], model="fake-esmc")


def fake_embed_fails(sequence: str) -> EmbeddingResult:
    return EmbeddingResult(status=QUERY_FAILED, model="fake-esmc", error="simulated network failure")


def failing_embed_fn(sequence: str) -> EmbeddingResult:
    raise AssertionError("embed_fn must not be called for ineligible or offline-cache-miss rows")


def fake_sae_ok(sequence: str) -> SaeResult:
    return SaeResult(status=SAE_EXTRACTED, vector=[0.1, 0.2, 0.3], model="fake-sae", layer=60)


# --------------------------------------------------------------------------
# Eligibility filtering
# --------------------------------------------------------------------------

def test_empty_input_produces_no_rows() -> None:
    results = embed_eligible_proteins([], {}, embed_fn=failing_embed_fn)
    assert results == []


def test_found_rows_are_skipped_as_not_applicable() -> None:
    rows = [make_row("fig|1.1.peg.1", "found")]
    results = embed_eligible_proteins(rows, {"fig|1.1.peg.1": "MKALIV"}, embed_fn=failing_embed_fn)

    assert len(results) == 1
    assert results[0].esmc_status == NOT_APPLICABLE_FOUND
    assert results[0].embedding == []


def test_query_failed_rows_are_recorded_but_never_embedded() -> None:
    """A failed PDB search is not confirmed evidence of "no homolog" (pitfall: never collapse
    a failure into an absence of hits). The embed_fn asserting on call proves this."""
    rows = [make_row("fig|1.1.peg.2", "query-failed")]
    results = embed_eligible_proteins(rows, {"fig|1.1.peg.2": "MKALIV"}, embed_fn=failing_embed_fn)

    assert len(results) == 1
    assert results[0].esmc_status == SKIPPED_QUERY_FAILED
    assert results[0].retrieval_status == "query-failed"
    assert results[0].embedding == []


def test_no_pdb_hit_tsv_mixes_no_hit_and_query_failed_and_both_are_handled(tmp_path: Path) -> None:
    """Mirrors the real pipeline: no_pdb_hit.tsv's `no_pdb_hit` flag is true for both
    retrieval_status values (score.py), so this module must filter on retrieval_status itself."""
    path = tmp_path / "no_pdb_hit.tsv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(make_row("x", "no-hit").keys()), delimiter="\t")
        writer.writeheader()
        writer.writerow(make_row("fig|1.1.peg.1", "no-hit"))
        writer.writerow(make_row("fig|1.1.peg.2", "query-failed"))

    rows = load_no_pdb_hit_rows(path)
    sequences = {"fig|1.1.peg.1": "MKALIV", "fig|1.1.peg.2": "MKALIV"}
    results = embed_eligible_proteins(rows, sequences, embed_fn=fake_embed_ok)

    by_id = {r.feature_id: r for r in results}
    assert by_id["fig|1.1.peg.1"].esmc_status == EMBEDDED
    assert by_id["fig|1.1.peg.2"].esmc_status == SKIPPED_QUERY_FAILED


def test_missing_sequence_is_recorded_not_silently_dropped() -> None:
    rows = [make_row("fig|1.1.peg.9", "no-hit")]
    results = embed_eligible_proteins(rows, {}, embed_fn=failing_embed_fn)

    assert len(results) == 1
    assert results[0].esmc_status == "skipped-no-sequence"
    assert "feature_id" in results[0].error or "sequence" in results[0].error


# --------------------------------------------------------------------------
# Successful embedding + metadata carried through
# --------------------------------------------------------------------------

def test_no_hit_row_is_embedded_and_specialty_flags_survive() -> None:
    rows = [make_row("fig|1.1.peg.3", "no-hit", virulence_amr="1.0", annotation_gap="1.0")]
    results = embed_eligible_proteins(rows, {"fig|1.1.peg.3": "MKALIVAAA"}, embed_fn=fake_embed_ok)

    assert len(results) == 1
    row = results[0]
    assert row.esmc_status == EMBEDDED
    assert row.embedding_dim == 3
    assert row.embedding == [9.0, 1.0, 2.0]
    assert row.specialty["virulence_amr"] == "1.0"
    assert row.specialty["annotation_gap"] == "1.0"
    assert row.sae_status == SAE_NOT_REQUESTED  # SAE off by default


def test_embedding_failure_is_recorded_as_query_failed_not_dropped() -> None:
    rows = [make_row("fig|1.1.peg.4", "no-hit")]
    results = embed_eligible_proteins(rows, {"fig|1.1.peg.4": "MKALIV"}, embed_fn=fake_embed_fails)

    assert len(results) == 1
    assert results[0].esmc_status == QUERY_FAILED
    assert "simulated network failure" in results[0].error


# --------------------------------------------------------------------------
# Caching + offline
# --------------------------------------------------------------------------

def test_offline_with_empty_cache_records_failure_and_never_calls_embed_fn(tmp_path: Path) -> None:
    rows = [make_row("fig|1.1.peg.5", "no-hit")]
    with JsonCache(tmp_path / "cache.sqlite") as cache:
        results = embed_eligible_proteins(
            rows, {"fig|1.1.peg.5": "MKALIV"}, embed_fn=failing_embed_fn, cache=cache, offline=True
        )

    assert len(results) == 1
    assert results[0].esmc_status == QUERY_FAILED
    assert "no cached ESMC embedding" in results[0].error


def test_offline_replays_a_previously_cached_embedding(tmp_path: Path) -> None:
    rows = [make_row("fig|1.1.peg.6", "no-hit")]
    sequences = {"fig|1.1.peg.6": "MKALIV"}
    cache_path = tmp_path / "cache.sqlite"

    with JsonCache(cache_path) as cache:
        online = embed_eligible_proteins(rows, sequences, embed_fn=fake_embed_ok, cache=cache, offline=False)
    assert online[0].esmc_status == EMBEDDED

    with JsonCache(cache_path) as cache:
        offline = embed_eligible_proteins(rows, sequences, embed_fn=failing_embed_fn, cache=cache, offline=True)
    assert offline[0].esmc_status == EMBEDDED
    assert offline[0].embedding == online[0].embedding


# --------------------------------------------------------------------------
# Optional SAE extraction
# --------------------------------------------------------------------------

def test_sae_is_off_by_default_and_on_when_sae_fn_given() -> None:
    rows = [make_row("fig|1.1.peg.7", "no-hit")]
    sequences = {"fig|1.1.peg.7": "MKALIV"}

    without_sae = embed_eligible_proteins(rows, sequences, embed_fn=fake_embed_ok)
    assert without_sae[0].sae_status == SAE_NOT_REQUESTED
    assert without_sae[0].sae == []

    with_sae = embed_eligible_proteins(rows, sequences, embed_fn=fake_embed_ok, sae_fn=fake_sae_ok)
    assert with_sae[0].sae_status == SAE_EXTRACTED
    assert with_sae[0].sae == [0.1, 0.2, 0.3]
    assert with_sae[0].sae_layer == 60


def test_sae_is_never_attempted_when_embedding_itself_failed() -> None:
    calls: list[str] = []

    def counting_sae_fn(sequence: str) -> SaeResult:
        calls.append(sequence)
        return fake_sae_ok(sequence)

    rows = [make_row("fig|1.1.peg.8", "no-hit")]
    embed_eligible_proteins(
        rows, {"fig|1.1.peg.8": "MKALIV"}, embed_fn=fake_embed_fails, sae_fn=counting_sae_fn
    )
    assert calls == []


# --------------------------------------------------------------------------
# G37 fixture integration: real 542-protein sequence join
# --------------------------------------------------------------------------

@pytest.mark.skipif(not GENOME_FIXTURE.exists(), reason="shared genome fixture not present")
def test_sequences_join_correctly_against_the_g37_fixture() -> None:
    from s2f.m2_triage.bvbrc_input import load_input

    bundle = load_input(GENOME_FIXTURE)
    sequence_by_id = {p.feature_id: p.sequence for p in bundle.proteins}

    # Build a synthetic no_pdb_hit.tsv-shaped input: one genuine no-hit, one query-failed,
    # one found, using real feature_ids and real sequences from the fixture.
    sample = bundle.proteins[:3]
    rows = [
        make_row(sample[0].feature_id, "no-hit", product=sample[0].product),
        make_row(sample[1].feature_id, "query-failed", product=sample[1].product),
        make_row(sample[2].feature_id, "found", product=sample[2].product),
    ]

    results = embed_eligible_proteins(rows, sequence_by_id, embed_fn=fake_embed_ok)
    by_id = {r.feature_id: r for r in results}

    assert by_id[sample[0].feature_id].esmc_status == EMBEDDED
    assert by_id[sample[0].feature_id].embedding_dim == 3
    assert by_id[sample[1].feature_id].esmc_status == SKIPPED_QUERY_FAILED
    assert by_id[sample[2].feature_id].esmc_status == NOT_APPLICABLE_FOUND

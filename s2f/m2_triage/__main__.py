"""M2 PDB-evidence triage CLI.

    python -m s2f.m2_triage --run runs/hs11286 [--limit N] [--workers 4] [--offline] [--dry-run]

Reads M1 output from ``<run>/m1/`` and writes ``<run>/m2_pdb/``. See
``docs/02a-m2-pdb-evidence.md`` for parameters, weights and the output contract.
"""

from __future__ import annotations

import argparse
import csv
import json
import time
from datetime import UTC, datetime
from pathlib import Path

from ..common.http import CachedJsonClient, JsonCache
from ..common.ids import IdMapper
from ..common.io import init_report, read_report, update_proteins, update_section
from .report_adapter import kg_section, protein_enrichment, seed_proteins
from .afdb import AlphaFoldClient
from .kg import KnowledgeGraphBuilder
from .bvbrc_input import load_input
from .pdb_evidence import (
    EVALUE_CUTOFF,
    GRAPHQL_URL,
    IDENTITY_CUTOFF,
    ROWS_PER_QUERY,
    SEARCH_URL,
    PdbEvidenceClient,
    apply_metadata,
)
from .score import (
    ANNOTATION_BONUS,
    WEIGHTS_VERSION,
    LIGAND_BONUS,
    MIN_EVALUE,
    MIN_QUERY_COVERAGE,
    PARTNER_BONUS,
    WEIGHTS,
    ProteinScore,
    rank_and_select,
    score_protein,
)

FIXTURE_DIR = Path(__file__).resolve().parents[2] / "fixtures" / "m2"

HIT_COLUMNS = [
    "feature_id", "entity_id", "entry_id", "identity", "query_coverage", "evalue", "bitscore",
    "experimental_method", "resolution", "description", "organism", "taxonomy_id",
    "uniprot_ids", "ligands", "metals", "partner_entities", "has_ec", "has_go", "hit_score",
]

PROTEIN_COLUMNS = [
    "rank", "feature_id", "product", "gene", "locus_tag", "pgfam", "triage_score",
    "pdb_evidence", "virulence_amr", "essential", "drug_target", "annotation_gap",
    "human_homolog_penalty", "retrieval_status", "best_entity_id", "best_identity",
    "best_coverage", "best_resolution", "best_method", "qualifying_hits", "distinct_entries",
    "has_ligand_in_entry", "holo_homolog", "metals_in_entry", "human_pdb_hit",
    "human_homolog_identity", "transporter", "metal_resistance", "no_pdb_hit",
    "pdb_hit_organism", "uniprot_of_best_hit", "selected", "reason", "error",
]

# Added only when --map-ids runs; see s2f/common/ids.py (issue #3) and afdb.py (issue #8).
XREF_COLUMNS = [
    "uniprot", "uniprot_route", "uniparc", "gene_name", "chembl", "xref_note",
    "afdb_status", "afdb_entry", "afdb_plddt", "afdb_confidence", "afdb_coverage",
    "afdb_usable", "afdb_reason", "afdb_cif_url",
]


def _write_tsv(path: Path, columns: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _protein_row(score: ProteinScore, protein_meta: dict[str, str]) -> dict[str, object]:
    top = score.best_hit
    components = score.components.as_dict()
    row: dict[str, object] = {
        "rank": score.rank,
        "feature_id": score.feature_id,
        "product": score.product,
        "gene": protein_meta.get("gene", ""),
        "locus_tag": protein_meta.get("locus_tag", ""),
        "pgfam": protein_meta.get("pgfam", ""),
        "triage_score": round(score.score, 6),
        "retrieval_status": score.retrieval_status,
        "best_entity_id": top.entity_id if top else "",
        "best_identity": round(top.identity, 4) if top else "",
        "best_coverage": round(top.query_coverage, 4) if top else "",
        "best_resolution": top.resolution if top and top.resolution is not None else "",
        "best_method": top.experimental_method if top else "",
        "qualifying_hits": score.qualifying_hits,
        "distinct_entries": score.distinct_entries,
        "selected": score.selected,
        "reason": score.reason,
        "error": score.error,
    }
    row.update({name: round(value, 6) for name, value in components.items()})
    row.update(score.flags)
    return row


def _hit_rows(feature_id: str, hits, scorer) -> list[dict[str, object]]:
    rows = []
    for hit in hits:
        rows.append(
            {
                "feature_id": feature_id,
                "entity_id": hit.entity_id,
                "entry_id": hit.entry_id,
                "identity": round(hit.identity, 4),
                "query_coverage": round(hit.query_coverage, 4),
                "evalue": hit.evalue,
                "bitscore": hit.bitscore,
                "experimental_method": hit.experimental_method,
                "resolution": hit.resolution if hit.resolution is not None else "",
                "description": hit.description,
                "organism": hit.organism,
                "taxonomy_id": hit.taxonomy_id if hit.taxonomy_id is not None else "",
                "uniprot_ids": ";".join(hit.uniprot_ids),
                "ligands": ";".join(hit.ligands),
                "metals": ";".join(hit.metals),
                "partner_entities": hit.partner_entities,
                "has_ec": hit.has_ec,
                "has_go": hit.has_go,
                "hit_score": round(scorer(hit), 6),
            }
        )
    return rows


def run(args: argparse.Namespace) -> int:
    from .score import hit_score  # local import keeps the module import list readable

    run_dir = Path(args.run)
    m1_dir = FIXTURE_DIR if args.dry_run else (Path(args.m1) if args.m1 else run_dir / "m1")
    out_dir = run_dir / "m2_pdb"
    out_dir.mkdir(parents=True, exist_ok=True)

    started = datetime.now(UTC)
    clock = time.monotonic()

    bundle = load_input(m1_dir)
    proteins = bundle.proteins[: args.limit] if args.limit else bundle.proteins
    groups = {}
    for protein in proteins:
        groups.setdefault(protein.sequence_hash, []).append(protein)
    sequences = {digest: members[0].sequence for digest, members in groups.items()}

    if args.cache:
        cache_path = Path(args.cache)
    elif args.dry_run:
        # The fixture cache is committed, which is what makes --dry-run work with no network.
        cache_path = FIXTURE_DIR / "cache.sqlite"
    else:
        cache_path = run_dir / "cache" / "m2_pdb.sqlite"
    offline = args.offline or args.dry_run
    with JsonCache(cache_path) as cache:
        client = CachedJsonClient(cache=cache, offline=offline)
        pdb = PdbEvidenceClient(client, rows=args.rows)

        results = pdb.search_sequences(sequences, workers=args.workers)

        all_hits = [hit for result in results.values() for hit in result.hits]
        metadata = pdb.fetch_entity_metadata(hit.entity_id for hit in all_hits)
        apply_metadata(all_hits, metadata)

        scores: list[ProteinScore] = []
        hit_rows: list[dict[str, object]] = []
        meta_by_id: dict[str, dict[str, str]] = {}
        for digest, members in groups.items():
            result = results.get(digest)
            for protein in members:
                status = result.status if result else "query-failed"
                error = result.error if result else "sequence was not searched"
                hits = list(result.hits) if result else []
                scores.append(score_protein(protein, hits, status, error=error if status == "query-failed" else ""))
                hit_rows.extend(_hit_rows(protein.feature_id, hits, hit_score))
                meta_by_id[protein.feature_id] = {
                    "gene": protein.gene,
                    "locus_tag": protein.locus_tag,
                    "pgfam": protein.pgfam,
                }

        xrefs_by_id: dict[str, object] = {}
        mapper: IdMapper | None = None
        if args.map_ids:
            # One mapper for the whole pipeline (pitfall #11): this module never maps its own IDs.
            mapper = IdMapper(client, taxon_id=args.taxon or None)
            by_feature = {p.feature_id: p for p in proteins}
            requests_ = [
                {
                    "feature_id": score.feature_id,
                    "sequence": by_feature[score.feature_id].sequence,
                    "locus_tag": by_feature[score.feature_id].locus_tag or None,
                    "pdb_hit_uniprot_ids": (
                        score.best_hit.uniprot_ids if score.best_hit is not None else ()
                    ),
                }
                for score in scores
            ]
            xrefs_by_id = mapper.map_many(requests_, workers=args.workers)

            # AlphaFold DB is keyed by UniProt accession, so it rides on the mapping step.
            # M3 must not re-fold what already exists (03-m3-fold.md).
            afdb = AlphaFoldClient(client)
            afdb_by_accession = afdb.lookup_many(
                (x.uniprot for x in xrefs_by_id.values() if x.uniprot), workers=args.workers
            )

        ranked = rank_and_select(scores, top_n=args.top)

        kg_summary: dict[str, object] = {"enabled": False}
        if args.kg:
            selected = [score for score in ranked if score.selected]
            builder = KnowledgeGraphBuilder(client, taxon_id=args.taxon or None)
            by_feature = {p.feature_id: p for p in proteins}
            kg_input = []
            for score in selected:
                protein = by_feature[score.feature_id]
                top = score.best_hit
                homologs = []
                if top is not None:
                    # The PDB hit's accessions are the homolog's, carrying the identity that
                    # justifies anything transferred from it (pitfall #4).
                    for accession in top.uniprot_ids:
                        homologs.append(
                            {
                                "accession": accession,
                                "identity": round(top.identity * 100, 1),
                                "coverage": round(top.query_coverage, 4),
                                "organism": top.organism,
                                "via": f"PDB sequence hit {top.entity_id}",
                                "name": top.description,
                            }
                        )
                kg_input.append(
                    {
                        "feature_id": score.feature_id,
                        "locus_tag": protein.locus_tag,
                        "product": protein.product,
                        "homologs": homologs,
                        "human_homolog_identity": protein.human_homolog_identity,
                    }
                )
            graph = builder.build(kg_input)
            kg_path = out_dir / "kg.json"
            kg_path.write_text(
                json.dumps(graph.to_dict(caps=builder.caps()), indent=2) + "\n", encoding="utf-8"
            )
            kg_summary = {
                "enabled": True,
                "proteins": len(kg_input),
                "path": str(kg_path),
                "caps": builder.caps(),
                "source_versions": graph.source_versions,
                "failures": len(graph.failures),
                **graph.summary(),
            }

        protein_rows = [_protein_row(score, meta_by_id.get(score.feature_id, {})) for score in ranked]
        columns = list(PROTEIN_COLUMNS)
        if args.map_ids:
            columns += XREF_COLUMNS
            for row in protein_rows:
                xrefs = xrefs_by_id.get(row["feature_id"])
                model = afdb_by_accession.get(xrefs.uniprot) if xrefs and xrefs.uniprot else None
                usable, reason = model.usable_for_docking() if model else (False, "not queried")
                row.update(
                    {
                        "uniprot": xrefs.uniprot or "" if xrefs else "",
                        "uniprot_route": xrefs.route if xrefs else "",
                        "uniparc": (xrefs.uniparc or "") if xrefs else "",
                        "gene_name": (xrefs.gene_name or "") if xrefs else "",
                        "chembl": ";".join(xrefs.chembl) if xrefs else "",
                        "xref_note": xrefs.note if xrefs else "",
                        "afdb_status": model.status if model else "not-queried",
                        "afdb_entry": model.entry_id if model else "",
                        "afdb_plddt": model.mean_plddt if model else "",
                        "afdb_confidence": model.confidence_band if model else "",
                        "afdb_coverage": model.coverage if model and model.coverage is not None else "",
                        "afdb_usable": usable,
                        "afdb_reason": reason,
                        "afdb_cif_url": model.cif_url if model else "",
                    }
                )
        _write_tsv(out_dir / "proteins.tsv", columns, protein_rows)
        _write_tsv(out_dir / "hits.tsv", HIT_COLUMNS, hit_rows)
        _write_tsv(out_dir / f"top{args.top}.tsv", columns, [r for r in protein_rows if r["selected"]])
        _write_tsv(out_dir / "no_pdb_hit.tsv", columns, [r for r in protein_rows if r["no_pdb_hit"]])

        status_counts: dict[str, int] = {}
        for score in ranked:
            status_counts[score.retrieval_status] = status_counts.get(score.retrieval_status, 0) + 1

        if args.report:
            # The shared contract (00-architecture.md). M1 owns proteins[]; M2 only enriches it,
            # so seed the list when M1 has not run for this run directory yet.
            init_report(run_dir, run_dir.name, input_file=str(m1_dir))
            existing = read_report(run_dir).get("proteins") or []
            if not existing:
                update_section(run_dir, "proteins", seed_proteins(proteins))
            update_proteins(
                run_dir,
                {
                    score.feature_id: protein_enrichment(
                        score,
                        xrefs=xrefs_by_id.get(score.feature_id) if args.map_ids else None,
                        model=(
                            afdb_by_accession.get(xrefs_by_id[score.feature_id].uniprot)
                            if args.map_ids
                            and xrefs_by_id.get(score.feature_id)
                            and xrefs_by_id[score.feature_id].uniprot
                            else None
                        ),
                        weights_version=WEIGHTS_VERSION,
                    )
                    for score in ranked
                },
            )
            if args.kg and kg_summary.get("enabled"):
                update_section(
                    run_dir, "kg", kg_section(graph.to_dict(caps=builder.caps()), path=str(kg_path))
                )

        manifest = {
            "module": "m2_pdb_evidence",
            "docs": "docs/02a-m2-pdb-evidence.md",
            "issues": [8, 12],
            "started_utc": started.isoformat(),
            "finished_utc": datetime.now(UTC).isoformat(),
            "elapsed_seconds": round(time.monotonic() - clock, 2),
            "input_dir": str(m1_dir),
            "output_dir": str(out_dir),
            "cache_path": str(cache_path),
            "offline": offline,
            "dry_run": args.dry_run,
            "parameters": {
                "search_url": SEARCH_URL,
                "graphql_url": GRAPHQL_URL,
                "evalue_cutoff": EVALUE_CUTOFF,
                "identity_cutoff": IDENTITY_CUTOFF,
                "rows_per_query": args.rows,
                "workers": args.workers,
                "limit": args.limit,
                "top_n": args.top,
                "score_gate": {"min_evalue": MIN_EVALUE, "min_query_coverage": MIN_QUERY_COVERAGE},
                "bonuses": {
                    "ligand": LIGAND_BONUS,
                    "partner": PARTNER_BONUS,
                    "annotation": ANNOTATION_BONUS,
                },
            },
            "weights": WEIGHTS,
            "counts": {
                "proteins_loaded": len(bundle.proteins),
                "proteins_scored": len(ranked),
                "unique_sequences": len(sequences),
                "hits_retained": len(hit_rows),
                "entities_with_metadata": len(metadata),
                "selected": sum(1 for r in protein_rows if r["selected"]),
                "no_pdb_hit": sum(1 for r in protein_rows if r["no_pdb_hit"]),
                "retrieval_status": status_counts,
                "specialty_rows_without_protein": len(bundle.specialty_without_protein),
                "table_rows_without_sequence": len(bundle.table_rows_without_sequence),
                "sequences_without_table_row": len(bundle.sequences_without_table_row),
            },
            "id_mapping": (
                {
                    "enabled": True,
                    "taxon_id": args.taxon or None,
                    "mapped": sum(1 for x in xrefs_by_id.values() if x.mapped),
                    "unmapped": sum(1 for x in xrefs_by_id.values() if not x.mapped),
                    "routes": {
                        route: sum(1 for x in xrefs_by_id.values() if x.route == route)
                        for route in sorted({x.route for x in xrefs_by_id.values()})
                    },
                    "alphafold": {
                        "accessions_queried": len(afdb_by_accession),
                        "status": {
                            status: sum(1 for m in afdb_by_accession.values() if m.status == status)
                            for status in sorted({m.status for m in afdb_by_accession.values()})
                        },
                        "confidence_band": {
                            band: sum(1 for m in afdb_by_accession.values() if m.confidence_band == band)
                            for band in sorted({m.confidence_band for m in afdb_by_accession.values() if m.found})
                        },
                        "usable_for_docking": sum(
                            1 for m in afdb_by_accession.values() if m.usable_for_docking()[0]
                        ),
                        "failures": len(afdb.failures),
                    },
                }
                if args.map_ids
                else {"enabled": False}
            ),
            "knowledge_graph": kg_summary,
            "failures": {
                "search": [
                    {"sequence_sha256": digest, "error": result.error}
                    for digest, result in results.items()
                    if result.status == "query-failed"
                ],
                "metadata": pdb.metadata_failures,
                "id_mapping": mapper.failures if mapper is not None else [],
            },
        }
        (out_dir / "run.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print(
        f"{len(ranked)} proteins scored from {len(sequences)} unique sequences; "
        f"{manifest['counts']['selected']} selected, {manifest['counts']['no_pdb_hit']} without a PDB hit."
    )
    print(f"Wrote {out_dir}")
    failed = len(manifest["failures"]["search"])
    if failed:
        print(f"Warning: {failed} sequence searches failed (recorded as query-failed in run.json).")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="python -m s2f.m2_triage", description=__doc__)
    parser.add_argument("--run", default="runs/dev", help="run directory (reads <run>/m1, writes <run>/m2_pdb)")
    parser.add_argument("--m1", default="", help="override the M1 input directory")
    parser.add_argument("--limit", type=int, default=0, help="score only the first N proteins")
    parser.add_argument("--top", type=int, default=50, help="how many proteins to select")
    parser.add_argument("--rows", type=int, default=ROWS_PER_QUERY, help="hits requested per search")
    parser.add_argument("--workers", type=int, default=4, help="concurrent RCSB searches")
    parser.add_argument("--cache", default="", help="override the SQLite cache path")
    parser.add_argument("--offline", action="store_true", help="replay from cache; never call the network")
    parser.add_argument("--dry-run", action="store_true", help="run from fixtures/m2 with no network")
    parser.add_argument(
        "--map-ids",
        action="store_true",
        help="resolve UniProt/UniParc/ChEMBL xrefs via s2f.common.ids (issue #3)",
    )
    parser.add_argument("--taxon", type=int, default=0, help="NCBI taxon ID of the genome, for id mapping and STRING")
    parser.add_argument(
        "--report",
        action="store_true",
        help="also write into runs/<id>/report.json through the shared contract (issue #2)",
    )
    parser.add_argument(
        "--kg",
        action="store_true",
        help="assemble the capped knowledge subgraph for selected proteins (issue #11)",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.dry_run and args.run == "runs/dev":
        args.run = "runs/dry-run"
    try:
        return run(args)
    except (OSError, ValueError) as exc:
        print(f"Error: {exc}")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

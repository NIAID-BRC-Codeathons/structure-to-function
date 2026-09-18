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
from datetime import datetime, timezone
from pathlib import Path

from ..common.http import CachedJsonClient, JsonCache
from ..common.ids import IdMapper
from ..common.io import init_report, read_report, update_proteins, update_section
from .report_adapter import kg_section, protein_enrichment, seed_proteins
from .afdb import AlphaFoldClient
from .human_homology import DEFAULT_MIN_COVERAGE, DiamondMissing, fetch_human_proteome
from .human_homology import search as human_homology_search
from .essentiality import DEFAULT_REFERENCE_LIMIT, resolve_reference_set
from .essentiality import search as essentiality_search
from .foldseek import (
    DEFAULT_DATABASES,
    DEFAULT_MAX_EVALUE,
    DEFAULT_MIN_COVERAGE as FOLDSEEK_MIN_COVERAGE,
    FoldseekClient,
    fetch_structure,
    summarize as foldseek_summarize,
)
from .kg import KnowledgeGraphBuilder
from .bvbrc_input import load_input
from .function import (
    FLAG_COLUMNS,
    InterProScanInstall,
    UniProtFunctionClient,
    annotate,
    discover_interproscan,
    parse_deeptmhmm,
    parse_eggnog,
    parse_interproscan,
    parse_psortb,
    parse_signalp6,
    describe_file,
    run_interproscan,
    write_terms_tsv,
    write_tool_fasta,
)
from .pdb_evidence import (
    EVALUE_CUTOFF,
    GRAPHQL_URL,
    IDENTITY_CUTOFF,
    ROWS_PER_QUERY,
    SEARCH_URL,
    PdbEvidenceClient,
    apply_metadata,
)
from . import taxonomy as taxonomy_mod
from .score import (
    ANNOTATION_BONUS,
    OPTIONAL_COMPONENTS,
    best_hit,
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
    "surface_bonus", "ligandable_homolog", "membrane_penalty", "human_homolog_penalty",
    "ligandable_basis",
    "amr_basis", "antibiotic_target_not_resistance",
    "surface_exposed_measured", "membrane_measured", "ligandable_measured", "retrieval_status", "best_entity_id", "best_identity",
    "best_coverage", "best_resolution", "best_method", "qualifying_hits", "distinct_entries",
    "has_ligand_in_entry", "holo_homolog", "metals_in_entry", "human_pdb_hit",
    "human_homolog_identity", "close_human_homolog", "human_homolog_source", "essential_source",
    "same_species_hit", "same_genus_hit",
    "foldseek_status", "foldseek_hit", "foldseek_description", "foldseek_evalue",
    "foldseek_identity", "foldseek_coverage", "foldseek_informative",
    "transporter", "metal_resistance", "no_pdb_hit", "pdb_search_failed",
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


def _failure_kinds(failures: list[dict[str, str]]) -> dict[str, int]:
    """Group failure messages so a pattern (rate limiting, timeouts) is visible at a glance."""
    kinds: dict[str, int] = {}
    for failure in failures:
        error = str(failure.get("error", ""))
        if "429" in error or "rate" in error.lower():
            kind = "rate-limited"
        elif "TIMEOUT" in error or "ended as" in error:
            kind = "ticket-not-complete"
        elif "timed out" in error.lower() or "timeout" in error.lower():
            kind = "client-timeout"
        elif "Connection" in error or "connection" in error:
            kind = "connection"
        else:
            kind = error[:60] or "unknown"
        kinds[kind] = kinds.get(kind, 0) + 1
    return dict(sorted(kinds.items(), key=lambda kv: -kv[1]))


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
    row.update(
        {
            "surface_exposed_measured": score.components_available.get("surface_bonus", False),
            "ligandable_measured": score.components_available.get("ligandable_homolog", False),
            "membrane_measured": score.components_available.get("membrane_penalty", False),
        }
    )
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


# Provider outputs the dry run reads, so --dry-run exercises the parsers and not only the
# built-in heuristic. Real runs pass explicit paths.
FIXTURE_ANNOTATION = {
    "eggnog": "eggnog.emapper.annotations",
    "deeptmhmm": "deeptmhmm.TMRs.gff3",
    "signalp": "signalp6_prediction_results.txt",
    "psortb": "psortb_terse.txt",
    "interproscan": "interproscan.tsv",
}

ANNOTATION_PARSERS = {
    "eggnog": parse_eggnog,
    "deeptmhmm": parse_deeptmhmm,
    "signalp": parse_signalp6,
    "psortb": parse_psortb,
    "interproscan": parse_interproscan,
}


def _protein_fasta(m1_dir: Path) -> Path | None:
    candidate = m1_dir / "proteins.faa"
    if candidate.exists():
        return candidate
    matches = sorted(m1_dir.glob("*.faa"))
    return matches[0] if matches else None


def _annotation_aliases(tool_fasta, xrefs_by_id: dict) -> dict[str, list[str]]:
    """Every identifier a provider might have echoed back for a protein.

    Two sources: the FASTA record that represented it (its own ID, or a duplicate's), and its
    UniProt accession when --map-ids resolved one.
    """
    aliases: dict[str, list[str]] = {}
    if tool_fasta is not None:
        for feature_id, representatives in tool_fasta.aliases.items():
            aliases.setdefault(feature_id, []).extend(representatives)
    for feature_id, xrefs in xrefs_by_id.items():
        if xrefs.uniprot:
            aliases.setdefault(feature_id, []).append(xrefs.uniprot)
    return aliases


def _load_annotation_providers(
    args: argparse.Namespace, m1_dir: Path, out_dir: Path, scan_fasta: Path | None = None
) -> tuple[dict[str, dict], object, list[str], dict[str, dict]]:
    """Resolve each provider's output file, parsing what exists and saying what did not.

    ``--interproscan auto`` looks for an installation (PATH, $INTERPROSCAN_HOME, the usual
    install directories) and runs it over the proteome when one is there; without one the
    module carries on with the other providers rather than failing.
    """
    notes: list[str] = []
    provider_files: dict[str, dict] = {}
    paths = {
        "eggnog": args.eggnog,
        "deeptmhmm": args.deeptmhmm,
        "signalp": args.signalp,
        "psortb": args.psortb,
        "interproscan": args.interproscan,
    }
    # Discovery executes `interproscan.sh --version`, so only look when asked to.
    install = (
        discover_interproscan([args.interproscan_path] if args.interproscan_path else [])
        if (paths["interproscan"] == "auto" or args.interproscan_path)
        else InterProScanInstall(error="not requested")
    )

    if paths["interproscan"] == "auto":
        # Prefer the cleaned, deduplicated FASTA this module wrote: InterProScan rejects stop
        # characters and its cost scales with unique sequences.
        fasta = scan_fasta if scan_fasta is not None and scan_fasta.exists() else _protein_fasta(m1_dir)
        if fasta is None:
            notes.append("interproscan auto: no protein FASTA to scan")
            paths["interproscan"] = ""
        elif install.available:
            notes.append(f"interproscan {install.version or 'unknown version'} at {install.path} ({install.found_via})")
            produced, error = run_interproscan(
                install,
                fasta,
                out_dir / "interproscan.tsv",
                applications=args.interproscan_appl,
                cpus=args.workers,
            )
            paths["interproscan"] = str(produced) if produced else ""
            if error:
                notes.append(f"interproscan run failed: {error}")
        else:
            notes.append(f"interproscan auto: {install.error}")
            paths["interproscan"] = ""

    if args.dry_run:
        for key, name in FIXTURE_ANNOTATION.items():
            if not paths[key]:
                candidate = FIXTURE_DIR / "annotation" / name
                if candidate.exists():
                    paths[key] = str(candidate)

    parsed: dict[str, dict] = {}
    for key, value in paths.items():
        if not value:
            continue
        path = Path(value)
        if not path.exists():
            notes.append(f"{key}: {path} does not exist, skipped")
            continue
        try:
            parsed[key] = ANNOTATION_PARSERS[key](path)
        except Exception as exc:  # noqa: BLE001 - an add-on must never lose the triage output
            notes.append(f"{key}: failed to parse {path} — {type(exc).__name__}: {exc}")
            continue
        # Recorded whether or not it parsed to anything: "this file, this digest, zero records"
        # is a more useful thing to find in a manifest than silence.
        provider_files[key] = describe_file(path)
        if not parsed[key]:
            # Parsed without error and produced nothing: almost always a format the parser did
            # not recognise. Silence here would look like "the tool found nothing".
            notes.append(f"{key}: {path} parsed to zero records — check the file format")
    return parsed, install, notes, provider_files


def run(args: argparse.Namespace) -> int:
    from .score import hit_score  # local import keeps the module import list readable

    run_dir = Path(args.run)
    m1_dir = FIXTURE_DIR if args.dry_run else (Path(args.m1) if args.m1 else run_dir / "m1")
    out_dir = run_dir / "m2_pdb"
    out_dir.mkdir(parents=True, exist_ok=True)

    started = datetime.now(timezone.utc)
    clock = time.monotonic()

    bundle = load_input(m1_dir)

    # Where the query organism came from, so an undetermined one is never mistaken for a
    # determination. The p3-CLI/web export carries no organism in its FASTA headers, and on
    # that layout the same-organism flags below are all False for want of anything to compare
    # against — which reads exactly like "every hit is a genuine cross-organism transfer".
    # Say so instead of implying it.
    # M1 determines the organism and writes it to report.json; read it rather than infer it
    # from FASTA header brackets (issue #55). --organism overrides, and the header path stays
    # for the standalone entry point where there is no genome section to read.
    query_taxonomy = taxonomy_mod.resolve(
        run_dir, override=args.organism, header_organism=bundle.organism
    )
    query_organism = query_taxonomy.scientific_name
    query_organism_source = query_taxonomy.source
    if not query_taxonomy.determined:
        print(
            "Warning: query organism undetermined — same-species and same-genus hits cannot "
            "be told apart from cross-organism transfers. Run M1 first, or pass --organism."
        )

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

        human_hits = {}
        human_summary: dict[str, object] = {"enabled": False}
        if args.human_homology:
            # BV-BRC precomputes Human Homolog rows for public genomes only, so a blinded run
            # has to compute them or the penalty silently never fires (issue #29).
            try:
                reference = fetch_human_proteome(client, Path(args.human_proteome))
                found = human_homology_search(
                    [(p.feature_id, p.sequence) for p in proteins],
                    reference,
                    min_coverage=args.human_min_coverage,
                    threads=args.workers,
                )
                human_hits = found.hits
                human_summary = {"enabled": True, **found.summary()}
            except (DiamondMissing, RuntimeError, OSError) as exc:
                human_summary = {"enabled": False, "error": str(exc)}
                print(f"Warning: human-homology search skipped — {exc}")

        essential_calls = {}
        essentiality_summary: dict[str, object] = {"enabled": False}
        if args.essentiality:
            # BV-BRC's FBA essentiality exists for public genomes only, so transfer it from
            # public relatives and record which one justified each call (issue #29).
            try:
                keyword, reference = resolve_reference_set(
                    client,
                    keyword=args.essentiality_keyword,
                    taxonomy=query_taxonomy,
                    limit=args.essentiality_limit,
                )
                if not reference:
                    print(
                        f"Warning: essentiality found no reference set (keyword {keyword!r}); "
                        "every protein is recorded as not-run, not as non-essential."
                    )
                found = essentiality_search(
                    [(p.feature_id, p.sequence) for p in proteins],
                    reference,
                    reference_query=keyword,
                    min_identity=args.essentiality_min_identity,
                    min_coverage=args.essentiality_min_coverage,
                    threads=args.workers,
                )
                essential_calls = found.calls
                essentiality_summary = {"enabled": True, **found.summary()}
            except (DiamondMissing, RuntimeError, OSError) as exc:
                essentiality_summary = {"enabled": False, "error": str(exc)}
                print(f"Warning: essentiality transfer skipped — {exc}")

        xrefs_by_id: dict[str, object] = {}
        afdb_by_accession: dict[str, object] = {}
        mapper: IdMapper | None = None
        if args.map_ids:
            # One mapper for the whole pipeline (pitfall #11): this module never maps its own IDs.
            mapper = IdMapper(client, taxon_id=args.taxon or None)
            # The homolog accessions come straight from the search results rather than from
            # ProteinScore: mapping has to finish before functional annotation, which has to
            # finish before scoring, so this step cannot depend on scores existing yet.
            best_hit_accessions: dict[str, tuple[str, ...]] = {}
            for digest, members in groups.items():
                result = results.get(digest)
                top = best_hit(result.hits) if result and result.hits else None
                accessions = tuple(top.uniprot_ids) if top is not None else ()
                for member in members:
                    best_hit_accessions[member.feature_id] = accessions
            requests_ = [
                {
                    "feature_id": protein.feature_id,
                    "sequence": protein.sequence,
                    "locus_tag": protein.locus_tag or None,
                    "pdb_hit_uniprot_ids": best_hit_accessions.get(protein.feature_id, ()),
                }
                for protein in proteins
            ]
            xrefs_by_id = mapper.map_many(requests_, workers=args.workers)

            # AlphaFold DB is keyed by UniProt accession, so it rides on the mapping step.
            # M3 must not re-fold what already exists (03-m3-fold.md).
            afdb = AlphaFoldClient(client)
            afdb_by_accession = afdb.lookup_many(
                (x.uniprot for x in xrefs_by_id.values() if x.uniprot), workers=args.workers
            )

        # Functional annotation (issue #10). Runs *before* scoring, because #12 now scores
        # surface_exposed and membrane from these flags; it used to run after, when it was
        # flags-only and could not move the selection. The weight change is logged in
        # docs/02a-m2-pdb-evidence.md, per pitfall #12.
        annotation_run = None
        tool_fasta = None
        if args.annotate:
            # Written before the providers run, because --interproscan auto scans it and
            # because the two-pass flow (run M2, run the tools, re-run M2 with their output)
            # needs a FASTA the tools will accept. See docs/02d-remote-runbook.md.
            tool_fasta = write_tool_fasta(out_dir / "annotate_all.faa", proteins)
            parsed, ipr_install, annotation_notes, provider_files = _load_annotation_providers(
                args, m1_dir, out_dir, scan_fasta=tool_fasta.path
            )
            uniprot_results: dict[str, object] = {}
            uniprot_client: UniProtFunctionClient | None = None
            if args.uniprot_function:
                if not args.map_ids:
                    annotation_notes.append("--uniprot-function needs --map-ids for accessions; skipped")
                else:
                    uniprot_client = UniProtFunctionClient(client)
                    uniprot_results = uniprot_client.lookup_many(
                        [(fid, x.uniprot) for fid, x in xrefs_by_id.items() if x.uniprot],
                        workers=args.workers,
                    )
            annotation_run = annotate(
                proteins,
                eggnog=parsed.get("eggnog"),
                deeptmhmm=parsed.get("deeptmhmm"),
                signalp=parsed.get("signalp"),
                psortb=parsed.get("psortb"),
                interproscan=parsed.get("interproscan"),
                uniprot=uniprot_results,
                use_heuristic=not args.no_heuristic,
                aliases=_annotation_aliases(tool_fasta, xrefs_by_id),
            )
            annotation_run.tool_fasta = tool_fasta
            annotation_run.provider_files = provider_files
            annotation_run.interproscan = ipr_install
            annotation_run.uniprot_failures = uniprot_client.failures if uniprot_client else []
            annotation_run.notes = annotation_notes

        scores: list[ProteinScore] = []
        hit_rows: list[dict[str, object]] = []
        meta_by_id: dict[str, dict[str, str]] = {}
        for digest, members in groups.items():
            result = results.get(digest)
            for protein in members:
                status = result.status if result else "query-failed"
                error = result.error if result else "sequence was not searched"
                hits = list(result.hits) if result else []
                human = human_hits.get(protein.feature_id)
                scores.append(
                    score_protein(
                        protein,
                        hits,
                        status,
                        error=error if status == "query-failed" else "",
                        human_identity=(human.identity if human is not None and human.counted else None),
                        human_homolog_source="diamond_human_proteome" if human is not None else "",
                        essential=(
                            essential_calls[protein.feature_id].essential
                            if protein.feature_id in essential_calls
                            else None
                        ),
                        essential_source=(
                            "ortholog_of_fba_essential" if protein.feature_id in essential_calls else ""
                        ),
                        annotation=(
                            annotation_run.annotations.get(protein.feature_id)
                            if annotation_run is not None
                            else None
                        ),
                        query_organism=query_organism,
                        query_taxonomy=query_taxonomy,
                        chembl_targets=(
                            xrefs_by_id[protein.feature_id].chembl
                            if xrefs_by_id.get(protein.feature_id)
                            else ()
                        ),
                    )
                )
                hit_rows.extend(_hit_rows(protein.feature_id, hits, hit_score))
                meta_by_id[protein.feature_id] = {
                    "gene": protein.gene,
                    "locus_tag": protein.locus_tag,
                    "pgfam": protein.pgfam,
                }

        foldseek_results: dict[str, object] = {}
        foldseek_summary: dict[str, object] = {"enabled": False}
        if args.foldseek:
            # Structure search only helps where a structure exists. The targets are proteins the
            # sequence search found nothing for but that do have a usable predicted model
            # (issue #9; the model comes from #8's AlphaFold lookup).
            fs = FoldseekClient(
                cache=cache,
                session=client.session,
                offline=offline,
                databases=tuple(args.foldseek_db),
            )
            structures_dir = run_dir / "structures"
            candidates = [
                score
                for score in scores
                if score.flags.get("no_pdb_hit")
                and afdb_by_accession.get(
                    (xrefs_by_id.get(score.feature_id).uniprot if xrefs_by_id.get(score.feature_id) else "") or ""
                )
            ]
            if args.foldseek_limit:
                candidates = candidates[: args.foldseek_limit]
            for score in candidates:
                xrefs = xrefs_by_id.get(score.feature_id)
                model = afdb_by_accession.get(xrefs.uniprot) if xrefs and xrefs.uniprot else None
                usable, reason = (
                    model.usable_for_structure_search() if model else (False, "no AlphaFold model")
                )
                if not model or not model.found or not usable:
                    # Recorded, never silently skipped: a protein we could not search is a
                    # different outcome from one we searched and found nothing for.
                    foldseek_results[score.feature_id] = fs.search_protein(
                        score.feature_id,
                        None,
                        structure_source="afdb",
                        structure_id=model.entry_id if model else "",
                    )
                    foldseek_results[score.feature_id].note = (
                        f"model not usable as a search query: {reason}"
                    )
                    continue
                try:
                    structure = fetch_structure(
                        client.session, model.cif_url, structures_dir / f"{model.entry_id}.cif"
                    ) if not offline else (structures_dir / f"{model.entry_id}.cif").read_bytes()
                except (OSError, Exception) as exc:  # noqa: BLE001 - recorded, never silent
                    foldseek_results[score.feature_id] = fs.search_protein(
                        score.feature_id, None, structure_source="afdb", structure_id=model.entry_id
                    )
                    fs.failures.append({"feature_id": score.feature_id, "error": f"structure fetch: {exc}"})
                    continue
                foldseek_results[score.feature_id] = fs.search_protein(
                    score.feature_id,
                    structure,
                    structure_source="afdb",
                    structure_id=model.entry_id,
                    max_evalue=args.foldseek_evalue,
                    min_coverage=args.foldseek_min_coverage,
                )
            foldseek_summary = {
                "enabled": True,
                "candidates": len(candidates),
                **foldseek_summarize(foldseek_results.values()),
                "failures": len(fs.failures),
                # The errors themselves, not just a count: a run that failed 121 searches needs
                # to say why, or the next person cannot tell rate limiting from a broken query.
                "failure_detail": fs.failures[:20],
                "failure_kinds": _failure_kinds(fs.failures),
                "rate_limited": fs.rate_limited,
                "note": (
                    "The public Foldseek server rate-limited this run. Searches that completed are "
                    "cached and will replay; the rest need the local binary with a downloaded "
                    "database, which is also the only route that yields a TM-score."
                    if fs.rate_limited or any("429" in f.get("error", "") for f in fs.failures)
                    else ""
                ),
            }

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
        if args.foldseek:
            for row in protein_rows:
                result = foldseek_results.get(row["feature_id"])
                hit = result.best_informative or result.best if result else None
                row.update({
                    "foldseek_status": result.status if result else "not-queried",
                    "foldseek_hit": hit.pdb_id if hit else "",
                    "foldseek_description": (hit.description[:80] if hit else ""),
                    "foldseek_evalue": hit.evalue if hit else "",
                    "foldseek_identity": hit.identity if hit else "",
                    "foldseek_coverage": round(hit.query_coverage, 4) if hit and hit.query_coverage else "",
                    "foldseek_informative": hit.informative if hit else "",
                })
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
        if annotation_run is not None:
            # The selected set is what the expensive tools should actually be run over:
            # DeepTMHMM and SignalP 6 over ~50 proteins is minutes, over a proteome is hours.
            selected_ids = {score.feature_id for score in ranked if score.selected}
            annotation_run.selected_fasta = write_tool_fasta(
                out_dir / "annotate_selected.faa",
                [protein for protein in proteins if protein.feature_id in selected_ids],
            )

        terms_written = 0
        if annotation_run is not None:
            columns += FLAG_COLUMNS
            for row in protein_rows:
                annotation = annotation_run.annotations.get(row["feature_id"])
                if annotation is not None:
                    row.update(annotation.flag_row())
            terms_written = write_terms_tsv(out_dir / "function_terms.tsv", annotation_run.annotations)

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
                        functional=(
                            annotation_run.annotations.get(score.feature_id)
                            if annotation_run is not None
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
            "finished_utc": datetime.now(timezone.utc).isoformat(),
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
            # Which components actually had evidence behind them. A component that measured
            # nothing contributed 0 to every protein, which is indistinguishable from "absent"
            # in the score alone — so a reader can tell a genome with no membrane proteins from
            # a run where nothing looked (pitfall #12).
            "components_available": {
                name: {
                    "measured_for": sum(1 for s_ in ranked if s_.components_available.get(name)),
                    "of": len(ranked),
                }
                for name in OPTIONAL_COMPONENTS
            },
            "query_organism": query_organism or None,
            "query_organism_source": query_organism_source,
            "query_taxonomy": query_taxonomy.as_dict(),
            # How much of the ranking rests on our own organism already being in the PDB. On a
            # novel genome this is 0; on the smoke-test genome it is most of the top of the list.
            "same_organism_hits": {
                "same_species_in_top_n": sum(
                    1 for s_ in ranked if s_.selected and s_.flags.get("same_species_hit")
                ),
                "same_species_total": sum(1 for s_ in ranked if s_.flags.get("same_species_hit")),
                "same_genus_in_top_n": sum(
                    1 for s_ in ranked if s_.selected and s_.flags.get("same_genus_hit")
                ),
                "same_genus_total": sum(1 for s_ in ranked if s_.flags.get("same_genus_hit")),
                # Zeros mean two different things, and only one of them is a measurement.
                "measured": query_taxonomy.determined,
                "note": (
                    "A hit against our own organism or genus is not a cross-organism transfer. "
                    "Weights calibrated on a genome with its own PDB structures will not carry "
                    "to a novel one."
                    if query_organism
                    else "Not measured: the query organism is unknown, so these counts are "
                    "absence of evidence, not evidence of absence. Pass --organism."
                ),
            },
            "components_nonzero": {
                name: sum(1 for s_ in ranked if s_.components.as_dict().get(name))
                for name in WEIGHTS
            },
            "counts": {
                "proteins_loaded": len(bundle.proteins),
                "proteins_scored": len(ranked),
                "unique_sequences": len(sequences),
                "hits_retained": len(hit_rows),
                "entities_with_metadata": len(metadata),
                "selected": sum(1 for r in protein_rows if r["selected"]),
                "disqualified_human_homolog": sum(1 for s in ranked if s.disqualified),
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
            "foldseek": foldseek_summary,
            "human_homology": human_summary,
            "essentiality": essentiality_summary,
            "knowledge_graph": kg_summary,
            "functional_annotation": (
                {"enabled": True, "terms_written": terms_written, **annotation_run.summary()}
                if annotation_run is not None
                else {"enabled": False}
            ),
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
    parser.add_argument(
        "--organism",
        default="",
        help="query organism name; read from the FASTA headers when not given",
    )
    parser.add_argument("--taxon", type=int, default=0, help="NCBI taxon ID of the genome, for id mapping and STRING")
    parser.add_argument(
        "--human-homology",
        action="store_true",
        help="compute human homologs with DIAMOND instead of trusting BV-BRC's precomputed rows (issue #29)",
    )
    parser.add_argument(
        "--human-proteome",
        default="data/reference/human_UP000005640_reviewed.fasta",
        help="reviewed human proteome FASTA; downloaded on first use",
    )
    parser.add_argument(
        "--human-min-coverage",
        type=float,
        default=DEFAULT_MIN_COVERAGE,
        help="minimum query coverage for a hit to count as a homolog",
    )
    parser.add_argument(
        "--foldseek",
        action="store_true",
        help="structure search for proteins the sequence search missed, via the Foldseek web API (issue #9)",
    )
    parser.add_argument(
        "--foldseek-db", action="append", default=[], help="Foldseek database (repeatable; default pdb100)"
    )
    parser.add_argument("--foldseek-limit", type=int, default=0, help="search only the first N candidates")
    parser.add_argument("--foldseek-evalue", type=float, default=DEFAULT_MAX_EVALUE)
    parser.add_argument("--foldseek-min-coverage", type=float, default=FOLDSEEK_MIN_COVERAGE)
    parser.add_argument(
        "--essentiality",
        action="store_true",
        help="transfer FBA essentiality from public relatives by orthology (issue #29)",
    )
    parser.add_argument(
        "--essentiality-keyword",
        default="",
        help="genus or species selecting the public relatives, e.g. Mycoplasma",
    )
    parser.add_argument("--essentiality-limit", type=int, default=DEFAULT_REFERENCE_LIMIT)
    parser.add_argument("--essentiality-min-identity", type=float, default=40.0)
    parser.add_argument("--essentiality-min-coverage", type=float, default=0.7)
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
    group = parser.add_argument_group(
        "functional annotation (issue #10)",
        "See docs/02c-m2-functional-annotation.md for how to produce each provider's output.",
    )
    group.add_argument(
        "--annotate",
        action="store_true",
        help="add function, localization and membrane flags to every protein",
    )
    group.add_argument("--eggnog", default="", help="eggNOG-mapper *.emapper.annotations file")
    group.add_argument("--deeptmhmm", default="", help="DeepTMHMM TMRs.gff3 file")
    group.add_argument("--signalp", default="", help="SignalP 6 prediction_results.txt file")
    group.add_argument("--psortb", default="", help="PSORTb output (terse or long format)")
    group.add_argument(
        "--interproscan",
        default="",
        help="InterProScan TSV/JSON output, or 'auto' to find a local install and run it",
    )
    group.add_argument(
        "--interproscan-path", default="", help="path to interproscan.sh, if discovery misses it"
    )
    group.add_argument(
        "--interproscan-appl", default="", help="value for InterProScan's -appl (default: all)"
    )
    group.add_argument(
        "--uniprot-function",
        action="store_true",
        help="fetch curated function and topology from UniProt (needs --map-ids)",
    )
    group.add_argument(
        "--no-heuristic",
        action="store_true",
        help="do not fall back to the built-in sequence heuristic; leave flags unset instead",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if not args.foldseek_db:
        args.foldseek_db = list(DEFAULT_DATABASES)
    if args.foldseek and not args.map_ids:
        # Foldseek searches with AlphaFold models, which are keyed by UniProt accession.
        print("Error: --foldseek needs --map-ids (structures are found via UniProt accessions).")
        return 2
    if args.dry_run and args.run == "runs/dev":
        args.run = "runs/dry-run"
    try:
        return run(args)
    except (OSError, ValueError) as exc:
        print(f"Error: {exc}")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

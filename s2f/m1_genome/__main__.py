"""M1: assembly in, annotated genome out (issue #5).

Two acquisition routes, one output contract. Both end at the same `report.json` sections
and the same `<run>/m1/` files, so M2 never needs to know which one produced a run.

**CGA route** — annotates *this* assembly with the BV-BRC Comprehensive Genome Analysis
service. Authoritative, needs an account and the `p3-*` CLI, takes hours:

    python -m s2f.m1_genome --run runs/<id> --contigs assembly.fna \\
        --ws-dir /USER@bvbrc/home/s2f

Order matters and is the point of this route: Minhash predicts the taxon *before* CGA runs,
because CGA's annotation quality depends on it and CGA cannot infer it. Annotating this
project's test genome with no taxon doubles the CDS count and halves mean protein length
(`docs/01a-cga-coverage.md`).

Offline re-parse of a job already retrieved:

    python -m s2f.m1_genome --run runs/<id> --from-cga-dir data/cga

**API route** — reads a genome BV-BRC has *already* annotated, from the public Data API.
No account, seconds not hours, but it characterises a reference genome for the organism
rather than the sample in hand, so it cannot support a claim about this isolate
(`docs/01b-m1-api-mode.md`):

    python -m s2f.m1_genome --run runs/<id> --from-bvbrc-api --contigs assembly.fna
    python -m s2f.m1_genome --run runs/<id> --from-bvbrc-api --species "Klebsiella pneumoniae"
    python -m s2f.m1_genome --run runs/<id> --from-bvbrc-api --genome-id 1125630.4

Either route can write the interactive HTML report with `--html`.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from ..common.http import (CachedJsonClient, HttpError, JsonCache,
                           OfflineCacheMiss)
from ..common.io import init_report, update_proteins, update_section
from . import cga
from .parse import load_cga, genome_section, proteins_section, run_fields, write_m1_dir
from .taxon import TaxonCallError, call_taxon


def _now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="python -m s2f.m1_genome", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--run", required=True, help="run directory, e.g. runs/mgen_G37")
    p.add_argument("--contigs", help="assembly FASTA (required unless --from-cga-dir)")
    p.add_argument("--ws-dir", help="BV-BRC workspace dir, e.g. /USER@bvbrc/home/s2f")
    p.add_argument("--out-name", default=None, help="CGA output name (default: the run id)")
    p.add_argument("--from-cga-dir", help="parse an already-retrieved CGA directory; no network")
    p.add_argument("--job-id", help="resume: skip submission and poll this job id")
    p.add_argument("--taxon-id", type=int, help="skip the Minhash call and use this taxon")
    p.add_argument("--genetic-code", type=int, help="translation table; required with --taxon-id")
    p.add_argument("--scientific-name", help="overrides the name from the taxon call")
    p.add_argument("--max-distance", type=float, default=0.2, help="Minhash cut-off (default 0.2)")
    p.add_argument("--max-hits", type=int, default=50)
    p.add_argument("--reference-scope", action="store_true",
                   help="search only reference/representative sketches (usually too narrow)")
    p.add_argument("--no-blind", action="store_true", help="keep original FASTA headers")
    p.add_argument("--poll-interval", type=int, default=30)
    p.add_argument("--poll-timeout", type=int, default=7200)
    p.add_argument("--allow-poor", action="store_true",
                   help="continue even if CGA reports Poor quality or sets a quality flag")
    p.add_argument("--dry-run", action="store_true", help="print the CGA submission and stop")

    api = p.add_argument_group(
        "BV-BRC Data API route",
        "characterises an already-annotated reference genome; no account needed")
    api.add_argument("--from-bvbrc-api", action="store_true",
                     help="use the Data API instead of the CGA service")
    api.add_argument("--species", help="species name, e.g. 'Klebsiella pneumoniae'")
    api.add_argument("--genome-id", help="explicit BV-BRC genome id, e.g. 1125630.4")
    api.add_argument("--annotation", choices=["PATRIC", "RefSeq"], default="PATRIC",
                     help="annotation source for CDS features (default PATRIC/RASTtk)")
    api.add_argument("--max-features", type=int, default=25000,
                     help="cap on CDS features to fetch (default 25000)")
    api.add_argument("--seq-cap", type=int, default=4000,
                     help="max protein sequences to fetch (default 4000)")
    api.add_argument("--tree-leaves", type=int, default=14,
                     help="max genomes on the gene-content phylogeny (default 14)")
    api.add_argument("--no-tree", action="store_true",
                     help="skip the gene-content phylogeny")
    api.add_argument("--cache", default=None,
                     help="HTTP cache path (default: <run>/m1/cache.sqlite)")
    api.add_argument("--offline", action="store_true",
                     help="replay from the cache only; fail on a miss")
    api.add_argument("--timeout", type=float, default=90.0,
                     help="per-request timeout in seconds (default 90)")
    api.add_argument("--retries", type=int, default=4,
                     help="attempts per request before giving up (default 4)")
    api.add_argument("--ca-bundle",
                     help="CA bundle for TLS verification, for hosts behind an "
                          "inspecting proxy that presents its own certificate")
    api.add_argument("--no-proxy", action="store_true",
                     help="ignore HTTP_PROXY/HTTPS_PROXY for BV-BRC requests; use when a "
                          "stale or broken proxy intermittently drops connections")

    out = p.add_argument_group("reporting")
    out.add_argument("--html", action="store_true",
                     help="write the interactive HTML report to <run>/m1/report.html")
    out.add_argument("--figures", action="store_true",
                     help="also write matplotlib figures (needs matplotlib)")
    return p


def _api_session(args: argparse.Namespace):
    """A requests session configured for this site's network, or None for the default.

    Both knobs exist because they were needed in practice, not speculatively: an
    inspecting proxy that presents its own certificate makes every BV-BRC call fail
    verification until it is given the bundle, and a stale environment proxy drops
    connections intermittently, which looks exactly like BV-BRC being down.

    There is no `--insecure`: disabling verification turns a loud, diagnosable TLS error
    into a silent trust of whatever answered, and `--ca-bundle` solves the real case.
    """
    if not (args.ca_bundle or args.no_proxy):
        return None
    import requests

    session = requests.Session()
    if args.ca_bundle:
        session.verify = args.ca_bundle
    if args.no_proxy:
        session.trust_env = False
    return session


def _rank_for(proteins: list[dict], features: list[dict],
              specialty_rows: list[dict]) -> dict[str, dict]:
    """Rank only the features that became proteins, so the ranks have no holes.

    `parse.proteins_section` keeps CDS rows only, but a parsed CGA directory also holds
    rRNA and tRNA features. Ranking all of them let a non-CDS row consume rank 1 — an
    rRNA annotated "haemolysin toxin secretion system" would outrank every real protein
    and then be dropped, leaving `proteins[]` with no rank 1 and a gap in the report's
    "all N proteins ranked" table.
    """
    known = {p["feature_id"] for p in proteins if p.get("feature_id")}
    from .priority import rank_proteins

    eligible = [f for f in features
                if (f.get("patric_id") or f.get("feature_id")) in known]
    return rank_proteins(eligible, specialty_rows)


def _write_priorities(run_dir: Path, proteins: list[dict],
                      ranked: dict[str, dict]) -> None:
    """Merge `m1_priority` onto `proteins[]`.

    Written through `update_proteins`, which requires the protein to exist already — so
    this runs after the `proteins` section and can never invent a feature_id. Note it
    lands on `m1_priority`, not `triage`: that key belongs to M2 and means something else
    (`priority.py`).
    """
    known = {p["feature_id"] for p in proteins if p.get("feature_id")}
    updates = {fid: {"m1_priority": priority}
               for fid, priority in ranked.items() if fid in known}
    if updates:
        update_proteins(run_dir, updates)


def _report_outputs(run_dir: Path, m1_dir: Path, *, tree: dict | None,
                    want_html: bool, want_figures: bool, seq_limit: int) -> None:
    """Write the reports and figures — all optional, none fatal."""
    from .report_html import collect_priorities, render_run
    from .report_md import render_run as render_markdown
    from ..common.io import read_report

    # The markdown summary is always written: it is small, diffable, greppable and
    # readable over SSH, and it is the run summary the original script printed.
    summary = render_markdown(run_dir)
    (m1_dir / "report.md").write_text(summary, encoding="utf-8")
    print(f"summary -> {m1_dir / 'report.md'}")

    if want_html:
        page = render_run(run_dir, seq_limit=seq_limit, tree=tree)
        target = m1_dir / "report.html"
        target.write_text(page, encoding="utf-8")
        print(f"interactive report -> {target} ({len(page) / 1024:.0f} KB)")

    if want_figures:
        from .figures import save_figures
        report = read_report(run_dir)
        counts = (report.get("genome") or {}).get("specialty_gene_counts") or {}
        made = save_figures(m1_dir / "figures",
                            priorities=collect_priorities(report),
                            specialty_counts=counts, tree=tree)
        if made:
            print("figures -> " + ", ".join(p.name for p in made))
        else:
            print("figures skipped (matplotlib not installed)", file=sys.stderr)


def _run_api_route(args: argparse.Namespace, run_dir: Path, m1_dir: Path,
                   started: float) -> int:
    """The Data API route: resolve a genome, collect, and write the same contract."""
    from .bvbrc_api import BvbrcApi, GenomeResolutionError, read_fasta_defline
    from .collect import (collect_all, fetch_sequences, genome_section as api_genome_section,
                          proteins_section as api_proteins_section,
                          run_fields as api_run_fields, taxonomy_tree_text,
                          write_m1_dir as api_write_m1_dir, write_section_csvs)
    from .pathogens import resolve_disease_profile

    if not (args.contigs or args.species or args.genome_id):
        print("--from-bvbrc-api needs --contigs, --species or --genome-id", file=sys.stderr)
        return 2

    cache_path = Path(args.cache) if args.cache else m1_dir / "cache.sqlite"
    client = CachedJsonClient(cache=JsonCache(cache_path), offline=args.offline,
                              session=_api_session(args),
                              timeout_seconds=args.timeout,
                              max_attempts=max(1, args.retries))
    api = BvbrcApi(client)

    try:
        print("resolving the genome in BV-BRC ...")
        genome, resolution = api.resolve_genome(
            genome_id=args.genome_id, species=args.species, fasta=args.contigs
        )
    except GenomeResolutionError as exc:
        print(f"genome resolution failed: {exc}", file=sys.stderr)
        return 3
    except OfflineCacheMiss as exc:
        print(f"offline: {exc}", file=sys.stderr)
        return 5

    genome_id = genome["genome_id"]
    print(f"  {genome.get('genome_name')} (genome_id {genome_id})")
    print(f"  {resolution}")

    print("collecting features / specialty genes / AMR / taxonomy / metadata ...")
    try:
        bundle = collect_all(api, genome, cap=args.max_features,
                             annotation=args.annotation)
        print(f"  {len(bundle['features'])} CDS, {len(bundle['specialty'])} specialty rows")
    except OfflineCacheMiss as exc:
        print(f"offline: {exc}", file=sys.stderr)
        return 5
    except HttpError as exc:
        print(f"BV-BRC request failed while collecting: {exc}", file=sys.stderr)
        return 6

    if args.contigs:
        _, n_seq, n_bases = read_fasta_defline(args.contigs)
        print(f"  input assembly: {n_seq} sequences, {n_bases} bp")

    # Rank before fetching sequences, so `--seq-cap` truncates by priority rather than by
    # position on the chromosome. Otherwise, on a genome with more CDS than the cap, the
    # top-ranked candidates can all fall outside the fetched window and never reach
    # `m1/proteins.faa`, which is M2's input.
    proteins = api_proteins_section(bundle)
    ranked = _rank_for(proteins, bundle["features"], bundle["specialty"])

    print("fetching protein sequences ...")
    try:
        sequences = fetch_sequences(api, bundle["features"], cap=args.seq_cap,
                                    priority=ranked)
        print(f"  {len(sequences)} sequences")
    except OfflineCacheMiss as exc:
        print(f"offline: {exc}", file=sys.stderr)
        return 5

    # --- the contract is written here, before anything optional -------------
    # Everything below this point is a figure or a convenience. An HTTP error while
    # building a decorative phylogeny used to leave `report.json` with no `genome` and no
    # `proteins` at all, which is the opposite of what tree.py promises.
    counts = api_write_m1_dir(bundle, sequences, m1_dir)
    update_section(run_dir, "genome", api_genome_section(genome, bundle, resolution))
    update_section(run_dir, "proteins", proteins)
    _write_priorities(run_dir, proteins, ranked)

    report = init_report(run_dir, run_dir.name)
    run_section = dict(report.get("run") or {})
    run_section.update(api_run_fields(genome, bundle, resolution))
    modules = dict(run_section.get("modules") or {})
    modules["m1"] = {"status": "ok", "route": "api", "finished_at": _now(),
                     "elapsed_seconds": round(time.time() - started, 1), **counts}
    run_section["modules"] = modules
    update_section(run_dir, "run", run_section)

    profile = resolve_disease_profile(genome)
    (m1_dir / "disease_profile.json").write_text(json.dumps(profile, indent=1, default=str))
    tables = write_section_csvs(bundle, m1_dir, priority=ranked)
    (m1_dir / "taxonomy_tree.txt").write_text(
        taxonomy_tree_text(bundle["taxonomy"]["lineage"], genome,
                           bundle["taxonomy"]["neighbors"]) + "\n", encoding="utf-8")
    print(f"tables -> {', '.join(tables)}, taxonomy_tree.txt")

    tree: dict | None = None
    if not args.no_tree:
        from .tree import build_gene_content_tree
        print("building the gene-content phylogeny ...")
        try:
            tree = build_gene_content_tree(
                api, genome, bundle["taxonomy"]["neighbors"], bundle["close_pathogens"],
                max_leaves=args.tree_leaves,
            )
        except OfflineCacheMiss as exc:
            tree = {"ok": False, "reason": f"offline cache miss ({exc})"}
        except HttpError as exc:
            tree = {"ok": False, "reason": f"BV-BRC request failed ({exc})"}
        if tree.get("ok"):
            from .report_html import render_tree_svg
            (m1_dir / "gene_content_tree.nwk").write_text(tree["newick"] + "\n")
            (m1_dir / "gene_content_tree.svg").write_text(render_tree_svg(tree),
                                                          encoding="utf-8")
            flagged = sum(1 for leaf in tree["leaves"] if leaf.get("is_pathogen"))
            print(f"  {tree['n_genomes']} genomes, {flagged} human-associated")
        else:
            print(f"  skipped: {tree.get('reason')}")

    if api.failures:
        print(f"\n{len(api.failures)} optional request(s) failed and were skipped:",
              file=sys.stderr)
        for failure in api.failures[:5]:
            print(f"  {failure.get('facet')}: {failure.get('error')}", file=sys.stderr)

    _report_outputs(run_dir, m1_dir, tree=tree, want_html=args.html,
                    want_figures=args.figures, seq_limit=args.seq_cap)

    print(f"\nwrote {len(proteins)} proteins, {counts['specialty_rows']} specialty rows "
          f"-> {m1_dir}")
    print("NOTE: the API route describes a reference genome for this organism, not the "
          "submitted assembly. Do not treat gene presence or absence here as a property "
          "of your isolate.")
    return 0


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    run_dir = Path(args.run)
    run_id = run_dir.name
    m1_dir = run_dir / "m1"
    m1_dir.mkdir(parents=True, exist_ok=True)
    started = time.time()

    init_report(run_dir, run_id, input_file=args.contigs, modules={})

    if args.from_bvbrc_api:
        if args.from_cga_dir:
            print("--from-bvbrc-api and --from-cga-dir are different routes; pick one",
                  file=sys.stderr)
            return 2
        return _run_api_route(args, run_dir, m1_dir, started)

    taxon_call: dict | None = None
    cga_dir = args.from_cga_dir

    if cga_dir is None:
        if not (args.contigs and args.ws_dir):
            print("--contigs and --ws-dir are required unless --from-cga-dir is given "
                  "(or use --from-bvbrc-api)", file=sys.stderr)
            return 2

        contigs = Path(args.contigs)
        if not args.no_blind:
            # Named for the run, and this is load-bearing. `cga.upload` shells out to
            # `p3-cp`, which uploads the file under its *local* basename and exits 0
            # **without overwriting** an existing workspace file. A constant basename
            # therefore means every run writes to the same workspace path, the second
            # run's upload silently no-ops, and Minhash calls the taxon of whichever
            # genome got there first — so CGA annotates with a taxon from a different
            # organism. Nothing downstream can detect that.
            #
            # Note `upload(name=...)` does NOT solve this: it only changes the path
            # string `upload` returns, not the name `p3-cp` writes, so passing it would
            # hand Minhash a path that does not exist.
            #
            # `test_blinded_contigs_are_named_for_the_run` pins this; it regressed once.
            contigs = cga.blind_contigs(contigs, m1_dir / f"{run_id}.contigs.fna")
            print(f"blinded contigs -> {contigs}")

        if args.taxon_id:
            if not args.genetic_code:
                print("--genetic-code is required with --taxon-id", file=sys.stderr)
                return 2
            taxon_call = {"called_by": "supplied", "taxon_id": args.taxon_id,
                          "genetic_code": args.genetic_code,
                          "taxon_name": args.scientific_name or str(args.taxon_id),
                          "called_rank": None, "hits": []}
        else:
            ws_fasta = cga.upload(contigs, args.ws_dir)
            print(f"uploaded -> {ws_fasta}\ncalling taxon via Minhash ...")
            try:
                taxon_call = call_taxon(ws_fasta, max_distance=args.max_distance,
                                        max_hits=args.max_hits,
                                        all_public=not args.reference_scope)
            except TaxonCallError as exc:
                print(f"taxon call failed: {exc}", file=sys.stderr)
                return 3
            print(f"  {taxon_call['n_hits']} hits; top {taxon_call['top_hit']['genome_id']} "
                  f"at d={taxon_call['top_hit']['distance']:.5f}")
            print(f"  called {taxon_call['called_rank']}: {taxon_call['taxon_name']} "
                  f"(taxon {taxon_call['taxon_id']}), genetic code {taxon_call['genetic_code']}")
        (m1_dir / "taxon_call.json").write_text(json.dumps(taxon_call, indent=1))

        out_name = args.out_name or run_id
        job_id = args.job_id
        if job_id is None:
            job_id = cga.submit(
                contigs, ws_dir=args.ws_dir, out_name=out_name,
                scientific_name=args.scientific_name or taxon_call["taxon_name"],
                taxon_id=taxon_call["taxon_id"], genetic_code=taxon_call["genetic_code"],
                label=run_id, dry_run=args.dry_run,
            )
            if args.dry_run:
                return 0
            print(f"submitted CGA job {job_id}")
        print("polling ...")
        print(f"  final status: {cga.wait(job_id, interval=args.poll_interval, timeout=args.poll_timeout)}")
        cga_dir = str(run_dir / "cga")
        cga.fetch(args.ws_dir, out_name, cga_dir)
        print(f"retrieved -> {cga_dir}")

    parsed = load_cga(cga_dir)
    proteins = proteins_section(parsed)
    counts = write_m1_dir(parsed, m1_dir)
    update_section(run_dir, "genome", genome_section(parsed, taxon_call))
    update_section(run_dir, "proteins", proteins)
    ranked = _rank_for(proteins, parsed.features, parsed.specialty)
    _write_priorities(run_dir, proteins, ranked)
    from .collect import write_ranked_csv
    write_ranked_csv(parsed.features, ranked, m1_dir)

    report = init_report(run_dir, run_id)
    run_section = dict(report.get("run") or {})
    run_section.update(run_fields(parsed))
    if taxon_call:
        run_section["taxon_call"] = {k: taxon_call.get(k) for k in
                                     ("called_by", "called_rank", "taxon_id", "taxon_name",
                                      "genetic_code", "n_hits")}
    modules = dict(run_section.get("modules") or {})
    modules["m1"] = {"status": "ok", "route": "cga", "finished_at": _now(),
                     "elapsed_seconds": round(time.time() - started, 1), **counts}
    run_section["modules"] = modules
    update_section(run_dir, "run", run_section)

    # Species-level metadata. CGA analyses an assembly BV-BRC has never seen — `2097.118`
    # is the id CGA minted for our own submission and has no public record — so growth,
    # isolation, host and disease cannot be looked up against our own genome id. They come
    # from a public record instead, and `metadata_provenance` records which one and how
    # close it is. Optional by design: the contract above is already on disk, and a network
    # failure here must not cost it. The catch is deliberately broad for that reason; the
    # exception type is printed so a real bug is still visible.
    try:
        from .bvbrc_api import BvbrcApi
        from .collect import collect_cga_metadata
        from .pathogens import resolve_disease_profile
        from ..common.http import CachedJsonClient, JsonCache
        cache_path = Path(args.cache) if getattr(args, "cache", "") else m1_dir / "cache.sqlite"
        api = BvbrcApi(CachedJsonClient(
            cache=JsonCache(cache_path), offline=args.offline, session=_api_session(args),
            timeout_seconds=args.timeout, max_attempts=max(1, args.retries)))
        section = dict((init_report(run_dir, run_id).get("genome")) or {})
        metadata = collect_cga_metadata(api, section)
        section.update(metadata)
        species = metadata["metadata_provenance"].get("species") or ""
        if species:
            section.setdefault("species", species)
            section.setdefault("genus", species.split(" ")[0])
        update_section(run_dir, "genome", section)
        profile = resolve_disease_profile({"species": species,
                                           "genus": section.get("genus") or "",
                                           "disease": metadata["disease"]})
        (m1_dir / "disease_profile.json").write_text(
            json.dumps(profile, indent=1, default=str), encoding="utf-8")
        p = metadata["metadata_provenance"]
        where = f" from {p['genome_id']} ({p['genome_name']})" if p["genome_id"] else ""
        print(f"species metadata: {p['basis']}{where}")
    except Exception as exc:
        print(f"species metadata skipped — {type(exc).__name__}: {exc}", file=sys.stderr)

    _report_outputs(run_dir, m1_dir, tree=None, want_html=args.html,
                    want_figures=args.figures, seq_limit=1200)

    print(f"\nwrote {len(proteins)} proteins, {counts['specialty_rows']} specialty rows -> {m1_dir}")
    print(f"genome {parsed.genome_id}, quality {parsed.quality or 'unknown'}")

    if parsed.is_poor:
        print(f"\nQUALITY GATE: {parsed.quality_reason()}", file=sys.stderr)
        print("A Poor annotation usually means the taxon or genetic code was wrong; "
              "the CDS count inflates and proteins fragment (docs/01a-cga-coverage.md).",
              file=sys.stderr)
        if not args.allow_poor:
            print("Refusing to hand this downstream. Re-run with --allow-poor to override.",
                  file=sys.stderr)
            return 4
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

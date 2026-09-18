# M1 — Genome analysis (BV-BRC CGA)

**Scope:** take an unknown bacterial assembly to genes, AMR calls, phylogeny and closest relatives, and write the `run`, `genome` and `proteins` sections.

> **Two routes.** This document is the CGA-service route, which annotates the submitted
> assembly and is authoritative. [01b-m1-api-mode.md](01b-m1-api-mode.md) is the BV-BRC
> Data API route (`--from-bvbrc-api`): no account, seconds instead of hours, but it
> describes a *reference genome for the organism* rather than your sample. Both write the
> same `report.json` sections and the same `<run>/m1/` files, so M2 cannot tell them apart.
> Reporting is shared too: `--html` writes the interactive report from either route, and
> both fill `proteins[].m1_priority`.

## Inputs

- Contig FASTA (blinded: headers stripped to `contig_1`, …).
- BV-BRC account; token from `p3-login`.

## Outputs

- `genome`: `taxonomy` (predicted taxon + how it was called), `closest_genomes[]` (genome_id, name, mash_distance, ani, snp_distance), `tree_newick`, `cga_job_id`.
- `proteins[]`: `feature_id`, `locus_tag`, `contig`, `start`, `end`, `strand`, `product`, `subsystems[]`, `specialty[]` (type: virulence|amr|drug_target|transporter, database, hit, identity, coverage).
- `proteins.faa` — all protein sequences.
- `proteins[].m1_priority`: pathogenesis priority score, rank, mechanism categories, a mechanism hypothesis and the full score breakdown. Distinct from M2's `triage` and expected to disagree with it — see [00a-data-contract.md](00a-data-contract.md).
- `run`: versions, job IDs, timestamps, seeds.
- `<run>/m1/report.md`: a plain-markdown summary — counts, the top-ranked candidates, the caveats. Built by `s2f/m1_genome/report_md.py`.
- `<run>/m1/proteins_ranked.csv`: the full ranking as a table.
- `<run>/m1/report.html` (with `--html`): a self-contained interactive report — searchable, sortable protein table with mechanism hypotheses, gene-content phylogeny with human-associated relatives flagged, mechanism → host effect → disease flow diagram, cited references. Built by `s2f/m1_genome/report_html.py`, which reads `report.json` so it works for either route and for runs made before it existed.
- `<run>/m1/figures/` (with `--figures`): publication PNGs, if matplotlib is installed.

## Tools and order

1. **Similar Genome Finder** (Mash) → predicted taxon and closest genomes. Run first: annotation needs a taxon at genus level or below.
2. **Comprehensive Genome Analysis** (`p3-submit-CGA`) → CDS, products, subsystems, specialty genes, tree, `FullGenomeReport.html`, JSON genome object. Record submit and finish times.
3. **skani** for ANI against the closest ~10 genomes — see [ANI to the closest genomes](#ani-to-the-closest-genomes) below. SNP distance is deferred to its own issue. BV-BRC Variation Analysis needs reads, so it does not apply.
4. Optional: BV-BRC Mobile Element Detection for plasmids/prophages; CheckM2 + QUAST for assembly sanity.

## Steps

1. `fetch_similar_genomes()` → cache raw response, parse to `closest_genomes`.
2. `submit_cga()` / `poll_cga()` → store job ID; poll politely; on failure fall through to step 5.
3. `parse_cga()` → `proteins[]`, `proteins.faa`, `genome.tree_newick`, subsystems, specialty genes with identity and coverage.
4. `distances.add_ani()` → ANI per closest genome (`--ani`). SNP distance is not implemented.
5. **Fallback path (build this early, it is also the fixture source):** pull the existing BV-BRC annotation for a known genome with the `p3-` CLI (`p3-get-genome-features`, specialty gene tables) and emit the same `proteins[]` shape. This unblocks M2–M6 before any CGA job finishes.

## ANI to the closest genomes

Issue #7. Similar Genome Finder returns a Mash distance and a shared-k-mer count for each
close genome. Neither is an identity — Mash distance approximates `1 - ANI` over a
1000-hash sketch, which is enough to call a taxon and not enough to print in a report. The
`--ani` step computes the identity itself:

```bash
python -m s2f.m1_genome --run runs/<id> --from-cga-dir data/<cga> \
    --contigs <the assembly that was submitted>.fna --ani
```

It downloads each closest genome's contigs from the BV-BRC **Data API** — plain HTTPS, no
`p3-` CLI and no token — caches them under `<run>/cache/` with a readable copy in
`<run>/m1/genomes/<genome_id>.fna`, and runs `skani dist` against the assembly. A second
run makes no network calls; `--offline` replays the cache and fails on a miss.

```
https://www.bv-brc.org/api/genome_sequence/?eq(genome_id,<id>)&limit(10000)&http_accept=application/dna+fasta
```

**Not `ftp.bvbrc.org`.** That host resolves (140.221.78.70) but refuses both port 80 and
port 443 — it serves FTP only. Checked from lambda0 on 2026-09-18: `curl` returns
`Connection refused` in under 3 ms, while `www.bv-brc.org/api` answers 200. `http_accept`
is passed as a query parameter rather than an `Accept` header so the URL recorded in
`ani.json` is pasteable and reproduces the exact bytes. For `243273.25` it returns 589,848
bytes in one record, byte-identical to `fixtures/genomes/mgen_G37/mgen_G37.fna` — which is
where that fixture came from.

Two Data API behaviours the code has to handle, and an FTP server would not have:

- **An unknown genome id answers `200` with an empty body**, not `404`. Whether a download
  worked is therefore decided from the body, never from the status code.
- **A truncated result looks exactly like a complete one.** The API pages at 25 rows and
  `limit()` did not visibly raise that in testing, so after each download the record count
  is checked against the `Content-Range` total from a one-row companion query. A short
  FASTA is refused with a recorded reason rather than used — half a reference genome
  yields a plausible, wrong ANI and fails nothing. A `Content-Range` of `items 0-0/*`
  means the total is unknown, which is not a mismatch and does not reject the download.

**It needs the assembly FASTA.** A retrieved CGA directory contains an annotated genome and
a tree but not the contigs that were submitted, so `--from-cga-dir` runs must pass
`--contigs` or `--ani-query`. There is no fallback to "whatever FASTA is lying around".

**skani, not fastANI.** skani estimates ANI by sparse chaining rather than fragment mapping:
two orders of magnitude faster, and — the part that matters — accurate on fragmented drafts,
where fastANI degrades. One static binary, no database:

```bash
curl -sSL -o ~/bin/skani https://github.com/bluenote-1577/skani/releases/download/latest/skani
chmod +x ~/bin/skani          # or: conda install -c bioconda skani / cargo install skani
```

### What lands in `report.json`

Each `genome.closest_genomes[]` row that was compared gains `ani`,
`ani_align_fraction_query`, `ani_align_fraction_ref`, `ani_source`, `ani_reference_url`,
`ani_retrieved_at` and `self_match`. `<run>/m1/ani.json` holds the skani version, the exact
argv, the parameters, per-genome download provenance and the reasons any genome has no
number. `run.tool_versions.skani` records the binary. The schema is unchanged — `genome` is
open, and `closest_genomes` is an untyped array.

**A row with no ANI keeps `ani: null` and a recorded reason.** skani drops a pair whose
aligned fraction is below `--min-af` (15% by default) instead of reporting a low identity,
because below that the estimate means nothing. The Mash distance is never converted into a
stand-in number. A genome BV-BRC would not serve is recorded the same way, with the HTTP
error.

**The self-match is kept and flagged.** A genome submitted blind is usually already in
BV-BRC, so the nearest hit is the assembly matching itself at ~100% ANI. It cannot be
detected by identifier — CGA mints a fresh genome id for the submission (`2097.118`), which
matches no public record — so it is called from the number: ANI ≥ `--ani-self-threshold`
(99.95 by default) over ≥ 90% of the query. `ani.json` counts the non-self rows separately
and warns when the *only* scored row is the assembly itself.

### Calibrated against known divergence

Run on the *M. genitalium* G37 fixture (580 kb) against copies of itself mutated at a
chosen substitution rate, skani 0.3.2, defaults:

| Substitutions introduced | ANI reported | Aligned fraction |
| --- | --- | --- |
| 0% (same file) | 100.00 | 99.99 |
| 0.05% | 99.95 | 99.99 |
| 0.5% | 99.45 | 99.92 |
| 2% | 97.92 | 99.49 |
| 12% | 87.64 | 76.62 |

Identity tracks the mutation rate closely down to the species boundary; past it the aligned
fraction falls away first, which is the signal that the pair is leaving the range where a
single ANI number is worth quoting. The whole six-genome comparison takes well under a
second; the download dominates the step.

### SNP distance is not implemented

Issue #7 also asked for a SNP distance per closest genome. It was deferred: it needs a
second aligner (Snippy `--ctgs`, Parsnp or MUMmer `dnadiff`), and it needs a defensible
answer to "how many SNPs" for pairs below ~95% ANI, where the question stops meaning
anything. `closest_genomes[]` rows keep `snp_distance: null`, and `ani.json` records the
deferral rather than leaving a silent hole. See the follow-up issue in
[issues.md](issues.md).

## Acceptance checks

- `--dry-run` produces a valid `genome` + `proteins` from fixtures, no network.
- Against the blinded known genome: correct species; CDS count within a few percent of RefSeq for the same assembly; known AMR and virulence genes present.
- Every specialty hit carries database, identity and coverage. No bare gene names.
- Re-running is idempotent: same input, same `report.json` apart from timestamps.

## Pitfalls

- Annotation quality depends on the taxon you pass; a genus-level or better taxon is required.
- CGA runtime is unknown to us — submit on day 1 and record the duration.
- The blinded genome will likely match itself; keep it but also report the closest non-identical genomes.
- Do not put sequences in `report.json`.

## Kickoff prompt

> Read the project docs `00-architecture.md` and `01-m1-genome.md`. Implement `s2f/m1_genome` as described: Similar Genome Finder lookup, CGA submit/poll/parse, ANI and SNP distance, plus the `p3-`CLI fallback path that emits the same `proteins[]` shape. Provide `--dry-run` against fixtures and the acceptance checks as tests. Use `common/http.py`, `common/ids.py` and `common.io.update_section`; do not change the schema — raise it if you think it needs changing.

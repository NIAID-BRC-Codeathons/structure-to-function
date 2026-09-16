# M1 — Genome analysis (BV-BRC CGA)

**Scope:** take an unknown bacterial assembly to genes, AMR calls, phylogeny and closest relatives, and write the `run`, `genome` and `proteins` sections.

## Inputs

- Contig FASTA (blinded: headers stripped to `contig_1`, …).
- BV-BRC account; token from `p3-login`.

## Outputs

- `genome`: `taxonomy` (predicted taxon + how it was called), `closest_genomes[]` (genome_id, name, mash_distance, ani, snp_distance), `tree_newick`, `cga_job_id`.
- `proteins[]`: `feature_id`, `locus_tag`, `contig`, `start`, `end`, `strand`, `product`, `subsystems[]`, `specialty[]` (type: virulence|amr|drug_target|transporter, database, hit, identity, coverage).
- `proteins.faa` — all protein sequences.
- `run`: versions, job IDs, timestamps, seeds.

## Tools and order

1. **Similar Genome Finder** (Mash) → predicted taxon and closest genomes. Run first: annotation needs a taxon at genus level or below.
2. **Comprehensive Genome Analysis** (`p3-submit-CGA`) → CDS, products, subsystems, specialty genes, tree, `FullGenomeReport.html`, JSON genome object. Record submit and finish times.
3. **skani** (or fastANI) for ANI against the closest ~10 genomes; **Snippy `--ctgs`** or **Parsnp** for SNP distance. BV-BRC Variation Analysis needs reads, so it does not apply.
4. Optional: BV-BRC Mobile Element Detection for plasmids/prophages; CheckM2 + QUAST for assembly sanity.

## Steps

1. `fetch_similar_genomes()` → cache raw response, parse to `closest_genomes`.
2. `submit_cga()` / `poll_cga()` → store job ID; poll politely; on failure fall through to step 5.
3. `parse_cga()` → `proteins[]`, `proteins.faa`, `genome.tree_newick`, subsystems, specialty genes with identity and coverage.
4. `compute_distances()` → ANI and SNP distance per closest genome.
5. **Fallback path (build this early, it is also the fixture source):** pull the existing BV-BRC annotation for a known genome with the `p3-` CLI (`p3-get-genome-features`, specialty gene tables) and emit the same `proteins[]` shape. This unblocks M2–M6 before any CGA job finishes.

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

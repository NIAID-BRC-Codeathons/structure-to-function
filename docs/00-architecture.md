# Project 9 — Architecture and Shared Conventions

Read this first in any coding thread. It is the contract between modules; everything else is module-local.

## Pipeline

1. **M1 genome** — BV-BRC Comprehensive Genome Analysis: genes/proteins, AMR, phylogeny, related genomes.
2. **M2 triage** — annotate proteins from lookups (PDB, AlphaFold DB, ESM Atlas, UniProt) and a knowledge-graph mapping to known disease/pathogen proteins; rank and select ~50.
3. **M3 fold** — reuse or predict structures for the selected proteins; define pockets.
4. **M4 ligands** — assemble and score candidate ligands per protein; dock the top ones.
5. **M5 disease** — virulence, growth/nutrient needs, non-drug mitigation strategies, RAG/agentic synthesis.
6. **M6 report** — one HTML report of ranked ligand candidates and supporting evidence.

Data flows one way. M2 may read M1; M4 may read M1–M3; nothing writes upstream.

## Repo layout (`NIAID-BRC-Codeathons/structure-to-function`)

```
s2f/
  common/        ids.py  schema.py  io.py  http.py  lit.py  argo.py  provenance.py
  m1_genome/     m2_triage/  m3_fold/  m4_ligands/  m5_disease/  m6_report/
fixtures/        report.fixture.json + small per-module fixtures
runs/<run_id>/   report.json, proteins.faa, structures/, docking/, cache/
tests/
```

- Python 3.11. One `requirements.txt` per module plus a root `requirements-common.txt`; no shared conda environment surgery during the codeathon.
- Every module is a CLI: `python -m s2f.m2_triage --run runs/<run_id> [--limit N] [--dry-run]`.
- `--dry-run` must work against `fixtures/` with no network. This is what lets six people build in parallel on day 1.

## The data contract

`runs/<run_id>/report.json` is the single state file. Each module owns exactly one top-level key and never edits another's:

| Key | Owner | Contents |
| --- | --- | --- |
| `run` | M1 | run_id, created_at, input file, tool versions, seeds |
| `genome` | M1 | taxonomy, closest_genomes[], tree_newick, cga_job_id |
| `proteins` | M1 writes the list; M2 writes `annotations`, `xrefs`, `flags`, `triage` | one record per CDS |
| `structures` | M3 | per selected protein: source, accession/job, confidence, pockets[] |
| `kg` | M2 | nodes[], edges[] of the assembled subgraph (or per-protein `kg_links`) |
| `ligands` | M4 | candidate ligands with score components |
| `docking` | M4 | poses, scores, controls |
| `analyses` | M3/M4/M5 | generic list: feature_id, method, parameters, result, confidence |
| `disease` | M5 | claims[] with evidence, non-drug strategies |
| `report` | M6 | rendered path, rendered_at, counts |

Rules:

- Protein amino-acid sequences live in `proteins.faa`, not in JSON. JSON holds IDs and a `faa_offset` if useful.
- Writes go through `common.io.update_section(run_dir, key, obj)`: read, replace one key, write temp file, atomic rename. Never hand-edit a live `report.json`.
- `schema.py` holds JSON Schema for each section plus `validate(section, obj)`. CI runs it; a module that fails validation is broken, not the schema. Implemented in `s2f/common/schema.py` and `s2f/common/io.py` — see [00a-data-contract.md](00a-data-contract.md) (issue #2).
- Schema changes are announced in the group channel and bumped in `schema_version`; fixtures update in the same commit.

## Identity

- Canonical protein ID = BV-BRC `feature_id`. Every record keys off it.
- `common/ids.py` is the only place that maps IDs: feature_id ↔ locus_tag ↔ UniProt ↔ PDB ↔ ChEMBL target ↔ gene name, using the UniProt ID-mapping API with a local cache. Everyone uses it; nobody writes their own mapping. Most cross-module bugs will come from this.

## Provenance and caching

- All outbound HTTP goes through `common/http.py`: on-disk cache under `runs/<run_id>/cache/`, retries with backoff, polite rate limiting, one user agent naming the project.
- Every derived record carries `source`, `source_version` (or release date), `retrieved_at`, `method`, `params`.
- Any LLM-written claim carries `evidence_feature_ids` and/or `evidence_pmcids`. Claims without evidence are dropped before M6.
- Docking and folding record seeds, model versions and parameters; runs must be repeatable.

## Credentials

Environment only, never in the repo: `ARGO_USER`, `RAGSTACK_API_KEY` (or `~/.patric_token` from `p3-login`). `.gitignore` covers `runs/`, `data/`, `.env`, tokens.

## Working in separate threads

Each module doc (`01-`…`06-`) carries: scope, inputs, outputs, tools/APIs, steps, acceptance checks, pitfalls, and a kickoff prompt to paste into a fresh conversation. A thread should read `00-architecture.md` plus its own module doc and nothing else.

Cross-module questions (schema, IDs, provenance) come back to this doc; if a thread needs a change here, it stops and raises it rather than editing its own copy of the contract.

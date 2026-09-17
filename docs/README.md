# Project 9 — Documentation

Structure-to-function pipeline: unknown bacterial genome → ranked therapeutic candidates with evidence.

NIAID-BRCs AI Codeathon 2.0 · Argonne National Laboratory · Sept 16–18, 2026.

## Start here

| Doc | Read it when |
| --- | --- |
| [00-architecture.md](00-architecture.md) | Always, before writing code. Data contract, repo layout, ID and provenance rules. |
| [00a-data-contract.md](00a-data-contract.md) | Writing or reading `report.json` — the schema, the section writer, the fixture. |
| [pipeline.md](pipeline.md) | Choosing a tool or database for a step. |
| [pitfalls.md](pitfalls.md) | Before trusting any result, and before designing a scoring or docking step. |
| [common-ids.md](common-ids.md) | Mapping any identifier — feature_id, UniProt, PDB, ChEMBL. Never write your own mapping. |
| [07-decisions-and-risks.md](07-decisions-and-risks.md) | Open decisions, thread dependencies, non-negotiables. |

## Modules

Each module doc carries scope, inputs, outputs, tools, steps, acceptance checks, pitfalls, and a kickoff prompt for a fresh coding thread.

| Module | Doc | Owns in `report.json` |
| --- | --- | --- |
| M1 genome | [01-m1-genome.md](01-m1-genome.md) | `run`, `genome`, `proteins` (list) |
| M2 triage | [02-m2-triage.md](02-m2-triage.md) | `proteins.*` enrichment, `kg` |
| M2a PDB evidence | [02a-m2-pdb-evidence.md](02a-m2-pdb-evidence.md) | query parameters, scoring weights, weight change log |
| M2b knowledge graph | [02b-m2-knowledge-graph.md](02b-m2-knowledge-graph.md) | source survey, edge types, caps, the homolog bridge |
| M2c functional annotation | [02c-m2-functional-annotation.md](02c-m2-functional-annotation.md) | source precedence, flag definitions, external-tool setup |
| M3 fold | [03-m3-fold.md](03-m3-fold.md) | `structures` |
| M4 ligands | [04-m4-ligands.md](04-m4-ligands.md) | `ligands`, `docking` |
| M5 disease | [05-m5-disease.md](05-m5-disease.md) | `disease` |
| M6 report | [06-m6-report.md](06-m6-report.md) | `report` |

## Ground rules

- One module owns one top-level key in `report.json`; no module edits another's.
- Every module has fixtures and a `--dry-run` that needs no network.
- Every external record carries source, version or access date, and retrieval time.
- No claim without evidence IDs. No clinical or treatment recommendations.
- Docking scores are rankings within one protein, never across proteins.

## Pipeline at a glance

```
contigs
  → M1  BV-BRC CGA: genes, AMR, phylogeny, relatives
  → M2  lookups (PDB / AFDB / ESM Atlas / UniProt) + knowledge graph → rank → ~50 proteins
  → M3  structures + pockets
  → M4  ligand candidates → score → dock (with controls)
  → M5  virulence, growth, nutrients, non-drug strategies (RAG + agentic)
  → M6  single self-contained HTML report
```

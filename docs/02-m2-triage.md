# M2 — Protein lookup, knowledge graph, and triage

**Scope:** enrich every protein from M1 with sequence, structure and knowledge-graph evidence, then rank and select ~50 for structural and ligand work.

## Inputs

`proteins[]` and `proteins.faa` from M1; `genome.taxonomy` for context.

## Outputs

- `proteins[].xrefs`: uniprot, pdb[], afdb, esm_atlas[], chembl_target[], gene_name.
- `proteins[].annotations[]`: each with source, hit, identity/coverage or TM-score, description.
- `proteins[].flags`: virulence, amr, secreted, membrane, human_homolog (with identity), essential_ortholog, structure_available.
- `proteins[].triage`: `score`, `components{}`, `rank`, `selected` (bool), `reason`.
- `kg`: `nodes[]`, `edges[]` (edge: source_id, target_id, type, evidence, provenance).

## Tools and APIs

| Need | Pick |
| --- | --- |
| ID mapping | **UniProt ID-mapping + REST** (via `common/ids.py`) |
| Sequence homologs | **RCSB Search API** (sequence search), DIAMOND vs Swiss-Prot locally, UniProt BLAST as backup |
| Experimental structures | **RCSB Data API**, PDBe REST — note ligands present in the entry |
| Predicted structures | **AlphaFold DB API** by UniProt accession |
| Structure search | **Foldseek** local binary against PDB + AFDB (+ ESM Atlas if available); Foldseek server REST API as fallback |
| ESM Metagenome Atlas | `esmatlas.com` fold/search endpoints — **test on day 1**, treat as bonus |
| Function | **eggNOG-mapper**, InterProScan (if time), PSORTb + SignalP 6, DeepTMHMM for membrane flag |
| Knowledge graph | **Open Targets GraphQL** (target–disease–drug), **STRING API**, DGIdb, HPIDB + PHI-base, IntAct/BioGRID, Reactome + KEGG REST, NCATS Translator (TRAPI) for federated queries |
| Graph assembly | **NetworkX** in-process; export nodes/edges into `kg` |

## Steps

1. `map_ids()` — every protein to UniProt where possible; cache aggressively.
2. `lookup_structures()` — PDB and AFDB hits; record whether a holo (ligand-bound) structure exists for the closest homolog.
3. `foldseek_search()` — for proteins with no informative sequence hit, especially hypotheticals.
4. `annotate_function()` — eggNOG, localization, membrane prediction.
5. `build_kg()` — pathogen protein → homolog (human or characterized pathogen) → disease/pathway → known ligands. Every edge keeps its source and the homology identity that justified it.
6. `score_and_select()` — weighted score from written-down components, e.g. virulence/AMR hit, essentiality ortholog, secreted or surface, ligandable homolog (known ligand exists), structure availability, **penalty for close human homolog**, penalty for predicted membrane protein. Emit rank, `selected` for the top ~50, and `reason` per protein. Keep the full ranked list, not just the winners.

## Acceptance checks

- Score components are recorded per protein and the selection is reproducible from them.
- Discarded proteins remain in `proteins[]` with `selected: false` and a reason.
- Spot-check: for the blinded genome, known virulence and AMR proteins rank in the top 50; a housekeeping protein with a close human homolog does not.
- `--dry-run` runs the whole scoring path from fixtures with no network.

## Pitfalls

- **This module decides the whole result.** Write criteria and weights into the doc before running them; do not tune them after seeing the answer without recording that you did.
- **Human-homolog paradox:** a close human homolog is useful for borrowing ligand scaffolds and bad for an antibacterial target. Keep both facts; treat selectivity as a flag, not a bonus.
- Entity resolution eats the time. Nothing here builds its own ID mapping.
- Open Targets is human-centric; it links through our homolog, not our protein. Say that in the edge type (`homolog_of`, then `target_of`).
- Do not attempt to ingest or host a full knowledge graph. Query APIs, keep our subgraph small.
- Keep `kg` small enough to ship in JSON; if it grows past a few thousand edges, write it to `runs/<id>/kg.json` and keep a summary in `report.json`.

## Kickoff prompt

> Read `00-architecture.md` and `02-m2-triage.md`. Implement `s2f/m2_triage`: ID mapping, PDB/AFDB lookups, Foldseek search, functional annotation, knowledge-graph assembly from Open Targets / STRING / HPIDB / PHI-base / ChEMBL, and the weighted triage score with the components listed. Keep every discarded protein with its reason. `--dry-run` must work from fixtures. Do not invent new ID mappings — extend `common/ids.py` if needed.

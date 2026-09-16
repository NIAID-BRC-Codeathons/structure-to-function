# M4 — Ligand triage, scoring, and docking

**Scope:** assemble candidate ligands per selected protein, score them transparently, dock the top ones, and report ranked candidates with controls.

## Inputs

`proteins[].triage/xrefs/flags`, `kg` (M2), `structures[]` with pockets (M3).

## Outputs

- `ligands[]`: ligand_id (InChIKey canonical), name, SMILES, source, approval status, `targets_via` (how it reached us: kg | chembl_homolog | literature | repurposing_set), `score`, `components{}`, `for_feature_id`, rank.
- `docking[]`: feature_id, ligand_id, pocket_id, engine, version, seed, exhaustiveness, `score`, pose path, `control_type` (none | positive | decoy), rescoring results.
- `analyses[]` entries for co-folding or affinity estimates.

## Tools and APIs

| Need | Pick |
| --- | --- |
| Approved / safe drugs | **ChEMBL API** (`max_phase=4`), **openFDA** label + approval endpoints (BioMCP wraps these), Broad Drug Repurposing Hub, ZINC in-stock subsets. DrugBank only interactively — its license restricts redistribution |
| Activity data | **ChEMBL** target bioactivities (pChEMBL), **BindingDB**, **PubChem PUG-REST** BioAssay |
| Knowledge-graph ligands | **Open Targets** known drugs for the mapped target, DGIdb |
| Literature candidates | **RAGStack** `/v1/retrieve` per target and per compound class; Europe PMC as backup |
| Chemistry | **RDKit**: canonical SMILES, InChIKey dedup, descriptors, property filters, PAINS flags |
| Docking | **AutoDock Vina** (breadth) or smina; gnina for CNN rescoring; **Boltz-2 co-folding + affinity** for the top ~10 pairs; DiffDock optional |

## Steps

1. `collect_ligands()` — union of the sources above per protein, each with provenance; dedup by InChIKey.
2. `score_ligands()` — weighted, written-down components: homology identity to the ligand's known target, measured potency, approval/phase status, property filters and PAINS, literature support, KG edge type. No learned model this week.
3. `cap_and_rank()` — hard cap per protein (start at 200 candidates in, top 20 to docking); record the cap in `run`.
4. `dock()` — Vina against the defined pocket only, fixed seed and exhaustiveness, batched. Save poses.
5. `controls()` — per protein, dock one known binder (positive) and ~20 decoys or random ligands (negative). Report the ranking of real candidates relative to those.
6. `refine_top()` — Boltz-2 co-folding with the ligand for the top handful; store its affinity estimate as an `analyses` entry, not as truth.

## Acceptance checks

- Every ligand row is traceable: source, retrieval date, and why it was proposed for that protein.
- Scores are only ever compared within one protein; the report must not rank across proteins by raw docking score.
- Controls exist for every docked protein, and the summary shows where real candidates fall against decoys.
- Re-running with the same seed reproduces the scores.

## Pitfalls

- **Docking scores are not affinities.** They are a ranking device, and a weak one on predicted apo structures. Every number ships with that caveat.
- **Blind docking is out.** Pocket from M3 or the protein is skipped.
- **"Approved" is not "safe for this."** No dose, route, spectrum or clinical claim anywhere in the output.
- **Compute:** 50 proteins × 200 ligands × Vina is a real job; benchmark one protein first, then decide the cap. Boltz-2 is minutes per pair — top candidates only.
- Cross-protein comparison is the most likely way this pipeline produces a confidently wrong answer.
- Watch for promiscuous compounds (detergents, dyes, PAINS) dominating the ranking; filter them before docking.

## Kickoff prompt

> Read `00-architecture.md` and `04-m4-ligands.md`. Implement `s2f/m4_ligands`: ligand collection from ChEMBL / openFDA / Open Targets / RAGStack / repurposing sets with provenance, RDKit dedup and filters, the weighted per-protein score, capped Vina docking into M3 pockets with fixed seeds, positive and decoy controls, and optional Boltz-2 refinement for the top pairs. Enforce within-protein-only comparisons in the output. `--dry-run` from fixtures, and benchmark one protein before any batch run.

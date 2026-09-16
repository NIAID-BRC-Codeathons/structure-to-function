# M3 — Structures and pockets

**Scope:** give every selected protein a usable structure and at least one defined pocket, or mark it undockable.

## Inputs

`proteins[]` where `triage.selected` is true, plus their `xrefs`.

## Outputs

- `structures[]`: feature_id, `source` (pdb | afdb | predicted), accession or job ID, path under `runs/<id>/structures/`, `confidence` (mean pLDDT, PAE summary, or experimental resolution), `holo_template` (PDB ID + ligand if one aligned), `pockets[]` (id, residues, score, method, `aligned_to_known_site` bool), `usable_for_docking` (bool + reason).
- `analyses[]` entries for anything computed here (normal modes, electrostatics).

## Tools

| Need | Pick |
| --- | --- |
| Existing structures | **Reuse PDB / AlphaFold DB from M2 first** — only fold what has nothing |
| Prediction | **BV-BRC Protein Structure Prediction** (Boltz-2, OpenFold 3, Chai-1, ESMFold, ESMFold2, AlphaFold 2 — web/CLI/API); ESMFold or Boltz-2 locally on ANL GPUs for throughput |
| Cleanup | PDBFixer (missing atoms, hydrogens), PDB2PQR (protonation), OpenBabel or Meeko (PDBQT for Vina) |
| Pockets | **P2Rank** (fast, ML-based) and/or fpocket; align the holo template with US-align or TM-align to transfer its ligand site |
| Optional biophysics | ProDy normal modes (flexibility), APBS for electrostatic surface — both cheap, both go in `analyses[]` |

## Steps

1. `collect_existing()` — copy or download PDB/AFDB structures; prefer experimental, prefer holo.
2. `predict_missing()` — one job per protein; record model, version, seed.
3. `quality_gate()` — mean pLDDT and per-region confidence; trim or mask low-confidence tails; if the whole model is poor, set `usable_for_docking: false`.
4. `find_pockets()` — P2Rank; then, where a holo template exists, superpose and mark the pocket that matches the known ligand site.
5. `prep_receptors()` — protonation and PDBQT for M4, saved per protein.

## Acceptance checks

- Every selected protein has a structure record or an explicit `usable_for_docking: false` with a reason.
- Pocket residues are written out, not just a score — M4 needs a box or residue list.
- 50 proteins complete within one working session; log wall-clock per protein.

## Pitfalls

- **Missing cofactors.** Predicted models have no Mg²⁺, Zn²⁺, heme or nucleotide. If the homolog's site needs one, the pocket is incomplete: mark it, and prefer aligning to a holo template.
- **Membrane proteins** fold and dock badly. Use the M2 membrane flag; deprioritize rather than silently including them.
- A pocket in a low-confidence or disordered region is not a pocket. Gate on local confidence, not just the global mean.
- Do not re-fold something AlphaFold DB already has.

## Kickoff prompt

> Read `00-architecture.md` and `03-m3-fold.md`. Implement `s2f/m3_fold`: collect existing structures, predict only the gaps via the BV-BRC structure service (or local ESMFold), apply the confidence gate, run P2Rank pockets, transfer holo-template sites by superposition, and prepare receptors for AutoDock Vina. Write the `structures` section exactly as specified, including `usable_for_docking` reasons. `--dry-run` from fixtures.

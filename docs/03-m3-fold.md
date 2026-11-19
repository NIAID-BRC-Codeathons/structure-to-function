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

## Pockets (issue #14)

The primary method needs no external tool and no network, because the evidence is already in
hand. M2 names the ligands and metals bound in every qualifying PDB entry, and M3 downloads
that entry's coordinates, so the site is read out of the file M3 already has: the residues of
*our* entity's chains within **4.5 Å** of the named component's atoms.

| Rule | Why |
| --- | --- |
| Only the components M2 named define a site | M2 already dropped waters, cryoprotectants and glycans (`classify_components`). M3 re-filters nothing by name. |
| At least 3 lining residues | Fewer is a surface brush, not a pocket. |
| Only our entity's chains line it | An entry can hold dozens of chains; a ligand bound to a neighbour is not our homolog's site. Chains come from `_struct_asym`. |
| A metal site is a site, flagged `cofactor_required` | pitfalls #2. A predicted model will not contain the Mg²⁺ or Zn²⁺, so M4 must be told the site is incomplete rather than dock into a hole. |
| A drug-like ligand outranks a metal | Both mark a real site; only the ligand shows the site takes a drug-sized molecule. |
| Reject a pocket whose own residues are poorly modelled | A pocket in a low-confidence region is not a pocket. Mean pLDDT over the pocket's residues, not the model's global mean. Only applied to AlphaFold models, where `B_iso_or_equiv` holds pLDDT; in an experimental entry the same column is a B-factor and no confidence is claimed. |

**The template M3 collects is the one it reads the site from.** M2 records the best
qualifying hit *and* the best hit that binds something (`flags.holo_template`); these are
often not the same entry. M3 prefers the holo one, because a site to transfer is worth more
than the last points of sequence identity (pitfalls #2), and records the swap and its
`identity_drop` in `template_preference`. On `runs/mgen_G37` this is the difference between
14 and 42 of 50 selected proteins having anything to transfer.

**No site is invented.** A site is only ever read from the entry that binds the component, so
no superposition is involved and `aligned_to_known_site` is true exactly when the collected
structure *is* that entry. A protein whose homolog binds nothing gets no pocket and a stated
reason — uncovered, not "no pocket exists", the distinction issue #90 draws. Those proteins
are what a detector such as P2Rank is for; `find_pockets` is where it plugs in, discovered on
`PATH` like M2's annotation providers rather than required.

## Steps

1. `collect_existing()` — copy or download PDB/AFDB structures; prefer experimental, prefer holo.
2. `predict_missing()` — one job per protein; record model, version, seed.
3. `quality_gate()` — mean pLDDT and per-region confidence; trim or mask low-confidence tails; if the whole model is poor, set `usable_for_docking: false`.
4. `find_pockets()` — read the holo template's own site (above). P2Rank for the rest, and superposition for a site held in an entry other than the collected one, are not built yet.
5. `prep_receptors()` — protonation and PDBQT for M4, saved per protein.

## Acceptance checks

- Every selected protein has a structure record or an explicit `usable_for_docking: false` with a reason.
- Pocket residues are written out, not just a score — M4 needs a box or residue list. Each
  residue carries chain and author numbering, because a bare residue number is ambiguous
  once an entry has more than one chain.
- A protein with no pocket says which case it fell into, in `pocket_search.reason`.
- 50 proteins complete within one working session; log wall-clock per protein.

## Pitfalls

- **Missing cofactors.** Predicted models have no Mg²⁺, Zn²⁺, heme or nucleotide. If the homolog's site needs one, the pocket is incomplete: mark it, and prefer aligning to a holo template.
- **Membrane proteins** fold and dock badly. Use the M2 membrane flag; deprioritize rather than silently including them.
- A pocket in a low-confidence or disordered region is not a pocket. Gate on local confidence, not just the global mean.
- Do not re-fold something AlphaFold DB already has.

## Kickoff prompt

> Read `00-architecture.md` and `03-m3-fold.md`. Implement `s2f/m3_fold`: collect existing structures, predict only the gaps via the BV-BRC structure service (or local ESMFold), apply the confidence gate, run P2Rank pockets, transfer holo-template sites by superposition, and prepare receptors for AutoDock Vina. Write the `structures` section exactly as specified, including `usable_for_docking` reasons. `--dry-run` from fixtures.

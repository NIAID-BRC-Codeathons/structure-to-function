# 15-minute presentation — outline

Issue [#21](https://github.com/NIAID-BRC-Codeathons/structure-to-function/issues/21).
NIAID-BRCs AI Codeathon 2.0, Argonne, 18 September 2026.

Source of truth for the deck. `slides.html` renders this; edit here first, then mirror.

**Rule for every number on a slide:** it comes from a tracked doc under `docs/`, never from
`runs/`, which is git-ignored (`.gitignore:2`) and therefore not reproducible for anyone reading
the repo. Where a figure is quoted, the doc it comes from is named in the speaker notes.

## Budget

| # | Slide | Time | Cumulative |
| --- | --- | --- | --- |
| 1 | Title | 0:30 | 0:30 |
| 2 | The problem | 0:45 | 1:15 |
| 3 | Six modules were scoped. Three exist. | 1:15 | 2:30 |
| 4 | M1 — genome in, proteins and evidence out | 1:00 | 3:30 |
| 5 | M2 — a transparent score | 1:45 | 5:15 |
| 6 | The weight we changed after seeing the data | 1:15 | 6:30 |
| 7 | **Our top hit was the genome matching itself** | 1:45 | 8:15 |
| 8 | **What we did about it** | 1:15 | 9:30 |
| 9 | M3 — reuse structures before predicting them | 1:00 | 10:30 |
| 10 | What we did not build | 0:45 | 11:15 |
| 11 | Demo | 1:30 | 12:45 |
| 12 | Next | 0:30 | 13:15 |

1:45 of slack against the 15:00 slot. If the demo misbehaves, cut slides 9 and 11 — the
argument survives without them.

## Slides

### 1. Title
Structure to function: from a bacterial genome to ranked, structure-backed drug targets.
Project 9. Name the team.

### 2. The problem
A newly sequenced pathogen genome gives thousands of proteins and almost no idea which ones
matter. A third to a quarter of them are annotated "hypothetical protein". The useful question
is not "what is every protein" but "which 50 are worth a structural and chemical look".

### 3. Pipeline, and what is actually built
Six modules were scoped. **Three exist: M1, M2, M3.** M4 (ligands), M5 (disease evidence) and
M6 (report) do not. Say this on slide 3, not in the last thirty seconds — the rest of the talk
is then about depth rather than coverage.

### 4. M1 — genome to proteins
BV-BRC Comprehensive Genome Analysis: genes, proteins, AMR and virulence rows, phylogeny.
Two genomes went through: *K. pneumoniae* HS11286, 5,523 proteins, 2,592 specialty rows; and
*M. genitalium* G37, 542 proteins. M1 writes an offline HTML report of its own.
*Notes: `docs/02a-m2-pdb-evidence.md:44`.*

### 5. M2 — triage scoring
Every protein gets a PDB sequence search, then a transparent weighted score. Nine components:

    0.40·pdb_evidence + 0.20·virulence_amr + 0.15·essential + 0.15·drug_target
  + 0.30·annotation_gap + 0.10·surface_bonus + 0.10·ligandable_homolog
  − 0.15·membrane_penalty − 0.25·(human_identity/100)

The point is that selection is reproducible from the recorded components alone — every
component, flag and reason is written out for all 542 proteins, not only the chosen 50.
*Notes: `docs/02a-m2-pdb-evidence.md:115`.*

### 6. The weight we changed after seeing the data
The first full HS11286 run selected **0 uncharacterized proteins out of 1,338**. A known drug
target collects 0.50 from virulence, essentiality and drug-target evidence; a hypothetical can
earn 0.10. That is backwards for a project about the unknown.

`annotation_gap` went 0.10 → 0.30, with the sensitivity recorded rather than tuned quietly:
0.20 → 6/50, 0.25 → 9/50, **0.30 → 10/50**, 0.35 → 15/50, 0.40 → 26/50.

Say plainly that this was changed *after* seeing a run, and that it is written down for exactly
that reason. *Notes: decision log, `docs/02a-m2-pdb-evidence.md:418`.*

### 7. What the data told us we had wrong — the centrepiece
A teammate noticed the top-ranked G37 protein was a hypothetical matched at **100% identity** to
PDB `1Q8C_1` — which is *Mycoplasmoides genitalium* P47. **The same organism.** The genome was
matching itself.

Measured across the whole run:

| | Top 50 | Whole genome |
| --- | --- | --- |
| Best hit is the same species | **30** | 34 |
| Best hit is the same genus | **47** | 108 |

88% of every same-species hit in the genome landed in our top 50. The ranking is not wrong — a
100% identity hit really is the strongest evidence available — but it is **easy** in a way a
novel genome will not be.

*Notes: `docs/02a-m2-pdb-evidence.md:310`.*

### 8. What we did about it
- **Flags, not penalties.** The score is unchanged; what changes is what the report may claim.
  Each protein says so in words, so a reader is told when an annotation was easy instead of
  being sold it as a discovery.
- **Do not calibrate the weights on this genome.** Written into the docs as a standing rule.

One detail that decided it: BV-BRC writes *Mycoplasma genitalium* where the PDB says
*Mycoplasmoides genitalium*. Comparing genus strings naively marks every self-match as foreign
— reporting the exact opposite of the truth, silently.

### 9. M3 — structures
Existing structures collected and quality-gated before anything is predicted: 153 AlphaFold
models for G37. Foldseek finds fold matches for proteins sequence search missed, at a median
best-hit identity of **16.6%** — below the 25% floor where sequence alignment stops working,
which is the whole reason structure search earns its place.
*Notes: `docs/02c-m2-structure-search.md:73`.*

### 10. What we did not build
No docking, no ligand scoring, no literature synthesis, no final report renderer. Nothing here
is experimentally validated and nothing implies a clinical claim. A short, flat slide.

### 11. Demo
Offline, from committed fixtures, no network and no credentials:

    python -m s2f.m2_triage --dry-run --run runs/demo --report

Show the written columns for a protein — score, components, flags, reason — and make the point
that the reason is a sentence a human can check, including "best structural hit is the same
species — not a cross-organism transfer".

### 12. Next
The three built modules are a working, reproducible front half. M4/M5/M6 are scoped and
documented, and the same-organism result is the thing to carry into whatever genome comes next.

## Rehearsal checklist

- [ ] Timed once, out loud, start to finish
- [ ] Demo run on the presenting machine, offline, network physically off
- [ ] Fallback ready: pre-rendered output or a recording, in case the live run fails
- [ ] Every number on a slide traced back to the doc named in its notes
- [ ] Slide 3 and slide 10 both delivered — coverage stated twice, on purpose

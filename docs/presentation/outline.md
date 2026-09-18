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
| 3 | **Flow chart** — six modules, one shared contract | 1:15 | 2:30 |
| 4 | **Input** — what actually goes in | 1:00 | 3:30 |
| 5 | M1 — APIs, filters, in/out | 1:00 | 4:30 |
| 6 | M2 — APIs, filters, in/out | 1:30 | 6:00 |
| 7 | M2 — the score, and the weight we changed | 1:30 | 7:30 |
| 8 | M3 — APIs, quality gate, in/out | 1:00 | 8:30 |
| 9 | **Output** — one row, fully shown | 1:15 | 9:45 |
| 10 | **Our top hit was the genome matching itself** | 1:30 | 11:15 |
| 11 | What we did about it | 1:00 | 12:15 |
| 12 | **Status** — what we have, what is left | 1:00 | 13:15 |
| 13 | Demo | 0:45 | 14:00 |
| 14 | Next | 0:30 | 14:30 |

0:30 of slack against the 15:00 slot, which is tight. If Q&A eats into it, cut slides 5 and 8
(M1 and M3 detail) for 12:30 — the flow chart on slide 3 already carries their shape, and the
argument survives.

## Slides

### 1. Title
Structure to function: from a bacterial genome to ranked, structure-backed drug targets.
Project 9. Name the team.

### 2. The problem
A newly sequenced pathogen genome gives thousands of proteins and almost no idea which ones
matter. A quarter to a third are annotated "hypothetical protein". The useful question is not
"what is every protein" but "which 50 are worth a structural and chemical look".

### 3. Flow chart — six modules, one shared contract
Inline SVG. M1/M2/M3 solid, M4/M5/M6 dashed and dimmed, with the built/not-built rule drawn
underneath and each module's files listed below it. Two things to say over it:

- **Data flows one way.** M2 may read M1; M4 may read M1-M3; nothing writes upstream.
- **`report.json` is the contract** — each module owns exactly one top-level key.

State here that three of six exist. Slide 12 repeats it; that is on purpose.

### 4. Input — what actually goes in
BV-BRC CGA output, three files joined on `patric_id`: `proteins.faa`, `genes_proteins.csv`,
`specialty_genes_all.csv`. Real rows on the slide, not a schema. HS11286 is 5,523 proteins and
2,592 specialty rows.

Worth one sentence: column names are matched through aliases, because BV-BRC's web export and
its `p3-` CLI disagree on spelling and both are in circulation.

### 5. M1 — genome
- **APIs:** CGA submit/poll/parse, `www.bv-brc.org/api`, Minhash Similar Genome Finder, `p3-` CLI fallback
- **Filters:** genome-quality gate (*poor* not passed downstream without a decision); codon-tree ingroup chosen from sequence, not taxonomy
- **In:** contigs or a CGA job directory - **Out:** the three `m1/` files, `report.json` keys `run`/`genome`/`proteins`, offline HTML report

### 6. M2 — triage
- **APIs:** RCSB `search/v2` + `data` GraphQL, AlphaFold DB, Foldseek web API, UniProt/UniParc, ChEMBL, STRING, DIAMOND locally
- **Filters:** search at e-value < 1e-3 and identity >= 25%, 25 rows; counts as evidence only at e-value < 1e-5 **and** coverage >= 50%; Foldseek 1e-4 / 40%; human homolog DIAMOND 1e-5 / 50% coverage, close at >= 40% identity; essential ortholog 40% identity / 70% coverage / 1e-10
- **In:** the three `m1/` files - **Out:** `proteins.tsv`, `top50.tsv`, `hits.tsv`, `kg.json`, `no_pdb_hit.tsv`, `run.json`

The two-cutoff design is the point worth making: a hit can be *retrieved* and still not *count*.

### 7. M2 — the score, and the weight we changed
Nine components, all recorded. Then the calibration story: the first full HS11286 run selected
**0 uncharacterized proteins out of 1,338**, because a known target collects 0.50 from
virulence/essentiality/drug-target where a hypothetical earns 0.10. `annotation_gap` went
0.10 -> 0.30, sensitivity recorded: 0.20->6, 0.25->9, **0.30->10**, 0.35->15, 0.40->26.

Say plainly that it was changed *after* seeing a run, and that this is why it is in the
decision log. *Notes: `docs/02a-m2-pdb-evidence.md:418`.*

### 8. M3 — fold
- **APIs:** `files.rcsb.org` for experimental structures, AlphaFold DB for predicted ones; 153 `.cif` collected for G37
- **Quality gate:** PDB at identity >= 25%, coverage >= 50%, resolution <= 3.5 A; AFDB at mean pLDDT >= 70 and coverage >= 80%; local pLDDT < 70 masked as low-confidence regions
- **In:** `report.json` proteins plus M2 xrefs - **Out:** `structures/`, and a pass/fail per structure with the evidence behind it

Reuse before predicting is the principle; the gate is what makes reuse safe.

### 9. Output — one row, fully shown
One real row of `top50.tsv` on screen: rank 1, `fig|243273.25.peg.503`, a hypothetical, score
0.65, PDB `1YDX_1` at 100% identity and 100% coverage, 2.3 A, and the sentence

> PDB 1YDX_1 identity 100%, coverage 100%, 2.30 A; uncharacterized product

50 columns, written for **all 542 proteins**, not only the chosen fifty. The score recomputes
from the recorded components alone.

### 10. Our top hit was the genome matching itself
A teammate noticed the top-ranked G37 protein was a hypothetical matched at **100% identity** to
PDB `1Q8C_1` - *Mycoplasmoides genitalium* P47. The same organism.

| | Top 50 | Whole genome |
| --- | --- | --- |
| Best hit is the same species | **30** | 34 |
| Best hit is the same genus | **47** | 108 |

88% of every same-species hit in the genome landed in our top 50. Not wrong - a 100% identity
hit is the strongest evidence there is - but **easy** in a way a novel genome will not be.

These are 30 distinct proteins, each with its own structure, at identities from 0.60 to 1.00
(median 0.84). Only 7 are at 100%. Do not let the room hear "30 self-matches at 100%".
*Notes: `docs/02a-m2-pdb-evidence.md:310`.*

### 11. What we did about it
- **Flags, not penalties.** The score is unchanged; what changes is what the report may claim.
- **Do not calibrate the weights on this genome.** Written into the docs as a standing rule.

One detail that decided it: BV-BRC writes *Mycoplasma genitalium* where the PDB says
*Mycoplasmoides genitalium*. Comparing genus strings naively marks every self-match as foreign
- reporting the exact opposite of the truth, silently.

### 12. Status - what we have, what is left
A table, module by module: M1 built, M2 built, M3 partial (collection and gate done, pockets
open), M4/M5/M6 not built. Then the limits in one line: no docking, no literature synthesis,
nothing experimentally validated, no clinical claim; scores are rankings not affinities and
predicted structures lack cofactors.

### 13. Demo
Offline, from committed fixtures, no network and no credentials:

    python -m s2f.m2_triage --dry-run --run runs/demo --report

Open one protein and read its score, components, flags and reason.

### 14. Next
The three built modules are a working, reproducible front half. M4/M5/M6 are scoped and
documented, and the same-organism result is the thing to carry into the next genome.

## Rehearsal checklist

- [ ] Timed once, out loud, start to finish
- [ ] Demo run on the presenting machine, offline, network physically off
- [ ] Fallback ready: pre-rendered output or a recording, in case the live run fails
- [ ] Every number on a slide traced back to the doc named in its notes
- [ ] Slide 3 and slide 12 both delivered — coverage stated twice, on purpose

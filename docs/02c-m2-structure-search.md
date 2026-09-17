# M2c — Foldseek structure search

Companion to [02-m2-triage.md](02-m2-triage.md) for M2 step 3 (`foldseek_search()`). Issue
[#9](https://github.com/NIAID-BRC-Codeathons/structure-to-function/issues/9). Implemented in
`s2f/m2_triage/foldseek.py`; run with
`python -m s2f.m2_triage --run runs/<id> --map-ids --foldseek`.

## Why, when we already search the PDB

The PDB route in [02a](02a-m2-pdb-evidence.md) is a **sequence** search, and its floor is real:
below roughly 20–25% identity, alignment cannot separate homology from chance. Structure survives
much longer than sequence, so a protein can share a fold, an active site and a function at 15%
identity and be invisible to every sequence method.

Measured on ten G37 proteins the sequence search found nothing for — **seven gained an
informative fold match**, every one at 10–22% identity, i.e. all below the 25% cutoff the
sequence search uses:

| Protein | e-value | Identity | Fold match |
| --- | --- | --- | --- |
| P47257 *hypothetical* | 5e-19 | 16.6% | ArgX/LysX amino-acid ligases, three consistent structures |
| Q49412 *putative esterase* | 1.6e-12 | 15.9% | bacterial esterases — confirms the annotation |
| Q49411 *uncharacterized MFS transporter* | 2.6e-08 | 14.5% | S1P transporter Spns2 — confirms the MFS fold |
| P47579 *putative Fe-S scaffold* | 2.1e-07 | 20.4% | SufU / Fe-S cluster assembly |
| P47595 *hypothetical* | 2.6e-05 | 17.7% | HU nucleoid-associated protein, DNA-bound |

The three that failed did so honestly: two had e-values too weak to trust, one matched a target
that is itself uncharacterized.

## Same fold is not same function

TIM barrels and Rossmann folds are scaffolds shared by enzymes doing unrelated chemistry, so a
fold match is a **lead**, not an assignment. Every hit carries its e-value, alignment coverage,
identity and the target's own description, and `informative_target` is false when the target's
name says nothing — a match to another "conserved hypothetical protein" is a recorded
non-result.

## Significance gate

`e-value ≤ 1e-4` and query coverage ≥ 0.4. The e-value threshold came from the pilot: a hit at
**e-value 1.89** named a real protein (a winged-helix DNA-binding protein) and meant nothing. A
summary that counts named targets without gating on significance will overstate the yield — mine
did, at first, reporting 10/10 where the honest number was 7/10.

## What gets searched

Foldseek needs a structure as the query, so the candidates are proteins that **have a structure
but no usable sequence hit**:

| G37 | Count |
| --- | --- |
| No qualifying PDB hit | 206 |
| → with an AlphaFold model usable as a search query | 113 |
| → no structure from any source | ~93 — needs a prediction first (M3), or the ESM Atlas |

**The search gate is deliberately looser than the docking gate.** `usable_for_docking()` requires
pLDDT ≥ 70 *and* ≥ 80% coverage, which a pocket needs; a search does not. A model at pLDDT 65
covering half the protein is a poor receptor and a perfectly good search query, and the alignment
span comes back per hit anyway. Using the docking gate here dropped 7 of 12 candidates before
they were ever searched; `usable_for_structure_search()` (pLDDT ≥ 50) searched 10 of 12.

## Yield, measured on G37

163 candidates (every protein with no qualifying PDB hit that has an AlphaFold model):

| Outcome | Count |
| --- | --- |
| Structural neighbour found | **22** — of which **21 informative**, 1 matching only another unknown |
| Searched, nothing above the gate | 29 |
| No model good enough to search with | 10 |
| Refused by the public server (see below) | 102 |

Median identity of the best hit: **16.6%** — far below the 25% floor of the sequence search, which
is the whole point.

## The public server will not do a whole genome

Measured, not assumed: batches of 10–25 searches succeed completely. A 163-search run is refused
with **HTTP 429** for most of it — 121 failures unpaced, still 102 with 2-second pacing and three
attempts. The 429 carries **no `Retry-After` and no rate-limit headers**, so the correct wait
cannot be computed, only guessed.

The client now backs off 30 s then 90 s on a 429 and **stops after five consecutive refusals**
with an actionable message rather than hammering a shared server. Searches that did complete are
cached and replay for free.

**So the boundary between the two routes is not speed, it is feasibility:**

| Route | Good for | Not good for |
| --- | --- | --- |
| Web API (default) | development, spot checks, teammates with no database | a genome-scale run — it will be refused |
| Local `foldseek` + database | the full genome, the final run, TM-scores | requires a multi-gigabyte download |

## Transport: the web API, not a local database

The Foldseek web server's ticket API (submit → poll → fetch). No multi-gigabyte download, so
teammates can run this as-is, and results are cached by structure checksum so a rerun costs
nothing and `--offline` replays with no network.

A local `foldseek` binary with a downloaded database would be faster for a final full run and is
**the only way to get a true TM-score** — the web API returns e-value, probability, identity and
alignment spans, but no TM-score. #9's definition of done asks for TM-scores, so either the local
route is used for the final run or that box is amended to match what the API provides. Recorded
rather than quietly substituted.

The response carries target C-alpha coordinates per hit (`tCa`), megabytes across a genome, so
only the fields we use are kept before caching.

## Status values

`found`, `no-hit`, `query-failed`, `not-queried-no-structure` — never collapsing a failure or an
absent structure into "no homolog exists".

Related: this work also fixed the upstream flag that fed it. `no_pdb_hit` was true for
`query-failed` rows too, so a network failure landed in the "no homolog" file and could have been
laundered into a structural claim. Found by @Ashita2619 while building [#43](https://github.com/NIAID-BRC-Codeathons/structure-to-function/issues/43)
on top of it; `pdb_search_failed` is now a separate flag.

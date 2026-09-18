# Evaluation — measuring what the pipeline's rankings actually do

Issues #49 and #12. Scope: measure the two scores that order proteins, report the numbers
before anyone touches a weight, and give the project a scoreboard it can quote.

    python -m s2f.evaluation --run runs/<run_id> --score m1_priority --top-k 50
    python -m s2f.evaluation --run runs/<run_id> --score triage
    python -m s2f.evaluation --dry-run          # committed fixture, no network, no credentials

Outputs land in `runs/<run_id>/eval/`: `metrics.json`, `truth_matches.tsv`, `summary.md`
(short enough to paste into an issue), and `run.json`. Three more appear as the pipeline
stages they depend on become available: `ablation.json`, `outcome.json`, `agreement.json`.
`--figures` writes `eval/figures/` (needs matplotlib; a missing figure never fails a run).

## Figures

One claim each, so any panel stands alone, plus a four-panel composite and the legends.

| file | claim |
| --- | --- |
| `seq_cap_loss` | M1's cap retains 72% of all CDS but only 42% of the hypothetical ones |
| `score_disagreement` | `m1_priority` rank against `triage` rank; 6 of 50 shared at the cut |
| `triage_ablation` | 5 of 9 triage components move the top 50; four move nothing |
| `outcome_lift` | precision 100% beside a 99.3% base rate — the label is not a test |
| `figure_composite` | all four as panels A–D at double-column width |
| `captions.md` | numbered legends, each self-contained |

The PDFs from the HS11286 run are committed under
[`figures/`](figures/) with their provenance and the command that regenerates them.

### `--figure-style publication` (default)

- **No title inside the artwork.** The legend lives in `captions.md`, where a journal wants
  it; a baked-in title has to be removed in production.
- **Real column widths**: 89 mm single, 120 mm for the ablation panel's long labels, 183 mm
  for the composite. Arial at 7 pt, which survives reduction.
- **PDF, SVG and 600 dpi PNG.** matplotlib defaults to `pdf.fonttype = 3` and
  `svg.fonttype = "path"`; Type 3 is rejected outright by Nature, Elsevier and IEEE, and
  outlined SVG text cannot be copy-edited or read by a screen reader. Publication mode sets
  fonttype 42 and live SVG text, and `tests/test_eval_figures.py` asserts both on the
  written files rather than trusting the setting.
- Mathtext is pinned to Arial too, or `$w$` and `$\rho$` fall back to DejaVu and the
  italics sit in a different typeface from every other glyph.

`--figure-style slide` puts one headline in each figure at larger type, for a talk.

### Two deliberate choices a reviewer will ask about

**The outcome panel uses a full 0–100 baseline.** Truncating it would magnify a 0.7-point
gap into a visible difference and argue the opposite of the panel's point: the two bars
being the same height *is* the finding.

**A component whose provider never ran is a hatched ghost reading "not measured",** not a
zero-length bar. "We did not measure this" and "we measured zero" are different claims and
a blank bar merges them.

Colour is two slots from a validated categorical palette, checked in both light and dark
for lightness band, chroma floor, CVD separation, normal-vision separation and surface
contrast. Identity never rests on colour alone; every series is direct-labelled.

## Four measurements, three of which need no curation

| what | needs | answers |
| --- | --- | --- |
| **calibration** | a curated truth set | does the ranking surface real virulence factors? |
| **outcome** | M3 | did the proteins it picked actually yield dockable structures? |
| **ablation** | M2 | which of M2's triage components decide the selection? |
| **agreement** | M1 and M2 | how far apart are the two scores, as the contract predicts? |

Only the first needs anyone to adjudicate biology. The other three score the pipeline
against its own recorded output, which is why they can be re-run on any genome without
curating a truth set for it first.

## Why this exists

The first live API run, *K. pneumoniae* HS11286, 5,523 CDS, put 14 proteins at the score
ceiling and four of those 14 were wrong: a bare `Outer membrane protein`, the periplasmic
chaperones Skp and PpiD, an IncF conjugative-transfer surface-exclusion protein, and a
T4SS inner-membrane protein. The `secretion` rule was matching envelope vocabulary rather
than secretion-system membership.

The fix is not nudged weights. Pitfall #12 is explicit: tuning after seeing results
invalidates the result. So the order is **measure, record, then change** — and the
recorded measurement has to survive someone asking how it was taken.

## This module owns no `report.json` key

`00-architecture.md` assigns every top-level key an owner and `common.schema.validate`
raises on a section it does not know, so writing an `evaluation` key would be exactly the
improvised schema change the contract forbids. The module reads a finished run and writes
its own directory, the way a reviewer would.

Adding an `evaluation` section is worth doing and is proposed as its own decision, not
smuggled in here. `tests/test_eval_dry_run.py::test_it_never_writes_a_report_section` is
what keeps the boundary honest.

## The truth set

`fixtures/eval/kpneumoniae_hs11286.truth.tsv` — 223 curated symbols: 115 positive,
97 negative, 11 excluded.

| column | meaning |
| --- | --- |
| `symbol` | gene symbol, the join key |
| `label` | `positive`, `negative` or `excluded` |
| `class` | capsule, siderophore, fimbriae, t6ss, urease, porin, lps, efflux, omp_biogenesis, general_secretion, plasmid_transfer, plasmid_maintenance, toxin_antitoxin, mobile_element |
| `basis` | why the call was made. Required — a row without one is refused at load |
| `note` | caveats, and the reason for every `excluded` row |

**Positives** are host-interaction factors: capsule export and regulation, O antigen, the
four siderophore systems plus the TonB energiser, type 1 / type 3 / common-pilus fimbriae,
type VI secretion, urease, the major porins, and the RND efflux pumps.

**Negatives** are the machinery that must *not* rank high and does under the current rules:
outer-membrane protein biogenesis (Skp, SurA, PpiD, the Bam and Lol complexes), the Sec and
Tat general export pathways, the IncF conjugative transfer region, plasmid replication and
partition, plasmid addiction toxin–antitoxin pairs, and mobile-element machinery.

**Excluded** is the third label, and it is the one that keeps the measurement honest. A
contested gene counted either way smuggles an unargued claim into the numerator. `traT` is
the worked example: issue #49 lists the IncF surface-exclusion protein among its false
positives, but TraT is also a well-described complement-resistance factor, so it is
recorded as contested and scored in neither direction. `degP`/`htrA`, `dsbA`, `lpp` and the
chromosomal `relBE`/`mazEF`/`higBA` systems are excluded for the same reason.

### Where the truth set is weakest

- **T4SS is context-dependent.** The `virB*` and `virD4` rows are negatives because on an
  IncF plasmid in *K. pneumoniae* they are conjugation machinery. In *Legionella*,
  *Bartonella* or *Brucella* a T4SS is a genuine host-directed effector translocator. This
  truth set is species-scoped and reusing it elsewhere would be wrong.
- **Efflux sits on a different axis.** `acr`, `oqx` and `eef` are AMR determinants, not
  host-interaction factors, and `WEIGHTS["amr"]` scores them +2. Issue #49 names efflux as
  part of the truth set so they are positives here, but a reader should know that a
  precision figure mixes "found a virulence factor" with "found a resistance determinant".
  Splitting the axes is the obvious next iteration.
- **Coverage is sparse.** 223 symbols against 5,523 CDS. See the precision definition below.

## What is measured

### Precision is over the labelled subset, not over K

Dividing by K would count every unadjudicated protein in the top 50 as a false positive,
which is false — most of them are simply proteins nobody has labelled. So

    precision = positives / (positives + negatives)   # both inside the top K
    recall    = distinct positive symbols in top K / positive symbols present in the genome

and `unlabelled` plus `labelled_coverage` are reported beside them. A precision of 0.8 over
10 labelled proteins is a different claim from 0.8 over 40, and the summary says which.

A symbol absent from the genome leaves the recall denominator rather than counting as a
miss. HS11286 is a ST11 carbapenem-resistant clinical isolate, not a hypervirulent strain,
so `rmpA`, the aerobactin locus and the salmochelin locus are expected to be missing.

### Ties

`ceiling_ties` and `cut_inside_tie_block` report what the ranking does at the cut. If the
top 50 ends inside a block of equal scores, the proteins above the line got there on the
length tie-break, not on evidence, and the precision figure inherits that arbitrariness.

### Baselines

| baseline | what it is | what beating it proves |
| --- | --- | --- |
| `keyword` | one regex (`PATHOGENESIS_RE`) over the product string, ties by length | that the nine weighted components earn their complexity |
| `length` | longest protein first | that the ranking is not a length proxy |
| `random` | seeded shuffle | the floor |

`keyword` is the one that matters. It is the cheapest thing anyone could have shipped
instead, and a score that does not beat it is not paying for itself.

### The leakage audit

`m1_priority` is computed from `gene + product + specialty classifications`. A truth set
built by reading product strings would score `MECHANISM_RULES` against itself and report
near-perfect precision however wrong the rules were. Two defences:

1. **The truth set is keyed on gene symbols from the gene-family literature**, not read off
   this genome's annotation text.
2. **Every top-K score is split** into evidence a curated database supplied (VFDB, Victors,
   CARD, NDARO, drug-target and essentiality tables) and evidence the product string
   supplied (`mechanism`, `surface`, `named_gene`). A positive that scored only on
   text-derived components is the regex agreeing with the string it read, and the summary
   names those proteins individually.

If `priority.py` ever emits a score reason the audit cannot attribute, the run prints a
warning and the summary says the split is incomplete, rather than quietly filing the
unknown reason as text-derived.
`tests/test_eval_metrics.py::test_every_reason_priority_can_emit_is_classified` fails first.

## First measurement — HS11286, 2026-09-17

Weights `275c864450ec`, recorded **before** any weight change. 5,523 scored proteins,
156 matched against the truth set, top 50.

| ranking | precision | recall | pos | neg | unlabelled |
| --- | --- | --- | --- | --- | --- |
| **`m1_priority`** | **84.2%** | **32.7%** | 16 | 3 | 30 |
| baseline `keyword` | 75.0% | 30.6% | 18 | 6 | 25 |
| baseline `length` | 75.0% | 10.2% | 6 | 2 | 42 |
| baseline `random` | 0.0% | 0.0% | 0 | 1 | 49 |

Precision is over the 19 labelled proteins in the top 50; recall over the 49 curated
positives present in this genome. Reproduce with

    python -m s2f.m1_genome --run runs/kp_hs11286 --from-bvbrc-api --genome-id 1125630.4
    python -m s2f.evaluation --run runs/kp_hs11286 --top-k 50 --include-product-matches

### What it says

**The ranking works, modestly.** It beats the one-regex `keyword` baseline on precision,
84.2% against 75.0%, by being more selective: `keyword` finds two more true positives and
twice as many false ones. Nine weighted components buy about nine points of precision over
a single regex, which is a defensible return but not a large one.

**It is not the regex agreeing with itself.** All 16 true positives in the top 50 have
curated-database evidence behind them and none scored on product text alone. That is the
result the leakage audit exists to produce, and it came out in the scorer's favour.

**The ties are worse than #49 reported.** The 14-way ceiling tie at score 6 is real, but
the cut at 50 lands inside a **74-protein block all tied at score 5**. Roughly half the
shortlist is chosen by the length tie-break rather than by evidence. Discriminating inside
that block matters more than fixing the ceiling.

**Tightening the `secretion` rule is necessary and not sufficient.** The three negatives in
the top 50 — SurA, Skp and VirB10 — each score 6, and in every case only +2 of that comes
from the mechanism rule. SurA and Skp take +3 from a curated VFDB virulence hit, VirB10
takes +3 from a drug-target table. Removing the envelope-vocabulary match drops them from
6 to 4 and leaves them in the top band. The curated hit is the larger term.

**A scorer bug, found on the way.** `MECHANISM_RULES["iron"]` contains `\btonb\b`, and a
hyphen is a word boundary, so "TonB-like" fires it. Two VirB10 proteins are consequently
tagged `iron` and told the reader they "likely scavenge host iron". Score impact is nil —
`mechanism` is a flat +2 and both already matched `secretion` — so this is a reporting
defect in `categories` and `mechanism_hypothesis`, not a ranking one. It is the same trap
as pitfall #27, inside the scorer rather than the harness.

**One of nine components is dead on this route.** BV-BRC returns no `gene` field for any of
HS11286's 5,523 CDS, confirmed against the live Data API, so `WEIGHTS["named_gene"]` can
never fire and the strong match tier has nothing to join on. The run says so rather than
reporting a silent zero.

### Caveats on this number

- Matched entirely through the **weak product tier**, because there are no gene symbols.
  Treat it as a first measurement, not a settled figure.
- Labelled coverage is 38%: 19 of the top 50 are adjudicated, 30 are not.
- Precision mixes the host-interaction and AMR axes, as described above.

## Second measurement — the triage score, HS11286, 2026-09-18

`python -m s2f.m2_triage --run runs/kp_hs11286 --top 300 --report` over 4,000 proteins (M1
caps sequence retrieval at 4,000), then `python -m s2f.m3_fold`, then the harness with
`--score triage`. Weights fingerprint **`4879c6253690`**.

> **Re-measured after #60.** The first pass ran against an eight-component triage score,
> fingerprint `f2713321fabd`. PR #60 added a ninth, `ligandable_homolog` at weight 0.10,
> **without bumping `WEIGHTS_VERSION`**, which is still `2026-09-16`. The fingerprint hashes
> the weight dict rather than trusting that constant, so the change was detected anyway —
> and it moved the result enough to overturn one of the conclusions below. Numbers taken
> under `f2713321fabd` do not compare to these.

**Caveat:** this M2 pass ran with no optional providers — no `--annotate`,
`--human-homology`, `--essentiality` or `--foldseek` — so three components are
under-exercised, and the harness reports that rather than scoring them as inert.

### The two scores select almost disjoint populations

| | |
| --- | --- |
| Spearman rho over the 4,000 both scored | **+0.348** |
| shared between the two top 50s | **6 of 50** (Jaccard 0.06) |

`m1_priority`'s top 50 is named virulence factors: FimD, fimbrial adhesins, TonB-dependent
receptors, phospholipase A1. `triage`'s top 50 is **33 of 50 hypothetical, putative or
uncharacterised** proteins with good structural evidence, against a 14% base rate in the
scored population.

`00a-data-contract.md` predicted disagreement. The magnitude at the top is near-total, and
it is the pipeline working as named: this project is a *hypothetical-protein* factory, and
`triage` is the score that acts on that brief. `m1_priority` ranks what we already know.

**The consequence matters more than the number.** M3 reads `triage.selected`, not
`m1_priority`. So the four false positives issue #49 opens with — SurA, Skp, VirB10 — sit
in M1's top band and **never reach M3 at all**. #49 is a defect in M1's report, not in the
pipeline's selection, and fixing it will not change one protein that gets folded or docked.

### Which components decide the selection

Leave-one-out at K=50, against the full score. The component list is read from `WEIGHTS`,
so #60's ninth component was swept without a code change.

| component | weight | fires on | moves out of top 50 | Δ precision vs M3 outcome |
| --- | ---: | ---: | ---: | ---: |
| `pdb_evidence` | +0.40 | 3056 | **34** | **−4.5 pp** |
| `annotation_gap` | +0.30 | 563 | **29** | 0 |
| `virulence_amr` | +0.20 | 250 | **23** | 0 |
| `ligandable_homolog` | +0.10 | 2190 | **17** | 0 |
| `drug_target` | +0.15 | 347 | 1 | 0 |
| `essential` | +0.15 | 126 | 0 | 0 |
| `surface_bonus` | +0.10 | **0** | 0 | 0 |
| `membrane_penalty` | −0.15 | **0** | 0 | 0 |
| `human_homolog_penalty` | −0.25 | 3 | 0 | 0 |

- **Five of nine components move the selection at all.** `drug_target` moves one protein;
  `essential` and `human_homolog_penalty` fire but change nothing at this cut.
- **`annotation_gap` is disproportionate to its reach.** It fires on 14% of the scored
  proteome (563 of 4,000) and is the second-largest mover. A +0.30 bonus for reading
  "hypothetical" is doing more work than its weight suggests.
- **`ligandable_homolog` is material on arrival**: it fires on 2,190 proteins, more than
  half the population, and moves 17 of the top 50 in its first measured run.
- **Two never fire.** `surface_bonus` and `membrane_penalty` are `OPTIONAL_COMPONENTS` whose
  providers did not run here. That is a data gap, not a finding about the design, and the
  run says so instead of reporting them as inert.

> **A conclusion this replaces.** Under the eight-component score, `annotation_gap` was the
> single largest mover (41 of 50, against `pdb_evidence`'s 33), and this document said so.
> With `ligandable_homolog` in play the ordering reverses and `pdb_evidence` leads. The
> claim was true of that scorer and is not true of this one. This is exactly the failure
> mode the fingerprint exists to catch, so it is recorded rather than quietly overwritten.

### The M3 outcome label cannot discriminate on this run

M3 collected structures for all 300 selected proteins: 0 required prediction, 0 failed, and
**298 of 300 passed the provisional quality gate**.

| | |
| --- | --- |
| precision @ 50 against `usable_for_docking` | 100.0% |
| base rate of `dockable` in the labelled set | 99.3% |
| **lift** | **+0.7 pp** |

A label that is 99.3% one class cannot separate a good ranking from a bad one, so the
headline 100% is not evidence that `triage` selects well. Two reasons, both actionable:

1. **The gate is not currently a test.** `provisional-2026-09-17` passes essentially
   everything handed to it. Its thresholds need tightening before it can validate anything.
2. **The measurement is partly circular by construction.** `pdb_evidence` is 0.40 of the
   triage score, and having a PDB hit is most of what makes a structure retrievable. The
   ablation is the non-circular part: only `pdb_evidence` has a nonzero Δ against this
   label (−4.5 pp), so the other eight contribute nothing to whether a pick turns out
   dockable.

### The finding underneath all three

**Neither instrument can currently evaluate `triage`**, and they fail in opposite
directions. The curated truth set cannot reach it — 0 of its top 50 are labelled, because
curated gene symbols are by definition characterised proteins and `triage` selects
uncharacterised ones. The outcome label cannot discriminate — 99.3% base rate.

That is a result about the *evaluability* of the pipeline's central decision, and it says
what a third instrument would have to be: a held-out set of proteins that were hypothetical
at annotation time and have since been characterised, so that "we would have picked this
one" can be checked against what it turned out to be. `gene-function-prediction` from the
2025 codeathon solved exactly this ground-truth design problem and is the place to start.

### The seam defect: M1's cap removes the population M2 exists to find

Not visible from inside either module, which is the argument for measuring across them.

`--seq-cap` defaults to 4,000 and truncates **by `m1_priority` rank** — deliberately, with
a code comment explaining that this beats truncating by position on the chromosome. That
reasoning is right. But `m1_priority` scores a hypothetical protein **0**: no gene name, no
mechanism keyword, no specialty hit. So hypothetical proteins sort to the bottom and are
cut first, and M2's `annotation_gap` then pays 0.30 for exactly that property.

| | |
| --- | --- |
| gap-like proteins in the genome | 1,333 of 5,523 |
| reached M2 | 560 |
| **never scored at all** | **773 (58%)** |
| overall drop rate | 28% |
| enrichment of the drop against gap-like proteins | **2.1×** |

The dropped set is exactly `m1_priority` ranks 4001–5523; every one scored 0 or −2.

Both policies are locally sensible and together they cancel. The project is named
*Hypothetical Protein and Binding-Evidence Factory*.

**Severity.** This does not fire on the chosen blinded test genome: USA300 is 2,741 CDS,
under the cap. It fires on *K. pneumoniae* at 5,523 — the genome the pipeline has actually
been developed against — and on *E. coli* and most Enterobacteriaceae. Chlamydiales, the
charter's named target at roughly 900 CDS, is also safe. Latent, but real, and it will bite
the first larger genome.

The fix is small and the choice is the team's: rank for the cap by something that does not
zero the hypotheticals, exempt gap-like products from truncation, or raise the default.
`rankings.population_stages` is what measures it, and `GAP_LIKE` there is a *population*
definition kept deliberately separate from M2's `annotation_gap` scoring component.

### The finding underneath all three

**Neither instrument can currently evaluate `triage`**, and they fail in opposite
directions. The curated truth set cannot reach it — 0 of its top 50 are labelled, because
curated gene symbols are by definition characterised proteins and `triage` selects
uncharacterised ones. The outcome label cannot discriminate — 99.3% base rate.

That is a result about the *evaluability* of the pipeline's central decision, and it says
what a third instrument would have to be: a held-out set of proteins that were hypothetical
at annotation time and have since been characterised, so that "we would have picked this
one" can be checked against what it turned out to be. `gene-function-prediction` from the
2025 codeathon solved exactly this ground-truth design problem and is the place to start.

## Matching, and the trap issue #49 names

> a naive regex over product text matches `ent` inside "ATP-dependent" and "TonB-dependent",
> which silently inflates any known virulence factor count. Anchor on gene symbols or word
> boundaries.

Word boundaries are not enough, and neither were the next two guards. Each was added only
after a live run produced the false match the previous ones let through — which is the
argument for having the harness at all.

| # | guard | what it lets through without it | cost |
| --- | --- | --- | --- |
| 1 | not part of a hyphenated compound | `\btonB\b` matches *TonB-dependent siderophore receptor*, because a hyphen is a word boundary | the trap #49 names |
| 2 | case-sensitive, over two spellings | case-insensitively `iroN` **is** the word `iron`, matching all 34 proteins that read "iron acquisition" | 34 false matches |
| 3 | not followed by a family word | `Transcriptional regulator YjdC, AcrR family` is a member of the AcrR family, not AcrR | 13 of 14 `acrR` matches |

Guard 2 needs care: `str.capitalize` would turn `iroN` into `Iron` and reintroduce the bug,
so `product_forms` upper-cases only the first character. Guard 3's word list is
`family`, `superfamily`, `domain`, `fold`, `motif`, `like`, `type`, and it rejects only the
occurrence it follows — `AcrB, RND family multidrug efflux transporter` still matches
`acrB`.

Two tiers, never merged:

| tier | rule | default |
| --- | --- | --- |
| `gene` | exact, case-insensitive, on `proteins[].gene` — a controlled symbol assigned by the annotation pipeline, so no guards needed | on |
| `product` | all three guards above, on `proteins[].product` | off, `--include-product-matches` |

Every matched row records which tier found it. After all three guards no symbol matches
more than four HS11286 proteins, and those are genuine paralogues (`hcp`, `fimD`).

## Reproducing a measurement

`metrics.json` and `run.json` carry the scoring weights, a 12-character
`weights_fingerprint`, the truth-set path and its SHA-256, the input report's SHA-256, the
git commit and whether the tree was dirty, the top-K cut, the baselines and the random
seed. A later run under a different fingerprint is measuring a different scorer, and the
two numbers do not compare — which is the point of recording it.

## Acceptance checks

- [x] `--dry-run` runs from `fixtures/eval/` with no network and no credentials
- [x] the dry-run fixture validates against `common.schema.validate_report`
- [x] no `report.json` section is written
- [x] precision states its denominator; unlabelled proteins never count against it
- [x] symbols absent from the genome leave the recall denominator
- [x] hyphenated compounds do not match a bare symbol, end to end
- [x] the weights fingerprint is recorded with every measurement

## Next, in order

Done: the first measurement of `m1_priority`, the extension to `triage`, and the outcome
and ablation passes above.

1. **Tighten the M3 quality gate.** `provisional-2026-09-17` passes 298 of 300, so nothing
   downstream is being filtered and no outcome measurement can mean anything until it is.
   Cheapest high-value fix on the board, and it unblocks every later validation.
2. **Decide what `annotation_gap` at 0.30 is for.** It is the dominant selector — 45 of the
   top 50 — and it is a bonus for the annotation being *absent*. That may be exactly right
   for a hypothetical-protein factory, in which case say so in the report and stop
   describing the output as a virulence ranking. If it is not intended, it is the first
   weight to revisit, and there is now a churn number to revisit it against.
3. **Discriminate inside the score-5 block.** 74 proteins tie there and the top 50 cuts
   through it, so half of M1's shortlist is ordered by protein length. Bigger than the
   14-way ceiling tie #49 opens with, and no weight change fixes it: the score is a small
   integer and needs a continuous component or an explicit secondary key.
4. **Then** tighten the `secretion` rule to require secretion-system membership rather than
   envelope vocabulary, and re-measure under a new fingerprint. Expect less than it looks:
   it removes +2 from SurA, Skp and VirB10, each of which keeps +3 from a curated database
   and stays in the top band — and none of them reach M3 anyway.
5. **Decide what a curated VFDB hit on SurA and Skp means.** They are genuinely listed as
   virulence-associated, because OMP biogenesis is required for virulence. If the truth set
   calls them negatives and VFDB calls them positives, the two are using "virulence factor"
   in different senses and the scorer cannot arbitrate.
6. Fix the `\btonb\b` hyphen match in `MECHANISM_RULES["iron"]` (pitfall #27), which
   mislabels VirB10 as iron acquisition in `categories` and in the hypothesis sentence.
7. Split the host-interaction and AMR axes so precision stops mixing them.
8. Re-run the ablation with the optional providers enabled (`--annotate`,
   `--human-homology`, `--essentiality`), so `surface_bonus` and `membrane_penalty` are
   measured rather than reported as never having fired.
9. Build the third instrument: proteins that were hypothetical at annotation time and have
   since been characterised, so `triage` can be scored on the population it actually
   selects. `gene-function-prediction` (2025) already solved this ground-truth design.
10. Propose the `evaluation` section key so the numbers can live in `report.json` and M6 can
    render them.

## Related

- [07-decisions-and-risks.md](07-decisions-and-risks.md) — standing risk 10, validation asymmetry:
  M1 and M2 can be scored against ground truth, M4 and M5 cannot be in three days
- [pitfalls.md](pitfalls.md) — #12, triage bias: tuning after seeing results invalidates the result
- [01-m1-genome.md](01-m1-genome.md) — what `m1_priority` is, and why it is not `triage`
- [00a-data-contract.md](00a-data-contract.md) — why `m1_priority` and `triage` disagree by design

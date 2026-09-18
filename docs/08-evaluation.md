# Evaluation — calibrating a ranking against a curated truth set

Issue #49. Scope: measure `proteins[].m1_priority` against ground truth, report the number
before anyone touches a weight, and give the project a scoreboard it can quote.

    python -m s2f.evaluation --run runs/<run_id> --top-k 50
    python -m s2f.evaluation --dry-run          # committed fixture, no network, no credentials

Outputs land in `runs/<run_id>/eval/`: `metrics.json`, `truth_matches.tsv`, `summary.md`
(short enough to paste into an issue), and `run.json`.

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

1. ~~Measure against a full live run and record the numbers before any weight change.~~
   Done, above.
2. **Discriminate inside the score-5 block.** 74 proteins tie there and the top 50 cuts
   through it, so half the shortlist is ordered by protein length. This is a bigger problem
   than the 14-way ceiling tie #49 opens with, and no weight change fixes it: the score is
   a small integer and needs either a continuous component or an explicit secondary key.
3. **Then** tighten the `secretion` rule to require secretion-system membership rather than
   envelope vocabulary, and re-measure under a new fingerprint. Expect it to help less than
   it looks: it removes +2 from SurA, Skp and VirB10, each of which keeps +3 from a curated
   database and stays in the top band.
4. **Decide what a curated VFDB hit on SurA and Skp means.** They are genuinely listed as
   virulence-associated, because OMP biogenesis is required for virulence. If the truth set
   calls them negatives and VFDB calls them positives, one of the two is using "virulence
   factor" in a sense the other does not, and the scorer cannot resolve it.
5. Fix the `\btonb\b` hyphen match in `MECHANISM_RULES["iron"]` (pitfall #27), which
   mislabels VirB10 as iron acquisition in `categories` and in the hypothesis sentence.
6. Split the host-interaction and AMR axes so precision stops mixing them.
7. Extend the harness to M2's `proteins[].triage`, which is a different question — structural
   tractability, not host interaction — and needs its own truth set.
8. Propose the `evaluation` section key so the numbers can live in `report.json` and M6 can
   render them.

## Related

- [07-decisions-and-risks.md](07-decisions-and-risks.md) — standing risk 10, validation asymmetry:
  M1 and M2 can be scored against ground truth, M4 and M5 cannot be in three days
- [pitfalls.md](pitfalls.md) — #12, triage bias: tuning after seeing results invalidates the result
- [01-m1-genome.md](01-m1-genome.md) — what `m1_priority` is, and why it is not `triage`
- [00a-data-contract.md](00a-data-contract.md) — why `m1_priority` and `triage` disagree by design

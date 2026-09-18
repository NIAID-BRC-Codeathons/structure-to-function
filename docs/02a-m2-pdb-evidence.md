# M2a — PDB evidence and interim triage weights

Companion to [02-m2-triage.md](02-m2-triage.md), covering M2 step 2 (`lookup_structures()`) and a
runnable first version of step 6 (`score_and_select()`). Issues
[#8](https://github.com/NIAID-BRC-Codeathons/structure-to-function/issues/8) and
[#12](https://github.com/NIAID-BRC-Codeathons/structure-to-function/issues/12).

`02-m2-triage.md` is the module contract. This doc holds the tunable specifics: query
parameters, scoring components, weights, and the change log that pitfall #12 requires. **The
weights below were written before the first full run.** Any later change gets a dated row in
[Weight change log](#weight-change-log) and is named in the commit message.

## Status

Standalone prototype for the `report.json` contract: `common/schema.py` and `common/io.py`
(issue #2) do not exist yet, so this module reads M1-shaped files directly and writes its own
TSV/JSON. When the schema lands, an adapter maps the output into `proteins[].xrefs`,
`annotations[]`, `flags` and `triage`.

Two shared pieces now live in `s2f/common/` rather than in this module:

- `common/http.py` (issue #4) — the cached, retrying HTTP client every outbound call uses.
  Provenance stamping, the other half of that issue, is still open.
- `common/ids.py` (issue #3) — identifier mapping. This module maps nothing itself
  (pitfall #11); `--map-ids` calls the shared helper. RCSB sequence search needs no UniProt
  accession, which is why the PDB half could run before the mapper existed.

## Input contract

A directory of M1 output, `runs/<run_id>/m1/`:

| File | Required | Contents |
| --- | --- | --- |
| `proteins.faa` | yes | protein FASTA, headers `>fig\|<genome>.peg.<n>  <product>` |
| `genes_proteins.csv` | yes | `patric_id`, `refseq_locus_tag`, `gene`, `product`, `aa_length`, `start`, `end`, `strand`, `plfam_id`, `pgfam_id` |
| `specialty_genes_all.csv` | no | `property`, `gene`, `product`, `source`, `classification`, `antibiotics_class`, `identity`, `query_coverage`, `patric_id` |

Canonical key: BV-BRC `patric_id` (the `feature_id` of `00-architecture.md`). Column names are
matched through aliases (`patric_id`/`BRC ID`/`feature_id`, `product`/`Product`, …) because
BV-BRC exports differ between the web UI and the `p3-` CLI.

`runs/` is git-ignored. Real M1 output and BV-BRC downloads stay local; this slot is a
placeholder until M1 (issues #5, #6) writes into it. Local development uses a
*K. pneumoniae* subsp. *pneumoniae* HS11286 export (5,523 proteins, 2,592 specialty rows), which
is **not committed**. Everything committed for tests lives in `fixtures/m2/`.

## PDB lookup

RCSB Search API, `service: sequence` (MMseqs2), `return_type: polymer_entity`:

| Parameter | Value | Why |
| --- | --- | --- |
| `evalue_cutoff` | `1e-3` | Permissive at retrieval; the scoring gate is stricter. |
| `identity_cutoff` | `0.25` | Below this, sequence-based function transfer is not defensible. |
| `rows` | 25 | Enough for a best hit plus context; the service reports thousands. |
| `results_verbosity` | `verbose` | Needed for `match_context`: identity, e-value, bitscore, query/subject ranges. |

Each distinct sequence is queried once (SHA-256 grouping), responses are cached on disk, and
`--offline` replays from cache with no network. Retrieval status per protein is recorded as
`found`, `no-hit` or `query-failed` — a failed query is never recorded as an absence of hits.

Entry and entity metadata come from batched RCSB GraphQL (100 entity IDs per call): experimental
method, resolution, entity description, EC numbers, Pfam/GO/InterPro annotations, non-polymer
components, other polymer entities in the entry, and source organism.

## Scoring

### Per hit

A hit counts only if `evalue ≤ 1e-5` **and** query coverage ≥ 0.5, where
`query_coverage = (query_end − query_beg + 1) / query_length`.

```
seq      = identity × query_coverage
quality  = 1.0   X-ray or EM, resolution ≤ 2.5 Å
           0.8   X-ray or EM, resolution ≤ 3.5 Å
           0.6   X-ray or EM, worse than 3.5 Å
           0.7   NMR
           0.6   anything else or unknown
bonus    = +0.10  entry has a non-additive ligand
           +0.05  entry has another polymer entity (possible partner)
           +0.05  entity carries an EC number or GO term
hit_score = seq × quality + bonus
```

`pdb_evidence` for a protein is its best `hit_score`, **normalized** by the best score
attainable (`1.0 + 0.10 + 0.05 + 0.05 = 1.20`), so it lands in 0–1 without clipping. Ties break
on the number of distinct PDB entries hit, then bitscore.

A ligand in an entry is **not** a binding measurement. The field is `has_ligand_in_entry`, never
"binder"; affinity data is M4's job. Crystallization additives are excluded:
HOH, GOL, EDO, PEG, PG4, PGE, SO4, PO4, ACT, DMS, MPD, TRS, IMD, FMT, ACY, CIT, MES, EPE, NA, K,
CL, MG, CA, ZN, MN, NI, CD, CU, IOD, BR.

Metals are on that list for ligand-bonus purposes, but a metal in a homolog's site matters to M3
(pitfall #2), so the metals present are still reported in `hits.tsv`.

### Interim triage score

Components come from PDB evidence plus the M1 specialty table. Components that need issue #10
(secreted/surface, membrane) are absent, not silently zero-weighted — they enter when that lands.

| Component | Weight | Definition |
| --- | --- | --- |
| `pdb_evidence` | **0.40** | above |
| `virulence_amr` | **0.20** | 1.0 if a `Virulence Factor`/`Virulance factor` (the export contains both spellings) or `Antibiotic Resistance` specialty row exists; else 0 |
| `essential` | **0.15** | 1.0 if an `Essential Gene` specialty row exists; else 0 |
| `drug_target` | **0.15** | 1.0 if a `Drug Target` specialty row exists (DrugBank, TTD); else 0 |
| `annotation_gap` | **0.30** | 1.0 if the product is hypothetical, uncharacterized, putative or a DUF; else 0 |
| `ligandable_homolog` | **0.10** | 1.0 when a real ligand is observed bound in the best hit's entry, or a ChEMBL target maps to the resolved accession; **unmeasured when there is no hit and no mapping** |
| `surface_bonus` | **0.10** | 1.0 if #10 calls the protein secreted or surface-exposed; 0 if it called it neither; **unmeasured when no provider covered it** |
| `membrane_penalty` | **−0.15** | 1.0 if #10 calls it a membrane protein — deprioritized, never excluded (pitfall #3) |
| `human_homolog` | **−0.25 × identity/100** | applied when a `Human Homolog` specialty row exists (identity is always populated in the export) |

```
triage_score = 0.40·pdb_evidence + 0.20·virulence_amr + 0.15·essential
             + 0.15·drug_target + 0.30·annotation_gap + 0.10·surface_bonus
             + 0.10·ligandable_homolog
             − 0.15·membrane_penalty − 0.25·(human_identity/100)
```

Rationale, so the numbers are arguable rather than arbitrary:

- `pdb_evidence` leads because this milestone is about structure-backed candidates, but at 0.40
  it cannot by itself outrank a virulent, essential protein with a weaker structure.
- `annotation_gap` is large enough for an uncharacterized protein to compete, because the project
  exists to interpret hypothetical proteins. It still cannot carry a protein on its own: with no
  PDB evidence, 0.30 leaves a protein far below the selection cutoff, so the ones that surface are
  uncharacterized **and** structurally supported.
- `human_homolog` is a penalty and never a bonus, per pitfall #4. Its usefulness for borrowing
  ligand scaffolds is kept as a flag for M4, not as score credit.

Flags carried but never scored: `has_ligand_in_entry`, `transporter`, `metal_resistance`,
`no_pdb_hit`, `human_homolog_identity`, `pdb_hit_organisms`.

### Selection

All proteins are ranked and kept. The top N (default 50, per issue #24) get `selected: true`.
Every protein carries `reason`, including discarded ones. Selection is reproducible from the
recorded components alone: `triage_score` is recomputed from `proteins.tsv` columns in a test.

## AlphaFold DB (issue #8)

`--map-ids` also looks up AlphaFold DB by the resolved UniProt accession, so M3 does not re-fold
what already exists (`03-m3-fold.md`). Two numbers are recorded, because "a model exists" is not
enough on its own:

- **Confidence** — mean pLDDT with AlphaFold's own bands (>90 very high, 70–90 confident,
  50–70 low, <50 very low). A pocket in a low-confidence region is not a pocket.
- **Coverage** — a prediction may cover only a fragment of the accession. SARS-CoV-2 orf1ab
  (7,096 aa) returns a high-confidence model for residues 1368–1493, which is 2% of the protein.
  `afdb_usable` therefore gates on pLDDT ≥ 70 **and** coverage ≥ 80%, and `afdb_reason` says
  which one failed.

API behaviour: 200 with a list = found, 404 = no model (e.g. titin), 400 = malformed accession.
A 404 is an answer, not a failure, and is recorded as `no-model`.

On G37: 476 accessions queried, all with models — 226 very high, 196 confident, 44 low, 10 very
low; **422 usable**. More useful for M3, of the 206 proteins with no PDB hit, **113 have a usable
predicted model**, leaving 93 with no structure from either source. Of the 147 uncharacterized
proteins with no PDB hit, 58 are covered this way.

A predicted model has no cofactors, metals or ligands (pitfall #2), so `holo_homolog` from the
PDB side remains the stronger signal for docking. **Experimental evidence outranks a prediction**
wherever both exist — the precedence order is in
[00a-data-contract.md](00a-data-contract.md). This is why AlphaFold contributes nothing to the
triage score: it records what M3 will have to work with, not how good a candidate the protein is.

## Human homology, computed rather than borrowed (issue #29)

`human_homolog_penalty` used to read BV-BRC's precomputed `Human Homolog` rows. Those exist **for
public genomes only**: a fresh CGA run on G37 returns 17 specialty rows where the public record
has 174. On a blinded genome the penalty would silently never fire — and it is the only thing
standing between the pipeline and recommending a target whose close human counterpart makes it a
toxicity risk (pitfall #3/#4).

`--human-homology` computes it: DIAMOND blastp against the reviewed human proteome (UniProt
UP000005640, 20,416 sequences), **`--very-sensitive`**, e-value ≤ 1e-5, query coverage ≥ 0.5.
1.5 s for 542 proteins.

**Use `--very-sensitive`, not the default.** Measured on G37:

| Mode | Time | Proteins with a hit |
| --- | --- | --- |
| default | 0.5 s | 107 |
| `--sensitive` | 2.5 s | 173 |
| **`--very-sensitive`** | **1.4 s** | **179** |
| `--ultra-sensitive` | 4.9 s | 180 |

Default mode misses 40% of them, including `rpsL` — a real 50.6% identity, 83-residue, gapless
alignment to `RT12_HUMAN`, confirmed with an independent pairwise aligner and matching BV-BRC's
own 50%. A missed human homolog is a target we fail to penalise, so recall matters more than the
one second.

**Spot check** (the DoD for #29): all five of the genome's known human homologs come back within
~1% of BV-BRC's identity — atpD 63.7 vs 63, tuf 55.3 vs 55, atpA 54.2 vs 52, rpsL 50.6 vs 50,
dnaK 50.1 vs 50.

**Coverage gate.** BV-BRC's rows carry an identity and no coverage at all, which would let 50%
identity over 20 residues penalise a target as hard as a full-length match. A hit counts only at
≥ 50% query coverage; 26 of 179 G37 hits fall below it and are recorded as hits that do not count.

**The penalty stays continuous** (identity/100 × 0.25), so a 22% homolog costs 0.055 and a 63% one
costs 0.16. A separate `close_human_homolog` flag marks ≥ 40% identity — that is what M6 should
state as a selectivity risk and M4 should read for scaffold borrowing, because "has a human
homolog at 22% identity" is not a risk worth reporting.

## Essentiality, transferred from relatives (issue #29)

Same problem as human homology: BV-BRC's `Essential Gene` rows come from flux-balance analysis on
**public** genomes. The public G37 record has 148; a fresh CGA run has none.

`--essentiality --essentiality-keyword <genus>` builds a reference set of FBA-essential proteins
from public relatives and transfers essentiality by orthology: DIAMOND `--very-sensitive`,
**≥ 40% identity and ≥ 70% query coverage**, the conventional bar for functional transfer between
bacteria. Every call records which relative genome and what identity justified it.

**Fetch the whole reference set.** Capping it costs recall, measured against G37's 148 known
essential genes:

| Reference proteins | Genomes | Calls | Agreement with the public record |
| --- | --- | --- | --- |
| 5,000 (capped) | 95 | 131 | 113 of 148 |
| **12,742 (all)** | **96** | **169** | **148 of 148** |

The first fetch takes ~95 s and is cached; the alignment itself is ~1 s. The 21 extra calls are
proteins with strong orthologs to essential genes in relatives that the public G37 record does
not list — plausible, and each one names its source so a reviewer can check it.

**What the call means, and does not.** FBA essentiality is a metabolic model's prediction that
deleting the gene stops growth *in silico*. Transferring it by homology adds a second inference on
top. Both are recorded, `evidence` stays `FBA`, and the annotation note says plainly that this is
not an experimental knockout. That matters when M6 writes the report.

## An "Antibiotic Resistance" row usually is not resistance

BV-BRC classifies most AMR specialty rows as **"antibiotic target in susceptible species"** —
gyrA and rpoB are what fluoroquinolones and rifamycins *hit*, not genes conferring resistance.
Scoring them as resistance evidence is the specific mistake flagged in issue #30, and the ranking
we published did exactly that.

`virulence_amr` now reads the classification:

| Evidence | Value | Example |
| --- | --- | --- |
| A classified resistance mechanism | **1.0** | efflux pump, antibiotic inactivation, cell-wall charge alteration |
| An AMR row with no classification recorded | **0.5** | the export omits it; that is uncertainty, not confirmation |
| Only "antibiotic target in susceptible species" | **0.0** | scored as `drug_target` instead, which is what it is |

Measured on *K. pneumoniae* HS11286, whose export carries classifications: of 130 proteins with an
AMR row, **95** are real resistance mechanisms, **19** are unclassified, and **16** were being
counted as resistance evidence when they are drug targets. Every protein records `amr_basis`, so
a reader sees which of the three applied.

G37's public export has no `classification` column at all, so all 15 of its AMR rows score 0.5.
The CGA output does carry it — 10 of its 14 AMR rows are "antibiotic target in susceptible
species" — so the correction takes effect as soon as M1's real output is used.

## Ligandability is tractability, not potency

`ligandable_homolog` is the last component from #12's list. It fires when a small molecule has
been **observed bound** to the best structural hit, or when compounds have been assayed against
the mapped accession (ChEMBL). It says the site is tractable — M4 needs a holo template and
pitfall #2 says prefer one — and it says nothing about affinity. A ligand in a PDB entry is not
a binding measurement, which is why additives, glycans and metals are excluded upstream: 6M0J's
only non-polymer components are a glycan, a metal and an additive, so the SARS-CoV-2 RBD scores
**0** here despite being a famous complex.

**The weight was chosen by measuring the trade it makes.** Ligand-bound homologs belong to
well-studied proteins, so this component pulls directly against `annotation_gap`, which exists to
keep hypotheticals in the list. On G37 (240 of 542 proteins have a ligand-bound homolog, but only
14 were in the top 50 before):

| Weight | Ligandable in top 50 | Displaced | Uncharacterized remaining |
| --- | --- | --- | --- |
| 0.05 | 20 | 6 | 30 |
| **0.10** | **27** | **13** | **25** |
| 0.15 | 39 | 25 | 15 |
| 0.20 | 45 | 31 | 9 |

0.10 roughly doubles the holo templates available to M4 while leaving half the selection
uncharacterized. At 0.15 the list turns over by half and the discovery half collapses.

**One known interaction.** A bound ligand counts twice, by 0.133 in total: 0.10 here, plus
0.033 inside `pdb_evidence`, whose per-hit `LIGAND_BONUS` raises the normalized hit score
(0.40 × 0.10/1.20). The two were written for different reasons — one rates how useful the
structure is, the other whether the site is tractable — but they are the same fact. It is pinned
in a test so it stays visible, and dropping `LIGAND_BONUS` now that a dedicated component exists
is the obvious simplification whenever the weights are next revisited.

## Measured or merely zero

`surface_bonus` and `membrane_penalty` depend on providers that may not have run. A component
with no provider behind it contributes 0 — the same number as "measured, and false" — so each
protein records `surface_exposed_measured` and `membrane_measured`, and `run.json` carries
`components_available` per component and `components_nonzero` for all of them.

Without this, a run where DeepTMHMM never ran looks exactly like a genome with no membrane
proteins. This is the same distinction #10 makes in its own flags, and the same class of bug
@Ashita2619 found in `no_pdb_hit`.

## The smoke-test genome matches itself

Raised by @cmmann21 on #12, and it changes how the G37 numbers should be read: *M. genitalium*
**has its own structures in the PDB**, so a 100%-identity "discovery" can be the organism
matching itself. That is not a cross-organism transfer, and nothing about it will reproduce on a
novel genome.

Measured on the full G37 run:

| | Top 50 | Whole genome |
| --- | --- | --- |
| Best hit is the **same species** | **30** | 34 |
| Best hit is the **same genus** (mostly *M. pneumoniae*) | **47** | 108 |

So 47 of the 50 selected proteins rest on a hit from the same genus, and 88% of every
same-species hit in the genome landed in the top 50. The ranking is not wrong — a 100% identity
hit is the strongest evidence there is — but it is **easy** in a way a blinded or novel genome
will not be.

**Where the organism comes from.** M1 determines it and writes it to `report.json`; M2 reads
that (`genome.taxonomy`) rather than guessing. The precedence is `--organism`, then M1's genome
section, then the FASTA headers, then `undetermined` — and `undetermined` is recorded as such,
because zero same-species hits and zero measurements are not the same claim.

The lineage settles two things spelling cannot:

- **The current genus.** BV-BRC writes *Mycoplasma genitalium*; NCBI's lineage and the PDB both
  say *Mycoplasmoides*. Reading the genus from `lineage_names` makes this a plain comparison
  instead of a stem-similarity test that has to separate *Mycoplasma*/*Mycoplasmoides* (a real
  renaming) from *Streptococcus*/*Streptomyces* (two unrelated genera).
- **Identity at a named rank.** `lineage_ids` carries the species taxon (2097) alongside the
  strain (243273), so a hit annotated at species level matches a strain-level query. ID equality
  never holds there, which is the case @cmmann21 raised on #12.

The string heuristics remain for the standalone entry point, where M2 runs against an `m1/` file
drop with no `report.json` to read.

`same_species_hit` and `same_genus_hit` are recorded per protein and summarised in `run.json`;
the protein's `reason` says "not a cross-organism transfer" in words, so M6 can explain why an
annotation was easy instead of presenting it as a discovery.

**They are flags, not penalties.** The evidence genuinely is the strongest available, so the
score is unchanged; what changes is what the report is allowed to claim. The operational rule
that follows: **do not calibrate `WEIGHTS` against this genome**, because the identity
distribution that produced them does not exist elsewhere.

### Knowing which organism we are

Both flags need the query organism, and it is only sometimes in the input. Two BV-BRC FASTA
layouts are in circulation: the G37 fixture carries `[Mycoplasma genitalium G37 | 243273.25]` on
every header, while the p3-CLI/web export used for HS11286 carries no organism at all. A
bracketed suffix is therefore accepted only when it reads like a binomial and covers at least
half the records — otherwise EC-name qualifiers become the answer, and on HS11286 exactly three
of 5,523 headers end in a bracket, all of them qualifiers such as `[decarboxylating]`.

When nothing qualifies, the run warns and records `query_organism_source: "undetermined"` with
`same_organism_hits.measured: false`. This matters because the flags fail *silently and in the
wrong direction*: an unknown query organism makes every hit compare as foreign, which is the
self-match artifact reported backwards. Pass `--organism` to settle it.

Genus comparison has to tolerate a rename — BV-BRC writes *Mycoplasma genitalium* where the PDB
and UniProt now say *Mycoplasmoides genitalium* — without merging distinct genera that share a
Greek root. Two names count as one genus when their shared prefix is at least 75% of the shorter
name: *Mycoplasma*/*Mycoplasmoides* is 0.90 and *Chlamydia*/*Chlamydophila* 0.78, while
*Streptococcus*/*Streptomyces* is 0.58 and *Enterococcus*/*Enterobacter* 0.50. A fixed-length
prefix does not separate these — seven characters is both *Chlamydia*/*Chlamydophila* and
*Streptococcus*/*Streptomyces*.

## Genome sensitivity of the weights

The weights were tuned on HS11286. They are **not automatically portable**, because the
components other than `pdb_evidence` come from BV-BRC's specialty table, whose coverage varies by
organism:

| Genome | Proteins | virulence/AMR | essential | drug target | uncharacterized | no PDB hit |
| --- | --- | --- | --- | --- | --- | --- |
| *K. pneumoniae* HS11286 | 5,523 | 6.2% | 3.7% | 6.5% | 24.2% | 33% |
| *M. genitalium* G37 (fixture) | 542 | 2.8% | 27.3% | 0.2% | 35.1% | 38% |

Same weights, very different selections: HS11286 gives 10 uncharacterized in the top 50, G37 gives
32. G37 is a minimal genome with almost no drug-target or virulence annotation, so characterized
proteins cannot collect the 0.15 + 0.20 those components carry, and structure-backed hypothetical
proteins rise instead. Neither result is wrong, but the balance is a property of the genome's
annotation, not only of the weights.

Practical rules:

- Read the composition line in `run.json` (`counts`) before trusting a top-50 from a new genome.
- Retune per genome only with a dated row in the change log below, never silently.
- A future fix, if this keeps biting: normalize each component by its prevalence in the genome, so
  a rare annotation counts for more than a common one. Not done here — it would make one run's
  scores incomparable with another's, which is its own trap.

## Outputs

Written to `runs/<run_id>/m2_pdb/`:

| File | Contents |
| --- | --- |
| `hits.tsv` | every retained hit: protein, entity, identity, coverage, e-value, bitscore, method, resolution, ligands, metals, organism, hit_score |
| `proteins.tsv` | every protein, all score components, best hit, rank, selected, reason, retrieval status |
| `top50.tsv` | the selected set |
| `no_pdb_hit.tsv` | proteins with no qualifying hit — the input to the Foldseek step (issue #9) |
| `run.json` | timestamps, parameters, weights, endpoints, counts, cache path, and every failed query with its error |

## Weight change log

| Date | Change | Reason |
| --- | --- | --- |
| 2026-09-16 | Initial weights, before any full run. | Written from the component list in `02-m2-triage.md` and what the HS11286 specialty export actually provides. |
| 2026-09-16 | `pdb_evidence` normalized by 1.20 instead of clipped at 1.0. Weights unchanged. | A 10-protein smoke run put two different-quality hits (raw 1.09 and 1.043) at an identical 0.40, because clipping discards the bonus range exactly where ranking matters. Found before any full run. |
| 2026-09-17 | Added `ligandable_homolog` (+0.10) — the last component on #12's list. | Fires on a ligand observed bound in the best hit's entry, or a ChEMBL target on the mapped accession. Weight chosen by measuring the trade against `annotation_gap`: at 0.10, holo templates in the top 50 go 14 → 27 while uncharacterized proteins stay at 25 of 50; at 0.15 the list turns over by half and hypotheticals drop to 15. Decision: project lead, 2026-09-17. Note the 0.033 overlap with `pdb_evidence`'s `LIGAND_BONUS`, documented above. |
| 2026-09-17 | Added `surface_bonus` (+0.10) and `membrane_penalty` (−0.15), fed by #10's localization flags. | #10 shipped the flags and deliberately left them unscored, noting the weights belonged here. Surface-exposed proteins are the accessible ones (and M5's vaccine/antibody angle); membrane proteins fold and dock badly, so pitfall #3 says deprioritize rather than exclude — hence a penalty a strong protein can still outweigh. On G37: 19 proteins surface-exposed, 91 membrane, 1 of each in the top 50. |
| 2026-09-17 | `virulence_amr` now reads the AMR classification instead of counting any AMR row. | Issue #30: most "Antibiotic Resistance" rows are "antibiotic target in susceptible species" — drug targets, not resistance genes. On HS11286, 16 of 130 such proteins were being scored as resistance evidence. Those now score `drug_target` instead, unclassified rows get 0.5, and `amr_basis` records which applied. Weight unchanged at 0.20. |
| 2026-09-17 | `essential` now comes from orthology to FBA-essential proteins in public relatives, not BV-BRC's precomputed rows. Weight unchanged at 0.15. | Same reason as the human-homolog change: those rows exist for public genomes only (issue #29). Recovers 148 of G37's 148 known essential genes, plus 21 more that each name their source relative and identity. |
| 2026-09-17 | `human_homolog_penalty` now comes from our own DIAMOND search, not BV-BRC's precomputed rows. Weight unchanged at −0.25. | Those rows exist for public genomes only, so on a blinded genome the penalty never fired (issue #29). The source change is larger than it sounds: 153 G37 proteins now carry a penalty where the public record listed 5. Ranking is barely affected — between penalising everything and penalising nothing, at most 5 of the top 50 change — but the recorded evidence is now ours and reproducible. |
| 2026-09-16 | `annotation_gap` 0.10 → 0.30. **Changed after seeing the first full run** (pitfall #12 — recorded here rather than left implicit). | The HS11286 run selected 0 uncharacterized proteins out of 1,338 (best rank 57), because a known target collects 0.50 from virulence/essential/drug-target while a hypothetical can earn 0.10. That contradicts the project's purpose. Sensitivity over the recorded components: 0.20 → 6/50, 0.25 → 9/50, **0.30 → 10/50**, 0.35 → 15/50, 0.40 → 26/50. 0.30 admits 10, every one with a PDB hit and 8 with a ligand-bound homolog, while keeping 40 characterized targets as the validation set. Decision: project lead, 2026-09-16. |
| 2026-09-18 | `annotation_gap` now graded from functional evidence instead of the product string. 1.0 no terms, 0.5 a domain match only, 0.0 an EC/KO/COG/GO assignment. Weight unchanged at 0.30. | `is_uncharacterized` reads the BV-BRC product, written before any provider runs, so the term never saw the annotation layer's output. On the lambda0 G37 run of 2026-09-18, 91 of the 178 proteins collecting a full gap had InterProScan terms, as did 17 of the 34 uncharacterized proteins in the top 50 — the term was paying full price for proteins whose function we had already recovered. Grading uses evidence we generate for any genome and needs no knowledge of the organism. |
| 2026-09-18 | `ligandable_homolog` withheld when the best hit's entry contains partner entities. Weight unchanged at 0.10. | The component fired on `top.ligands`, which is entry-level. A ribosomal protein's best hit is a ribosome, and ribosome entries routinely carry antibiotics that survive the additive/glycan/metal filter, so every subunit inherited credit for a ligand bound to the rRNA or another chain. Measured on the 2026-09-18 lambda0 run: enabling the component correctly lifted DNA gyrase A 34 → 9, ribose-5-phosphate isomerase 21 → 3 and phosphoglycerate mutase 32 → 7, but also put L33p and L16p at ranks 13-14 on borrowed ligands. Attributing a ligand to a chain needs per-entity contacts we do not fetch, so the credit is withheld and `ligandable_basis` records why. |
| 2026-09-18 | Human homolog at or above 40% identity now **disqualifies** a protein from selection outright, rather than only penalising it. `human_homolog_penalty` stays at −0.25 for ranking. | Selectivity is a requirement, not a preference. As a penalty it could be outvoted: on the 2026-09-18 lambda0 run `essential` (+0.15) beat a 50.6%-identity penalty (−0.127) and returned SSU ribosomal protein S12p to the top 50, whose human mitochondrial counterpart MRPS12 is an established aminoglycoside liability. Disqualified proteins keep their rank, score and reason so the report can state what was removed and why, and the slot passes to the next eligible protein so the caller still gets `top_n`. Count recorded as `counts.disqualified_human_homolog`. |

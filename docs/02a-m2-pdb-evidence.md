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
| `human_homolog` | **−0.25 × identity/100** | applied when a `Human Homolog` specialty row exists (identity is always populated in the export) |

```
triage_score = 0.40·pdb_evidence + 0.20·virulence_amr + 0.15·essential
             + 0.15·drug_target + 0.30·annotation_gap − 0.25·(human_identity/100)
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
| 2026-09-16 | `annotation_gap` 0.10 → 0.30. **Changed after seeing the first full run** (pitfall #12 — recorded here rather than left implicit). | The HS11286 run selected 0 uncharacterized proteins out of 1,338 (best rank 57), because a known target collects 0.50 from virulence/essential/drug-target while a hypothetical can earn 0.10. That contradicts the project's purpose. Sensitivity over the recorded components: 0.20 → 6/50, 0.25 → 9/50, **0.30 → 10/50**, 0.35 → 15/50, 0.40 → 26/50. 0.30 admits 10, every one with a PDB hit and 8 with a ligand-bound homolog, while keeping 40 characterized targets as the validation set. Decision: project lead, 2026-09-16. |

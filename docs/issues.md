# Task issues

Source of truth for `scripts/create_issues.py`. One `##` heading per issue; the `labels:` line is consumed by the script, everything after it is the issue body. Issues are intentionally unassigned — claim one by commenting on it.

## common: report.json schema, fixtures, and section writer
labels: module:common, priority:p0, type:code
**Goal:** the data contract everyone else builds against. Blocking for all other work.

**Docs:** `docs/00-architecture.md`

**Scope**
- JSON Schema per `report.json` section, plus `validate(section, obj)`.
- `common/io.py`: `update_section(run_dir, key, obj)` — read, replace one key, temp file, atomic rename.
- `fixtures/report.fixture.json`: a small but complete example with 5 proteins, 1 structure, 3 ligands, 2 docking rows, 1 disease claim.
- Run-directory layout and `run` manifest fields (versions, seeds, caps).

**Definition of done**
- [ ] `validate()` accepts the fixture and rejects a mangled copy
- [ ] Two modules can write different sections concurrently without clobbering
- [ ] Fixture committed and referenced by every other module's `--dry-run`

**Watch out for:** schema drift. Changes bump `schema_version` and update fixtures in the same commit.

## common: identifier mapping helper
labels: module:common, priority:p0, type:code
**Goal:** one place that maps BV-BRC `feature_id` ↔ locus_tag ↔ UniProt ↔ PDB ↔ ChEMBL target ↔ gene name.

**Docs:** `docs/00-architecture.md`, `docs/pitfalls.md` (#11)

**Scope**
- UniProt ID-mapping + REST lookups, on-disk cache, batch submission.
- Stable `xrefs` dict written onto each protein record.
- Clear behaviour when no mapping exists (explicit `null`, never a guess).

**Definition of done**
- [ ] Round-trips a known genome's feature IDs to UniProt accessions
- [ ] Cached: a second run makes no network calls
- [ ] No other module contains its own mapping code

**Watch out for:** this is where most cross-module bugs will come from. Prefer explicit failure over a plausible-looking wrong mapping.

## common: cached HTTP client and provenance stamping
labels: module:common, priority:p0, type:code
**Goal:** every external call goes through one client that caches, retries, rate-limits and records provenance.

**Docs:** `docs/00-architecture.md`

**Scope**
- `common/http.py`: on-disk cache under `runs/<id>/cache/`, backoff, polite rate limits, project user agent.
- `common/provenance.py`: stamp `source`, `source_version` or release date, `retrieved_at`, `method`, `params` onto derived records.
- Helper for recording tool versions and seeds into the `run` manifest.

**Definition of done**
- [ ] Re-running any module offline works from cache
- [ ] Every record produced in the repo carries provenance fields
- [ ] Rate-limit and 5xx handling tested against one real API

## M1: Similar Genome Finder and CGA submit/poll/parse
labels: module:m1, priority:p0, type:code
**Goal:** assembly in, `genome` + `proteins` + `proteins.faa` out.

**Docs:** `docs/01-m1-genome.md`

**Scope**
- Similar Genome Finder first (predicted taxon, closest genomes) — CGA annotation needs a genus-level-or-better taxon.
- `p3-submit-CGA`, poll, retrieve; record job ID and wall-clock duration.
- Parse CDS, products, subsystems, specialty genes (with database, identity, coverage), tree Newick.

**Definition of done**
- [ ] Valid `genome` and `proteins` sections from a real run
- [ ] Specialty hits carry database, identity and coverage — no bare gene names
- [ ] Sequences in `proteins.faa`, not in JSON
- [ ] CGA runtime recorded in the issue as a comment for the team

## M1: p3-CLI fallback annotation path
labels: module:m1, priority:p0, type:code, good-first-task
**Goal:** produce the same `proteins[]` shape from BV-BRC's existing annotation for a known genome, so nobody waits on a CGA job. Doubles as the fixture source.

**Docs:** `docs/01-m1-genome.md`, `docs/pitfalls.md` (#15)

**Scope**
- `p3-get-genome-features` plus specialty-gene tables for a chosen genome ID.
- Emit `proteins[]` identical in shape to the CGA parser's output.
- Flag in the `run` manifest which path produced the data.

**Definition of done**
- [ ] Output validates against the same schema as the CGA path
- [ ] Runs in under a couple of minutes
- [ ] Used to generate `fixtures/report.fixture.json`

## M1: ANI and SNP distance to closest genomes
labels: module:m1, priority:p1, type:code, good-first-task
**Goal:** numeric distances for the closest ~10 genomes.

**Docs:** `docs/01-m1-genome.md`

**Scope**
- skani (or fastANI) for ANI; Snippy `--ctgs` or Parsnp for SNP distance.
- Download reference genomes via the `p3-` CLI; cache them.
- Write into `genome.closest_genomes[]`.

**Definition of done**
- [ ] ANI and SNP distance populated per closest genome
- [ ] Self-match (the blinded genome matching itself) reported but not the only row
- [ ] Runs on a laptop

## M2: sequence and structure lookups
labels: module:m2, priority:p0, type:code
**Goal:** populate `xrefs` and `annotations` from PDB, AlphaFold DB and UniProt.

**Docs:** `docs/02-m2-triage.md`

**Scope**
- RCSB Search API sequence search; RCSB/PDBe Data API for entry details, including whether a ligand is bound.
- AlphaFold DB by UniProt accession.
- Record identity, coverage and e-value on every hit.

**Definition of done**
- [ ] Every protein has `xrefs` populated or explicitly empty
- [ ] Holo (ligand-bound) homolog structures are flagged for M3
- [ ] `--dry-run` works from fixtures

## M2: Foldseek structure search (and ESM Atlas probe)
labels: module:m2, priority:p1, type:code
**Goal:** find relatives by shape for proteins whose sequence matches nothing.

**Docs:** `docs/02-m2-triage.md`, `docs/pitfalls.md` (#14)

**Scope**
- Foldseek local binary against PDB + AlphaFold DB; server REST API as fallback.
- Record TM-score, coverage and the hit's annotation.
- Day-1 probe of the ESM Metagenome Atlas endpoints: availability, length caps, rate limits. Report findings in a comment; treat as bonus.

**Definition of done**
- [ ] Top hits with TM-scores written into `annotations`
- [ ] Hypothetical proteins get structural neighbours where they exist
- [ ] ESM Atlas verdict documented (usable / not usable, with limits)

## M2: functional annotation layer
labels: module:m2, priority:p1, type:code
**Goal:** function, localization and membrane flags for every protein.

**Docs:** `docs/02-m2-triage.md`

**Scope**
- eggNOG-mapper (COG, GO, KEGG KO, EC); InterProScan if time allows.
- PSORTb and SignalP 6 for surface/secreted; DeepTMHMM for the membrane flag.
- Write into `annotations[]` and `flags`.

**Definition of done**
- [ ] `flags.secreted` and `flags.membrane` set for every protein
- [ ] Runtime for a full proteome measured and recorded

## M2: knowledge graph assembly
labels: module:m2, priority:p1, type:code
**Goal:** a small provenance-carrying subgraph linking our proteins to known disease/pathogen proteins, pathways and ligands.

**Docs:** `docs/02-m2-triage.md`, `docs/pitfalls.md` (#4, #10)

**Scope**
- Open Targets GraphQL, STRING, DGIdb, HPIDB, PHI-base, IntAct/BioGRID, Reactome/KEGG; NCATS Translator optional.
- Edges typed explicitly (`homolog_of`, `target_of`, `interacts_with`, …) with the homology identity that justified them.
- NetworkX assembly; write `kg` (or `kg.json` plus a summary if it grows large).

**Definition of done**
- [ ] Every edge carries source, retrieval date and supporting identity
- [ ] Human-homolog edges record identity so M2 scoring can penalize them
- [ ] No attempt to ingest or host a full KG

## M2: triage scoring and selection
labels: module:m2, priority:p0, type:code
**Goal:** rank all proteins and select ~50, reproducibly.

**Docs:** `docs/02-m2-triage.md`, `docs/pitfalls.md` (#12)

**Scope**
- Weighted score with written-down components: virulence/AMR hit, essentiality ortholog, secreted/surface, ligandable homolog, structure availability; penalties for close human homolog and predicted membrane protein.
- Write `triage.score`, `components`, `rank`, `selected`, `reason` for every protein — including the discarded ones.

**Definition of done**
- [ ] Selection reproducible from recorded components alone
- [ ] Weights committed before the first real run; later changes noted in the commit message
- [ ] Sanity check: known virulence/AMR proteins rank high; a conserved housekeeping protein with a close human homolog does not

## M3: structure collection, prediction, confidence gate
labels: module:m3, priority:p0, type:code
**Goal:** a usable structure for every selected protein, or an explicit "not dockable".

**Docs:** `docs/03-m3-fold.md`

**Scope**
- Reuse PDB/AFDB hits from M2 first; prefer experimental, prefer holo.
- Predict the remainder via the BV-BRC Protein Structure Prediction service (or local ESMFold); record model, version, seed.
- Confidence gate on pLDDT/PAE; trim or mask low-confidence regions; set `usable_for_docking` with a reason.

**Definition of done**
- [ ] 50 proteins processed in one working session, with per-protein wall-clock logged
- [ ] Nothing re-folded that AlphaFold DB already has
- [ ] `structures` section validates

## M3: pockets and receptor preparation
labels: module:m3, priority:p0, type:code
**Goal:** at least one defined pocket per dockable protein, plus prepared receptors for docking.

**Docs:** `docs/03-m3-fold.md`, `docs/pitfalls.md` (#1, #2)

**Scope**
- P2Rank (and/or fpocket); superpose a holo template with US-align/TM-align to transfer its known ligand site and mark `aligned_to_known_site`.
- PDBFixer, PDB2PQR, Meeko/OpenBabel to produce docking-ready receptors.
- Flag pockets whose homolog site requires a cofactor the model lacks.

**Definition of done**
- [ ] Pocket residue lists (not just scores) written for M4
- [ ] Cofactor-dependent sites flagged
- [ ] Pockets in low-confidence regions rejected

## M4: ligand candidate collection
labels: module:m4, priority:p1, type:code
**Goal:** per-protein candidate ligand sets with provenance.

**Docs:** `docs/04-m4-ligands.md`

**Scope**
- ChEMBL (`max_phase=4` for approved), openFDA labels/approvals, Broad Drug Repurposing Hub, ZINC in-stock subsets; Open Targets and DGIdb for KG-derived ligands; BindingDB and PubChem BioAssay for activity data.
- RAGStack literature pass for candidates not in the databases.
- RDKit canonicalization and InChIKey dedup; record `targets_via` for every ligand.

**Definition of done**
- [ ] Every ligand traceable to source, retrieval date and why it was proposed for that protein
- [ ] DrugBank used interactively only, never shipped
- [ ] Counts per source reported

## M4: ligand scoring and caps
labels: module:m4, priority:p1, type:code
**Goal:** a transparent per-protein ranking that decides what gets docked.

**Docs:** `docs/04-m4-ligands.md`, `docs/pitfalls.md` (#12, #16)

**Scope**
- Weighted components: homology identity to the ligand's known target, measured potency (pChEMBL), approval status, property filters and PAINS flags, literature support, KG edge type. No learned model this week.
- Caps recorded in the `run` manifest; start at 200 candidates in, top 20 docked.

**Definition of done**
- [ ] `score` plus `components` on every ligand row
- [ ] PAINS and promiscuous compounds filtered before docking
- [ ] Weights committed before the first real run

## M4: docking with controls
labels: module:m4, priority:p1, type:code
**Goal:** docked poses for the top ligands per protein, with controls that make the ranking interpretable.

**Docs:** `docs/04-m4-ligands.md`, `docs/pitfalls.md` (#1, #5, #13)

**Scope**
- AutoDock Vina (or smina) into M3 pockets only — no blind docking. Fixed seed and exhaustiveness, recorded.
- Per protein: one known binder (positive) and ~20 decoys or random ligands (negative).
- gnina rescoring of top poses if time allows. Save poses under `runs/<id>/docking/`.
- Benchmark one protein before any batch run and report the timing.

**Definition of done**
- [ ] Controls present for every docked protein
- [ ] Output structurally prevents cross-protein score comparison (scores carry protein context)
- [ ] Same seed reproduces the same scores

## M4: Boltz-2 co-folding refinement for top pairs
labels: module:m4, priority:p2, type:code
**Goal:** a second, independent read on the top handful of protein–ligand pairs.

**Docs:** `docs/04-m4-ligands.md`, `docs/pipeline.md`

**Scope**
- BV-BRC Protein Structure Prediction with a ligand (CCD code or SMILES) for the top ~10 pairs.
- Store the affinity estimate as an `analyses` entry, explicitly not as truth.
- Compare rank order against Vina; note disagreements.

**Definition of done**
- [ ] Top pairs have a co-folded model and an affinity estimate
- [ ] Runtime per pair recorded
- [ ] Disagreements with Vina surfaced in the report

## M5: deterministic evidence gathering
labels: module:m5, priority:p1, type:code
**Goal:** the factual base for disease analysis, with no LLM involved.

**Docs:** `docs/05-m5-disease.md`

**Scope**
- Virulence and AMR from M1 specialty genes; PHI-base phenotypes; secretion systems (TXSScan).
- BacDive and MediaDive APIs for culture conditions of the closest relatives; GapMind for amino-acid/vitamin gaps; subsystems for pathway completeness.
- Write `disease.growth` with source type per field (strain record vs. genomic inference).

**Definition of done**
- [ ] Every growth or nutrient value names its source and source type
- [ ] Re-runs are identical
- [ ] Clear distinction between "measured in a relative" and "inferred from our genome"

## M5: RAG synthesis and claim verification
labels: module:m5, priority:p1, type:code
**Goal:** a mechanism narrative and mitigation strategies where every claim cites evidence.

**Docs:** `docs/05-m5-disease.md`, `docs/pitfalls.md` (#7, #21)

**Scope**
- RAGStack `/v1/retrieve` per gene and organism; cache passages with PMCIDs.
- Argo synthesis pass (haiku while prototyping, opus for final) with tools wrapping RAGStack, ChEMBL, UniProt, Open Targets. Non-streaming path for tool calling; explicit `max_tokens`; hard iteration cap.
- Separate verification pass that checks each claim against its cited evidence and drops or downgrades unsupported ones; keep the dropped list.
- Mitigation strategies, each tied to a specific finding.

**Definition of done**
- [ ] Zero claims without feature IDs or PMCIDs in the final section
- [ ] Dropped-claim count reported
- [ ] Model, prompt version and iteration cap recorded in the run manifest

## M6: report renderer
labels: module:m6, priority:p0, type:code
**Goal:** one self-contained HTML report — the deliverable the judges see.

**Docs:** `docs/06-m6-report.md`, `docs/pitfalls.md` (#5, #9)

**Scope**
- Jinja2 templates for the seven sections; graceful "not run" per missing section.
- 3Dmol.js receptor/pocket/pose views; RDKit ligand SVGs; tree rendering.
- Candidates grouped by protein, never ranked across proteins. Score breakdown and provenance per row.
- Fixed limitations block: docking scores are rankings not affinities, predicted structures lack cofactors, no experimental validation, no clinical implication.

**Definition of done**
- [ ] Renders from `fixtures/report.fixture.json` on day 1, no network, no other module implemented
- [ ] Renders from a partial `report.json`
- [ ] Nothing computed in templates
- [ ] Opens offline in a browser

## CI: schema validation and tests
labels: module:common, priority:p1, type:ops, good-first-task
**Goal:** stop schema drift and broken fixtures from reaching other threads.

**Scope**
- GitHub Actions: run `validate()` against fixtures, run unit tests, lint.
- Fail the build when a module's output does not validate.
- Optional: check that every module's `--dry-run` completes.

**Definition of done**
- [ ] CI runs on pull requests and on pushes to `main`
- [ ] A deliberately broken fixture fails the build

## Ops: day-1 access smoke tests
labels: priority:p0, type:ops, good-first-task
**Goal:** find the broken credential or dead endpoint before it costs someone an afternoon.

**Scope**
- BV-BRC: `p3-login`, a trivial service submission, CLI queries.
- RAGStack: `GET /v1/stats/tenants?counts=false` with both auth routes (API key and `~/.patric_token`).
- Argo: `examples/00-smoke-test.sh`, a tool-calling round trip on the non-streaming path.
- Public APIs: UniProt, RCSB, AlphaFold DB, ChEMBL, openFDA, Open Targets, STRING, BacDive, MediaDive — one call each, note registration requirements.
- ESM Atlas: availability, length cap, rate limit.

**Definition of done**
- [ ] A table of endpoint → works / needs key / broken, posted as a comment
- [ ] Anything broken raised as its own issue

## Decision: test genome, triage count, ligand caps
labels: priority:p0, type:decision
**Goal:** close the open parameters before the pipeline runs on real data.

**Docs:** `docs/07-decisions-and-risks.md`

**To decide**
- Test genome: *S. aureus* USA300 (default), *K. pneumoniae* HS11286, or *P. aeruginosa* PAO1 — blinded either way.
- Triage target: 20 / 50 (default) / 100 proteins, with the pipeline proven on 10 first.
- Ligand caps: candidates in per protein (default 200) and number docked (default 20).
- Docking engine and rescoring: Vina + gnina (default).
- Argo models: haiku for prototyping, opus for final (default).

**Definition of done**
- [ ] Decisions recorded in this issue and reflected in `docs/07-decisions-and-risks.md`
- [ ] Defaults taken silently if nobody objects by the first checkpoint

## M2: essentiality and human-homolog evidence
labels: module:m2, priority:p0, type:code
**Goal:** produce the two triage components that no upstream service gives us.

**Docs:** `docs/01a-cga-coverage.md`, `docs/02-m2-triage.md`, `docs/pitfalls.md` (#3)

A fresh CGA run returns 17 specialty-gene rows. The public BV-BRC record for the
same genome has 174, and the extra 157 are mostly `Essential Gene` (148, from
flux-balance analysis) and `Human Homolog` (5, from BLASTP). Those are
precomputed for public genomes only, so a blinded genome gets neither — yet
#12's score has an essentiality component and a human-homolog penalty.

**Scope**
- Human homology: BLASTP or DIAMOND of every protein against the human
  proteome (UniProt UP000005640). Record identity, coverage and e-value per hit.
- Essentiality: orthologs against an essential-gene set (DEG, OGEE, or the
  FBA essential genes on close BV-BRC public relatives). Record which relative
  and what identity justified the call.
- Write both into `annotations[]` with provenance, in the shape `triage`
  expects. Explicit null when there is no hit; never a silent zero.

**Definition of done**
- [ ] Every protein has a human-homolog record or an explicit "no hit"
- [ ] Essentiality calls name their source and the identity behind them
- [ ] Spot check against the public record for 243273.25: the 5 known human
      homologs (atpD, tuf, atpA, rpsL, dnaK) come back at comparable identity
- [ ] Runtime for a full proteome measured and recorded

**Watch out for:** the human-homolog paradox (pitfall #3) - a close human
homolog both helps find ligands and disqualifies the target. Record the
identity, do not collapse it to a boolean.

## M1: normalize specialty-gene evidence, and run AMRFinderPlus and RGI locally
labels: module:m1, priority:p1, type:code
**Goal:** make every specialty-gene row carry evidence a reader can judge, and
stop AMR calls from being read as resistance when they are not.

**Docs:** `docs/01a-cga-coverage.md`, `docs/01-m1-genome.md`

Two findings from the first CGA run. First, `K-mer Search` rows have null
`identity`, `query_coverage`, `subject_coverage` and `e_value`; their confidence
lives per feature in `annotation.genome` as `quality.hit_count`,
`quality.weighted_hit_count` and `quality.priority`. Only `DIAMOND` rows carry
identity and coverage. Second, 10 of the 14 AMR rows are classified "antibiotic
target in susceptible species", 2 alter cell wall charge, 1 confers resistance
via absence - none is an acquired resistance gene.

**Scope**
- Carry `classification` and `evidence` onto every specialty record, not just
  `property`, and expose them to M2 scoring and the M6 report.
- For k-mer rows, carry the per-feature k-mer confidence instead of faking an
  identity. A null identity stays null.
- Run AMRFinderPlus and RGI locally on the proteins, since CGA's copies returned
  header-only files, and merge their hits with identity and coverage.
- Decide and document what "AMR hit" means for the triage score, given that most
  rows are conserved drug targets.

**Definition of done**
- [ ] No specialty record loses its classification on the way into `report.json`
- [ ] Every record carries either identity plus coverage, or a k-mer confidence,
      and which one it is is explicit
- [ ] AMRFinderPlus and RGI produce non-empty output for at least one genome
- [ ] #5's "database, identity and coverage" check updated to match reality

**Watch out for:** scoring gyrA and rpoB as resistance evidence. In this genome
they are simply the drug targets of fluoroquinolones and rifamycins.

## Ops: find the Similar Genome Finder access path
labels: module:m1, priority:p0, type:ops, good-first-task
**Goal:** settle how we call Similar Genome Finder, since CGA does not give us
closest genomes.

**Docs:** `docs/01a-cga-coverage.md`, `docs/01-m1-genome.md`

`close_genomes` in the CGA output is an empty list. The codon tree names 10
relative genome IDs with branch lengths, but no Mash distance and no identity.
M1 needs a predicted taxon before CGA can annotate a blinded genome, and
`genome.closest_genomes[]` needs real distances. No `p3-` command for the
service has been identified; `p3-closest-seqs` exists but has not been shown to
be the same Mash service.

**Scope**
- Establish the access path: a `p3-` command, the service API, or the website.
- Confirm what it returns (genome IDs, Mash distance, predicted taxon) and at
  what taxonomic resolution.
- Write the finding into `docs/01-m1-genome.md` and note it in this issue.

**Definition of done**
- [ ] One documented, repeatable call that returns closest genomes for a FASTA
- [ ] Output fields mapped onto `genome.closest_genomes[]`
- [ ] If no programmatic path exists, that is recorded and the manual fallback
      is written down

## M1: phylogeny gate on the codon tree
labels: module:m1, priority:p1, type:code
**Goal:** refuse to publish a tree built from too few genes.

**Docs:** `docs/01a-cga-coverage.md`, `docs/06-m6-report.md`

`codontree.genesPerGenome.txt` from the first run reports
`Filtered_SingleCopy = 5` for every genome in the ingroup: the tree came from
five single-copy genes. That is fine for a smoke test and not defensible in the
report.

**Scope**
- Parse `codontree.genesPerGenome.txt` and `codontree.homologAlignmentStats.txt`
  and record the gene count in the `genome` section.
- Set a minimum gene count below which the tree is marked not publishable, and
  say so in the report rather than hiding it.
- If the count is too low, fall back to the BV-BRC Phylogenetic Tree service
  with a wider ingroup.

**Definition of done**
- [ ] Gene count and alignment stats recorded per run
- [ ] A tree under the threshold is flagged in `report.json` and in M6
- [ ] Threshold written down with the reason for the number chosen

## Ops: verify CGA on the real test genome before relying on virulence and AMR phenotype
labels: module:m1, priority:p0, type:ops
**Goal:** find out which CGA analyses only work for well-characterized pathogens,
before the pipeline depends on them.

**Docs:** `docs/01a-cga-coverage.md`, `docs/pipeline.md`

Several analyses ran on *M. genitalium* G37 and produced nothing:
`genome_amr.json` is `[]` despite two `models.spcAb.*` AMR classification models
executing, cgMLST called 0 percent of loci, and the AMRFinderPlus and RGI tables
are header-only. Zero virulence rows came back, and `specialty-blast.txt` shows
only TCDB and TTD were searched, so we do not know whether VFDB and Victors run
at all. `pipeline.md` assumes virulence factors come from CGA specialty genes,
and M5 assumes AMR phenotype data exists.

**Scope**
- Run CGA on *S. aureus* USA300, the decided test genome, blinded.
- Record which of these populate: virulence rows (VFDB, Victors), `genome_amr`
  MIC and SIR predictions, MLST and cgMLST, AMRFinderPlus, RGI.
- Record the runtime for a 2.8 Mb genome against the 164.5 s for 580 kb.
- Post the results as a comment and update `docs/pipeline.md` where it is wrong.

**Definition of done**
- [ ] A table of analysis to populated/empty for USA300, posted in this issue
- [ ] `pipeline.md` corrected wherever it assumes data CGA does not produce
- [ ] M5 told whether AMR phenotype data will exist for the real genome

## M1: define behaviour when Similar Genome Finder returns no usable hit
labels: module:m1, priority:p2, type:code
**Goal:** decide what M1 does when Minhash finds nothing close enough to give a
taxon, instead of falling through to a domain-level guess.

**Docs:** `docs/01a-cga-coverage.md`, `docs/01-m1-genome.md`

Deferred on purpose: the current test genome always resolves, so this does not
block the pipeline. It matters for a genuinely novel organism.

Our test resolved only because G37 is itself a BV-BRC reference genome: with
scope limited to reference and representative sketches and `max_distance` 0.2,
Minhash returned exactly one hit, 243273.25 at distance 0 with 1000/1000 shared
k-mers. A novel organism in that configuration returns an empty list, and the
floor run shows what annotating without a taxon costs: CDS count doubles, mean
protein length falls from 347 to 142 aa, quality drops to Poor.

**Scope**
- Retry policy when the list is empty: widen `max_distance`, then drop to
  full-scope search (`include_reference` and `include_representative` both 0).
- Decide the fallback taxon and genetic code when still empty, and record that
  the taxon was a guess.
- Refuse to proceed, or mark the whole run low confidence, rather than emitting
  a Poor-quality annotation that looks normal downstream.

**Definition of done**
- [ ] Documented ladder of retries and the final fallback
- [ ] `run` manifest records predicted taxon, its Mash distance and whether a
      fallback was used
- [ ] A synthetic no-hit case is tested (e.g. an absurd `max_distance`)

**Watch out for:** an empty result and a self-match look equally "successful" to
a caller that only checks for errors.

## M2: annotation_gap is priced against the product string, not against the evidence

labels: module:m2, priority:p1, type:decision

**Goal:** decide what `annotation_gap` should be worth now that the annotation
layer can tell us whether the gap is real.

**Docs:** `docs/02a-m2-pdb-evidence.md` (weight change log), `docs/02c-m2-functional-annotation.md`

`score.py` sets `annotation_gap=1.0 if protein.is_uncharacterized else 0.0`, and
`is_uncharacterized` is a regex over the BV-BRC product string, written before any
provider runs. The term therefore cannot know that InterProScan has since named
the protein.

Measured on the lambda0 G37 run, 2026-09-18 (530 proteins, InterProScan `-appl Pfam`):

- 178 proteins collect the full 0.30.
- **91 of those 178 already have InterProScan terms** — 21 with an EC/KO/COG/GO
  assignment, 70 with a domain match and no functional assignment.
- 34 of the top 50 are `hypothetical protein`; 17 of those 34 have terms.

At 0.30, `annotation_gap` is the second-largest positive weight, equal to
`essential` + `drug_target` combined, so this is not a rounding error: proteins
whose function we have already recovered are being paid in full for not having
it.

A fix was written and reverted on `m2/triage-scoring-fixes`. Both shapes break
`test_m2_function.test_annotation_now_feeds_the_triage_score`, which reconstructs
the whole triage score from `surface_bonus` and `membrane_penalty` and so forbids
*any* annotation-driven component from carrying weight. That invariant is worth
keeping or worth changing, but it is a decision, not a patch.

**Scope**
- Decide whether annotation may reach the triage score through more than two
  components, and update the invariant deliberately if so.
- If yes: either grade `annotation_gap` from functional terms, or add an
  offsetting component. A graded definition that worked in testing was 1.0 no
  terms, 0.5 a domain match only, 0.0 an EC/KO/COG/GO assignment.
- Re-price 0.30 either way. The sweep in the change log (0.20 to 6/50, 0.30 to
  10/50, 0.40 to 26/50) was run on HS11286 against the *old*, looser definition,
  which fired on roughly twice as many proteins. Same number, different quantity.

**Definition of done**
- [ ] The invariant is either reaffirmed in the test docstring or changed on
      purpose, with the decision recorded
- [ ] `annotation_gap` reflects whether function is actually unknown
- [ ] A dated row in the weight change log, per pitfall #12
- [ ] Top-50 composition before and after recorded on the same genome

**Watch out for:** re-pricing the weight and redefining the term in one commit
makes the effect of neither measurable.

## M2: essentiality reference set can include the query's own species

labels: module:m2, priority:p1, type:code

**Goal:** record, per run, how much of the transferred essentiality came from
relatives that are the same species as the query.

**Docs:** `docs/02a-m2-pdb-evidence.md`, `docs/09-pitfalls.md`

`resolve_reference_set` now selects FBA-essential relatives by BV-BRC genome id
from M1's ingroup, which is right: names drift between databases and ids do not.
Measured 2026-09-18 against the `sp_gene` core — `eq(taxon_id,2097)` 0 rows,
`eq(genome_id,2097.118)` 0 rows, every name in this genome's current NCBI lineage
0 rows (BV-BRC still indexes the superseded genus *Mycoplasma*), and
`in(genome_id,(243273.25,272634.6))` 299 rows.

But `243273.25` **is G37**. On the last run, 158 essential calls were transferred
from 1,707 reference proteins across 10 genomes, and at least one of those genomes
is the query's own species. On a novel organism there would be no such relative,
so the transfer would be strictly harder than the number suggests.

The PDB path already solves exactly this: `same_organism_hits` records
`same_species_in_top_n`, `same_species_total` and a note saying weights calibrated
on a genome with its own structures will not carry to a novel one. Essentiality
needs the same treatment.

**Scope**
- Record how many reference genomes and how many transferred calls are
  same-species and same-genus as the query, using `QueryTaxonomy.same_species_as`
  / `same_genus_as`.
- Surface it in `run.json` next to `essentiality`, and make it available to M6.
- Decide whether same-species relatives should be excluded, kept and flagged, or
  reported both ways. Flagging is the smaller change and loses nothing.

**Definition of done**
- [ ] `essentiality` in `run.json` carries same-species and same-genus counts
- [ ] A run on a genome with a same-species relative and one without are
      distinguishable from the manifest alone

**Watch out for:** excluding same-species relatives silently would make the G37
numbers drop with no visible explanation.

## M2: counts.no_pdb_hit does not count a failed search as a miss

labels: module:m2, priority:p2, type:code

**Goal:** stop the run summary reporting zero misses when every search failed.

**Docs:** `docs/02a-m2-pdb-evidence.md`

Offline with a cold cache, all 530 sequence searches fail and stdout says:

```
530 proteins scored from 530 unique sequences; 50 selected, 0 without a PDB hit.
Warning: 530 sequence searches failed (recorded as query-failed in run.json).
```

`counts.no_pdb_hit` is 0 when the true figure is 530. `query-failed` lands in
neither the hit bucket nor the miss bucket, so the two lines contradict each
other. The per-protein records are honest — every one carries
`triage.retrieval_status: "query-failed"` and a matching reason — and live runs
reconcile exactly (193 = 163 `no-hit` + 30 below the retention threshold, and
530 - 337 = 193). So this is a defect in the summary line and one counter, not
in the data.

**Scope**
- Count `query-failed` separately and print it, or fold it into `no_pdb_hit`
  with the failure count stated alongside.
- Make the stdout summary agree with `counts` in every case.

**Definition of done**
- [ ] A fully failed offline run reports 530, not 0
- [ ] A test covers the all-failed case

**Watch out for:** `no_pdb_hit.tsv` is built from the same flag and has the same
blind spot.

## Runner: --dry-run reports a report.json it never writes

labels: module:common, priority:p2, type:code

**Goal:** stop `--dry-run` printing a success line for a file that does not exist.

`python -m s2f.run --run runs/x --from-cga-dir data/<cga> --dry-run` prints the
two subcommands it would run, then prints:

```
done. report: runs/x/report.json
```

`find runs/x -type f | wc -l` returns **0**. Verified on lambda0, 2026-09-18.

Cosmetic, but it is a success claim for a nonexistent artefact, and a caller that
checks for the line rather than the file would be misled.

**Scope**
- Print what a dry run actually did, and name no output path it did not write.

**Definition of done**
- [ ] A dry run's output names no file it did not create
- [ ] A test asserts a dry run writes nothing

## M1: the genome section carries no contig count or genome length

labels: module:m1, priority:p2, type:code

**Goal:** carry the genome's size into the contract, since CGA already supplies it.

**Docs:** `docs/00a-data-contract.md`, `docs/01-m1-genome.md`

`report.json`'s `genome` holds `genome_id`, `taxon_id`, `taxonomy`, `quality`,
`annotation_route`, `cga_job_id`, `closest_genomes`, `tree_ingroup` and
`tree_newick` — but neither the contig count nor the genome length.

`data/<cga>/.annotation/load_files/genome.json` supplies both:
`contigs: 1`, `genome_length: 580076` for the G37 test genome. A report about a
genome that never states the genome's size is incomplete, and M6 has nothing to
quote.

**Scope**
- Carry `contigs` and `genome_length` from the CGA load files into `genome`.
- Do the same on the BV-BRC API route so both routes stay symmetric.

**Definition of done**
- [ ] `genome.contigs` and `genome.genome_length` populated on both routes
- [ ] Schema updated and validated

**Watch out for:** the API route reports these under different field names.

## M2: annotation confidence is reported per protein, not per flag

labels: module:m2, priority:p1, type:code

**Goal:** stop a single confidence histogram implying that every flag on every
protein came from a real predictor.

**Docs:** `docs/02c-m2-functional-annotation.md`

After a full provider pass on lambda0 (2026-09-18, DeepTMHMM over all 530 plus
InterProScan), `run.json` reports:

```
providers            {'deeptmhmm': 530, 'interproscan': 435, 'heuristic': 530}
confidence           {'predicted': 530}
membrane_flag_source {'deeptmhmm': 530}
heuristic_only       0
```

`confidence` is computed per protein, from the sources behind the two required
flags. But `lipoprotein` is **not** one of them, and DeepTMHMM does not predict
lipoprotein status at all — it predicts topology. Scored against DeepTMHMM over
530 proteins, `lipoprotein` gives P 1.00, R 1.00, MCC 1.00, which looks like a
perfect heuristic and actually means *the flag never changed because no provider
covers it*. There is a `membrane_flag_source` but no `lipoprotein_source`.

For contrast, the flags that are covered: `membrane` MCC 0.95 (P 0.98, R 0.94),
`signal_peptide` MCC 0.41 (P 0.42, R 0.47), `secreted` MCC 0.36 (P 0.46, R 0.32).
The last two are poor enough that a reader needs to know which tier each flag came
from.

**Scope**
- Give every flag a source, or report confidence per flag rather than per protein.
- Make M6 able to say which flags are predicted and which are still heuristic.

**Definition of done**
- [ ] `lipoprotein` carries a source, or is reported as uncovered
- [ ] A run where one flag is predicted and another is not is distinguishable
      from the manifest alone

**Watch out for:** `{'predicted': 530}` is the number most likely to be quoted in
the paper. It currently overstates what was measured.

## Ops: the lambda runbook describes a different node than the one we run on

labels: module:common, priority:p2, type:ops

**Goal:** make `docs/02d-lambda-runbook.md` match the machine people actually use.

Checked against lambda0 on 2026-09-18. Every row below is wrong in the doc:

| runbook says | lambda0 has |
| --- | --- |
| node lambda13 | lambda0 |
| `scripts/setup_lambda.sh` | `scripts/setup_env.sh` |
| 173 tests | 390 passed, 4 skipped |
| InterProScan at `/nfs/lambda_stor_01/homes/cmann/software/interproscan-5.78-109.0` | `/homes/cmann/software/interproscan-5.78-109.0` |
| `INTERPROSCAN_HOME` exported in `~/.bashrc` | unset; discovery finds it by directory search |

Measured runtimes worth adding, none of which the doc currently has:

| Step | Route | Sequences | Wall |
| --- | --- | --- | --- |
| Pass one (RCSB + UniProt + AFDB) | network, serial | 530 | 23 m 47 s |
| InterProScan `-appl Pfam` | local | 530 | 1 m 5 s |
| DeepTMHMM | BioLib cloud | 530 | 9 m 8 s |
| Pass two ingest | `--offline` | 530 | 1.8 s |

Pass one is ~99.5% network wait (7.5 s user time, ~2.7 s per unique sequence,
serial) and is the pipeline's bottleneck. DeepTMHMM in the cloud is close enough
to the local V100 figure (6 m 45 s) that installing it locally is not worth it
for single-genome work.

**Scope**
- Rewrite against lambda0, or parameterise the node as
  `docs/02d-remote-runbook.md` already does.
- Add the runtime table and note that InterProScan is fast here only because
  *M. genitalium* sits in EBI's pre-calculated lookup.

**Definition of done**
- [ ] Every command in the doc runs as written on lambda0
- [ ] Pass-one runtime recorded, closing that part of #10's definition of done

## Runner: --limit means different things to different stages

labels: module:common, priority:p1, type:code

**Goal:** stop one flag silently truncating a run.

**Docs:** `docs/00-architecture.md`

`s2f/run.py` forwards `--limit` to both M2 and M3, where it means different things: to M2
it caps how many proteins are *scored*, to M3 how many structures are *fetched*. Passing
it to cap M3's downloads therefore caps M2's whole analysis.

Observed on lambda0, 2026-09-18. `--limit 10` intended for structure collection produced:

```
10 proteins scored from 10 unique sequences; 10 selected, 3 without a PDB hit.
```

on a 530-protein proteome. The run completed, wrote a valid `report.json`, and the only
sign was one line in the stage output. A reader of that report would see a ten-protein
triage presented exactly like a complete one.

**Scope**
- Give each stage its own flag — `--limit` keeps its M2 meaning, M3 gets
  `--structure-limit` or similar — or namespace them (`--m2-limit`, `--m3-limit`).
- Whichever is chosen, `run.py` should refuse to pass a flag to a stage that did not ask
  for it, rather than forwarding by coincidence of name.
- Record the caps that were in force in the `run` manifest, so a truncated run is visible
  from the report rather than only from the console.

**Definition of done**
- [ ] Capping structure collection leaves protein scoring untouched
- [ ] `run` records the per-stage caps actually applied
- [ ] A test pins that a cap meant for one stage does not reach another

**Watch out for:** a truncated run that still validates is the dangerous shape here. The
schema cannot tell 10 proteins from 530.

## Runner: --human-homology and --essentiality cannot be reached from the pipeline

labels: module:common, priority:p1, type:code

**Goal:** let a single `python -m s2f.run` invocation produce the same ranking as running
M2 by hand.

**Docs:** `docs/00-architecture.md`, `docs/02a-m2-pdb-evidence.md`

`run.py` forwards `--annotate`, `--map-ids`, `--kg`, `--offline`, the five provider paths,
and now `--organism` and `--taxon`. It does **not** forward `--human-homology`,
`--human-proteome`, `--essentiality` or the four `--essentiality-*` tuning flags.

Those drive two scoring components:

| component | weight | reachable from `s2f.run`? |
| --- | --- | --- |
| `human_homolog_penalty` | −0.25 | no |
| `essential` | +0.15 | no |

That is 0.40 of the weight mass, and `human_homolog_penalty` is the largest negative
weight in the model. A pipeline run therefore produces a ranking with the selectivity
filter switched off and no way to switch it on, and the disqualification added in
`m2/triage-scoring-fixes` never fires. Every full run on 2026-09-18 needed a second,
hand-written `m2_triage` invocation afterwards, and then a re-render of the report.

**Scope**
- Forward the flags, or give the runner a single `--full-evidence` switch that turns on
  the components that need external data and records that it did.
- Make the absence visible: if a weighted component had no provider behind it, the run
  manifest already records `components_available`; the runner should say so on the console
  too, rather than producing a quietly different ranking.

**Definition of done**
- [ ] One `s2f.run` invocation reproduces a hand-run M2 with the same flags
- [ ] A run where a weighted component could not be computed says so at the end
- [ ] A test pins that every M2 flag the runner claims to support actually arrives

**Watch out for:** the two rankings look equally finished. Nothing downstream can tell that
0.40 of the weight mass was absent.

## Report: structure and ligand viewers, once there is something to view

labels: module:m6, priority:p2, type:code

**Goal:** the two items of #66's scope that could not be built while the modules feeding
them did not exist.

**Docs:** `docs/06-m6-report.md`, `docs/pitfalls.md` (#2, #5)

#66 shipped and closed on its definition of done: renders from the fixture, from a partial
report and with no network, nothing computed in the template, opens offline. Two scope
items were deliberately left, because a viewer with nothing to display is a stub rather
than a feature:

- **3Dmol.js receptor / pocket / pose views.** M3 collects existing PDB and AlphaFold
  structures but defines no pockets yet, and poses need M4, which is not written. A
  self-contained report also means 3Dmol.js has to be vendored into the file rather than
  loaded from a CDN — `tests/test_m6_report.py::test_the_file_has_no_external_references`
  and the full report's equivalent both pin that, and they should keep passing.
- **RDKit ligand SVGs.** Needs `ligands` from M4. Depictions are generated at render time
  from SMILES and inlined; RDKit is a build-time dependency of the report, not a runtime
  one for the reader.

**Scope**
- Vendor 3Dmol.js (or an equivalent small viewer) inline; render receptor, pocket and pose
  for each candidate that has a structure.
- Inline RDKit-generated SVG per ligand.
- Keep both behind the same "not run" treatment every other section has, so a run without
  M3 or M4 renders exactly as it does today.

**Definition of done**
- [ ] A report with structures shows them; one without renders unchanged
- [ ] The no-external-references tests still pass with the viewer vendored
- [ ] Ligand depictions are inline SVG, never a link or a runtime fetch
- [ ] File size with a full run recorded, since a vendored viewer is not small

**Watch out for:** a predicted structure has no cofactors and a docking pose is a ranking,
not an affinity. The limitations block already says both; a 3D view makes them look far
more authoritative than they are, so the caveat has to sit next to the viewer, not only at
the bottom of the page.

## M1: SNP distance to the closest genomes

labels: module:m1, priority:p2, type:code

**Goal:** a SNP count per close genome, to sit beside the ANI that #7 landed.

**Docs:** `docs/01-m1-genome.md`, `docs/pipeline.md`

Issue #7 asked for ANI *and* SNP distance. ANI shipped (`--ani`, skani); SNP distance did
not, and `genome.closest_genomes[].snp_distance` is still `null` on every row —
`<run>/m1/ani.json` records the deferral so the gap is stated rather than silent.

Two things were left open, and both need deciding before code:

**Which aligner.** Snippy `--ctgs` and Parsnp are what `01-m1-genome.md` originally named.
Snippy is genuinely pairwise — one run per reference, which is the shape
`closest_genomes[]` wants — but it shreds the contigs into synthetic reads and pulls in the
bwa/samtools/freebayes stack, and its counts are reference-biased. Parsnp aligns all ~11
genomes at once and gives every pair from one core alignment, but the core genome collapses
as divergence grows, so the distant rows come back empty anyway. MUMmer4 `dnadiff`
(`nucmer` → `delta-filter` → `show-snps`) is the third option: pairwise, reports
`TotalSNPs` and aligned bases directly, installs from a release tarball with no conda, and
is the lightest of the three. It is not what the original issue named, which is the only
argument against it.

**What the number means below the species boundary.** A SNP count between genomes at 87%
ANI over 76% of their length is not a distance anyone should quote — the alignable
fraction, not the substitutions, is what differs. Whatever lands should refuse to emit a
number below a stated ANI floor (95% is the conventional species line), the same way
`--ani` leaves `ani: null` and records the reason when skani drops a pair below `--min-af`.

**Scope**
- Pick the aligner, write down why in `docs/01-m1-genome.md`.
- Reuse the reference genomes `--ani` already downloaded into `<run>/m1/genomes/` — they
  are cached, and re-fetching them would be the only slow part of this.
- Write `snp_distance`, plus the aligned length the count is over, onto each row; a count
  with no denominator is not comparable between pairs.
- Record the tool, version and argv in `<run>/m1/` beside `ani.json`, as `--ani` does.
- Refuse to emit a number below the ANI floor, and record the reason instead.

**Definition of done**
- [ ] `snp_distance` populated for every close genome above the ANI floor, `null` with a
      recorded reason below it
- [ ] Every count carries the aligned length it was measured over
- [ ] Tool, version and command in the run record; a second run makes no network calls
- [ ] Calibrated the way `--ani` was: run against copies of a fixture genome mutated at a
      known substitution rate, and check the recovered count against the number injected

**Watch out for:** counting SNPs against a reference whose own assembly is poor.
`243273.27` in the *M. genitalium* hit list is `genome_quality: Poor` — its SNP distance
would be measuring that assembly's errors as much as real divergence.

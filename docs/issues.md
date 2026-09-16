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

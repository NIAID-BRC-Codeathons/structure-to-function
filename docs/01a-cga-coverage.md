# M1a — What the CGA run covers, and what it does not

Written after the first real Comprehensive Genome Analysis run, 2026-09-17.
Source data: `fixtures/cga/mgen_G37/` (trimmed) and `data/cga/` (full, gitignored).

## The run

| | |
| --- | --- |
| Input | `fixtures/genomes/mgen_G37/mgen_G37.fna`, header blinded to `contig_1` |
| Genome | *M. genitalium* G37, 580,076 bp, one contig |
| Output genome ID | 243273.147 |
| Genetic code | 4 (confirmed in the annotation output) |
| Runtime | **164.5 s** (`elapsed_time` in the job's `annotation` file, job 23587516) |
| Result | 530 CDS, 35 tRNA, 3 rRNA, 16 repeat regions; quality Good |

The known taxon was passed directly, so this run is not blinded.

## Covered by CGA

| Downstream need | Issue | Coverage | Source file |
| --- | --- | --- | --- |
| CDS, products, coordinates, protein FASTA | #5 | Full. prodigal + glimmer3; annotated by `kmer_search` and `annotate_proteins_similarity` | `load_files/genome_feature.json`, `annotation.gff`, `annotation.feature_protein.fasta` |
| Subsystems | #5, #19 | Full. 85 subsystems, 341 role bindings, 8 superclasses including STRESS RESPONSE, DEFENSE, VIRULENCE (10 subsystems, 31 genes) | `load_files/subsystem.json`, `annotated.genome` |
| Phylogenetic tree | #5 | Present but weak: 11 genomes from only 5 filtered single-copy genes | `codontree_tree.nwk`, `codontree.genesPerGenome.txt` |
| Job ID and wall-clock runtime | #5 | Full | `annotation` |
| Tool versions, parameters, timestamps | #4 | Full. 24 tools, e.g. AMRFinder db 2025-07-16.1 / sw 4.0.3, MLST 2.23.0 | `analysis_events` in `annotated.genome` |
| EC, GO, pathways, protein families | #10 | Partial. 161 EC, 127 GO, 113 features across 42 pathways, PGFAM 524, PLFAM 521. No COG, no KEGG KO | `load_files/pathway.json`, features |
| AMR calls | #12, #19 | Partial. 14 rows with drug lists, PMIDs and a classification; no identity or coverage on the 13 k-mer rows | `load_files/sp_gene.json`, `quality.amr_genes` |
| Drug target and transporter hits | #12 | Full for what it finds: 1 TTD, 2 TCDB, with identity and coverage | `specialty-blast.txt` |
| Assembly QC | #5 (optional) | Full. Completeness 100, contamination 0, coarse/fine consistency 99.9/99.7, N50/L50/N70/N90 | `quality.json` |
| BV-BRC access smoke test | #23 | Login and a real service submission both proven | — |

## Not covered by CGA

| Gap | Owned by |
| --- | --- |
| Closest genomes, predicted taxon. `close_genomes` is **empty**; the tree ingroup gives 10 genome IDs and branch lengths but no Mash distance | #5, and the new access-path issue |
| ANI to the closest genomes | #7 — done, `--ani` (skani); see [01-m1-genome.md](01-m1-genome.md#ani-to-the-closest-genomes) |
| SNP distance to the closest genomes | deferred out of #7 to its own issue |
| Essentiality and human homology. Confirmed absent from all 43 output files | new issue |
| Localization and membrane flags (`flags.secreted`, `flags.membrane`) | #10 |
| Structures, pockets, ligands, docking | #8, #9, #13-#18 |
| Knowledge graph | #11 |
| Growth and nutrients (BacDive, MediaDive, GapMind), PHI-base, secretion systems | #19 |
| RAG synthesis and claim verification | #20 |
| Mobile element detection (optional M1 item) | not run |

## Findings

1. **Specialty genes from a fresh CGA run are far fewer than from a public BV-BRC
   genome record.** 17 rows here, versus 174 for the public G37 record
   (243273.25), which adds 148 FBA essential genes and 5 human homologs. Those
   come from BV-BRC's precomputed analysis of public genomes, not from CGA, so
   the blinded pipeline cannot rely on them.

2. **Most "Antibiotic Resistance" rows are not resistance determinants.** By
   `amr_gene_summary`: 10 are "antibiotic target in susceptible species" (gyrA,
   gyrB, rpoB, rpoC, S10p, S12p, EF-G, EF-Tu, folA, Iso-tRNA), 2 "protein
   altering cell wall charge conferring antibiotic resistance" (PgsA, GdpD), and
   1 "gene conferring resistance via absence" (gidB). Anything scoring AMR must
   read `classification`, not just `property`.

3. **Confidence fields depend on the evidence type.** `K-mer Search` rows have
   null `identity`, `query_coverage`, `subject_coverage` and `e_value`;
   confidence lives per feature in `annotation.genome` as `quality.hit_count`,
   `quality.weighted_hit_count` and `quality.priority`. `DIAMOND` rows carry
   identity and coverage. Issue #5's "database, identity and coverage" check
   needs rewording or a separate AMRFinderPlus/RGI run.

4. **Several analyses ran but produced nothing for this species.**
   `p3x-compute-amr-classification` ran two `models.spcAb.*` models yet
   `genome_amr.json` is `[]`; cgMLST called 0 percent of loci;
   `specialty-amrfinder.txt` and `specialty-rgi.txt` are header-only. Whether
   these populate for a well-characterized pathogen is untested.

5. **Virulence detection is unverified.** `pipeline.md` expects virulence factors
   from CGA specialty genes via VFDB and Victors. This run returned zero
   virulence rows, and `specialty-blast.txt` shows only TCDB and TTD were
   searched.

6. **The codon tree used 5 genes.** `codontree.genesPerGenome.txt` reports
   `Filtered_SingleCopy = 5` for every genome in the ingroup. Fine as a smoke
   test, not defensible in the report.

7. **The 164-second runtime weakens the case for the p3-CLI fallback path (#6).**
   That issue exists so nobody waits on a CGA job. A 580 kb genome finished in
   under three minutes; USA300 at 2.8 Mb will be slower but likely still
   minutes. #6 remains useful as an offline fixture source, but is no longer the
   critical unblock and should drop from p0.

## Output inventory

43 files, 16 MB. The pieces worth parsing:

- `.annotation/load_files/*.json` - the indexed tables: `genome`, `genome_feature`,
  `sp_gene`, `subsystem`, `pathway`, `feature_sequence`, `genome_sequence`,
  `genome_typing`; `genome_amr` and `taxonomy` came back empty.
- `annotated.genome` / `.annotation/annotation.genome` - the full genome object:
  features with k-mer confidence, `analysis_events`, `quality`, `subsystems`,
  `close_genomes` (empty), `typing`, `job_metadata`.
- `annotation` - job record: parameters as submitted, `start_time`, `end_time`,
  `elapsed_time`, and the 31 output files.
- `.annotation/quality.json` - quality metrics, `amr_genes`, summaries,
  `problematic_roles_report`.
- `codontree*` - tree, ingroup, per-genome gene counts, alignment stats.
- Flat formats we do not need to parse: `annotation.gb`, `.embl`, `.gff`,
  `.txt`, `.xls`, `.tar.gz`, `FullGenomeReport.html`, `circos.png/svg`.

Do not parse `genome_sequence.json` or `feature_sequence.json` into
`report.json`: sequences belong in `proteins.faa`.

## Taxon resolution experiment (2026-09-17)

Same contigs, submitted twice: once with the exact taxon and genetic code, once
as a blinded unknown where the only fact available is "it is a bacterium". The
floor run used taxon 2 (Bacteria) and genetic code 11, which is what BV-BRC's
own taxonomy record for taxon 2 returns.

| | Exact: taxon 243273, code 4 | Floor: taxon 2, code 11 |
| --- | --- | --- |
| Genome ID | 243273.147 | 2.119822 |
| Runtime | 164.5 s | 156.4 s |
| CDS | 530 | **1070** |
| Mean protein length | 346.5 aa | **142.2 aa** |
| Median protein length | 284 aa | 103 aa |
| Total coding aa | 183,623 | **152,170** |
| Proteins under 100 aa | 44 (8%) | **517 (48%)** |
| Hypothetical | 139 (34.2%) | 333 (41.1%) |
| Genome quality | Good | **Poor** |
| Quality flags | none | **Abnormal CDS ratio** |
| Fine consistency | 99.7 | 92.1 |
| PLFAM assignments | 521 | **none** |
| PGFAM assignments | 524 | 1007 |
| Specialty genes | 17 | 38 |
| Subsystem rows / distinct | 341 / 85 | 571 / 84 |
| Pathway rows / distinct | 226 / 42 | 452 / 41 |
| Tree ingroup | same 10 genomes | same 10 genomes |

The mechanism is UGA read as a stop instead of tryptophan, so genes fragment:

- `rpoB` is one 1390 aa protein in the exact run; in the floor run it is three
  pieces (89, 358, 898 aa).
- `ileS` is one 895 aa protein; in the floor run, five pieces (71, 45, 81, 90,
  121 aa and more).
- `gyrA` survives intact at 836 aa in both, so the damage is uneven.

Downstream consequences:

1. **Fragment inflation looks like more data, not less.** Specialty genes go
   from 17 to 38 because one gene is counted several times: gyrB appears 5
   times, Iso-tRNA 5, EF-G, GdpD, rpoB and rpoC 3 each. Subsystem and pathway
   row counts roughly double while distinct counts stay flat. Any triage score
   that counts hits rather than distinct genes will be badly skewed.
2. **PLFAM assignment disappears entirely**, since local families are defined
   per genus and there is no genus. PGFAM assignment still works.
3. **CGA does flag it**: quality Poor, `Abnormal CDS ratio`, fine consistency
   92.1. M1 should treat those as a hard gate, not a note.
4. **The codon tree ingroup is unaffected** - the same 10 genomes, including
   243273.25, in both runs. Ingroup selection is sequence-based, not taxon-based,
   so it is a usable relative-finder even when the taxon is wrong.

Conclusion: taxon and genetic code must come from Minhash before CGA runs, and
the M1 parser should refuse a run whose genome quality is Poor or whose quality
flags are non-empty.

## USA300 versus M. genitalium (2026-09-17, issue #33)

Everything that came back empty for *M. genitalium* populates for *S. aureus*
USA300_FPR3757, so the analyses run for every genome — the smoke-test organism
simply has little to find. Blinded contigs, taxon called by Minhash at species
level (1280), genetic code 11, genome 1280.69637.

| Analysis | USA300 (2.92 Mb) | M. genitalium (580 kb) |
| --- | --- | --- |
| CDS | 2767 | 530 |
| CGA compute | **682 s** | 252 s |
| Specialty rows, total | 375 | 17 |
| Virulence: VFDB | 89 | none |
| Virulence: Victors | 35 | none |
| AMR: PATRIC k-mer | 41 | 14 |
| AMR: CARD via RGI (homolog + variant models) | 24 | none |
| AMR: NDARO via AMRFinderPlus | 10 | none |
| AMR phenotype (`genome_amr`, MIC/SIR) | **15** | none |
| Transporter: TCDB | 140 | 2 |
| Drug target: TTD | 27 | 1 |
| Metal resistance: BacMet | 9 | none |
| cgMLST, loci called | **96.35%** | 0% |
| MLST sequence type | **none called** | none called |
| Subsystem rows | 1417 | 340 |
| Pathway rows | 1287 | 225 |

Consequences:

1. **`pipeline.md` is correct about virulence and AMR sources**, but their yield is
   organism-dependent. A poorly-characterized organism returns nothing from VFDB,
   Victors, CARD or AMRFinderPlus, and that is a real absence of evidence, not a
   broken pipeline. M2 and M5 must distinguish "no hits" from "not run" — the
   specialty summary in `quality.json` says which analyses executed.

2. **AMR phenotype data will exist for the real test genome, with caveats.**
   `genome_amr` carries 15 rows for USA300 and none for *M. genitalium*: 3 MIC
   predictions (daptomycin 1.0, oxacillin 4.0, vancomycin 2.0 mg/L) and 12 SIR
   calls. Every row is XGBoost model output - `evidence: Computational Method`,
   `vendor: BVBRC`, `computational_method_version: 20250225` - so these are
   predictions, not laboratory measurements, and M5 must present them as such.

   | Call | Antibiotics |
   | --- | --- |
   | Resistant | ciprofloxacin, clindamycin, erythromycin, cefoxitin, methicillin, oxacillin, penicillin, tetracycline, daptomycin |
   | Susceptible | fusidic acid, `geamycin`, `co_trimoxazole` |

   Three handling requirements for M5:

   - **Gate on `computational_method_performance`.** Each row carries its own F1 or
     W1 score with a confidence interval. Most are 0.91-0.99; daptomycin SIR is
     **F1 0.3, CI[-0.26, 0.85]**, an interval spanning zero, and should be dropped
     rather than reported.
   - **MIC and SIR can contradict each other.** Daptomycin is called Resistant by
     the SIR model while the MIC model returns 1.0 mg/L, the susceptible
     breakpoint. A documented reconciliation rule is needed; here the
     low-confidence SIR row is the one to discard.
   - **Normalize antibiotic names before joining to ChEMBL or openFDA.** This run
     contains `geamycin` (gentamicin) and `co_trimoxazole`. Same class of problem
     as BV-BRC's `Virulance factor` spelling, which its own data also carries.

   The remaining calls are consistent with MRSA: methicillin, oxacillin, cefoxitin
   and penicillin resistant, the oxacillin MIC of 4.0 mg/L agreeing with its SIR
   call, and vancomycin susceptible at 2.0 mg/L.

3. **Human homolog is still missing.** BV-BRC's public record for 451515.3 has 648
   specialty rows including 21 human-homolog rows and 34 DrugBank rows; a fresh CGA
   run gives 375 and neither category. So the finding above holds even for a
   well-characterized pathogen: M2 has to compute human homology itself. (DrugBank
   is not shippable anyway.)

4. **MLST is unreliable.** `p3x-compute-mlst` ran and called no sequence type, even
   though cgMLST called 96.35% of loci and BV-BRC's public record types this genome
   as ST8. Do not depend on MLST.

5. **Runtime scales sub-linearly with CDS count**: 5.2x the CDS for 2.7x the compute.
   Queue time still dominates wall-clock.

`scripts/cga_coverage_report.py` regenerates this table from two retrieved CGA
directories.

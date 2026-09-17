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
| ANI and SNP distance | #7 |
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

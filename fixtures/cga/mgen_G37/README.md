# CGA output fixture: M. genitalium G37

Trimmed output of one BV-BRC Comprehensive Genome Analysis run, so the M1
parser can be built and tested without submitting a job.

## Run

- Input: `fixtures/genomes/mgen_G37/mgen_G37.fna`, header blinded to `contig_1`
- Command: `p3-submit-CGA --contigs ... --scientific-name "Mycoplasma genitalium"
  --taxonomy-id 243273 --code 4 --domain Bacteria --label smoke`
- Submitted 2026-09-17 by cmann@bvbrc; output genome ID 243273.147
- Runtime: 164.5 s (`elapsed_time` in the job's `annotation` file; job 23587516,
  pear.cels.anl.gov). Job stages ran 1789655083-1789655222 epoch.
- Genetic code 4 confirmed in the annotation output. Mycoplasma reads UGA as
  tryptophan; code 11 would truncate many CDS.
- The known taxon was passed directly, so this run is NOT blinded. The real M1
  path predicts the taxon with Similar Genome Finder first.

## Files

- `genome.json` - genome-level record: 530 CDS, 35 tRNA, 3 rRNA, quality Good,
  coarse/fine consistency 99.9/99.7, 34.2 percent hypothetical
- `sp_gene.json` - all 17 specialty-gene rows (14 antibiotic resistance,
  2 transporter, 1 drug target)
- `genome_feature.sample.json` - 5 CDS covering the cases the parser must
  handle: gyrA (AMR), gpmI (TTD drug target), an amino-acid permease
  (transporter), a hypothetical protein, and ribosomal protein S12p
- `quality.summary.json` - the summary blocks of `.annotation/quality.json`
- `amr_genes.json` - the `amr_genes` table and `amr_gene_summary` grouping,
  which carry the AMR classification per gene
- `annotation_genome.feature_sample.json` - the same 5 features as they appear
  in `annotation.genome`, with k-mer confidence (`quality.hit_count`,
  `weighted_hit_count`, `priority`), annotation events, GO/EC terms and
  PGFAM/PLFAM assignments; protein sequences stripped
- `specialty-blast.txt` - the 3 DIAMOND specialty hits with identity and
  coverage
- `codontree_tree.nwk`, `tree_ingroup.txt` - phylogenetic tree and the 10
  genomes in it

Full output (16 MB) is not committed. It lives under `data/cga/`, which is
gitignored, and can be regenerated with the command above.

## Findings that affect the pipeline

1. CGA specialty genes do NOT include Essential Gene or Human Homolog
   properties. BV-BRC's public genome record for 243273.25 carries 174 rows
   including 148 FBA essential genes and 5 human homologs; a fresh CGA run
   gives 17. M2 triage must compute essentiality and human homology itself.

2. Most of the 14 "Antibiotic Resistance" rows are not resistance
   determinants. By `amr_gene_summary`: 10 are "antibiotic target in
   susceptible species" (gyrA, gyrB, rpoB, rpoC, S10p, S12p, EF-G, EF-Tu,
   folA, Iso-tRNA), 2 "protein altering cell wall charge conferring antibiotic
   resistance" (PgsA, GdpD), 1 "gene conferring resistance via absence"
   (gidB). Triage must read `classification`, not just `property`, or it will
   score conserved drug targets as resistance evidence.

3. Confidence fields differ by evidence type:
   - `K-mer Search` rows (14, all PATRIC): `identity`, `query_coverage`,
     `subject_coverage` and `e_value` are all null. Confidence lives per
     feature in `annotation.genome` as `quality.hit_count`,
     `quality.weighted_hit_count` and `quality.priority`, plus `pmid` on the
     sp_gene row.
   - `DIAMOND` rows (3: 2 TCDB, 1 TTD): identity and both coverages present,
     also in `specialty-blast.txt`.
   Issue #5's "database, identity and coverage" check therefore needs either
   rewording to accept k-mer confidence, or AMRFinderPlus/RGI run separately.

4. Empty or unpopulated outputs: `genome_amr.json` and `taxonomy.json` are
   `[]`, `specialty-amrfinder.txt` and `specialty-rgi.txt` are header-only,
   cgMLST called 0 percent of loci, and `similarity_associations` is empty on
   the features checked.

5. The tree contains 243273.25, the public G37 genome, at ~0 distance from our
   run: the expected self-match. Report closest non-identical genomes too.

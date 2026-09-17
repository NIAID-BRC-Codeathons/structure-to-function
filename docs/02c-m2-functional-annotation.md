# M2c — Functional annotation, localization and membrane flags

Companion to [02-m2-triage.md](02-m2-triage.md), covering M2 step 4 (`annotate_function()`).
Issue [#10](https://github.com/NIAID-BRC-Codeathons/structure-to-function/issues/10).
Implemented in `s2f/m2_triage/function.py`.

The acceptance check is "`flags.secreted` and `flags.membrane` set for every protein". The
second requirement is the one that shapes the design: **the flag has to carry the evidence that
produced it.** Pitfall #3 deprioritizes membrane proteins explicitly rather than silently, and a
deprioritization that rests on a hydropathy guess is not the same claim as one resting on
DeepTMHMM. Both are recorded, with the source named, so M3 and M6 can tell them apart.

## Sources and precedence

Nothing here installs a database. Whichever of these is present gets merged; the rest are
absent, not false.

| Rank | Source | How it arrives | Confidence |
| --- | --- | --- | --- |
| 60 | UniProt, experimental evidence (ECO:0000269, ECO:0007744) | REST, keyed by the accession `ids.py` resolved | `experimental` |
| 50 | DeepTMHMM | `TMRs.gff3` you generated | `predicted` |
| 50 | SignalP 6 | `prediction_results.txt` | `predicted` |
| 45 | PSORTb | `-o terse` or long output | `predicted` |
| 40 | InterProScan | TSV or JSON; run locally when an install is found | `predicted` |
| 30 | UniProt, non-experimental | same REST call, by similarity / sequence analysis | `inferred` |
| 20 | eggNOG-mapper | `*.emapper.annotations` | `inferred` |
| 10 | built-in sequence heuristic | computed in-process, no network, no database | `heuristic` |

This is the same principle [00a-data-contract.md](00a-data-contract.md) applies to structures —
experimental evidence outranks a prediction wherever both exist — carried over to annotation.
The difference worth knowing is that a bacterial UniProt entry is very often itself a
prediction, so the evidence code, not the database, decides where it sits.

Two rules make the table mean something:

1. **A provider that said nothing sets nothing.** A missing DeepTMHMM run does not become
   `membrane: false`. Each flag records which source decided it; `--no-heuristic` leaves the
   flag's source empty and writes an **empty cell**, not `False`, so "unknown" survives the trip
   to TSV.
2. **The highest-ranked source that answered wins**, and its evidence string goes into
   `*_evidence`. A UniProt entry carrying experimental evidence outranks a predictor; a
   predictor outranks UniProt's own "by similarity" line, because a bacterial UniProt entry is
   very often itself a prediction.

Names are ranked separately (`DESCRIPTION_RANK`: UniProt, then eggNOG, then InterProScan).
How good a source is at predicting membrane topology says nothing about how good its protein
name is.

## Flag definitions

| Flag | Definition |
| --- | --- |
| `tm_helices` | count of transmembrane segments (helices, plus beta-barrel strands where the source reports them) |
| `membrane` | `tm_helices > 0`, or a localization of cytoplasmic/outer membrane |
| `signal_peptide` | an N-terminal secretion signal |
| `lipoprotein` | lipobox / SignalP LIPO or TATLIPO / UniProt lipidation feature |
| `secreted` | localization is extracellular, **or** a signal peptide with no TM segment and no lipid anchor |
| `surface_exposed` | secreted, or outer membrane / cell wall, or a lipoprotein |
| `localization` | one of: cytoplasm, cytoplasmic membrane, periplasm, outer membrane, cell wall, extracellular, or empty |

A lipoprotein has a signal peptide but stays tethered to a membrane, so it is **not** secreted
while still being surface-exposed. Triage wants both facts; collapsing them loses the
distinction that matters for an antibody or a surface-acting compound.

## The built-in heuristic

Every protein must get a flag, and the proteins triage cares about most — hypotheticals with no
UniProt accession and no ortholog — are exactly the ones the external tools reach last. So the
module carries a dependency-free fallback and labels it `heuristic`:

- **Transmembrane:** Kyte-Doolittle hydropathy, 19-residue window, mean ≥ 1.6, sustained across
  ≥ 9 consecutive windows. Spans are reported as the *centres* of the passing windows, so the
  numbers are residue positions and two helices across a short loop stay two spans.
- **Signal peptide:** von Heijne's three parts, all required — a K or R in residues 1–12; an
  h-region of 8–20 residues with no charged residue, no N/Q and no proline, mean KD ≥ 2.2,
  starting in residues 3–18; and a small residue at −1 and −3 of a cleavage site 14–45 residues
  in and at least 3 residues past the h-region (the c-region).
- **Lipoprotein:** lipobox `[LVI][ASTVIG][GAS]C` with the cysteine in residues 12–45 *and* at
  the end of a called h-region. The motif alone, without that context, fires on several percent
  of ordinary sequences.
- A signal peptide and an N-terminal TM helix are the same shape on a hydropathy plot. The
  heuristic attributes its N-terminal segment to a signal peptide when it called one; and if a
  higher-ranked source (SignalP, DeepTMHMM, UniProt) calls a signal peptide the heuristic
  missed, the merge takes one helix back off the heuristic's count. Both moves are written into
  `membrane_evidence`, so the ambiguity stays visible instead of being resolved silently.

**The thresholds were calibrated, not chosen by eye.** Both the window rule and the lipobox
were first written loose enough that ~48% of random sequences of *E. coli* amino-acid
composition came out as membrane proteins and ~4% as lipoproteins. The values above were picked
by sweeping them against that random-sequence null and against the two real proteomes.

**Known limitations.** It misses signal peptides with a neutral n-region; it does not
distinguish Sec from Tat; it has no Gram-positive/negative model; beta-barrel outer-membrane
proteins are invisible to a hydropathy scan. It is a ranking aid, not a substitute for
DeepTMHMM. `run.json` reports `heuristic_only`, so any run can say how much of its proteome
rests on it.

### Measured behaviour

*M. genitalium* G37, 542 proteins, heuristic only, no network (`fixtures/genomes/mgen_G37`),
against 2,000 random 300-residue sequences of *E. coli* amino-acid composition as the null:

| | M. genitalium G37 | random null |
| --- | --- | --- |
| Runtime | **53 ms for 542 proteins** (98 µs/protein), Python 3.11, single core | |
| membrane | 91 (16.8%) | 5.5% |
| signal peptide | 32 (5.9%) | 3.2% |
| lipoprotein | 8 (1.5%) | 0.0% |
| secreted | 11 (2.0%) | |

Sanity check against M1 product names, which the heuristic never sees:

| Product matches | n | membrane | signal | lipobox |
| --- | --- | --- | --- | --- |
| permease | 20 | 18 | 4 | 0 |
| ribosomal protein | 52 | 0 | 0 | 0 |
| synthetase / ligase | 28 | 0 | 0 | 0 |
| lipoprotein | 19 | 3 | 7 | 5 |
| hypothetical | 161 | 45 | 16 | 0 |

Permeases are almost all membrane; ribosomal proteins and aminoacyl-tRNA synthetases none at
all. That is the shape it should have. It is **not a benchmark** — nothing here was scored
against a labelled set, the random null is a floor rather than a false-positive rate, and the
19 "lipoprotein"-named proteins show the recall problem plainly: 7 get a signal peptide and 5 a
lipobox. Mycoplasma lipoprotein signal peptides are unusual, but the honest reading is that
this heuristic finds roughly a third of them. Run DeepTMHMM and SignalP over the selected set
before anything depends on these flags.

## Running it

```bash
# heuristic only, no network, no external tools
python -m s2f.m2_triage --run runs/<id> --annotate

# with tool outputs you generated elsewhere
python -m s2f.m2_triage --run runs/<id> --annotate \
    --eggnog     out/proteins.emapper.annotations \
    --deeptmhmm  out/TMRs.gff3 \
    --signalp    out/prediction_results.txt \
    --psortb     out/psortb_terse.txt

# find a local InterProScan install and run it over the proteome
python -m s2f.m2_triage --run runs/<id> --annotate --interproscan auto

# curated annotation for the proteins that mapped to UniProt
python -m s2f.m2_triage --run runs/<id> --annotate --map-ids --uniprot-function
```

`--dry-run --annotate` reads the fixture provider outputs in `fixtures/m2/annotation/`, so the
parsers are exercised offline rather than only the heuristic.

`--interproscan auto` looks, in order, at `--interproscan-path`, `$INTERPROSCAN_HOME`,
`$INTERPROSCAN`, `PATH`, then `/opt`, `/usr/local`, `/software`, `/apps`, `/share/apps` and the
matching directories under `$HOME`. It records the path and version it found in `run.json`
(pitfall #19). No install found is a note, not a failure: the run continues on the other
providers. An existing `interproscan.tsv` in the output directory is reused rather than
re-scanned, because a proteome-scale scan is hours.

Through the pipeline runner, which always passes `--report`:

```bash
python -m s2f.run --run runs/<id> --from-cga-dir data/cga --annotate --interproscan auto
```

Passing any provider path to the runner implies `--annotate`, so `--psortb out.txt` alone does
what it looks like it does.

## Outputs

### `report.json` (the shared contract)

With `--report`, the flags go into `proteins[].flags` and the evidence into
`proteins[].annotations[]`, through `report_adapter.py` (see
[00a-data-contract.md](00a-data-contract.md)):

| Field | Contents |
| --- | --- |
| `flags.membrane`, `flags.secreted` | the two required flags — **`null` when no source called them**, which is not `false` |
| `flags.localization`, `tm_helices`, `signal_peptide`, `lipoprotein`, `surface_exposed` | the rest of the topology picture |
| `flags.functional_evidence` | `membrane_source`, `membrane_evidence`, `signal_source`, `signal_evidence`, `localization_source`, `confidence`, and every source that contributed |
| `annotations[]` | one row per source, each with its terms, its evidence string, its confidence and `retrieved_at` — merged flags live in `flags`, but `source` is what keeps the evidence apart |

Without `--annotate`, every one of those fields stays `null`: the adapter reserved those slots
for this issue and a protein nobody annotated must not reach M3 looking like a confirmed
cytoplasmic one.

### Side files

Added to `runs/<run_id>/m2_pdb/`:

| File | Contents |
| --- | --- |
| `proteins.tsv` | 23 new columns: the flags, the terms, and `*_source` / `*_evidence` / `annotation_confidence` for each flag |
| `function_terms.tsv` | long format, one row per (protein, term, source): GO, EC, KEGG KO, COG, COG category, Pfam, InterPro, gene name |
| `annotate_all.faa` | cleaned, deduplicated FASTA of every protein — the input for the external tools |
| `annotate_selected.faa` | the same for the selected ~50, for the tools that are slow per sequence |
| `run.json` | `functional_annotation`: per-provider counts, unmatched provider rows, flag totals, localization histogram, which source set each membrane flag, `heuristic_only`, the FASTA stats and the InterProScan install record |

`run.json`'s `functional_annotation.notes` carries anything that went wrong: a file that did not
exist, a file that failed to parse, and — the quiet one — a file that parsed without error and
produced **zero records**, which almost always means an unrecognised format rather than a tool
that found nothing. A broken provider file is a note, never a failed run: the triage output is
written regardless.

**This step does not touch the triage score.** `WEIGHTS` in `score.py` is unchanged and the
ranking is bit-identical with and without `--annotate` (there is a test that asserts this).
Feeding `surface_exposed` and `membrane` into the score is issue #12's change, and it gets its
own dated row in the [weight change log](02a-m2-pdb-evidence.md#weight-change-log).

## Identifier matching

Provider outputs are keyed by whatever FASTA header they were given. Each protein is looked up
by `feature_id` first, then by its UniProt accession (when `--map-ids` ran), then by locus tag
and gene name. Rows that match nothing are counted in `run.json` under `unmatched_ids` rather
than dropped in silence — a non-zero count there usually means the FASTA fed to the tool had
different headers, which is the failure mode worth catching early.

Keep `fig|` IDs in the FASTA you give these tools. Some of them mangle headers; if a tool
insists on short names, write the mapping out and rename its output before ingesting.

---

# Setting up the external tools

None of these are required — the module runs without them.

The commands below are the documented shape of each tool's CLI, not a transcript of a verified
install: versions, container tags and flag names move. Check the tool's own instructions before
a first install, and fix anything here that turns out to be stale.

In rough order of value per hour spent: DeepTMHMM and SignalP 6 (they set the two required flags properly), then eggNOG-mapper
(function terms for everything, including hypotheticals), then PSORTb, then InterProScan.

## InterProScan

Already installed on lambda13 at
`/nfs/lambda_stor_01/homes/cmann/software/interproscan-5.78-109.0` (version 5.78-109.0,
installed 2026-09-03 for the GenSLM-ESM homology work), with `INTERPROSCAN_HOME` exported in
`~/.bashrc` — so `--interproscan auto` finds it there without being told. Member database
versions recorded with that install: Gene3D 4.3.0, FunFam 4.3.0, Pfam 38.2, SUPERFAMILY 1.75.

Elsewhere, to check from a shell:

```bash
which interproscan.sh || ls -d /opt/interproscan* /software/interproscan* 2>/dev/null
interproscan.sh --version
export INTERPROSCAN_HOME=/path/to/interproscan-5.x   # if discovery misses it
```

Fresh install (needs Java 11+ and tens of GB for the member databases). Take the current
version and its exact tarball name from the release directory listing at
<https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/> — the shape is:

```bash
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/<version>/interproscan-<version>-64-bit.tar.gz
tar -pxzf interproscan-<version>-64-bit.tar.gz
cd interproscan-<version> && python3 setup.py -f interproscan.properties   # indexes the HMMs
```

The version numbers below are deliberately left as placeholders: check the install
instructions shipped with the release rather than trusting a line copied from a doc.

Run it directly if you would rather not go through the module:

```bash
interproscan.sh -i proteins.faa -f TSV -o interproscan.tsv -goterms --cpu 8
python -m s2f.m2_triage --run runs/<id> --annotate --interproscan interproscan.tsv
```

`--interproscan auto` scans the cleaned, deduplicated FASTA this module writes rather than
`m1/proteins.faa`: InterProScan rejects `*` stop characters and its cost scales with *unique*
sequences. Duplicate groups are collapsed to one record and the result is fanned back out to
every member (the same reason the GenSLM-ESM homology run deduplicated by MD5 across 2,195
genomes). The three CGA proteomes currently in `data/` happen to contain no stops and no
duplicates, so this is insurance for the next genome rather than a fix for these.
Expect roughly 1–5 s per protein per core with the default analysis set; `-appl Pfam,TIGRFAM`
cuts that by an order of magnitude when only domain terms are wanted. The EBI web service
(`https://www.ebi.ac.uk/Tools/services/rest/iprscan5`) is a fallback for a handful of
sequences, not for a proteome — it is one job per sequence and requires an email parameter.

## DeepTMHMM

No local install needed; it runs on BioLib.

```bash
pip install pybiolib
biolib run DTU/DeepTMHMM --fasta proteins.faa   # writes biolib_results/TMRs.gff3
python -m s2f.m2_triage --run runs/<id> --annotate --deeptmhmm biolib_results/TMRs.gff3
```

It needs outbound network and a BioLib account for larger jobs, and it is slow on a full
proteome — a few thousand sequences is an hour-ish. Running it only over the selected ~50 is a
reasonable compromise: the rest keep their heuristic flag, and `run.json` says which is which.

## SignalP 6.0

Free for academic use, downloaded after accepting the licence at
<https://services.healthtech.dtu.dk/services/SignalP-6.0/> (choose the "fast" model unless you
need the slow one).

```bash
pip install signalp-6-package/    # the tarball is a pip-installable package
signalp6 --fastafile proteins.faa --output_dir signalp_out \
         --organism other --format none --mode fast
python -m s2f.m2_triage --run runs/<id> --annotate \
       --signalp signalp_out/prediction_results.txt
```

`--organism other` is the bacterial setting. `--format none` skips the per-protein plots, which
is most of the runtime.

## eggNOG-mapper

The database is ~50 GB; the web server at <http://eggnog-mapper.embl.de> takes a FASTA and
emails a result, which is the faster route for one proteome.

```bash
conda install -c bioconda eggnog-mapper
download_eggnog_data.py -y                     # ~50 GB into $EGGNOG_DATA_DIR
emapper.py -i proteins.faa -o proteins --cpu 16 --itype proteins
python -m s2f.m2_triage --run runs/<id> --annotate \
       --eggnog proteins.emapper.annotations
```

The parser reads columns by name from the `#query` header line, so a different eggNOG-mapper
release is fine as long as it still writes that line.

## PSORTb 3.0

Easiest as a container; pick the Gram stain that matches the genome (M1's taxonomy call).

```bash
# image name and tag from https://hub.docker.com/r/brinkmanlab/psortb_commandline
docker run --rm -v "$PWD":/data brinkmanlab/psortb_commandline:<tag> \
       -i /data/proteins.faa -r /data -n -o terse     # -n negative, -p positive, -a archaea
python -m s2f.m2_triage --run runs/<id> --annotate --psortb psortb_terse.txt
```

A web submission at <https://www.psort.org/psortb/> works for a single proteome too. Both the
terse TSV and the long report parse.

Gram stain matters: running the Gram-negative model on a Gram-positive genome invents a
periplasm. If M1's taxonomy is uncertain, say so in the run notes rather than picking one.

## UniProt

No setup — it rides on `--map-ids`, which resolves accessions through `common/ids.py`
(pitfall #11: nothing here maps its own identifiers). Coverage is whatever the mapping achieved;
for a novel genome that can be most of the proteome or very little of it, and the count is in
`run.json`.

## Related docs

- [02-m2-triage.md](02-m2-triage.md) — the module contract
- [00a-data-contract.md](00a-data-contract.md) — `report.json`, the schema and the section writer
- [02d-lambda-runbook.md](02d-lambda-runbook.md) — installing and running the providers on Lambda
- [02a-m2-pdb-evidence.md](02a-m2-pdb-evidence.md) — triage weights and their change log
- [02b-m2-knowledge-graph.md](02b-m2-knowledge-graph.md) — the knowledge subgraph
- [pitfalls.md](pitfalls.md) — #3 membrane proteins, #12 triage bias, #19 reproducibility

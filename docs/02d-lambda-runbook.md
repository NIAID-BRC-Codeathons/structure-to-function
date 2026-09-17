# M2d — Running the annotation layer for real, on Lambda

Companion to [02c-m2-functional-annotation.md](02c-m2-functional-annotation.md). That doc says
what the layer does and how it decides; this one is the operational sequence for getting real
provider output instead of the built-in heuristic.

**Why this doc exists.** Issue #10 landed the code and every parser, and closed on the merge of
PR #37. But no provider has ever been run: on any real genome today, 100% of `flags.membrane`
and `flags.secreted` come from the built-in sequence heuristic, and `run.json` says so
(`functional_annotation.heuristic_only`). Turning that into `predicted` confidence is what the
steps below do.

## What has and has not been verified

| | |
| --- | --- |
| M1 on real CGA output (`data/cga_sgf`, 530 proteins) | **verified**, 0.2 s |
| M2 `--annotate --report` on that M1 output | **verified** offline, 0.6 s, `proteins` section validates |
| Flags set for every protein | **verified** — 530/530, all `heuristic` |
| Tool-ready FASTA emission | **verified** — `annotate_all.faa` 530 records, `annotate_selected.faa` 50 |
| Two-pass ingest on real feature IDs | **verified with stand-in outputs** (see below) |
| Any external provider actually run | **not done** |
| Any provider's runtime on a real proteome | **not measured** |

The two-pass ingest was exercised by generating DeepTMHMM `TMRs.gff3` and SignalP 6
`prediction_results.txt` in their real formats, carrying the 50 real `fig|2097.118.peg.*`
identifiers from `annotate_selected.faa`, and running pass two over them. Result:

```
providers              {'deeptmhmm': 50, 'signalp6': 50, 'heuristic': 530}
unmatched_ids          {}
membrane_flag_source   {'deeptmhmm': 50, 'heuristic': 480}
confidence             {'heuristic': 480, 'predicted': 50}
heuristic_only         480
```

All 50 selected proteins moved from `heuristic` to `predicted`, nothing went unmatched, and the
other 480 kept their heuristic flags. **Those were synthetic files, not predictions** — the
numbers say the plumbing joins up on real identifiers, and say nothing about any protein.

The three CGA proteomes in `data/` (`cga`, `cga_blind`, `cga_sgf`) were checked for the two
things that break these tools: they contain **no `*` stop characters and no duplicate
sequences**. So the cleaning and deduplication in `write_tool_fasta` is insurance for the next
genome, not a fix for these.

## Node facts

From the GenSLM-ESM environment work on the same node (recorded 2026-09-03):

| | |
| --- | --- |
| Node | lambda13, Ubuntu 22.04 |
| Python | 3.13.2 |
| Java | 11.0.32 (OpenJDK) — above the InterProScan minimum |
| Writable | `/nfs/lambda_stor_01/homes/<user>/` only, **not** the `lambda_stor_01` root |
| InterProScan | `/nfs/lambda_stor_01/homes/cmann/software/interproscan-5.78-109.0`, `INTERPROSCAN_HOME` in `~/.bashrc` |
| InterProScan member DBs | Gene3D 4.3.0 · FunFam 4.3.0 · Pfam 38.2 · SUPERFAMILY 1.75 |
| Network | the InterProScan pre-calculated lookup service works from this node, so outbound HTTPS is available — RCSB, UniProt and AlphaFold DB should be reachable too |

## Three machines, and which one runs what

Getting this wrong wastes an afternoon, so it is worth stating plainly.

| Machine | Role |
| --- | --- |
| **Mac** (`~/10 Coding Workspace/structure-to-function`) | the git working copy. Patches are applied here, branches are pushed from here, CGA data lives here |
| **GitHub** | how code reaches lambda. Nothing else does |
| **lambda13** | where the providers run. `git clone` / `git pull` only — never `git am` |

**Code reaches lambda through GitHub, not by copying files.** A patch file is a Mac-side step:
apply it on the Mac, push the branch, then pull on lambda. Never `scp` a patch to lambda and
apply it there — lambda then sits on a commit nobody else has, and the next `git pull` conflicts
with it.

**Data does not go through GitHub.** `data/` and `runs/` are gitignored, so CGA directories go
Mac → lambda by rsync, and results come back the same way.

## 1. Get the code onto lambda

On the **Mac**, land whatever is outstanding and push it:

```bash
cd "/Users/cmann/10 Coding Workspace/structure-to-function"
git switch main && git pull
git switch -c <branch>                  # skip if the work is already on main
git am --3way <patch>                   # Mac-side only
python -m pytest -q
git push -u origin <branch>
```

On **lambda13**, pull it:

```bash
ssh lambda13
cd /nfs/lambda_stor_01/homes/$USER
git clone https://github.com/NIAID-BRC-Codeathons/structure-to-function.git   # first time
cd structure-to-function
git fetch origin && git checkout <branch> && git pull                          # or: git checkout main && git pull

bash scripts/setup_lambda.sh
```

The script picks the newest Python ≥ 3.10 it can find, builds `.venv/`, installs
`requirements-common.txt` + `requirements-m2.txt` + pytest, imports every module, runs the test
suite, and then prints which annotation providers exist on the node. It stops at the first
failure rather than reporting success over a broken install, and re-running reuses the venv.

Expected tail:

```
imports ok
173 passed
annotation providers on this node:
  interproscan  FOUND  /nfs/.../interproscan-5.78-109.0/interproscan.sh (5.78-109.0, via $INTERPROSCAN_HOME)
  deeptmhmm     no     pip install pybiolib, then: biolib run DTU/DeepTMHMM --fasta <faa>
  ...
```

If `interproscan` says `no`, `$INTERPROSCAN_HOME` is not exported in a non-interactive shell —
`export INTERPROSCAN_HOME=/nfs/lambda_stor_01/homes/cmann/software/interproscan-5.78-109.0` or
pass `--interproscan-path`.

Then get a CGA directory onto the node. Run this on the **Mac** — `data/` is gitignored, so it
never travels through GitHub:

```bash
cd "/Users/cmann/10 Coding Workspace/structure-to-function"
rsync -av --exclude __pycache__ \
  data/cga_sgf \
  lambda13:/nfs/lambda_stor_01/homes/$USER/structure-to-function/data/
```

The CGA layout has a hidden `.annotation/` directory that M1 needs; rsync of the parent
directory includes it, but if you copy by hand, copy that too. Verify on lambda before running
anything:

```bash
ls data/cga_sgf/annotation data/cga_sgf/.annotation/load_files/genome_feature.json
```

Results come back the same way when you want them on the Mac:

```bash
# on the Mac
rsync -av lambda13:/nfs/lambda_stor_01/homes/$USER/structure-to-function/runs/sgf runs/
```

## 2. Pass one — pipeline runs, FASTAs come out  *(on lambda)*

```bash
source .venv/bin/activate
python -m s2f.run --run runs/sgf --from-cga-dir data/cga_sgf --annotate --map-ids
```

This is the full M1 → M2 path with network: taxon resolution, RCSB sequence search per unique
sequence, UniProt/ChEMBL id mapping, AlphaFold DB lookups, and the annotation layer on the
heuristic tier. It writes:

| File | Use |
| --- | --- |
| `runs/sgf/report.json` | the contract — `proteins[].flags` / `.annotations[]` |
| `runs/sgf/m2_pdb/annotate_all.faa` | **cleaned, deduplicated FASTA of every protein** — feed this to the tools |
| `runs/sgf/m2_pdb/annotate_selected.faa` | the top ~50 only — feed this to the slow tools |
| `runs/sgf/m2_pdb/run.json` | `functional_annotation`: provider counts, `heuristic_only`, FASTA stats, InterProScan install record |

Sanity-check before spending tool time:

```bash
python -c "
import json; m=json.load(open('runs/sgf/m2_pdb/run.json'))['functional_annotation']
print('annotated', m['proteins_annotated'], 'heuristic_only', m['heuristic_only'])
print('fasta', m['tool_fasta'])
print('notes', m['notes'])"
```

`tool_fasta.id_truncation_warning`, if present, matters: SignalP 6 and some other tools truncate
FASTA identifiers, and a truncated ID matches no `feature_id`. It does not corrupt anything —
those rows land in `run.json`'s `unmatched_ids` — but it means that provider contributed
nothing. Check `unmatched_ids` after every provider run.

## 3. Run the providers  *(on lambda)*

Ranked by value per hour spent. Nothing here is required; each one that runs upgrades some
flags from `heuristic` to `predicted`.

### 3a. DeepTMHMM and SignalP 6 — the two required flags

These are what make `flags.membrane` and `flags.secreted` mean something. Run them over
`annotate_selected.faa` (50 proteins) first; the proteome can wait.

```bash
pip install pybiolib
biolib run DTU/DeepTMHMM --fasta runs/sgf/m2_pdb/annotate_selected.faa
# writes biolib_results/TMRs.gff3
```

SignalP 6 needs a one-time download after accepting the academic licence at
<https://services.healthtech.dtu.dk/services/SignalP-6.0/> (take the "fast" model):

```bash
pip install ./signalp-6-package/          # the tarball is a pip-installable package
signalp6 --fastafile runs/sgf/m2_pdb/annotate_selected.faa \
         --output_dir runs/sgf/signalp --organism other --format none --mode fast
```

`--organism other` is the bacterial setting; `--format none` skips the per-protein plots, which
are most of the runtime.

### 3b. InterProScan — function terms, already installed

```bash
python -m s2f.m2_triage --run runs/sgf --annotate --report --offline \
    --interproscan auto --interproscan-appl Pfam
```

`auto` finds the install, scans `annotate_all.faa` (cleaned and deduplicated — InterProScan
rejects `*` and its cost scales with unique sequences), writes
`runs/sgf/m2_pdb/interproscan.tsv`, and records path and version in `run.json` (pitfall #19). An
existing `interproscan.tsv` is reused rather than re-scanned, so a re-run is cheap.

**Start with `-appl Pfam`.** The install has 18 analyses enabled and the default is all of
them. Measured on this node during the GenSLM-ESM work, with four member databases and the
lookup service: 60 proteins in 64 s, and 89 s with `-dp` (lookup disabled, everything local).
Extrapolating to 530 proteins is roughly 10 minutes — but that was four databases, not
eighteen, so budget accordingly and widen `-appl` only once one pass has completed.

### 3c. eggNOG-mapper — best function coverage, heaviest setup

The database is ~50 GB. The web server at <http://eggnog-mapper.embl.de> takes
`annotate_all.faa` and emails the result, which is faster than installing for one proteome.
Locally:

```bash
conda install -c bioconda eggnog-mapper
download_eggnog_data.py -y
emapper.py -i runs/sgf/m2_pdb/annotate_all.faa -o runs/sgf/eggnog --cpu 16 --itype proteins
```

### 3d. PSORTb — localization

```bash
docker run --rm -v "$PWD":/data brinkmanlab/psortb_commandline:<tag> \
    -i /data/runs/sgf/m2_pdb/annotate_all.faa -r /data/runs/sgf -n -o terse
```

`-n` Gram-negative, `-p` Gram-positive, `-a` archaea. **Pick from M1's taxon call, not by
default:** running the Gram-negative model on a Gram-positive genome invents a periplasm.
*M. genitalium* has no outer membrane at all, so for the current test genomes PSORTb is the
least informative of the four.

## 4. Pass two — ingest what the tools produced  *(on lambda)*

```bash
python -m s2f.m2_triage --run runs/sgf --annotate --report --offline \
    --deeptmhmm   biolib_results/TMRs.gff3 \
    --signalp     runs/sgf/signalp/prediction_results.txt \
    --interproscan runs/sgf/m2_pdb/interproscan.tsv \
    --eggnog      runs/sgf/eggnog.emapper.annotations \
    --psortb      runs/sgf/psortb_terse.txt
```

`--offline` replays the PDB/UniProt work from pass one's cache, so this is seconds rather than
minutes. Pass only the files that exist; a missing one is a note, not a failure.

Then check what actually landed:

```bash
python -c "
import json; m=json.load(open('runs/sgf/m2_pdb/run.json'))['functional_annotation']
print('providers      ', m['providers'])
print('unmatched_ids  ', m['unmatched_ids'])
print('membrane source', m['membrane_flag_source'])
print('confidence     ', m['confidence'])
print('heuristic_only ', m['heuristic_only'])
print('notes          ', m['notes'])"
```

What each line tells you:

- **`providers`** — how many proteins each source covered. A provider you passed that shows 0
  did not match; look at `unmatched_ids`.
- **`unmatched_ids`** — rows whose identifier matched no protein. Non-zero almost always means
  the tool rewrote or truncated the FASTA headers.
- **`membrane_flag_source`** — which source decided each membrane flag. This is the number that
  says whether the run is still resting on the heuristic.
- **`confidence`** — the histogram to quote in the report: `experimental` / `predicted` /
  `inferred` / `heuristic`.
- **`heuristic_only`** — proteins where the heuristic still set a required flag. Driving this
  toward zero for the selected set is the point of the whole exercise.

## 5. Record it

Definition of done on issue #10 asks for a measured full-proteome runtime, which is still
outstanding. When a provider run completes, put its wall-clock time and the
`membrane_flag_source` histogram in the follow-up issue and in
[02c](02c-m2-functional-annotation.md#measured-behaviour), next to the heuristic's numbers.

## Failure modes worth knowing before you start

| Symptom | Cause | What to do |
| --- | --- | --- |
| `interproscan auto: no interproscan.sh found` | `$INTERPROSCAN_HOME` not set in a non-interactive shell | export it, or `--interproscan-path /path/to/interproscan.sh` |
| Provider passed, `providers` shows 0 for it | headers rewritten or truncated by the tool | check `unmatched_ids`; rename the tool's output IDs back to `feature_id` |
| `parsed to zero records — check the file format` in `notes` | the parser did not recognise the file | confirm you passed the right file (`TMRs.gff3`, `prediction_results.txt`, `*.emapper.annotations`, terse PSORTb, InterProScan TSV/JSON) |
| InterProScan rejects the input | stop characters | already handled — but confirm you passed `annotate_all.faa`, not `m1/proteins.faa` |
| `pip install jsonschema` fails | Python < 3.8 | `scripts/setup_lambda.sh` refuses to build below 3.10 for exactly this reason |
| M2 writes no `report.json` | `--report` missing | the runner always passes it; a direct `python -m s2f.m2_triage` call does not |
| `git pull` on lambda conflicts | a patch was applied on lambda instead of on the Mac | `git reset --hard origin/<branch>`; apply patches on the Mac and push |
| M1 fails with a missing-file error | the hidden `.annotation/` directory did not come across | re-rsync the CGA parent directory, then check `ls data/<cga>/.annotation/load_files/` |

## Related docs

- [02c-m2-functional-annotation.md](02c-m2-functional-annotation.md) — precedence, flag definitions, calibration
- [00a-data-contract.md](00a-data-contract.md) — `report.json` and the schema
- [pitfalls.md](pitfalls.md) — #3 membrane proteins, #19 reproducibility

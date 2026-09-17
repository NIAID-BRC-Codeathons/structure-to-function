# M2d — Running the annotation providers on a remote Linux node

Companion to [02c-m2-functional-annotation.md](02c-m2-functional-annotation.md). That doc says
what the annotation layer does and how it decides; this one is the operational sequence for
getting real provider output instead of the built-in heuristic, on a compute node with more
cores and more installed software than a laptop.

Nothing here is site-specific. Set two variables for your own node and the commands below work
as written:

```bash
REMOTE=user@host.example.org                     # ssh/rsync target
PROJECT=/path/on/remote/structure-to-function    # where the repo lives there
RUN=mygenome                                     # run id
CGA=mycga                                        # CGA directory name under data/
```

**Why this doc exists.** Issue #10 landed the code and every parser, and closed on the merge of
PR #37. But no provider has ever been run: on any real genome today, 100% of `flags.membrane`
and `flags.secreted` come from the built-in sequence heuristic, and `run.json` says so
(`functional_annotation.heuristic_only`). Turning that into `predicted` confidence is what the
steps below do.

## What has and has not been verified

Everything below ran on a real CGA directory (*M. genitalium*, 530 proteins) on a V100 node.

| | |
| --- | --- |
| M1 on a real CGA directory | **verified**, 0.2 s |
| M2 `--annotate --report` on that M1 output | **verified**, `proteins` section validates |
| Flags set for every protein | **verified** — 530/530 |
| Tool-ready FASTA emission | **verified** — `annotate_all.faa` 530 records, `annotate_selected.faa` 50 |
| Pass one against the network | **verified** — 337/530 proteins have a PDB hit, id mapping ran |
| InterProScan via `--interproscan auto` | **verified** — 435/530 proteins got Pfam/InterPro terms |
| DeepTMHMM, BioLib cloud | **verified** — 50 sequences in 76 s (1.32 seq/s) |
| Heuristic `membrane` vs DeepTMHMM | **verified** — precision 0.98, recall 0.94, MCC 0.95 over 530 proteins |
| Heuristic `signal_peptide` vs DeepTMHMM | **verified, and poor** — precision 0.42, recall 0.47, MCC 0.41 |
| DeepTMHMM, local install | **verified** on a V100 with `torch==1.5.0+cu101` — 530 sequences in **6m45s** |
| Two-pass ingest, real provider output | **verified** — 530/530 flags from a real predictor, `heuristic_only: 0` |
| Provider runtime on a full proteome | **measured** — see below |
| eggNOG-mapper, PSORTb, SignalP 6 | **not run** |

Pass two over real DeepTMHMM (all 530) and InterProScan output:

```
providers              {'deeptmhmm': 530, 'interproscan': 435, 'heuristic': 530}
unmatched_ids          {}
membrane_flag_source   {'deeptmhmm': 530}
confidence             {'predicted': 530}
heuristic_only         0
membrane 99 · secreted 19 · signal_peptide 36 · lipoprotein 11 · with_function_terms 435
```

**That is what a finished annotation run looks like**: every flag set by a real predictor, no
protein left on the fallback, nothing unmatched. The heuristic still ran for all 530 — it simply
lost to DeepTMHMM everywhere, which is the precedence table doing its job.

Two findings worth carrying forward. **The heuristic's two flags are of very different quality.**
Over the same 530 proteins its membrane call tracks DeepTMHMM closely — precision 0.98, recall
0.94, MCC 0.95, against 81.3% accuracy for an always-negative baseline — while its signal-peptide
call does not: precision 0.42, recall 0.47, MCC 0.41, and an accuracy *below* the always-negative
baseline. Treat a heuristic-tier `membrane` as usable and a heuristic-tier `signal_peptide`, or
the `secreted` derived from it, as not. Full numbers and caveats in
[02c](02c-m2-functional-annotation.md#checked-against-deeptmhmm). **And DeepTMHMM found helices
in only 2 of the 50 triage-selected proteins** against 99 across the proteome — the selection
criteria showing through, since selecting for PDB evidence and drug-target status selects against
membrane proteins, which rarely crystallize. A membrane penalty in the triage score would have
very little to do among the selected set.

### Measured runtimes

DeepTMHMM over the full 530-protein proteome, locally on one V100-SXM2-32GB:

| Step | Time | Rate |
| --- | --- | --- |
| 1. Load transformer model | — | one-off |
| 2. Generate ESM embeddings | 36 s | 14.5 seq/s |
| 3. Predict topologies | 5 m 45 s | 1.53 seq/s |
| **Total** | **6 m 45 s** | |

Step 3 is four-fifths of the cost. `predict.py` hardcodes `chunk_size = 1` at module level, so
topologies are predicted one sequence at a time even though the batching machinery is there —
`chunk_with_constraints` exists to isolate sequences longer than `max_length_for_batching = 5000`
from a batch, and the input is pre-sorted longest-first, which is how you prepare for batching.
Raising that constant is an operator-side change to third-party code, and the package is licensed
CC BY-NC-SA 4.0, so shipping an adapted copy would put a NonCommercial ShareAlike obligation on
this repository. Not worth it for a single genome: seven minutes is not the bottleneck. If you
are scanning many proteomes, ask BioLib to expose a `--batch-size` flag rather than forking.

For scale: the same 50-protein subset took 76 s through the BioLib cloud (1.32 seq/s), so a local
V100 is roughly at parity on topology and about 3x faster on embeddings.

**InterProScan finished far faster than its own benchmarks suggest** because *M. genitalium* is a
model organism: most of its proteins are already in EBI's pre-calculated match lookup, so the
scan is a lookup rather than an HMM run. A novel genome will not be.

The CGA proteomes tested so far contain **no `*` stop characters and no duplicate sequences**,
so the cleaning and deduplication in `write_tool_fasta` is insurance for the next genome rather
than a fix for these.

## What the node needs

| | |
| --- | --- |
| Python | **3.10 or newer**. `jsonschema>=4.20` will not install below 3.8, and `scripts/setup_env.sh` refuses to build below 3.10 |
| Java | 11 or newer, and only if InterProScan is to run |
| Network | outbound HTTPS for pass one (RCSB, UniProt, AlphaFold DB) and for InterProScan's pre-calculated lookup service. Without it, run everything `--offline` on the heuristic tier |
| Disk | the repo plus `runs/`; InterProScan itself is tens of GB if you install it fresh |

Versions are per node, so record what the node actually has rather than assuming: `setup_env.sh`
prints the Python it picked and the providers it found, and `run.json` records the InterProScan
path and version for the run (pitfall #19). Member database versions belong in the run notes
too — they change between InterProScan releases and the annotations are not comparable across
them.

## Three places, and which one does what

Getting this wrong wastes an afternoon, so it is worth stating plainly.

| Where | Role |
| --- | --- |
| **Your working copy** (laptop or workstation) | git working copy. Patches are applied here, branches are pushed from here, CGA data starts here |
| **GitHub** | how code reaches the node. Nothing else does |
| **The remote node** | where the providers run. `git clone` / `git pull` only — never `git am` |

**Code reaches the node through GitHub, not by copying files.** A patch file is a working-copy
step: apply it there, push the branch, then pull on the node. Never `scp` a patch to the node
and apply it there — the node then sits on a commit nobody else has, and the next `git pull`
conflicts with it.

**Data does not go through GitHub.** `data/` and `runs/` are gitignored, so CGA directories go
working copy → node by rsync, and results come back the same way.

## 1. Get the code onto the node

In your working copy, land whatever is outstanding and push it:

```bash
git switch main && git pull
git switch -c <branch>                  # skip if the work is already on main
git am --3way <patch>                   # working-copy side only
python -m pytest -q
git push -u origin <branch>
```

On the node, pull it:

```bash
ssh $REMOTE
cd "$PROJECT"
git clone https://github.com/NIAID-BRC-Codeathons/structure-to-function.git .   # first time only
git fetch origin
git branch -r | grep <branch>           # confirm the push landed before checking it out
git checkout <branch> && git pull

bash scripts/setup_env.sh
```

`pathspec '<branch>' did not match` means the branch is not on GitHub yet: the push above has
not happened or did not succeed. `scripts/setup_env.sh: No such file or directory` right after
it is the same cause, not a second problem — the script only exists on that branch.

The script picks the newest Python >= 3.10 it can find, builds `.venv/`, installs
`requirements-common.txt` + `requirements-m2.txt` + pytest, imports every module, runs the test
suite, and then prints which annotation providers exist on the node. It stops at the first
failure rather than reporting success over a broken install, and re-running reuses the venv.

Expected tail:

```
imports ok
173 passed
annotation providers on this node:
  interproscan  FOUND  /path/to/interproscan-<version>/interproscan.sh (<version>, via <route>)
  deeptmhmm     no     pip install pybiolib, then: biolib run DTU/DeepTMHMM --fasta <faa>
  ...
```

Missing providers are not an error — the layer falls back to the heuristic and records
`heuristic_only`. Discovery looks at `--interproscan-path`, `$INTERPROSCAN_HOME`,
`$INTERPROSCAN`, `PATH`, then `/opt`, `/usr/local`, `/software`, `/apps`, `/share/apps` and the
matching directories under `$HOME`. Two things worth knowing: `$INTERPROSCAN_HOME` set in
`~/.bashrc` is often **not** set in the non-interactive shell the script runs in, so the
directory search is what usually finds an existing install; and on a cluster `$HOME` and the
project directory can be different filesystems, so an install under `$HOME` is found while the
project lives elsewhere. If discovery misses it, pass `--interproscan-path`.

Then get a CGA directory onto the node. Run this from your working copy — `data/` is gitignored,
so it never travels through GitHub:

```bash
rsync -av --exclude __pycache__ data/$CGA "$REMOTE:$PROJECT/data/"
```

The CGA layout has a hidden `.annotation/` directory that M1 needs; rsync of the parent
directory includes it, but if you copy by hand, copy that too. Verify on the node before running
anything:

```bash
ls data/$CGA/annotation data/$CGA/.annotation/load_files/genome_feature.json
```

Results come back the same way when you want them locally:

```bash
rsync -av "$REMOTE:$PROJECT/runs/$RUN" runs/
```

## 2. Pass one — pipeline runs, FASTAs come out  *(on the remote node)*

```bash
source .venv/bin/activate
python -m s2f.run --run runs/$RUN --from-cga-dir data/$CGA --annotate --map-ids
```

This is the full M1 → M2 path with network: taxon resolution, RCSB sequence search per unique
sequence, UniProt/ChEMBL id mapping, AlphaFold DB lookups, and the annotation layer on the
heuristic tier. It writes:

| File | Use |
| --- | --- |
| `runs/$RUN/report.json` | the contract — `proteins[].flags` / `.annotations[]` |
| `runs/$RUN/m2_pdb/annotate_all.faa` | **cleaned, deduplicated FASTA of every protein** — feed this to the tools |
| `runs/$RUN/m2_pdb/annotate_selected.faa` | the top ~50 only — feed this to the slow tools |
| `runs/$RUN/m2_pdb/run.json` | `functional_annotation`: provider counts, `heuristic_only`, FASTA stats, InterProScan install record |

Sanity-check before spending tool time:

```bash
python -c "
import json; m=json.load(open('runs/$RUN/m2_pdb/run.json'))['functional_annotation']
print('annotated', m['proteins_annotated'], 'heuristic_only', m['heuristic_only'])
print('fasta', m['tool_fasta'])
print('notes', m['notes'])"
```

`tool_fasta.id_truncation_warning`, if present, matters: SignalP 6 and some other tools truncate
FASTA identifiers, and a truncated ID matches no `feature_id`. It does not corrupt anything —
those rows land in `run.json`'s `unmatched_ids` — but it means that provider contributed
nothing. Check `unmatched_ids` after every provider run.

## 3. Run the providers  *(on the remote node)*

Ranked by value per hour spent. Nothing here is required; each one that runs upgrades some
flags from `heuristic` to `predicted`.

### 3a. DeepTMHMM and SignalP 6 — the two required flags

These are what make `flags.membrane` and `flags.secreted` mean something. Run them over
`annotate_selected.faa` (50 proteins) first; the proteome can wait.

**Cloud (quickest to a first result).** Runs on BioLib's compute:

```bash
pip install pybiolib
cp runs/$RUN/m2_pdb/annotate_selected.faa selected.faa   # see note below
biolib run DTU/DeepTMHMM --fasta selected.faa
# writes biolib_results/TMRs.gff3
```

Copy the FASTA into the working directory first: passing a nested relative path fails in the
cloud job with `FileNotFoundError` on a file that plainly exists locally.

Note that this sends your sequences to a third party. For a published genome that is harmless;
for a blinded or embargoed one it is a decision somebody should make deliberately. The local
install below avoids the question entirely.

**Local.** See [Installing DeepTMHMM locally](#installing-deeptmhmm-locally) — the output is the
same `TMRs.gff3`, so nothing downstream changes.

SignalP 6 needs a one-time download after accepting the academic licence at
<https://services.healthtech.dtu.dk/services/SignalP-6.0/> (take the "fast" model):

```bash
pip install ./signalp-6-package/          # the tarball is a pip-installable package
signalp6 --fastafile runs/$RUN/m2_pdb/annotate_selected.faa \
         --output_dir runs/$RUN/signalp --organism other --format none --mode fast
```

`--organism other` is the bacterial setting; `--format none` skips the per-protein plots, which
are most of the runtime.

### 3b. InterProScan — function terms, already installed

```bash
python -m s2f.m2_triage --run runs/$RUN --annotate --report --offline \
    --interproscan auto --interproscan-appl Pfam
```

`auto` finds the install, scans `annotate_all.faa` (cleaned and deduplicated — InterProScan
rejects `*` and its cost scales with unique sequences), writes
`runs/$RUN/m2_pdb/interproscan.tsv`, and records path and version in `run.json` (pitfall #19). An
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
emapper.py -i runs/$RUN/m2_pdb/annotate_all.faa -o runs/$RUN/eggnog --cpu 16 --itype proteins
```

### 3d. PSORTb — localization

```bash
docker run --rm -v "$PWD":/data brinkmanlab/psortb_commandline:<tag> \
    -i /data/runs/$RUN/m2_pdb/annotate_all.faa -r /data/runs/$RUN -n -o terse
```

`-n` Gram-negative, `-p` Gram-positive, `-a` archaea. **Pick from M1's taxon call, not by
default:** running the Gram-negative model on a Gram-positive genome invents a periplasm.
*M. genitalium* has no outer membrane at all, so for the current test genomes PSORTb is the
least informative of the four.

## Installing DeepTMHMM locally

Worth doing: it removes the cloud round-trip, has no sequence limit, and lets the whole pipeline
run without sending anything to a third party.

The package is **request-only** — the app page at <https://dtu.biolib.com/DeepTMHMM> links a
short form under "Running DeepTMHMM Locally", and access has been granted immediately in
practice. Free for academic use; commercial users running it on their own servers need a licence
from BioLib. Roughly 1.7 GB zipped, 2.8 GB unpacked: `predict.py`, five cross-validation models,
and an embedded ESM-1b state dict. Keep it outside the repo.

```bash
cd <software dir>
unzip -q DeepTMHMM-Academic-License-v1.0.zip
cd DeepTMHMM-Academic-License-v1.0
```

### It needs its own Python 3.8 environment

The pinned dependencies are from 2022 and will not resolve on a modern Python. Build a separate
environment and keep it separate — the pipeline venv stays as `setup_env.sh` made it.

```bash
conda create -y -n deeptmhmm python=3.8
deactivate                     # LEAVE THE PIPELINE VENV FIRST
conda activate deeptmhmm
python -V                      # must say 3.8.x
which python                   # must be under .../envs/deeptmhmm/bin
```

**An active venv shadows a conda environment.** `conda activate` sets `CONDA_PREFIX` but the
venv's `bin` still wins on `PATH`, so `pip` silently installs into the venv instead. The symptom
is `cp311` wheels in the pip output and a torch resolution failure; the damage is packages in
the wrong environment. Check `which python` before installing anything.

### Installing, around the torch pin

```bash
python -m pip install wheel Cython==0.29.37 pkgconfig==1.5.5
python -m pip install torch==<build> -f https://download.pytorch.org/whl/torch_stable.html
grep -v '^torch==' requirements.txt > requirements-notorch.txt
python -m pip install -r requirements-notorch.txt
python -c "import torch, esm, h5py, Bio; print(torch.__version__, torch.cuda.is_available())"
```

`requirements.txt` pins `torch==1.5.0+cu92`, a local version that exists only on PyTorch's own
index and not on PyPI, so installing the requirements file as-is always fails. Install torch
first from the PyTorch index, then the rest with that line stripped.

**Choosing `<build>` decides whether you get the GPU.** torch 1.5.0 ships exactly three:
`1.5.0+cpu`, `1.5.0+cu92`, `1.5.0+cu101` — there is no cu102 at this version.

| Your GPU | Build | Why |
| --- | --- | --- |
| Volta (V100, sm_70) or older | `1.5.0+cu101` | CUDA 10.1 supports up to sm_70; needs driver ≥ 418.39 (`+cu92` needs only ≥ 396) |
| Ampere (A100, sm_80) or newer | `1.5.0+cpu` | no CUDA build at this version reaches sm_80. A newer torch (≤ 2.4.1 on Python 3.8) might work, but is untested against this package's 2022-era code |
| None | `1.5.0+cpu` | `predict.py` branches on `torch.cuda.is_available()`, so the CPU build makes that deterministic |

Confirm with `nvidia-smi --query-gpu=name,driver_version --format=csv` before choosing.

### Running it

```bash
time python3 predict.py --fasta sample.fasta --output-dir result1        # their sample first
time python3 predict.py --fasta <project>/runs/$RUN/m2_pdb/annotate_all.faa \
                        --output-dir <project>/runs/$RUN/deeptmhmm_local
```

**`--output-dir` must not already exist** — `predict.py` errors rather than overwriting, so use a
fresh name per run. Output is `TMRs.gff3`, `predicted_topologies.3line` and
`deeptmhmm_results.md`; the gff3 is identical in format to the cloud's, so `--deeptmhmm` ingests
it unchanged.

The two environments meet only through that file: conda `deeptmhmm` runs the predictor, the
pipeline venv runs the pipeline. Do not merge them — torch 1.5 and its 2022 pins in the pipeline
venv would be much harder to undo than to keep apart.

## 4. Pass two — ingest what the tools produced  *(on the remote node)*

```bash
python -m s2f.m2_triage --run runs/$RUN --annotate --report --offline \
    --deeptmhmm   biolib_results/TMRs.gff3 \
    --signalp     runs/$RUN/signalp/prediction_results.txt \
    --interproscan runs/$RUN/m2_pdb/interproscan.tsv \
    --eggnog      runs/$RUN/eggnog.emapper.annotations \
    --psortb      runs/$RUN/psortb_terse.txt
```

`--offline` replays the PDB/UniProt work from pass one's cache, so this is seconds rather than
minutes. Pass only the files that exist; a missing one is a note, not a failure.

Then check what actually landed:

```bash
python -c "
import json; m=json.load(open('runs/$RUN/m2_pdb/run.json'))['functional_annotation']
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

For each run, four things belong in the run notes, because none of them can be reconstructed
afterwards from the outputs alone:

- **Wall-clock per provider**, and which hardware. "DeepTMHMM 6m45s" means nothing without
  "530 sequences, one V100".
- **Tool and database versions.** `run.json` captures the InterProScan path and version
  automatically; the member database versions and the DeepTMHMM package version do not capture
  themselves.
- **The `confidence` and `membrane_flag_source` histograms.** These are what let the report say
  how much of the annotation is a real prediction and how much is the fallback heuristic.
- **Whether any sequences left the building.** The BioLib cloud path sends them to a third
  party; the local install does not. For a blinded or embargoed genome that distinction belongs
  in the methods section, not in someone's memory.

## Failure modes worth knowing before you start

| Symptom | Cause | What to do |
| --- | --- | --- |
| `interproscan auto: no interproscan.sh found` | `$INTERPROSCAN_HOME` not set in a non-interactive shell | export it, or `--interproscan-path /path/to/interproscan.sh` |
| Provider passed, `providers` shows 0 for it | headers rewritten or truncated by the tool | check `unmatched_ids`; rename the tool's output IDs back to `feature_id` |
| `parsed to zero records — check the file format` in `notes` | the parser did not recognise the file | confirm you passed the right file (`TMRs.gff3`, `prediction_results.txt`, `*.emapper.annotations`, terse PSORTb, InterProScan TSV/JSON) |
| InterProScan rejects the input | stop characters | already handled — but confirm you passed `annotate_all.faa`, not `m1/proteins.faa` |
| `pip install jsonschema` fails | Python < 3.8 | `scripts/setup_env.sh` refuses to build below 3.10 for exactly this reason |
| M2 writes no `report.json` | `--report` missing | the runner always passes it; a direct `python -m s2f.m2_triage` call does not |
| `git pull` on the node conflicts | a patch was applied on the node instead of on the workstation | `git reset --hard origin/<branch>`; apply patches where you hold the working copy, then push |
| M1 fails with a missing-file error | the hidden `.annotation/` directory did not come across | re-rsync the CGA parent directory, then check `ls data/<cga>/.annotation/load_files/` |

## Related docs

- [02c-m2-functional-annotation.md](02c-m2-functional-annotation.md) — precedence, flag definitions, calibration
- [00a-data-contract.md](00a-data-contract.md) — `report.json` and the schema
- [pitfalls.md](pitfalls.md) — #3 membrane proteins, #19 reproducibility

# Structure-to-Function

**NIAID-BRCs AI Codeathon 2.0** · September 16-18, 2026 · Argonne National Laboratory

Takes a bacterial assembly and produces ranked, evidence-carrying functional
annotations: genes and specialty genes from BV-BRC, structural and binding
evidence, and a report in which every claim names its source.

Project page: https://niaid-brc-codeathons.github.io/projects/structure-to-function/

## Status

| Stage | Module | State |
| --- | --- | --- |
| M1 genome | `s2f.m1_genome` | taxon call, CGA submit/poll/parse, quality gate, BV-BRC Data API route, priority ranking, interactive report |
| M2 triage | `s2f.m2_triage` | PDB evidence, AlphaFold, knowledge graph, triage scoring |
| M3 fold | - | not started |
| M4 ligands | - | not started |
| M5 disease | - | not started |
| M6 report | - | not started |
| shared | `s2f.common` | `report.json` schema, atomic section writer, HTTP cache, id mapping |

`python -m s2f.run` runs the stages that exist, in order, into one
`runs/<id>/report.json`.

## Requirements

- **Python 3.10 or newer**, the floor `scripts/setup_env.sh` enforces. `jsonschema>=4.20`
  will not install below 3.8; the setup script stops at 3.10 rather than build an
  environment nobody has tested. Nothing in `s2f/` uses syntax or stdlib newer than that.
- `curl`, for the BV-BRC services that reject some HTTP client user agents.
- Optional, for M1: `scipy` and `numpy` for the gene-content phylogeny, `matplotlib`
  for publication figures. Without them the run succeeds and says what it skipped.
- For M1's CGA route only: the **BV-BRC command-line tools**
  (<https://www.bv-brc.org/docs/cli_tutorial/index.html>) and a BV-BRC account.
  Everything else runs offline from fixtures.

## Install

```bash
git clone https://github.com/NIAID-BRC-Codeathons/structure-to-function.git
cd structure-to-function
python3 -m venv .venv && source .venv/bin/activate   # python3 must be 3.10 or newer
pip install -r requirements-common.txt -r requirements-m1.txt -r requirements-m2.txt
```

### Notes

The code assumes nothing about the host beyond the requirements above: no shared
filesystem, no module system, no site-specific paths. If `python3` is older than 3.10, a
`module load python` or a conda environment provides a newer one. `scripts/setup_env.sh`
does the whole install: it picks the newest Python >= 3.10 on `PATH`, builds `.venv/`,
installs the requirements, imports every module and runs the suite, stopping at the first
failure rather than reporting success over a broken environment.

For M2's optional annotation providers (InterProScan, DeepTMHMM, SignalP, PSORTb,
eggNOG-mapper), see `docs/02c-m2-functional-annotation.md`. None are required, and
discovery is by `$INTERPROSCAN_HOME`, `PATH` and the usual install prefixes rather
than any hard-coded location.

## Check it works, without credentials or network

```bash
# M1 parses a trimmed real CGA output committed as a fixture
python -m s2f.m1_genome --run runs/demo --from-cga-dir fixtures/m1/cga_sample --html
# -> runs/demo/m1/report.html, openable offline

# M1 without a BV-BRC account, straight from the public Data API (needs network).
# NOTE: this characterises a reference genome for the organism, not your assembly
# (docs/01b-m1-api-mode.md). Re-run with --offline afterwards to replay from cache.
python -m s2f.m1_genome --run runs/kp --from-bvbrc-api --genome-id 1125630.4 --html

# M2 from its own fixtures, in its own run directory
python -m s2f.m2_triage --run runs/m2demo --dry-run --report

pytest -q
```

Use a **separate run directory** for `--dry-run`: M2's fixtures are a synthetic
genome, and merging them into a run M1 populated fails on unknown feature IDs (by
design - `feature_id` is the join key everywhere).

## Run the pipeline

Your workspace path is built from your BV-BRC account, and the suffix differs
between accounts (`@bvbrc` on some, `@patricbrc.org` on others). Find yours, then
create the directory the run will write to:

```bash
p3-login YOUR_BVBRC_USERNAME          # writes ~/.patric_token
un=$(tr '|' '\n' < ~/.patric_token | sed -n 's/^un=//p')   # your account, token not printed
p3-ls "/$un/home"                     # confirm the prefix
p3-mkdir "/$un/home/s2f" "/$un/home/s2f/uploads"
```

Then:

```bash
python -m s2f.run \
    --run runs/mgen_G37 \
    --contigs fixtures/genomes/mgen_G37/mgen_G37.fna \
    --ws-dir "/$un/home/s2f" \
    --limit 50
```

What M1 does, in order:

1. **Blinds the contig headers** to `contig_1`, `contig_2`, ...
2. **Uploads** the FASTA to your BV-BRC workspace.
3. **Calls the taxon** with Similar Genome Finder (the `Minhash` JSON-RPC service),
   choosing the rank the Mash distance justifies: species at <= 0.05, genus at
   <= 0.20, and refusing to guess past that. The genetic code comes from the called
   taxon's own BV-BRC record.
4. **Submits CGA** with that taxon and code, polls, and retrieves the output.
5. **Parses** it into `proteins[]`, `genome`, the `run` manifest, and the
   `runs/<id>/m1/` tables M2 reads.
6. **Gates on quality**: a `Poor` genome or any quality flag stops the run unless
   you pass `--allow-poor`.

Step 3 is not optional ceremony. Annotating this project's test genome with no
taxon (so genetic code 11 instead of 4) doubled the CDS count, halved mean protein
length and produced a `Poor` genome — `docs/01a-cga-coverage.md` has the numbers.

Useful flags: `--from-cga-dir DIR` reuses a job you already retrieved, `--job-id ID`
resumes polling a submitted job, `--only m1` / `--skip m2` select stages,
`--map-ids` and `--kg` turn on M2's id resolution and knowledge graph.

`--ani` turns the Mash distances in `genome.closest_genomes[]` into real identities: it
downloads each close genome from BV-BRC over HTTPS and runs skani against the assembly.
It needs the assembly FASTA — a retrieved CGA directory does not contain it — so pass
`--contigs` or `--ani-query` alongside `--from-cga-dir`. See
[docs/01-m1-genome.md](docs/01-m1-genome.md#ani-to-the-closest-genomes).

## Output

```
runs/<id>/
  report.json          every module's section; the single state file
  m1/
    proteins.faa       one record per CDS, headers are fig| feature IDs
    genes_proteins.csv  specialty_genes_all.csv  taxon_call.json
    ani.json           with --ani: skani version, argv, params, per-genome provenance
    genomes/           with --ani: the close genomes' contigs, named <genome_id>.fna
  cga/                 the retrieved CGA output
  cache/               downloads and API responses; what --offline replays
  m2_pdb/              proteins.tsv, hits.tsv, top50.tsv, no_pdb_hit.tsv, run.json
```

`report.json` is the contract between modules: each owns exactly one top-level key
and never edits another's. Writes go through `s2f.common.io.update_section`, which
is locked and atomic, so two modules can write different sections at once.
`s2f.common.schema.validate` checks a section before it is written.

## Documentation

`docs/00-architecture.md` is the contract; read it first in any coding thread.
Module docs are `01-` through `06-`, `01a-cga-coverage.md` records what CGA does and
does not give us, `pitfalls.md` collects the ways this goes wrong, and
`07-decisions-and-risks.md` holds the settled parameters.

Tasks live in GitHub issues, generated from `docs/issues.md` by
`scripts/create_issues.py` (idempotent; skips titles that already exist). Issues are
unassigned on purpose - comment to claim one.

## Development

```bash
pytest -q                    # 115 tests, no network required
```

Schema changes bump `SCHEMA_VERSION` and update the fixtures in the same commit.

---

# Original project pitch

> The proposal below is what the organizing team circulated to seed the project.
> It is kept for the record; scope and methods have moved since.

## Goal (proposed)

Automate functional interpretation of poorly characterized pathogen proteins by combining sequence, structure, binding databases, and literature evidence.

## Three-Day MVP (proposed)

Analyze 50–100 hypothetical proteins from a selected pathogen group such as *Chlamydiales*. The agent should retrieve sequences from BV-BRC, run similarity and structure tools, identify structural neighbors, search PDB/BindingDB/ChEMBL evidence, extract supporting statements from papers, and generate ranked functional annotations.

AutoPDB components can create a provenance-aware training dataset of protein–ligand or protein–protein interactions.

## Model and evaluation (proposed)

Train or calibrate an embedding-based function classifier or evidence-ranking model. Compare recommendations with curated annotations or held-out known proteins.

## Leads

- Jiahui Chen
- Moises Gualapuro

Team assignments are still being finalized. Participants can review their project, and request a reassignment, in the participant spreadsheet circulated by the organizing team.

## Setup

Python 3.10+. Python dependencies:

```bash
pip install -r requirements-common.txt   # jsonschema, requests — needed by s2f/common
pip install -r requirements-m2.txt       # adds networkx for the M2 knowledge subgraph
```

### External tools

Not pip-installable, so record anything you add here with its version and what needs it.

| Tool | Version | Install | Needed by |
| --- | --- | --- | --- |
| DIAMOND | 2.2.7 | `brew install diamond` (or `conda install -c bioconda diamond`) | M2 human-homology search against the human proteome (#29). 542 × 20,600 proteins in seconds, where BLASTP takes minutes. |
| skani | 0.3.2 | `conda install -c bioconda skani`, `cargo install skani`, or the single static binary from [the releases page](https://github.com/bluenote-1577/skani/releases/latest) | M1 `--ani`: ANI between the assembly and its closest genomes (#7). Six 580 kb genomes in under a second; accurate on fragmented drafts, where fastANI is not. |

Everything else the pipeline uses is a web API reached through `s2f/common/http.py`, which caches
responses, so a rerun costs nothing and `--offline` replays a run with no network at all.

Tests need no network and no credentials: `python -m pytest`.

## Working here

This repository is the team's working space for the codeathon — code, notebooks, data pointers, and notes. Replace this README with the real thing once the charter is written. Team members get access through the [NIAID-BRC-Codeathons](https://github.com/NIAID-BRC-Codeathons) organization; accept the invitation if you have not already.

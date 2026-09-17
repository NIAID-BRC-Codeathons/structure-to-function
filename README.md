# Structure-to-Function: Hypothetical Protein and Binding-Evidence Factory

**NIAID-BRCs AI Codeathon 2.0** · September 16–18, 2026 · Argonne National Laboratory

Automating functional interpretation of poorly characterized pathogen proteins by combining sequence, structure, binding databases, and literature evidence.

Project page: https://niaid-brc-codeathons.github.io/projects/structure-to-function/

---

> **This is a draft pitch, not a plan.**
>
> What follows is a one-slide proposal from the organizing team. It exists
> to seed a team, not to constrain one. Scope, methods, target organism,
> and success criteria are all still open — expect them to change
> substantially. Turning this into a real plan is the team's first job, and
> it lands in the project charter due August 28, 2026.

---

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

Python 3.11+. Python dependencies:

```bash
pip install -r requirements-common.txt   # jsonschema, requests — needed by s2f/common
pip install -r requirements-m2.txt       # adds networkx for the M2 knowledge subgraph
```

### External tools

Not pip-installable, so record anything you add here with its version and what needs it.

| Tool | Version | Install | Needed by |
| --- | --- | --- | --- |
| DIAMOND | 2.2.7 | `brew install diamond` (or `conda install -c bioconda diamond`) | M2 human-homology search against the human proteome (#29). 542 × 20,600 proteins in seconds, where BLASTP takes minutes. |

Everything else the pipeline uses is a web API reached through `s2f/common/http.py`, which caches
responses, so a rerun costs nothing and `--offline` replays a run with no network at all.

Tests need no network and no credentials: `python -m pytest`.

## Working here

This repository is the team's working space for the codeathon — code, notebooks, data pointers, and notes. Replace this README with the real thing once the charter is written. Team members get access through the [NIAID-BRC-Codeathons](https://github.com/NIAID-BRC-Codeathons) organization; accept the invitation if you have not already.

# The `report.json` contract — how to read and write it

Implements the data contract in [00-architecture.md](00-architecture.md). Issue
[#2](https://github.com/NIAID-BRC-Codeathons/structure-to-function/issues/2).
Code: `s2f/common/schema.py` (shapes) and `s2f/common/io.py` (writing).

## Writing a section

```python
from s2f.common.io import init_report, update_section, update_proteins

init_report(run_dir, run_id, input_file=str(contigs))   # once, writes `run`
update_section(run_dir, "structures", structures)        # your own key only
update_proteins(run_dir, {feature_id: {"flags": {...}}}) # enrichment onto M1's list
```

`update_section` validates first and writes nothing if validation fails, so a bad section can
never reach the file. Never write `report.json` directly.

## Why not `json.dump`

Two modules finishing at the same moment both read the file, each replaces its own key, and the
second writer saves a copy of what it read *before* the first writer's change — silently deleting
a section. So a write is: take an exclusive `flock`, read, replace one key, write a temp file in
the same directory, `os.replace` (atomic on POSIX).

This is measured, not assumed. With six processes writing different sections at the same
instant: **all six survive with the lock; only one survives without it** (3 runs each way). A
test reproduces it with real subprocesses.

## Sections and owners

| Key | Owner | Required fields |
| --- | --- | --- |
| `run` | M1 | `run_id`, `created_at` |
| `genome` | M1 | — |
| `proteins` | M1 writes the list; M2 enriches | `feature_id` per record |
| `proteins[].m1_priority` | M1 | — (see below) |
| `proteins[].triage` | M2 | `score`, `components`, `selected`, `reason` |
| `structures` | M3 | `feature_id`, `source` (`pdb`/`afdb`/`predicted`/`esm_atlas`) |
| `kg` | M2 | every edge: `source_id`, `target_id`, `type`, `provenance` |
| `ligands` | M4 | `ligand_id`, `provenance` |
| `docking` | M4 | `feature_id`, `ligand_id` |
| `analyses` | M3/M4/M5 | `method` |
| `disease` | M5 | `claim` per claim |
| `report` | M6 | — |

## Permissive now, tighter later

Each section requires only what makes a record **identifiable and attributable**; everything else
is optional and unknown properties are allowed. Ship a field before its shape is settled — the
schema will not block you. The shapes tighten once M3–M6 exist and we can see what they really
produce. Tightening bumps `SCHEMA_VERSION` and updates the fixtures in the same commit
(pitfall #18).

Three things are enforced from day one, because they are what makes the report defensible:

1. **`feature_id` is the canonical key** everywhere (pitfall #11).
2. **Anything retrieved externally carries `source` and `retrieved_at`.**
3. **Anything transferred from another protein carries the identity that justified it** (pitfall #4).

## Two scores, two owners, one protein

`proteins[].m1_priority` and `proteins[].triage` are different keys on purpose and neither
module may write the other's.

| | `m1_priority` (M1) | `triage` (M2) |
| --- | --- | --- |
| Ranks by | host-interaction evidence in the annotation | structural tractability |
| Available | before any structure is looked at | after PDB and AlphaFold lookups |
| Weights | virulence +3, drug target +3, essential/AMR/mechanism +2, surface/named/transporter +1, human homolog −2 | `pdb_evidence` 0.40, AlphaFold 0.00, … |
| Shape | `score`, `rank`, `categories`, `mechanism_hypothesis`, `breakdown`, `selected`, `weights` | `score`, `components`, `selected`, `reason` |

**They are expected to disagree.** A virulence factor with no solved homolog ranks high in
M1 and low in M2; that is the pipeline working, not a bug. Collapsing them into one
`score` would either discard the annotation evidence M1 has and M2 does not, or overwrite
the structural evidence M3 gates on.

`m1_priority` carries its full `breakdown` — the points and the sentence that earned them —
and the `weights` it ran with, so M2 can re-rank with its own weights instead of inheriting
M1's, and so a reader can audit a ranking rather than trusting a number. Code:
`s2f/m1_genome/priority.py`.

## Two shapes worth arguing about

**`annotations[]` requires `source`.** M2 sequence hits, Foldseek (#9) and eggNOG (#10) all write
here. Without `source`, a TM-score row and a sequence-identity row are indistinguishable.

**Structure evidence has two independent axes**, kept apart in `flags`:

- `experimental_homolog` — PDB hit: identity, coverage, resolution, `holo` (has a real ligand).
- `predicted_model` — AlphaFold: `mean_plddt`, `confidence_band`, **covered residue range** and
  `coverage`, `usable_for_docking` with its `reason`.

A single `structure_available` boolean would erase the distinction M3 gates on. From the G37 run:
422 usable predicted models, but only 14 of the top 50 have a ligand-bound experimental homolog —
and a predicted model has no cofactors or ligands at all (pitfall #2).

**The two axes are not equal, and the order is deliberate: experimental evidence outranks a
prediction.** Where both exist, use `experimental_homolog`; fall back to `predicted_model` only
when there is no usable experimental one. Consumers should follow this order rather than deciding
per module:

1. **Experimental homolog, holo** (`experimental_homolog.holo == true`) — real coordinates with a
   real ligand in the site. Best available: the site is observed, not inferred, and M3 can
   transfer it by superposition.
2. **Experimental homolog, apo** — real coordinates, no ligand. The fold and the pocket geometry
   are observed; what binds there is not.
3. **Predicted model, usable** (`predicted_model.usable_for_docking == true`) — a plausible fold
   with **no cofactors, metals or ligands** (pitfall #2). Docking into it is the weakest link in
   the whole pipeline, so pair it with a holo template where one exists and flag it where none
   does.
4. **Nothing usable** — say so, and let it fail the `usable_for_docking` gate with a reason rather
   than docking into something unusable.

The triage score already reflects this: `pdb_evidence` carries weight 0.40 while the AlphaFold
result carries **none**. A protein does not rank higher for having a prediction; the prediction
only tells M3 whether there is something to work with once the protein is already selected.

## The fixture

`fixtures/report.fixture.json` — 5 proteins, 1 structure, 3 ligands, 2 docking rows, 1 disease
claim. The `run`, `genome`, `proteins` and `kg` sections are **trimmed from a real G37 run**; the
M3–M6 sections are illustrative until those modules exist and can replace them with real output.

Build against it before M1 has produced anything: that is what `--dry-run` is for.

## Open: which G37 is canonical?

Two fixtures on `main` describe the same organism with **different canonical keys**:

| Fixture | Genome | CDS | Feature IDs |
| --- | --- | --- | --- |
| `fixtures/genomes/mgen_G37/` (#25) | 243273.**25** | 542 | `fig\|243273.25.peg.N` |
| `fixtures/cga/mgen_G37/` (#28) | 243273.**147** | 530 | `fig\|243273.147.peg.N` |

`feature_id` is what every record joins on, so these two cannot both be the reference. The
M2 fixture above uses 243273.25 because that is what M2 has run on; the real CGA output is
243273.147. **M1 should decide which one the pipeline uses**, and the fixture should be rebuilt
from it. This is exactly the drift the contract exists to prevent, so it is worth settling early.

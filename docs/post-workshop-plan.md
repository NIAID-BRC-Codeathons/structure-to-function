# Post-workshop plan

**29 September 2026.** The codeathon ran 16–18 September. Five of the eight participants are
continuing, most of them available from **mid-October**.

| Continuing | GitHub | Worked on |
| --- | --- | --- |
| **Jiahui Chen** (lead) | `Jiahuic` | M2 triage, `common/` contract, M6 report |
| Carla Mann | `cmmann21` | M1, M2 functional annotation, most of the documentation |
| Alexander Partin | `adpartin` | M3 structure collection |
| Ping Hu | `PingPnG` | M1 scripts, the docking and pocket prototype |
| Ashita Jawali | `Ashita2619` | ESMC embeddings |

Not continuing: Moises Gualapuro (co-lead), MD Shakhaowat Hossain, Heng Ma, Archit Vasan.
**There is no longer a co-lead**, so sequencing decisions land in one place.

Five people working evenings is the constraint every proposal below is measured against.

This is a proposal, not a decision. It exists so the restart argues about a written thing.

## Where the project actually is

M1 (genome), M2 (triage) and M3 (structure collection) run end to end: a full pass over
*M. genitalium* G37 takes about 5½ minutes — 542 CDS in, 50 proteins selected, 50 structures
collected, `report.json` validating against the schema, and a self-contained HTML report out.
418 tests pass with no network.

**M4 (ligands and docking) and M5 (disease evidence) were never built.** M4 is the
*binding evidence* in "Hypothetical Protein and Binding-Evidence Factory", so the project's
headline deliverable is half-built. Ping's `bacteria_structure_docking.py` does much of M4 and
M3's pocket work, but it sits outside the pipeline with its own ranking and its own output
directories.

Three pull requests have been open since 18 September: **#62** (evaluation harness, ~7,900
lines), **#76** (declares M6's `jinja2` dependency), **#77** (M1 species metadata). 36 issues are
open.

### What we have not shown

The pipeline has never run on a genome whose **relatives** rather than **itself** are in the PDB.
On G37, 47 of the 50 selected proteins rest on a same-genus structural hit and 30 on the same
species. The triage weights are uncalibrated (#61), and two circularity defects surfaced late:
the self-match artifact (#12) and an essentiality reference set that can include the query
genome itself (#68).

Everything below assumes that gap is the most important thing to close, in any direction.

## The shape of the code, and why generalizing is cheaper than it looks

About 2,500 lines have **no BV-BRC coupling at all**:

| Portable today | Lines |
| --- | --- |
| `m2_triage/foldseek.py`, `pdb_evidence.py`, `afdb.py` | 972 |
| `m3_fold/` — structure collection and quality gates | 923 |
| `common/io.py`, `common/schema.py`, `common/http.py` | ~500 |
| `m6_report/`, `full_report/` | ~460 |

BV-BRC coupling is concentrated in **M1** (24 references in `collect.py` alone, plus
`bvbrc_api.py`, `parse.py`, `priority.py`, `taxon.py`) and in **three M2 files**:
`essentiality.py`, `bvbrc_input.py`, `taxonomy.py`.

So the project is really two things stacked:

1. a **generic protein-evidence factory** — given protein sequences, find structures, homologs,
   ligandability and knowledge-graph links, and report all of it with provenance; and
2. an **antibacterial target-selection policy** — the weights that say human homology is a
   penalty, essentiality is a bonus, and membrane proteins are deprioritized.

Generalizing the pipeline means making (2) pluggable. It is not a rewrite of (1).

## Direction A — finish what exists

Sparse-time work, suitable for everyone. In rough order of value:

- [ ] **Decide the three stale PRs.** #76 is a dependency declaration and should just land; #77
      is M1's, and Carla is continuing. **#62 is the awkward one**: it is the only *measurement*
      work anyone did, it is ~7,900 lines, and its author is not continuing — so merging it
      means owning it. Its most important finding is now filed separately as **#85** (M1's
      `--seq-cap` removes 58% of hypothetical proteins before M2 sees them), so the finding
      survives whichever way the PR goes. Closing it with thanks is a reasonable outcome;
      leaving it open indefinitely is not.
- [ ] **The defect backlog** — #67, #68, #69, #72, #80, #81. Each is one sitting. #68 matters
      most: essentiality transferred from a reference set that can contain the query genome is
      circular, and it is in a number we have quoted.
- [ ] **#14, pockets and receptor preparation** (p0, unstarted). Ping's script already produces
      prepared receptors, pocket boxes and Boltz-2 inputs. Adapting it is cheaper than starting.
- [ ] **#22, CI.** 418 tests exist and nothing runs them automatically. Everything merged so far
      was verified by whoever remembered to run `pytest`.
- [ ] **Close #21** (the presentation) — it is finished and historical.

## Direction B — the biggest missing piece

**M4: ligand candidates, scoring and docking** (#15, #16, #17). This is the one gap that changes
what the pipeline *claims to be*, and it needs sustained attention rather than spare evenings.

**It is also the best-staffed gap we have**, which was not obvious until the team shrank:
Ping's `bacteria_structure_docking.py` already produces prepared receptors, pocket boxes and
Boltz-2 inputs; Alex holds #18 (Boltz-2 co-folding) and built M3's structure collection, which
is what M4 consumes; and Ashita holds #15 (ligand candidate collection). Three of the five
continuing people are already pointed at this.

The work is less "build M4" than "bring Ping's prototype inside the contract" — it currently has
its own ranking and its own output directories, and reads none of `report.json`.

M5 (#19, #20) is the same shape but one step further out; it can wait.

## Direction C — generalize, starting with a document

**Write the seam down before refactoring anything.** One document that answers: where does the
generic evidence factory end and the domain policy begin, and what does a domain have to supply?

Use **breast cancer** as the test case, because it inverts the policy in ways that expose whether
the seam is real:

| Component | Bacterial pathogen | Human cancer target |
| --- | --- | --- |
| `human_homolog_penalty` (−0.25) | avoid human similarity — toxicity risk | **the target is human**; becomes a paralog-selectivity question instead |
| `essential` | BV-BRC flux-balance prediction | **DepMap CRISPR screens** — stronger evidence than we have now |
| Open Targets | human-centric, reachable only through a homolog | **native**; this is what it is built for |
| `surface_bonus` | vaccine and antibody angle | antibody and ADC targets — still applies |
| ChEMBL ligandability | applies | applies, with far more data |
| Input | genome assembly → BV-BRC annotation | a gene or protein list from expression or mutation data |

The parts that were *weakest* for bacteria are *strongest* here, which is the best evidence that
the split is genuine rather than wishful.

**Food-science peptides are a poorer fit, and worth saying so early.** Short peptides (5–20
residues) defeat the evidence modules specifically: PDB sequence search at that length returns
noise against our 1e-5 gate, AlphaFold and ESMFold are unreliable, and Foldseek needs a fold that
does not exist. What transfers is the *architecture* — the report contract, the provenance
discipline, the transparent weighted score, the rule that "not run" is not "zero" — not the
evidence modules. That is still worth something, but it is a different kind of reuse.

**CEPI belongs here too.** A teammate noted that services their group is developing could feed
several pipeline steps. CEPI (Coalition for Epidemic Preparedness Innovations) funds vaccine and
countermeasure work, including epitope prediction and immunogen design, and was among the
collaborators on the AlphaFold Database's viral protein complex release of 24 September 2026.
Their public work is largely viral and vaccine-facing while ours has only ever run on bacteria,
so the fit needs establishing rather than assuming. It is the same question as the seam: can a
third party's evidence service be plugged in without touching the core?

## Sequencing

A and C are not in competition. The refactor that makes the policy layer pluggable is also what
lets us run a **null-policy control** on the bacterial case — which is exactly what calibration
(#61) needs. Doing a slice of C serves A.

What should not happen is building B and C at once on an unvalidated core.

## Open questions for the restart

1. Do we finish the bacterial case first, or is a second domain more motivating for people
   working evenings? With five people this is a real trade, not a both-at-once.
2. **M4**: does Ping's prototype come inside the contract, and who drives it — Ping, Alex, or
   the two together? This is the highest-value thing on the board and the only one with three
   people already adjacent to it.
3. What do we do with #62? Merging ~7,900 lines whose author has left means owning them.
   The key finding is preserved in #85 either way.
4. Which CEPI services, concretely, and do they have a callable API?
5. Do we need a genome that is **not** in the PDB before any of this means anything? Nothing
   we have measured so far survives that question unchanged.

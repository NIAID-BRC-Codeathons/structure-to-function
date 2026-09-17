# Pitfalls and Mitigations

Known failure modes for this pipeline, with the mitigation each one requires. Every item here is enforced somewhere in code or in the report — the "where" column says which module owns it.

## Scientific validity

| # | Pitfall | Mitigation | Owner |
| --- | --- | --- | --- |
| 1 | **Docking on predicted apo structures is the weakest link.** Blind docking over a whole model produces numbers that look quantitative and mean little. | Defined pockets only (P2Rank, or a holo template's site transferred by superposition). Rank *within* a protein, never across proteins. Include positive and decoy controls for every docked protein. | M3, M4 |
| 2 | **Predicted structures lack cofactors and metals.** Many bacterial targets need Mg²⁺, Zn²⁺, heme or a nucleotide in the site; without them the pocket is wrong. | Prefer targets whose closest structural neighbour has a holo structure to align to. Flag models whose homolog site requires a cofactor. | M3 |
| 3 | **Membrane proteins** are common among surface and efflux targets and unreliable for both folding and docking. | Membrane flag from DeepTMHMM in triage; deprioritize explicitly rather than silently including them. | M2, M3 |
| 4 | **The human-homolog paradox.** A close human homolog is useful for borrowing ligand scaffolds and disqualifying for an antibacterial target (toxicity). | Record homology identity on every knowledge-graph edge; treat a close human homolog as a scoring penalty and a stated selectivity risk. | M2, M6 |
| 5 | **Docking scores are not affinities.** | Present as a ranking device with its caveat attached to the number, not only in a limitations section. Boltz-2 affinity estimates are recorded as analyses, not truth. | M4, M6 |
| 6 | **"Approved drug" is not "safe for this."** Dose, route, tissue penetration and spectrum all matter. | Candidates are framed as computational hypotheses for wet-lab triage. No clinical or treatment claims anywhere in the output. | M4, M5, M6 |
| 7 | **Unsupported LLM claims.** Fluent mechanism narratives with no basis are the default failure mode. | Every claim carries feature IDs and/or PMCIDs; a separate verification pass drops or downgrades unsupported claims and the dropped count is reported. | M5, M6 |
| 8 | **Inference presented as measurement.** BacDive describes *relatives*, not our isolate; pathway gaps suggest but do not establish auxotrophy. | Label each growth or nutrient statement with its source type (strain record vs. genomic inference). | M5 |
| 9 | **Validation asymmetry.** Steps 1–2 can be scored against a blinded known genome; steps 4–5 cannot be validated in three days. | Say so plainly in the report rather than implying validation. | M6 |

## Engineering and scope

| # | Pitfall | Mitigation | Owner |
| --- | --- | --- | --- |
| 10 | **Knowledge-graph scope creep.** Building or ingesting a full KG will consume the whole codeathon. | Query existing APIs, assemble a small subgraph with provenance, cap it. If edges exceed a few thousand, write `kg.json` separately and keep a summary in `report.json`. | M2 |
| 11 | **Entity resolution is where the time actually goes** (locus tag ↔ UniProt ↔ gene name ↔ ChEMBL target). | One shared ID-mapping helper (`common/ids.py`) built in the first hour; no module writes its own mapping. | common |
| 12 | **Triage bias.** Both triage steps decide the outcome; tuning them after seeing results invalidates the result. | Write criteria and weights down before running. Keep discarded proteins and ligands with reasons. Record any post-hoc change to the weights. | M2, M4 |
| 13 | **Compute budget.** 50 folds is cheap; 50 proteins × hundreds of ligands is not, and Boltz-2 co-folding is minutes per pair. | Benchmark one protein before any batch run. Set the ligand cap (start: 200 in, top 20 docked) and record it in the run manifest. | M3, M4 |
| 14 | **ESM Atlas dependency.** Endpoints and limits have changed before, and metagenomic hits are largely uncharacterized. | Test on day 1; treat as a bonus path, never a required step. | M2 |
| 15 | **CGA runtime is unknown to us.** | Submit on day 1 and record the duration. The `p3-`CLI fallback path (existing BV-BRC annotation for a known genome) exists for this and doubles as the fixture source. | M1 |
| 16 | **Promiscuous compounds** (detergents, dyes, PAINS) dominate docking rankings. | RDKit property and PAINS filters before docking. | M4 |
| 17 | **Cross-module blocking.** Five people waiting on one module means nothing ships. | Fixtures plus `--dry-run` for every module, shipped before real data exists. | common, all |
| 18 | **Schema drift.** Improvised changes to `report.json` break other threads silently. | JSON Schema validation in CI; schema changes are announced and bump `schema_version`, with fixtures updated in the same commit. | common |
| 19 | **Reproducibility.** | Pin database versions and access dates, docking seeds and exhaustiveness, model versions and prompt versions. | all |
| 20 | **Licensing.** DrugBank, CARD and DisGeNET restrict redistribution. | Interactive use only; prefer ChEMBL, PubChem, openFDA, UniProt, Open Targets and STRING for anything shipped. | M4, M5 |
| 21 | **Runaway agent loops** on Argo are expensive and attributed to a named account. | Hard iteration cap, explicit `max_tokens`, prototype on `claudehaiku45`. | M5 |
| 22 | **BV-BRC output is not guaranteed UTF-8.** A real USA300 `sp_gene.json` carried a 0xa0 byte (Latin-1 non-breaking space), and a strict `read_text(encoding="utf-8")` turns one stray codepoint in a metadata field into a hard failure of the whole module. | Every reader of a CGA directory goes through `parse._read_text`, which falls back utf-8 → utf-8-sig → cp1252 → latin-1 before replacing. | M1 |
| 23 | **`p3-cp` exits 0 without overwriting.** Two runs uploading the same basename leave the first genome in the workspace, so Minhash calls *its* taxon and CGA annotates with a taxon from the wrong organism — which #15 and `01a-cga-coverage.md` measure as roughly doubling the CDS count and halving mean protein length. Nothing downstream can detect it. | Two guards, both needed: `cga.upload` verifies the workspace copy by size (detects a stale one), and M1 names the blinded contigs for the run id (prevents the collision). `upload` takes no `name` parameter — `p3-cp` writes under the *local* basename, so a differing `name` returns a path with no file at it. `tests/test_m1_upload_naming.py` covers the second guard and the `name` trap. | M1 |

## The single most likely way this produces a wrong headline

Comparing docking scores across different proteins and calling the highest one the best candidate. Enforced against in M4 (scores carry their protein context) and M6 (the report groups candidates by protein and never ranks across them).

## Related docs

- [pipeline.md](pipeline.md) — tools and APIs per step
- [00-architecture.md](00-architecture.md) — data contract and shared conventions
- [07-decisions-and-risks.md](07-decisions-and-risks.md) — open decisions and the parallel-thread plan

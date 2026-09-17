# Open decisions and standing risks

## Decisions (recorded 2026-09-16, issue #24)

| Decision | Chosen | Note |
| --- | --- | --- |
| Smoke-test genome | *M. genitalium* G37, BV-BRC 243273.25, 580,076 bp, 542 CDS | `fixtures/genomes/mgen_G37/`. Pipeline shakedown and fixtures only, never a reported result. Human pathogen with AMR, human homologs and a TTD drug target; no virulence-factor hits, so that triage component is untested by it. |
| Test genome | *S. aureus* USA300 (MRSA), blinded | Default taken. Headers stripped, name never in prompts or filenames. |
| Triage target count | 50 proteins | Default taken. Pipeline proven on 10 first. |
| Ligand cap per protein | 200 candidates in, top 20 docked | Default taken. Caps recorded in the `run` manifest. |
| Docking engine | Vina, with gnina rescoring on the top poses | Default taken. Fixed seed and exhaustiveness, both recorded. |
| Argo model | haiku for prototyping, opus for final | Default taken. Cap `max_tokens`, cap loop iterations. |

Options considered for each are in the issue. Changes from here are noted in the
commit message that makes them.

### Still open

| Decision | Options | Status |
| --- | --- | --- |
| Module owners | one per M1-M6 | Claimed: M2 sequence/structure lookups and triage scoring (#8, #12, Jiahuic), M2 Foldseek (#9, Ashita2619). M1, M3, M4, M5, M6 have no owner. |

## Standing risks

1. **Cross-protein score comparison** — the most likely way this produces a confidently wrong headline. Enforced in M4 and M6.
2. **Predicted apo structures without cofactors** — pockets may be incomplete. Prefer holo templates; flag the rest.
3. **Human-homolog paradox** — a close human homolog helps find ligands and disqualifies the target. Record both.
4. **Knowledge-graph scope** — query existing APIs, assemble a small subgraph, never ingest a full KG. ID mapping is the real work.
5. **ESM Atlas availability** — endpoints and limits have changed before; test day 1, treat as bonus.
6. **Unsupported LLM claims** — the verify pass in M5 and the evidence rule everywhere.
7. **CGA runtime unknown** — M1's `p3-`CLI fallback path exists for this; build it early, it doubles as the fixture source.
8. **Licensing** — ChEMBL, PubChem, openFDA, UniProt, Open Targets, STRING are shippable. DrugBank, CARD, DisGeNET are not; use interactively only.
9. **Compute budget** — benchmark one protein before any batch fold or dock; Boltz-2 is minutes per pair.
10. **Validation asymmetry** — M1 and M2 can be scored against the blinded genome's ground truth; M4 and M5 cannot be validated in three days. The report says so plainly.

## Non-negotiables

- No claim without evidence IDs.
- No clinical or treatment recommendation.
- Every external record carries source, version/date and retrieval time.
- Fixtures and `--dry-run` for every module, so nobody is blocked on anyone else.
- Schema changes are announced, not improvised.

## Parallel-thread plan

| Thread | Module | Depends on | Can start |
| --- | --- | --- | --- |
| 1 | `common/` + schema + fixtures | — | first, and quickly |
| 2 | M1 genome | common | immediately after |
| 3 | M2 triage | common + fixtures | immediately (fixtures) |
| 4 | M3 fold | common + fixtures | immediately (fixtures) |
| 5 | M4 ligands | common + fixtures | immediately (fixtures) |
| 6 | M5 disease | common + fixtures | immediately (fixtures) |
| 7 | M6 report | common + fixtures | immediately (fixtures) |

Thread 1 is the bottleneck for everyone, so it ships the schema and fixtures before anything else — ideally in the first hour.

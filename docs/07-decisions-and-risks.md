# Open decisions and standing risks

## Decisions needed (owner: team lead, today)

| Decision | Options | Default if nobody decides |
| --- | --- | --- |
| Test genome | *S. aureus* USA300 (MRSA), *K. pneumoniae* HS11286, *P. aeruginosa* PAO1 | USA300, blinded (headers stripped, name never in prompts or filenames) |
| Triage target count | 20 / 50 / 100 proteins | 50, with the pipeline proven on 10 first |
| Ligand cap per protein | candidates in / docked | 200 in, top 20 docked |
| Docking engine | Vina, smina, gnina rescoring | Vina + gnina rescoring on the top poses |
| Module owners | one per M1–M6 | assign before any code is written; every module needs a name |
| Argo model | haiku for prototyping, opus for final | as stated; cap `max_tokens`, cap loop iterations |

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

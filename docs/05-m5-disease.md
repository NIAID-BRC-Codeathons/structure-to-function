# M5 — Disease analysis and mitigation strategies

**Scope:** explain how this organism causes disease, what it needs to grow, and which non-drug strategies are plausible — every statement tied to evidence.

## Inputs

`proteins[]` with flags and annotations (M1/M2), `kg` (M2), `ligands`/`docking` summary (M4), `genome.closest_genomes` (M1).

## Outputs

- `disease.claims[]`: text, `kind` (virulence | growth | nutrient | mechanism | mitigation), `evidence_feature_ids[]`, `evidence_pmcids[]`, `confidence`, `contradicts[]` if sources disagree.
- `disease.strategies[]`: strategy, rationale, evidence, feasibility note.
- `disease.growth`: temperature, oxygen, pH, media candidates, auxotrophy predictions — each with source.

## Tools and APIs

| Need | Pick |
| --- | --- |
| Virulence basis | **BV-BRC specialty genes** (VFDB, Victors), PHI-base for phenotype evidence, secretion systems from TXSScan |
| Growth conditions | **BacDive API** (temperature, oxygen, pH, media for close relatives), **MediaDive API** (recipes) |
| Nutrient needs | **GapMind** (amino-acid and vitamin biosynthesis gaps), BV-BRC subsystems for pathway completeness; cobrapy FBA on a gapseq model if anyone has one ready |
| Literature | **RAGStack** `/v1/retrieve` (passages + PMCIDs), Europe PMC as backup |
| Non-drug strategies | Auxotrophy → nutrient restriction; surface/secreted proteins → vaccine or antibody targets (BepiPred, DiscoTope); anti-virulence and biofilm targets; phage candidates (PhagesDB, INPHARED); host-directed options via Open Targets |
| Synthesis | **Argo** — prototype on `claudehaiku45`, final on `claudeopus5`; tool calling on the non-streaming path only; `max_tokens` ≤ ~21k; hard iteration cap |

## Steps

1. `gather_evidence()` — deterministic lookups first: specialty genes, BacDive/MediaDive, GapMind, subsystems. No LLM involved.
2. `retrieve_literature()` — RAGStack per gene and per organism; cache passages with PMCIDs.
3. `synthesize()` — agentic pass over the gathered evidence with thin tools (RAGStack, ChEMBL, UniProt, Open Targets). The prompt requires every claim to cite feature IDs or PMCIDs.
4. `verify()` — second pass, separate call: check each claim against its cited evidence and drop or downgrade unsupported ones. Keep the dropped list for the report appendix.
5. `strategies()` — propose mitigations, each tied to a specific finding (an auxotrophy, a surface protein, an AMR profile).

## Acceptance checks

- Zero claims without evidence in the final `disease` section; the count of dropped claims is reported.
- Growth statements name their source (BacDive strain record vs. inference from pathway gaps) and never present inference as measurement.
- Deterministic lookups produce the same result on re-run; LLM passes record model, prompt version and temperature.

## Pitfalls

- The LLM will happily write a fluent mechanism narrative with no support. The verify pass is not optional.
- BacDive describes *relatives*, not our isolate. Phrase accordingly.
- Growth and nutrient claims cannot be tested here; they are inferences and must be labeled.
- No clinical or treatment recommendations — this is a research hypothesis document.
- Keep the agent loop bounded; a runaway Argo loop is expensive and gets noticed.

## Kickoff prompt

> Read `00-architecture.md` and `05-m5-disease.md`. Implement `s2f/m5_disease`: deterministic evidence gathering (specialty genes, BacDive, MediaDive, GapMind, subsystems), RAGStack retrieval, an Argo synthesis pass whose claims must cite feature IDs or PMCIDs, a separate verification pass that drops unsupported claims, and the mitigation-strategy generator. Record model, prompt version and dropped claims. `--dry-run` from fixtures with a stubbed LLM.

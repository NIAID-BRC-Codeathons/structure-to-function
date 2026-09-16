# Pipeline Reference — Tools and APIs

Project 9, NIAID-BRCs AI Codeathon 2.0 (Argonne, Sept 16–18, 2026).

This document is the **what to use** reference for each pipeline step. For implementation contracts (data schema, module boundaries, acceptance checks) see the module docs listed in [Related docs](#related-docs); those are authoritative when the two disagree.

**Bold** marks the choice recommended for the codeathon. Everything listed is open or free with registration unless noted otherwise.

## Contents

- [1. Genome analysis](#1-genome-analysis)
- [2. Protein triage](#2-protein-triage)
- [3. Structure prediction](#3-structure-prediction)
- [4. Ligand triage and docking](#4-ligand-triage-and-docking)
- [5. Disease analysis](#5-disease-analysis)
- [6. Report](#6-report)
- [Access and licensing](#access-and-licensing)
- [Related docs](#related-docs)

## 1. Genome analysis

Input: unknown bacterial genome assembly (contig FASTA). Output: coding genes and proteins, AMR calls, phylogeny, closest relatives.

| Output | Source | Notes |
| --- | --- | --- |
| Coding genes / proteins | **BV-BRC Comprehensive Genome Analysis (RASTtk)** | CDS, products, subsystems, protein FASTA, JSON genome object |
| AMR | **CGA specialty genes** (CARD, NCBI reference gene catalog) | Offline cross-check: AMRFinderPlus, or RGI against CARD |
| Phylogenetics | **CGA tree** (Newick, NEXUS) | Alternative: BV-BRC Phylogenetic Tree service (codon trees) |
| Related species / genomes | **BV-BRC Similar Genome Finder** (Mash/MinHash) | Run *before* CGA — annotation needs a taxon at genus level or below |
| Identity to relatives | **skani** (or fastANI) | CGA reports Mash distance, not ANI |
| SNP distance | **Snippy `--ctgs`** or **Parsnp** | BV-BRC Variation Analysis requires reads, so it does not apply to an assembly |
| Viral input (if scope widens) | BV-BRC annotation with VIGOR4 or LowVan; Subspecies Classification service | Bacterial gene callers are wrong for viruses |

Access: `p3-login`, then `p3-submit-CGA`, `p3-get-genome-features` and related BV-BRC CLI commands.

## 2. Protein triage

Input: proteins from step 1. Output: annotated, ranked proteins; ~50 selected for structural and ligand work.

### Lookups

| Need | Tools / APIs |
| --- | --- |
| ID mapping | **UniProt ID-mapping + REST API** — one shared helper, never per-module mappings |
| Sequence homologs | **RCSB Search API** (MMseqs2-based sequence search), DIAMOND against Swiss-Prot locally |
| Experimental structures | **RCSB Data API**, PDBe REST — record whether the entry has a bound ligand |
| Predicted structures | **AlphaFold DB API** by UniProt accession |
| Structure search | **Foldseek** local binary against PDB + AlphaFold DB (+ ESM Atlas); the Foldseek web server also exposes a ticket-based REST API |
| ESM Metagenome Atlas | `esmatlas.com` fold and search endpoints — verify availability, length caps and rate limits on day 1; hits are largely uncharacterized |
| Functional annotation | **eggNOG-mapper** (COG, GO, KEGG KO, EC), InterProScan (heavier), PSORTb + SignalP 6 (surface / secreted), DeepTMHMM (membrane flag) |

### Knowledge graph

| Need | Tools / APIs |
| --- | --- |
| Target–disease–drug edges | **Open Targets Platform GraphQL API** (human-centric; reached through our homolog) |
| Functional networks | **STRING API** (includes bacterial proteomes) |
| Drug–gene edges | DGIdb API, Pharos / TCRD |
| Host–pathogen edges | HPIDB, PHI-base |
| Protein interactions | IntAct, BioGRID REST |
| Pathways | Reactome REST, KEGG REST |
| Federated queries | NCATS Biomedical Data Translator (ARAX / TRAPI), Monarch API |
| Assembly | **NetworkX** in-process; export nodes and edges with provenance |

### Selection

Weighted, written-down score. Suggested components: virulence or AMR hit; essentiality by orthology (DEG / Geptop-style); surface or secreted; ligandable homolog (a known ligand exists); structure obtainable; **penalty** for a close human homolog; **penalty** for predicted membrane protein. Keep the full ranked list with reasons, not only the selected set.

## 3. Structure prediction

Input: selected proteins. Output: a usable structure and at least one defined pocket per protein, or an explicit "not dockable".

| Need | Tools |
| --- | --- |
| Existing structures | **Reuse PDB and AlphaFold DB hits from step 2 first** — fold only what has nothing usable |
| Prediction | **BV-BRC Protein Structure Prediction** (Boltz-2, OpenFold 3, Chai-1, ESMFold, ESMFold2, AlphaFold 2; web, CLI, API); ESMFold or Boltz-2 locally on ANL GPUs for throughput |
| Complexes | Same service with two chains (Boltz-2, Chai-1); AlphaFold-Multimer long term |
| Receptor prep | PDBFixer, PDB2PQR (protonation), Meeko or OpenBabel (PDBQT) |
| Pockets | **P2Rank** and/or fpocket; US-align or TM-align to transfer a holo template's ligand site |
| Confidence | pLDDT per residue, PAE for single chains, ipTM for complexes |
| Optional biophysics | ProDy normal modes (flexibility), APBS + PDB2PQR (electrostatic surface) |

## 4. Ligand triage and docking

Input: selected proteins with pockets, plus knowledge-graph edges. Output: ranked ligand candidates per protein with docking poses and controls.

| Sub-step | Tools / APIs |
| --- | --- |
| Approved / usable drugs | **ChEMBL API** (`max_phase=4`), **openFDA** label and approval endpoints, Broad Drug Repurposing Hub, ZINC in-stock subsets. DrugBank interactively only — its license restricts redistribution |
| Activity data | **ChEMBL** bioactivities (pChEMBL), **BindingDB**, **PubChem PUG-REST** BioAssay |
| Knowledge-graph ligands | **Open Targets** known drugs for the mapped target, DGIdb |
| Literature candidates | **RAGStack** `/v1/retrieve` per target and compound class; Europe PMC as backup |
| Chemistry handling | **RDKit** — canonical SMILES, InChIKey dedup, descriptors, property filters, PAINS flags |
| Per-protein score | **Transparent weighted score**, not a learned model: homology identity to the ligand's known target, measured potency, approval status, property filters, literature support, knowledge-graph edge type |
| Docking | **AutoDock Vina** (or smina) for breadth, gnina for CNN rescoring, **Boltz-2 co-folding with affinity estimate** for the top ~10 pairs, DiffDock as an alternative pose generator |
| Controls | Per protein: one known binder (positive) plus decoys or random ligands (negative) |

Compare scores **within one protein only**. See [pitfalls.md](pitfalls.md).

## 5. Disease analysis

Input: everything above. Output: mechanism, growth and nutrient needs, non-drug mitigation strategies — each tied to evidence.

| Sub-step | Tools / APIs |
| --- | --- |
| Virulence basis | **CGA specialty genes** (VFDB, Victors), PHI-base phenotypes, TXSScan / MacSyFinder for secretion systems |
| Growth conditions | **BacDive API** (temperature, oxygen, pH, media for close relatives), **MediaDive API** (recipes) |
| Nutrient needs | **GapMind** (amino-acid and vitamin biosynthesis gaps), BV-BRC subsystems for pathway completeness; cobrapy FBA on a gapseq or ModelSEED model (long term) |
| Literature | **RAGStack** `/v1/retrieve`; Europe PMC as backup |
| Non-drug strategies | Auxotrophy → nutrient restriction; surface and secreted proteins → vaccine or antibody targets (BepiPred, DiscoTope); anti-virulence and biofilm targets; phage candidates (PhagesDB, INPHARED); host-directed options via Open Targets |
| Essentiality (best targets) | Orthology to known essential genes; BV-BRC Tn-Seq Analysis service where data exist |
| Expression | BV-BRC RNA-Seq Analysis service, or Expression Import for published datasets |
| Synthesis | **Argo** — `claudehaiku45` while prototyping, `claudeopus5` for final runs. Tool calling works only on the non-streaming path; pass `max_tokens` ≤ ~21,000; cap agent loop iterations |

Evidence rule: every generated claim carries feature IDs and/or PMCIDs, and a separate verification pass drops unsupported claims.

## 6. Report

| Need | Tools |
| --- | --- |
| Data contract | One `report.json` per run — see [00-architecture.md](00-architecture.md) |
| Report | **Python + Jinja2 → single self-contained HTML file** (no server, no external calls at view time) |
| Structures and poses | **3Dmol.js** inline — receptor, pocket, top pose |
| Ligand depictions | RDKit to SVG |
| Trees | Phylocanvas, or a static image from ETE3 |
| Reproducibility | Record database versions and access dates, model versions, docking seeds and exhaustiveness, ligand caps |

## Access and licensing

| Resource | Note |
| --- | --- |
| BV-BRC services and CLI | Accounts for all team members; `p3-login` token also authenticates RAGStack |
| RAGStack | `X-API-Key` header, or a raw BV-BRC token in `Authorization` (no Bearer prefix); never both |
| Argo | Argonne collaborator username is the credential; no API key |
| BacDive / MediaDive | Free registration for API access |
| ChEMBL, PubChem, openFDA, UniProt, Open Targets, STRING, PDB, AlphaFold DB | Open; cite versions and access dates |
| DrugBank, CARD, DisGeNET | Redistribution restricted — interactive use only, check terms before shipping derived data |
| Bakta, GTDB-Tk, InterProScan | Large reference databases; not laptop-friendly without planning |

Credentials live in environment variables (`ARGO_USER`, `RAGSTACK_API_KEY`) or `~/.patric_token`, never in the repository.

## Related docs

| Doc | Contents |
| --- | --- |
| [00-architecture.md](00-architecture.md) | Data contract, repo layout, ID and provenance rules |
| [01-m1-genome.md](01-m1-genome.md) | Module 1 — genome analysis |
| [02-m2-triage.md](02-m2-triage.md) | Module 2 — protein triage and knowledge graph |
| [03-m3-fold.md](03-m3-fold.md) | Module 3 — structures and pockets |
| [04-m4-ligands.md](04-m4-ligands.md) | Module 4 — ligands and docking |
| [05-m5-disease.md](05-m5-disease.md) | Module 5 — disease analysis |
| [06-m6-report.md](06-m6-report.md) | Module 6 — report |
| [07-decisions-and-risks.md](07-decisions-and-risks.md) | Open decisions, thread plan, non-negotiables |
| [pitfalls.md](pitfalls.md) | Known failure modes and mitigations |

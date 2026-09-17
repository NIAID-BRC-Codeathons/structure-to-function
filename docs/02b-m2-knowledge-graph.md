# M2b — Knowledge subgraph

Companion to [02-m2-triage.md](02-m2-triage.md) for M2 step 5 (`build_kg()`). Issue
[#11](https://github.com/NIAID-BRC-Codeathons/structure-to-function/issues/11). Implemented in
`s2f/m2_triage/kg.py`; run with `python -m s2f.m2_triage --run runs/<id> --kg --taxon <id>`.

**Nothing is ingested or hosted** (pitfall #10). Public APIs are queried for the selected
proteins only, results are capped, and what we keep is the subgraph plus the evidence for it.

## Source survey (checked 2026-09-17, against our two genomes)

| Source | Usable for our bacteria? | Evidence |
| --- | --- | --- |
| **STRING** | **Yes, directly** | Accepts BV-BRC locus tags as-is: `gyrA` → `1125630.KPHS_37060`, `MG_401` → `243273.MG_401`. No mapping needed. |
| **ChEMBL** | Through a homolog | Only 25 targets for *K. pneumoniae*, but they are KPC and SHV-1. *E. coli* gyrase B has 339 measured activities — the value is in the homolog. |
| KEGG | Yes, not yet wired | `kpn:KPN_02467` resolves; needs a taxon → KEGG organism-code map. |
| PHI-base | Site up, not yet wired | Bacterial-specific phenotypes. |
| IntAct | API up, thin | Returned nothing for the protein tried. |
| Open Targets | **Human only** | Searching "gyrase" returns TOP2A. Reachable only through a human homolog. |
| HPIDB | **Unreachable** | Connection refused. Not a blocker. |

## The homolog bridge

Drug and disease knowledge sits on characterized proteins in *other* organisms, so most edges
are reached through a homolog. That hop stays explicit and carries the identity that justifies
it — never collapsed into "our protein is a drug target" (pitfall #4).

A real chain from the G37 run:

```
brc:fig|243273.25.peg.405  Isoleucyl-tRNA synthetase
  --homolog_of-->   uniprot:P41972   identity 37.6%, coverage 99.9%, via PDB hit 1FFY_2
  --represented_by--> chembl:CHEMBL1982  (Isoleucine--tRNA ligase, S. aureus)
  <--target_of--   compound:CHEMBL538163  IC50 300 nM, pChEMBL 6.52
```

Isoleucyl-tRNA synthetase is a genuine antibacterial target (mupirocin's), **and** 37.6% identity
is weak. Both facts travel together on the edge; a reader who sees only the compound has been
misled.

## Edge types

| Type | From → To | Evidence carried |
| --- | --- | --- |
| `interacts_with` | our protein → STRING partner | combined/experimental/database/textmining scores, with a note that a STRING score is aggregated evidence, not a measured interaction |
| `homolog_of` | our protein → UniProt protein | percent identity, coverage, organism, and how it was found |
| `represented_by` | UniProt accession → ChEMBL target | the accession, verified to be a component of that target |
| `target_of` | compound → ChEMBL target | assay type, value, units, pChEMBL, assay ID |

Human homologs from BV-BRC's specialty table get a `homolog_of` edge to a node that is
explicitly **unaccessioned** (BV-BRC gives a name and identity, not an accession). The identity
is recorded because M2 scoring penalizes it and M6 must state the selectivity risk.

## Caps

Recorded in `kg.json` so a reader knows the graph is bounded by choice, not by what happened to
exist: 10 STRING partners per protein at score ≥ 700, 5 ChEMBL targets per accession, 10
compounds per target, and only the selected proteins.

G37, top 50: **270 nodes, 336 edges** — 238 `interacts_with`, 51 `homolog_of`, 41 `target_of`,
6 `represented_by`. Small enough to ship inside `report.json` if wanted; it is written to
`kg.json` with a summary in `run.json`.

## A trap worth knowing

**ChEMBL silently ignores unknown filter names.** `target_organism__icontains=Klebsiella
pneumoniae` returned 18,552 targets — exactly the unfiltered count — full of human proteins like
*Maltase-glucoamylase*. The working filter is `organism__icontains` (25 results), and for the
homolog bridge, `target_components__accession`, which does filter correctly.

Because a wrong filter looks like a successful query, every ChEMBL target we keep is re-checked:
the accession we asked about must appear in the target's own component list, or the target is
dropped and the mismatch recorded in `failures`. A test covers this.

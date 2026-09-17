# `common/ids.py` — identifier mapping

Issue [#3](https://github.com/NIAID-BRC-Codeathons/structure-to-function/issues/3). The single
place that maps BV-BRC `feature_id` to UniProt, UniParc, RefSeq, PDB, ChEMBL and gene name.
No module writes its own mapping (pitfall #11), because entity resolution is where this project
will lose time and where cross-module bugs come from.

## Use it

```python
from s2f.common.http import CachedJsonClient, JsonCache
from s2f.common.ids import IdMapper

with JsonCache("runs/<id>/cache/ids.sqlite") as cache:
    mapper = IdMapper(CachedJsonClient(cache=cache), taxon_id=243273)
    xrefs = mapper.map_protein(feature_id, sequence, locus_tag="MG_401")
```

`map_many([...], workers=4)` takes a list of those keyword sets and returns `{feature_id: Xrefs}`.
In M2: `python -m s2f.m2_triage --run runs/<id> --map-ids --taxon 243273`.

## Routes, strongest first

`Xrefs.route` always records which one produced the mapping, so a weak mapping is never mistaken
for a strong one.

| Route | How | Notes |
| --- | --- | --- |
| `proteome_checksum` | One bulk fetch of the genome's taxon (`uniprotkb/stream`), then an exact CRC64 match in memory | The fast path. Entries come back with the xref fields already attached, so a matched protein needs **no** further request. |
| `proteome_locus` | Same bulk index, matched on ordered locus name | For proteins whose sequence differs from UniProt's (annotation drift) but whose locus is unambiguous. |
| `sequence_crc64` | The same checksum, looked up in UniParc | For proteins the taxon index does not cover. UniProt's checksum is CRC64-ISO reflected (poly `0xD800000000000000`) — verified against `sequence.crc64` in a test, since the common ECMA-182 variant gives a different value. |
| `patric_crossref` | UniParc indexes BV-BRC/PATRIC feature IDs: `database:PATRIC AND dbid:"fig\|..."` | Catches proteins whose sequence has drifted from the UniProt record. Quote the ID and scope by database — an unquoted feature ID matches loosely. |
| `locus_tag` | UniProt search by locus tag within the taxon, **verified** against the entry's ordered locus names | UniProt writes `MG401` where BV-BRC writes `MG_401`, so the comparison ignores punctuation. UniProt's `gene:` field does not match ordered locus names, so this is a free-text search plus verification — the verification is what makes it safe. |
| `pdb_hit` | UniProt accessions the caller already holds from PDB sequence hits | Recorded in `related_uniprot`, **never** in `uniprot`: a homolog structure's accession is not our protein's identity. |
| `unmapped` | nothing matched | `uniprot` is `None` and `note` says why. An explicit gap beats a plausible wrong answer. |

## What you get

`Xrefs`: `uniprot`, `uniprot_entry_name`, `uniparc`, `gene_name`, `protein_name`, `taxon_id`,
`refseq[]`, `pdb[]`, `chembl[]`, `related_uniprot[]`, `route`, `note`, and `mapped`.
`as_dict()` serializes it with `mapped` included.

Every response is cached on disk, so a second run makes no network calls and `offline=True`
replays a run with none at all.

## Cost

Pass `taxon_id` — it is the difference between one request and thousands.

| | Requests | Measured |
| --- | --- | --- |
| With `taxon_id` | 1 bulk fetch, plus fallbacks only for proteins it misses | *M. genitalium* G37: 542 proteins in **1.8 min**, 88% mapped |
| Without | 2–4 per protein | the same genome was on track for **~100 min** |

The bulk fetch is one request regardless of size: 483 entries for G37 in 3 s, 5,728 for
*K. pneumoniae* HS11286 in 131 s. Everything is cached, so a second run makes no network calls;
point modules at the same SQLite file to share it.

## What "unmapped" means

On G37, 66 of 542 proteins do not map. Their sequences **are** in UniParc, but the UniProtKB
entries are no longer active — BV-BRC's annotation includes ORFs that UniProt has retired or
never had. That is a real gap in the reference data, not a lookup bug, and it is recorded as
`unmapped` with the reason in `note` rather than resolved to a near-neighbour.

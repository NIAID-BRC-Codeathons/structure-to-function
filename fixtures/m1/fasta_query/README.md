Input for M1's `--from-fasta` route tests (#89). Three header styles, three kinds of protein:

| Record | Header style | What it exercises |
| --- | --- | --- |
| `gyrA` | plain id + product + `[organism]` bracket | a characterised bacterial protein |
| `tr|TEST00001|TEST1_KLEPN` | UniProt-style with `OS=` tail; **made-up test accession** | an uncharacterised protein, and a trailing `*` stop |
| `sp|P0DTC2|SPIKE_SARS2` | UniProt-style | a viral protein |

The sequences are **short fragments for testing**, not full-length entries, so exact-sequence
UniProt mapping (`--map-ids`) is not expected to match them. Neither the bracket nor the `OS=`
organism may reach `m1/proteins.faa`.

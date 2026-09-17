# `fixtures/m1/bvbrc_api/` — canned BV-BRC Data API responses

`klebsiella_hs11286.json` is one bundle of trimmed **real** BV-BRC Data API responses for
*Klebsiella pneumoniae* subsp. *pneumoniae* HS11286 (`genome_id 1125630.4`), the genome the
API route was developed against. 13 KB, so it can live in git; the run it came from was
~4 MB, which cannot.

It backs `tests/test_m1_bvbrc_api.py`, which serves these payloads through a fake
`requests` session. No test in this repo touches the network.

## Where it deviates from the live API, and why that matters

BV-BRC omits unset fields from a record rather than returning them as null, so a real
`genome` response carries only the fields that have a value: 68 for this genome, and
`genetic_code` is **not** among them.

This fixture's `genome` record nonetheless carries `cell_shape`, `collection_year`,
`disease`, `gram_stain`, `host_name`, `isolation_country`, `oxygen_requirement` and
`publication`, which the live record for `1125630.4` does not return. They are kept
because they exercise the "this genome" branch of `collect_growth` and its neighbours --
but note the consequence, because it is the opposite of what you want: the
*species-facet fallback*, which is the branch that actually runs against the live API for
every one of those fields, is the branch this fixture does **not** cover.

`genetic_code` was moved rather than kept. It has no fallback, so with it present in the
`genome` record the tests passed while every live run wrote `genetic_code: null`. It now
sits in the `taxonomy` record, which is where BV-BRC really returns it and where the CGA
route reads it from (`taxon.py`).

## Keys

| Key | Serves | Core |
| --- | --- | --- |
| `genome` | the resolved genome record | `genome` |
| `features` | 9 CDS features | `genome_feature` |
| `specialty` | 7 specialty-gene rows (4 VFDB virulence, 3 AMR) | `sp_gene` |
| `amr_phenotypes` | 8 laboratory/computational AMR records | `genome_amr` |
| `neighbors` | 3 reference/representative genus neighbours | `genome` |
| `taxonomy` | the lineage with ranks, and `genetic_code` | `taxonomy` |
| `subsystems` | 4 metabolic subsystems | `subsystem` |
| `close_pathogens_facet` | 6 human-associated *Klebsiella* species | `genome` facet |
| `sequences` | one stub sequence per feature md5 | `feature_sequence` |

## Why these nine features

The slice is chosen so every branch of `priority.score_protein` fires on real annotation
text rather than on invented strings:

* `peg.997` / `peg.998` — CFA/I fimbrial usher and adhesin: **adhesion**, VFDB-classified
  `Adherence`, so the regex and the classification path agree.
* `peg.1937` / `peg.3382` — aerobactin and yersiniabactin receptors: **iron**, and
  `peg.3382` also matches `outermembrane`, exercising the surface bonus.
* `peg.1591` — BasR/PmrA two-component regulator: **regulation**, and an AMR row whose
  `classification` is the "regulator modulating expression of antibiotic resistance
  genes" string, which is deliberately *not* a mechanism keyword.
* `peg.5465` — Tet(G) MFS efflux pump: AMR with no mechanism category.
* `peg.994` — EcpR transcriptional regulator: a VF with no product-text mechanism.
* `peg.5518` / `peg.5521` — hypothetical proteins: the zero-score, no-hypothesis path.

`aa_sequence_md5` values are sequential stubs, not real md5s, and the sequences are a
repeating stub string. The tests assert on plumbing (which feature got which sequence),
never on sequence content, so a real proteome would add megabytes and no coverage.

## Refreshing it

Re-run the API route against the same genome and re-trim, keeping the nine features above
so the scoring assertions stay meaningful:

    python -m s2f.m1_genome --run runs/kp --from-bvbrc-api --genome-id 1125630.4

BV-BRC recurates, so `patric_id`s and counts can move. If they do, the fixture and the
tests change in the same commit (pitfall #18).

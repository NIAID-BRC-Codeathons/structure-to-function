# M1 — the BV-BRC Data API route

Implements a second acquisition route for M1 alongside the CGA service in
[01-m1-genome.md](01-m1-genome.md). Code: `s2f/m1_genome/bvbrc_api.py` (client and genome
resolution) and `s2f/m1_genome/collect.py` (collectors and the contract adapters).

    python -m s2f.m1_genome --run runs/kp --from-bvbrc-api --contigs assembly.fna
    python -m s2f.m1_genome --run runs/kp --from-bvbrc-api --species "Klebsiella pneumoniae"
    python -m s2f.m1_genome --run runs/kp --from-bvbrc-api --genome-id 1125630.4

This is the fallback path step 5 of [01-m1-genome.md](01-m1-genome.md) asked for — "pull the
existing BV-BRC annotation for a known genome and emit the same `proteins[]` shape" — built
on the Data REST API rather than the `p3-` CLI, so it needs no account and no CLI install.

## The one thing that must not be misread

**The API route does not annotate your assembly.** It resolves the organism to a genome
BV-BRC has *already* annotated — normally a reference or representative genome — and
describes that. So:

* an isolate-specific gene (a plasmid-borne carbapenemase, a prophage) is **invisible**;
* a gene the reference carries and your isolate does not is a **false positive**;
* AMR phenotypes and isolation metadata belong to *that* genome, not to your sample.

It is the right tool for orienting a run, building the report, and unblocking M2–M6 before
any CGA job finishes. It cannot support a claim about the sample in hand. `genome.resolution`
records exactly which genome was chosen and why, because a species name resolves to a
different reference as BV-BRC is recurated, and `genome.annotation_route` is `api` so nothing
downstream has to guess. The HTML report prints the caveat at the top of the page.

| | CGA route | API route |
| --- | --- | --- |
| Describes | your assembly | a reference genome for the organism |
| Needs | account + `p3-` CLI | nothing |
| Runtime | hours | seconds |
| Per-feature subsystems | yes | no (`subsystems` is `[]`) |
| Codon tree | yes (`tree_newick`) | no (gene-content tree instead) |
| Isolate metadata | as submitted | the reference genome's |
| Laboratory AMR phenotypes | no | yes, where BV-BRC holds them |

## One contract, two routes

`collect.proteins_section`, `collect.genome_section`, `collect.run_fields` and
`collect.write_m1_dir` deliberately mirror their `parse.py` counterparts, and
`write_m1_dir` imports `parse.FEATURE_COLUMNS` and `parse.SPECIALTY_COLUMNS` rather than
keeping its own copy. A test asserts that an API protein record has exactly the keys a CGA
one has, and another loads an API run through M2's own `bvbrc_input.load_input`. Two
differently-shaped `proteins[]` depending on which route ran is precisely the drift
[00a-data-contract.md](00a-data-contract.md) exists to prevent.

API-only metadata is additive and namespaced under `genome` — `growth`, `isolation`,
`nutrition`, `close_human_pathogens`, `amr_phenotypes`, `specialty_gene_counts` — never a
redefinition of a field `parse.py` already owns.

## Query correctness

Three things about the Data API will silently corrupt a run if you get them wrong.

**Pin the annotation source.** BV-BRC holds both PATRIC/RASTtk and RefSeq calls for many
genomes. A `genome_feature` query without `eq(annotation,PATRIC)` returns both, so every
gene appears twice with a different `patric_id` and every count in the report doubles.

**Quote multi-word RQL values.** `eq(species,Klebsiella pneumoniae)` is not an exact match;
`eq(species,"Klebsiella pneumoniae")` is. `bvbrc_api.rql_value` quotes then percent-encodes.

**RQL goes in the raw query string.** Handing it to `requests` as a parameter mapping
re-encodes the parentheses and commas, and BV-BRC then reads the whole expression as a
field name. `BvbrcApi._url` builds the URL by hand for this reason.

## Paging, caching and offline replay

Every request goes through `s2f.common.http.CachedJsonClient` (issue #4), so runs are
cached, retried and replayable. The API needs two things a body-only accessor cannot give:
a `Range` **request** header to select a page, and the `Content-Range` **response** header,
which is the only place the result-set total appears. `CachedJsonClient.get_envelope`
returns `{"status", "headers", "body"}` and caches the whole envelope — replaying a paged
fetch offline needs the totals as much as the rows. Request headers are part of the cache
key, since a different `Range` is a different response.

`--offline` replays a previous run from `<run>/m1/cache.sqlite` and exits 5 on a miss.
`Authorization` and `Accept-Encoding` are rejected as `get_envelope` headers: neither varies
the resource, and caching a credentialled response under a shared key would leak it between
runs. The route reads public data and sends no credentials.

Metadata facets are best-effort. A failed facet is appended to `BvbrcApi.failures` and
reported at the end of the run rather than raised — a missing growth-condition facet should
degrade the report, not fail the analysis. Findings are never best-effort this way.

## What a run writes

Beyond the contract (`report.json` and the three files M2 reads), the API route writes the
per-section tables, because a spreadsheet of virulence factors is what a microbiologist
actually opens:

    m1/amr_genes.csv                       m1/growth_conditions.csv
    m1/virulence_factors.csv               m1/isolation_genome.csv
    m1/pathogenesis_host_invasion.csv      m1/isolation_species_distribution.csv
    m1/close_human_pathogens.csv           m1/nutrition_biosynthesis.csv
    m1/proteins_ranked.csv                 m1/amr_phenotypes.csv
    m1/taxonomy_tree.txt                   m1/disease_profile.json
    m1/gene_content_tree.nwk + .svg        m1/report.md
    m1/report.html      (--html)           m1/figures/*.png   (--figures)

All of them are views of `report.json`, never a second contract. `pathogenesis_host_invasion.csv`
is the one that carries evidence the ranking does not: `priority.PATHOGENESIS_RE` is
deliberately broader than the mechanism rules, so it catches a "colonization factor antigen"
or a "pilin subunit" that no curated database flagged and no mechanism regex matched. Its
`evidence` column distinguishes a VFDB hit from a keyword match, because those are very
different claims.

**Order of operations matters here.** The contract sections are written *before* the
phylogeny, the reports and the figures. Everything after that point is a figure or a
convenience, and a transient HTTP error while building one must not leave a run with no
`genome` and no `proteins`.

## Inferences, labelled as such

Two sections are derived, not measured, and both say so in the data as well as the report.

`genome.nutrition` infers requirement from encoded metabolic subsystems: a genome that
encodes a biosynthesis pathway need not be fed that metabolite, one that lacks it probably
must be. It reports counts and carries `basis: "inferred from encoded metabolic subsystems;
not a growth assay"`, because a missing subsystem *annotation* and a genuinely absent
*pathway* are indistinguishable here.

`genome.close_human_pathogens` counts BV-BRC genomes with `host_name=Homo sapiens` in the
same genus, falling back to family. A species listed there means "BV-BRC holds human
isolates of it", not "it is a validated human pathogen"; `basis` records the query.

Growth conditions are per-row attributed: `source` is either `this genome` or
`species-typical (N genomes)`. "This genome is anaerobic" and "most genomes of this species
are anaerobic" are different claims and the report must not present the second as the first.

## Acceptance checks

* An API `proteins[]` record has exactly the keys a CGA one has.
* `report.json` from an API run passes `validate_report`.
* M2's `bvbrc_input.load_input` reads `<run>/m1/` from an API run unchanged.
* `genome.annotation_route == "api"` and `genome.resolution` is non-empty.
* `--offline` replays a cached run with no outbound request, and exits 5 on a miss.
* The whole suite passes with scipy, numpy and matplotlib absent.

`tests/test_m1_bvbrc_api.py`, `tests/test_m1_tree_and_cli.py`.

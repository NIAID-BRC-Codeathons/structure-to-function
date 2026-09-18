"""Collect the M1 sections from the BV-BRC Data API and shape them like a CGA run.

The point of this module is that the *output* is indistinguishable from `parse.py`'s:
`genome_section`, `proteins_section`, `run_fields` and the three files in `<run>/m1/` come
out in the same shapes, so M2 reads an API run with no knowledge that it was one and
`report.json` stays a single contract (`docs/00a-data-contract.md`).

Where a CGA run has no equivalent the field is additive and namespaced under `genome`
(`growth`, `isolation`, `nutrition`, `close_human_pathogens`), never a redefinition of a
field `parse.py` already owns.

What the API gives that CGA does not: curated isolate and growth metadata, laboratory AMR
phenotypes, and species-level context. What it does not give: an annotation of the
assembly in hand. See `bvbrc_api` for why that distinction is load-bearing.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

from ..common.http import HttpError
from .bvbrc_api import (GENOME_URL, BvbrcApi, nonempty_facets, rql_value, top_nonempty)
from .parse import FEATURE_COLUMNS, SPECIALTY_COLUMNS, _as_float

# `cell_arrangement` and `ph_range` are deliberately absent: neither is in BV-BRC's
# `genome` Solr schema, so faceting on them returns HTTP 400 `undefined field` for every
# genome, on every run. Confirmed against the live API 2026-09-17:
#   "msg": "... undefined field: \"cell_arrangement\"", "code": 400
# Re-adding one costs a guaranteed-dead request and a Solr stack trace in the console.
GROWTH_FIELDS = [
    ("gram_stain", "Gram stain"),
    ("cell_shape", "Cell shape"),
    ("motility", "Motility"),
    ("sporulation", "Sporulation"),
    ("optimal_temperature", "Optimal temperature (C)"),
    ("temperature_range", "Temperature range"),
    ("oxygen_requirement", "Oxygen requirement"),
    ("salinity", "Salinity"),
    ("habitat", "Habitat"),
]
ISOLATION_FIELDS = [
    ("isolation_country", "Isolation country"),
    ("geographic_location", "Geographic location"),
    ("geographic_group", "Geographic group"),
    ("isolation_source", "Isolation source"),
    ("isolation_site", "Isolation site"),
    ("host_name", "Host"),
    ("host_common_name", "Host (common)"),
    ("host_health", "Host health"),
    ("body_sample_site", "Body sample site"),
    ("disease", "Disease"),
    ("collection_date", "Collection date"),
    ("collection_year", "Collection year"),
    ("latitude", "Latitude"),
    ("longitude", "Longitude"),
]

FEATURE_SELECT = (
    "select(patric_id,refseq_locus_tag,gene,product,aa_length,start,end,strand,"
    "plfam_id,pgfam_id,aa_sequence_md5)"
)
SPECIALTY_SELECT = (
    "select(patric_id,gene,product,property,source,property_source,classification,"
    "antibiotics_class,function,identity,query_coverage,subject_coverage,evidence,"
    "source_id,same_species,same_genus)"
)
AMR_SELECT = (
    "select(antibiotic,resistant_phenotype,measurement,measurement_unit,"
    "laboratory_typing_method,laboratory_typing_platform,evidence,computational_method)"
)

AMR_PHENOTYPE_COLUMNS = [
    "antibiotic", "resistant_phenotype", "measurement", "measurement_unit",
    "laboratory_typing_method", "evidence", "computational_method",
]

#: Per-section CSVs. `report.json` is the contract and holds all of this, but a
#: spreadsheet of virulence factors is what a microbiologist actually opens, so the
#: tables are written too rather than making everyone parse JSON.
SECTION_CSVS: dict[str, list[str]] = {
    "amr_genes.csv": ["gene", "product", "function", "source", "classification",
                      "antibiotics_class", "identity", "query_coverage", "patric_id"],
    "virulence_factors.csv": ["gene", "product", "function", "source", "classification",
                              "identity", "query_coverage", "source_id", "patric_id"],
    "pathogenesis_host_invasion.csv": ["patric_id", "gene", "product", "evidence",
                                       "source"],
    "close_human_pathogens.csv": ["species", "human_associated_genomes"],
    "growth_conditions.csv": ["property", "value", "source"],
    "isolation_genome.csv": ["property", "value"],
    "isolation_species_distribution.csv": ["field", "value", "genomes"],
    "nutrition_biosynthesis.csv": ["category", "class", "subsystem_name"],
    "proteins_ranked.csv": ["rank", "score", "patric_id", "gene", "locus_tag", "product",
                            "aa_length", "categories", "specialty_types", "selected",
                            "mechanism_hypothesis"],
}


def _as_int(value: Any) -> int | None:
    """BV-BRC writes coordinates as strings; the schema wants integers or null."""
    if value is None or value == "":
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _present(value: Any) -> bool:
    return value not in (None, "", [])


# --- individual collectors --------------------------------------------------
def collect_features(api: BvbrcApi, genome_id: str, *, cap: int = 25000,
                     annotation: str = "PATRIC") -> list[dict[str, Any]]:
    """CDS features from **one** annotation source.

    Pinning `annotation` is not a preference, it is a correctness requirement: BV-BRC
    holds both PATRIC/RASTtk and RefSeq calls for many genomes, and querying without it
    returns both, so every gene appears twice with different `patric_id`s and every count
    in the report doubles.
    """
    return api.query_all(
        "genome_feature",
        f"and(eq(genome_id,{rql_value(genome_id)}),eq(feature_type,CDS),"
        f"eq(annotation,{rql_value(annotation)}))&{FEATURE_SELECT}&sort(+start)",
        cap=cap,
    )


def collect_specialty(api: BvbrcApi, genome_id: str, *,
                      cap: int = 25000) -> tuple[list[dict[str, Any]], dict[str, list]]:
    """All specialty-gene rows, plus the same rows bucketed by `property`."""
    rows = api.query_all(
        "sp_gene",
        f"eq(genome_id,{rql_value(genome_id)})&{SPECIALTY_SELECT}&sort(+property)",
        cap=cap,
    )
    buckets: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        buckets.setdefault(row.get("property") or "Unknown", []).append(row)
    return rows, buckets


def collect_amr_phenotypes(api: BvbrcApi, genome_id: str) -> list[dict[str, Any]]:
    """Laboratory AMR phenotypes. Empty for most genomes, and that is not an error."""
    return api.query_all(
        "genome_amr",
        f"eq(genome_id,{rql_value(genome_id)})&{AMR_SELECT}&sort(+antibiotic)",
        cap=20000,
    )


def collect_taxonomy(api: BvbrcApi, genome: dict[str, Any]) -> dict[str, Any]:
    """Lineage with ranks, plus reference/representative neighbours in the same genus."""
    taxon_id = genome.get("taxon_id")
    names = genome.get("taxon_lineage_names")
    ids = genome.get("taxon_lineage_ids")
    ranks = None
    # The genome core does not return `genetic_code` at all -- BV-BRC omits unset fields
    # rather than returning null, so the key is simply absent from a real genome record
    # (verified live against 1125630.4, 68 fields, no genetic_code). The taxonomy core
    # does carry it, which is where the CGA route reads it too (taxon.py).
    genetic_code = None
    if taxon_id:
        hits = api.query(
            "taxonomy",
            f"eq(taxon_id,{rql_value(taxon_id)})&select(taxon_id,taxon_name,taxon_rank,"
            "genetic_code,lineage_names,lineage_ranks,lineage_ids)&limit(1)",
        )
        if hits:
            names = hits[0].get("lineage_names", names)
            ranks = hits[0].get("lineage_ranks")
            ids = hits[0].get("lineage_ids", ids)
            genetic_code = hits[0].get("genetic_code")

    names = names or []
    ranks = ranks or [""] * len(names)
    ids = ids or [""] * len(names)
    lineage = [{"rank": r, "name": n, "taxon_id": i}
               for r, n, i in zip(ranks, names, ids)]

    neighbors: list[dict[str, Any]] = []
    genus = genome.get("genus")
    if genus:
        neighbors = api.query(
            "genome",
            f"and(eq(genus,{rql_value(genus)}),or(eq(reference_genome,Reference),"
            "eq(reference_genome,Representative)))"
            "&select(genome_id,genome_name,species,reference_genome)&sort(+species)&limit(15)",
        )
    return {"taxon_id": taxon_id, "lineage": lineage, "neighbors": neighbors,
            "genetic_code": genetic_code}


def collect_growth(api: BvbrcApi, genome: dict[str, Any], *,
                   own_label: str = "this genome") -> list[dict[str, Any]]:
    """Growth conditions for this genome, falling back to the species-typical value.

    Each row carries its own `source`, because "this genome is anaerobic" and "most
    genomes of this species are anaerobic" are different claims and the report must not
    present the second as the first.
    """
    species = genome.get("species")
    rows: list[dict[str, Any]] = []
    for key, label in GROWTH_FIELDS:
        value = genome.get(key)
        if _present(value):
            rows.append({"property": label, "value": value, "source": own_label})
            continue
        if not species:
            continue
        value, count = top_nonempty(
            api.facet("genome", f"eq(species,{rql_value(species)})", key, limit=8)
        )
        if value is not None:
            rows.append({"property": label, "value": value,
                         "source": f"species-typical ({count} genomes)"})
    return rows


def collect_isolation(api: BvbrcApi, genome: dict[str, Any], *,
                      own_label: str = "this genome") -> dict[str, Any]:
    """This isolate's provenance, plus the species-wide distribution as context."""
    own = [{"property": label, "value": genome[key], "source": own_label}
           for key, label in ISOLATION_FIELDS if _present(genome.get(key))]
    distribution: dict[str, list[list[Any]]] = {}
    species = genome.get("species")
    if species:
        for key, label in [("isolation_country", "Isolation country"),
                           ("isolation_source", "Isolation source"),
                           ("host_name", "Host"),
                           ("geographic_group", "Geographic group")]:
            items = nonempty_facets(
                api.facet("genome", f"eq(species,{rql_value(species)})", key, limit=8)
            )
            if items:
                # Stored as lists, not tuples: this goes through JSON, which has no tuple,
                # so writing tuples would make the in-memory and round-tripped shapes
                # differ and only the round-tripped one is ever read back.
                distribution[label] = [[value, count] for value, count in items[:6]]
    return {"genome": own, "species_distribution": distribution}


def collect_nutrition(api: BvbrcApi, genome: dict[str, Any]) -> dict[str, Any]:
    """Nutritional requirement *inferred* from encoded metabolic subsystems.

    This is a proxy, not an assay. The logic is the standard auxotrophy inference: a
    genome that encodes a biosynthesis pathway need not be fed that metabolite, and one
    that lacks it probably must be. Counts are reported rather than a verdict, because
    a missing subsystem annotation and a genuinely absent pathway look identical here.
    """
    genome_id = genome["genome_id"]
    species = genome.get("species")
    superclasses = dict(api.facet("subsystem", f"eq(genome_id,{rql_value(genome_id)})",
                                  "superclass", limit=30))
    classes = dict(api.facet(
        "subsystem",
        f"and(eq(genome_id,{rql_value(genome_id)}),eq(superclass,Metabolism))",
        "class", limit=40,
    ))
    # Bucketed client-side rather than by a `class` query: BV-BRC class names for amino
    # acid and cofactor metabolism are not stable enough to match on server-side.
    rows = api.query_all(
        "subsystem",
        f"and(eq(genome_id,{rql_value(genome_id)}),eq(superclass,Metabolism))"
        "&select(class,subclass,subsystem_name)&sort(+class)",
        cap=20000,
    )
    amino_acid: dict[str, str] = {}
    cofactor: dict[str, str] = {}
    for row in rows:
        cls = (row.get("class") or "").lower()
        name = row.get("subsystem_name")
        if not name:
            continue
        if "amino acid" in cls:
            amino_acid[name] = row.get("class") or ""
        elif "cofactor" in cls or "vitamin" in cls:
            cofactor[name] = row.get("class") or ""

    biosynthesis = (
        [{"category": "Amino acid biosynthesis", "class": c, "subsystem_name": n}
         for n, c in sorted(amino_acid.items())]
        + [{"category": "Cofactor/vitamin biosynthesis", "class": c, "subsystem_name": n}
           for n, c in sorted(cofactor.items())]
    )

    oxygen = genome.get("oxygen_requirement")
    if not oxygen and species:
        value, _ = top_nonempty(
            api.facet("genome", f"eq(species,{rql_value(species)})",
                      "oxygen_requirement", limit=8)
        )
        oxygen = f"{value} (species-typical)" if value else None

    return {
        "oxygen_requirement": oxygen,
        "metabolism_superclasses": superclasses,
        "metabolism_classes": classes,
        "amino_acid_biosynthesis_count": len(amino_acid),
        "cofactor_vitamin_biosynthesis_count": len(cofactor),
        "biosynthesis_subsystems": biosynthesis,
        "basis": "inferred from encoded metabolic subsystems; not a growth assay",
    }


#: Everything a public genome record can contribute to the metadata sections.
_METADATA_SELECT = sorted(
    {key for key, _ in GROWTH_FIELDS}
    | {key for key, _ in ISOLATION_FIELDS}
    | {"genome_id", "genome_name", "species", "genus", "disease"}
)


def resolve_metadata_genome(
    api: BvbrcApi,
    *,
    genome_id: str = "",
    species: str = "",
    closest: Any = (),
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Which public BV-BRC record the species-level metadata may be taken from.

    A CGA run analyses an assembly BV-BRC has never seen. `2097.118` is the id CGA minted
    for our own submission, and `eq(genome_id,2097.118)` returns nothing, so growth,
    isolation, host and disease — which are held against *public* records — have to come
    from somewhere else: this genome if it happens to be public, otherwise the nearest
    relative the taxon call found, otherwise species-wide facets alone.

    The basis is returned rather than folded away, because "this isolate came from a human
    in Australia" and "a relative of this isolate did" are different claims and the report
    must not print the second as the first.
    """
    select = f"&select({','.join(_METADATA_SELECT)})&limit(1)"

    def fetch(gid: str) -> dict[str, Any]:
        rows = api.query("genome", f"eq(genome_id,{rql_value(gid)}){select}")
        return rows[0] if rows else {}

    if genome_id:
        record = fetch(genome_id)
        if record:
            return record, {"basis": "this-genome", "genome_id": str(genome_id),
                            "genome_name": record.get("genome_name") or "",
                            "mash_distance": 0.0, "ani": 100.0}
    for hit in closest or ():
        gid = str((hit or {}).get("genome_id") or "")
        if not gid:
            continue
        record = fetch(gid)
        if record:
            return record, {"basis": "relative", "genome_id": gid,
                            "genome_name": record.get("genome_name") or hit.get("name") or "",
                            "mash_distance": hit.get("mash_distance"),
                            "ani": hit.get("ani")}
    return {}, {"basis": "species" if species else "none", "genome_id": "",
                "genome_name": "", "mash_distance": None, "ani": None}


_METADATA_NOTE = {
    "this-genome": ("This assembly has its own public BV-BRC record, so the isolate "
                    "metadata below describes this genome."),
    "relative": ("This assembly is not in BV-BRC. The isolate metadata below belongs to "
                 "the nearest public relative, not to the analysed genome."),
    "species": ("This assembly is not in BV-BRC and no public relative was identified, so "
                "only species-wide distributions are shown."),
    "none": "No species was resolved, so no metadata could be looked up.",
}


def metadata_candidates(genome: dict[str, Any]) -> list[dict[str, Any]]:
    """Public genomes to try as a metadata donor, nearest first.

    Two sources, because only one of them is ever populated on a given run.
    ``closest_genomes`` comes from the Minhash taxon call and carries distances, but is
    empty whenever the taxon was supplied or a CGA directory was parsed directly.
    ``tree_ingroup`` comes from CGA's own codon tree and is a bare list of genome ids with
    no distances — but it is present exactly when the other is not, and its first entry is
    the nearest genome CGA placed next to ours.

    Reading only ``closest_genomes`` is why a run against a pre-fetched CGA directory fell
    through to species-wide facets and reported "Aerobic" for an organism whose nearest
    relative's record says "Facultative".
    """
    seen: set[str] = set()
    out: list[dict[str, Any]] = []
    for hit in genome.get("closest_genomes") or []:
        gid = str((hit or {}).get("genome_id") or "")
        if gid and gid not in seen:
            seen.add(gid)
            out.append(dict(hit))
    for gid in genome.get("tree_ingroup") or []:
        gid = str(gid or "")
        if gid and gid not in seen:
            seen.add(gid)
            # No distance: the codon tree gives order, not a metric. Left null rather than
            # invented, and filed as its own issue.
            out.append({"genome_id": gid, "name": "", "mash_distance": None, "ani": None})
    return out


def collect_cga_metadata(api: BvbrcApi, genome: dict[str, Any]) -> dict[str, Any]:
    """Species-level metadata for a CGA run, with an honest account of its origin."""
    taxonomy = genome.get("taxonomy") or {}
    species = (taxonomy.get("scientific_name") or "").strip()
    record, provenance = resolve_metadata_genome(
        api,
        genome_id=str(genome.get("genome_id") or ""),
        species=species,
        closest=metadata_candidates(genome),
    )
    basis = provenance["basis"]
    provenance["species"] = species
    provenance["note"] = _METADATA_NOTE[basis]

    donor = dict(record)
    donor.setdefault("species", species)
    own_label = {"this-genome": "this genome",
                 "relative": f"relative {provenance['genome_id']}",
                 "species": "species-typical",
                 "none": "species-typical"}[basis]

    disease = donor.get("disease") or []
    if isinstance(disease, str):
        disease = [disease]

    return {
        "metadata_provenance": provenance,
        "growth": collect_growth(api, donor, own_label=own_label) if species else [],
        "isolation": (collect_isolation(api, donor, own_label=own_label) if species else {}),
        "amr_phenotypes": (collect_amr_phenotypes(api, provenance["genome_id"])
                           if provenance["genome_id"] else []),
        "disease": list(disease),
    }


def collect_close_pathogens(api: BvbrcApi, genome: dict[str, Any]) -> dict[str, Any]:
    """Related species that BV-BRC records isolated from humans.

    Human-host *genome count* is the evidence, so a species appearing here means "BV-BRC
    holds human isolates of it", not "it is a validated human pathogen". The genus is
    tried before the family because a family-level answer is usually too broad to mean
    anything.
    """
    species = genome.get("species")
    for level, field in (("genus", "genus"), ("family", "family")):
        value = genome.get(field)
        if not value:
            continue
        facet = api.facet(
            "genome",
            f"and(eq({field},{rql_value(value)}),eq(host_name,{rql_value('Homo sapiens')}))",
            "species", mincount=1, limit=25,
        )
        rows = [{"species": s, "human_associated_genomes": c}
                for s, c in nonempty_facets(facet) if s != species]
        if rows:
            return {"level": level, "query_species": species, "species": rows[:15],
                    "basis": "count of BV-BRC genomes with host_name=Homo sapiens"}
    return {"level": None, "query_species": species, "species": [],
            "basis": "count of BV-BRC genomes with host_name=Homo sapiens"}


def fetch_sequences(api: BvbrcApi, features: list[dict[str, Any]], *, cap: int = 4000,
                    priority: dict[str, dict[str, Any]] | None = None) -> dict[str, str]:
    """`{feature_id: aa_sequence}` for up to `cap` features, de-duplicated by md5.

    Identical proteins share one md5 in BV-BRC, so fetching by md5 rather than by feature
    collapses paralogue families into a single request each.

    `cap` truncates by **priority** when a ranking is supplied, not by genome coordinate.
    `collect_features` sorts by `+start`, so a coordinate-ordered cap on a genome with
    more CDS than `cap` would drop whichever proteins happen to sit late on the
    chromosome — possibly every top-ranked candidate — from `m1/proteins.faa`, which is
    M2's input. Selected proteins come first, then by rank.

    A failed chunk costs its sequences, not the run: sequences are an output, and the
    contract in `report.json` and the two CSVs does not depend on them.
    """
    candidates = [f for f in features if f.get("aa_sequence_md5") and f.get("patric_id")]
    if priority:
        def order(feature: dict[str, Any]) -> tuple[int, int]:
            entry = priority.get(feature["patric_id"]) or {}
            return (0 if entry.get("selected") else 1, entry.get("rank") or 10**9)
        candidates.sort(key=order)
    wanted = candidates[:cap]

    unique: list[str] = []
    seen: set[str] = set()
    for feature in wanted:
        md5 = feature["aa_sequence_md5"]
        if md5 not in seen:
            seen.add(md5)
            unique.append(md5)

    by_md5: dict[str, str] = {}
    for start in range(0, len(unique), 200):
        chunk = unique[start:start + 200]
        rql = "in(md5,(" + ",".join(chunk) + "))&select(md5,sequence)"
        try:
            records = api.query_all("feature_sequence", rql, cap=len(chunk) + 5)
        except HttpError as exc:
            api.failures.append({"facet": f"feature_sequence[{start}:{start + len(chunk)}]",
                                 "error": str(exc)})
            continue
        for record in records:
            if record.get("md5"):
                by_md5[record["md5"]] = record.get("sequence", "")

    return {f["patric_id"]: by_md5[f["aa_sequence_md5"]]
            for f in wanted if by_md5.get(f["aa_sequence_md5"])}


# --- assembly into report.json shapes ---------------------------------------
def collect_all(api: BvbrcApi, genome: dict[str, Any], *, cap: int = 25000,
                annotation: str = "PATRIC") -> dict[str, Any]:
    """Run every collector for one genome and return the raw section bundle."""
    from .priority import pathogenesis_proteins

    genome_id = genome["genome_id"]
    features = collect_features(api, genome_id, cap=cap, annotation=annotation)
    specialty, buckets = collect_specialty(api, genome_id, cap=cap)
    virulence = [row for prop, rows in buckets.items()
                 if "virul" in prop.lower() for row in rows]
    return {
        # The broader keyword tier: proteins whose product reads like pathogenesis but
        # which no curated database flagged. Weaker evidence than a VFDB hit and labelled
        # as such, but it is the only tier that catches an unannotated adhesin.
        "pathogenesis": pathogenesis_proteins(features, virulence),
        "features": features,
        "specialty": specialty,
        "specialty_buckets": buckets,
        "specialty_counts": {k: len(v) for k, v in sorted(buckets.items())},
        "virulence": virulence,
        "amr_genes": buckets.get("Antibiotic Resistance", []),
        "amr_phenotypes": collect_amr_phenotypes(api, genome_id),
        "taxonomy": collect_taxonomy(api, genome),
        "growth": collect_growth(api, genome),
        "isolation": collect_isolation(api, genome),
        "nutrition": collect_nutrition(api, genome),
        "close_pathogens": collect_close_pathogens(api, genome),
        "annotation": annotation,
    }


def _subsystems_placeholder() -> list[dict[str, Any]]:
    """The API route has no per-feature subsystem bindings, so the key stays empty.

    `parse.py` fills `subsystems` from CGA's role bindings. The Data API exposes
    subsystems per genome, not per feature, so there is nothing to put here. Emitting
    `[]` keeps the protein record the same shape for consumers; inventing bindings from
    the genome-level list would be fabricated evidence.
    """
    return []


def proteins_section(bundle: dict[str, Any]) -> list[dict[str, Any]]:
    """`proteins[]` for report.json, shaped exactly like `parse.proteins_section`."""
    specialty_by_feature: dict[str, list[dict[str, Any]]] = {}
    for row in bundle["specialty"]:
        fid = row.get("patric_id")
        if not fid:
            continue
        antibiotics = row.get("antibiotics_class") or []
        # Field for field with `parse._specialty_by_feature`, including the numeric
        # normalization: DIAMOND rows carry `identity` as a string and AMRFinderPlus as
        # a float, and downstream comparisons need one type. Reusing `parse._as_float`
        # rather than repeating it is what stops the two routes drifting apart.
        specialty_by_feature.setdefault(fid, []).append({
            "type": row.get("property"),
            "classification": row.get("classification"),
            "database": row.get("source"),
            "evidence": row.get("evidence"),
            "gene": row.get("gene"),
            "hit": row.get("source_id") or row.get("product"),
            "identity": _as_float(row.get("identity")),
            "coverage": _as_float(row.get("query_coverage")),
            "subject_coverage": _as_float(row.get("subject_coverage")),
            "e_value": _as_float(row.get("e_value")),
            "same_species": row.get("same_species"),
            "same_genus": row.get("same_genus"),
            "antibiotics": antibiotics if isinstance(antibiotics, list) else [antibiotics],
            "pmid": row.get("pmid") or [],
        })

    proteins = []
    for row in bundle["features"]:
        fid = row.get("patric_id")
        if not fid:
            continue
        proteins.append({
            "feature_id": fid,
            "locus_tag": row.get("refseq_locus_tag"),
            "gene": row.get("gene"),
            "product": row.get("product"),
            "contig": row.get("sequence_id") or row.get("accession"),
            "start": _as_int(row.get("start")),
            "end": _as_int(row.get("end")),
            "strand": row.get("strand"),
            "aa_length": _as_int(row.get("aa_length")),
            "subsystems": _subsystems_placeholder(),
            "specialty": specialty_by_feature.get(fid, []),
            "families": {"pgfam": row.get("pgfam_id"), "plfam": row.get("plfam_id")},
            "kmer_confidence": None,
        })
    return proteins


def genome_section(genome: dict[str, Any], bundle: dict[str, Any],
                   resolution: str) -> dict[str, Any]:
    """`genome` for report.json, shaped like `parse.genome_section` plus API-only keys.

    `annotation_route` is the field that stops the two M1 routes being confused. Anything
    reading `genome` can tell whether the features describe the submitted assembly (`cga`)
    or a reference genome for the same organism (`api`), which changes what a gene's
    presence or absence is allowed to mean.
    """
    taxonomy = bundle["taxonomy"]
    closest = [
        {
            "genome_id": n.get("genome_id"),
            "name": n.get("genome_name"),
            "mash_distance": None,
            "pvalue": None,
            "shared_kmers": None,
            "ani": None,
            "snp_distance": None,
            "reference_genome": n.get("reference_genome"),
        }
        for n in taxonomy.get("neighbors") or []
    ]
    lineage_names = [node["name"] for node in taxonomy.get("lineage") or []]
    return {
        "genome_id": str(genome.get("genome_id") or ""),
        "taxon_id": int(genome["taxon_id"]) if genome.get("taxon_id") else None,
        "taxonomy": {
            "scientific_name": genome.get("genome_name"),
            "lineage_names": lineage_names,
            "lineage_ranked": taxonomy.get("lineage"),
            "genetic_code": taxonomy.get("genetic_code") or genome.get("genetic_code"),
            "called_by": "bvbrc_api",
            "called_rank": "species" if genome.get("species") else None,
            "top_hit_distance": None,
        },
        "closest_genomes": closest,
        "tree_ingroup": [],
        "tree_newick": None,
        "cga_job_id": None,
        # --- API-route provenance and metadata (no CGA equivalent) ---
        "annotation_route": "api",
        "annotation_source": bundle.get("annotation", "PATRIC"),
        "resolution": resolution,
        "bvbrc_url": GENOME_URL.format(gid=genome.get("genome_id")),
        "species": genome.get("species"),
        "genus": genome.get("genus"),
        "family": genome.get("family"),
        "order": genome.get("order"),
        "gc_content": genome.get("gc_content"),
        "disease": _as_list(genome.get("disease")),
        # Carried through so `pathogens.build_reference_list` can cite the primary genome
        # paper. Without it here the PubMed citations were unreachable, because the
        # report builds its reference list from `report.json`, not from the API record.
        "publication": genome.get("publication"),
        "growth": bundle.get("growth", []),
        "isolation": bundle.get("isolation", {}),
        "nutrition": bundle.get("nutrition", {}),
        "close_human_pathogens": bundle.get("close_pathogens", {}),
        "amr_phenotypes": bundle.get("amr_phenotypes", []),
        "specialty_gene_counts": bundle.get("specialty_counts", {}),
        "quality": {
            "genome_quality": genome.get("genome_quality") or "",
            "genome_quality_flags": [],
            "coarse_consistency": genome.get("coarse_consistency"),
            "fine_consistency": genome.get("fine_consistency"),
            "cds_ratio": genome.get("cds_ratio"),
            "hypothetical_cds_ratio": genome.get("hypothetical_cds_ratio"),
            "contigs": genome.get("contigs"),
            "genome_length": genome.get("genome_length"),
            "genome_status": genome.get("genome_status"),
            "patric_cds": genome.get("patric_cds"),
            "feature_summary": None,
            "protein_summary": None,
        },
    }


def _as_list(value: Any) -> list[Any]:
    if value is None or value == "":
        return []
    return value if isinstance(value, list) else [value]


def run_fields(genome: dict[str, Any], bundle: dict[str, Any],
               resolution: str) -> dict[str, Any]:
    """`run` manifest fields for an API-route run."""
    return {
        "bvbrc_api": {
            "genome_id": str(genome.get("genome_id") or ""),
            "resolution": resolution,
            "annotation_source": bundle.get("annotation", "PATRIC"),
            "features": len(bundle.get("features") or []),
            "specialty_rows": len(bundle.get("specialty") or []),
        },
        "tool_versions": {"bvbrc_data_api": "v1 (unversioned REST endpoint)"},
    }


def write_m1_dir(bundle: dict[str, Any], sequences: dict[str, str],
                 m1_dir: str | Path) -> dict[str, int]:
    """Write the three files M2 reads, with the column names `parse.py` defines.

    The columns come from `parse.FEATURE_COLUMNS` / `parse.SPECIALTY_COLUMNS` rather than
    a local copy, so the two M1 routes cannot drift into writing different headers for
    the same file (pitfall: a second contract is worse than no contract).
    """
    m1_dir = Path(m1_dir)
    m1_dir.mkdir(parents=True, exist_ok=True)

    written = 0
    with (m1_dir / "proteins.faa").open("w", encoding="utf-8") as handle:
        for row in bundle["features"]:
            fid = row.get("patric_id")
            sequence = sequences.get(fid or "")
            if not fid or not sequence:
                continue
            handle.write(f">{fid} {row.get('product') or ''}\n")
            for i in range(0, len(sequence), 60):
                handle.write(sequence[i:i + 60] + "\n")
            written += 1

    feature_rows = 0
    with (m1_dir / "genes_proteins.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FEATURE_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        for row in bundle["features"]:
            if not row.get("patric_id"):
                continue
            writer.writerow({k: row.get(k, "") for k in FEATURE_COLUMNS})
            feature_rows += 1

    specialty_rows = 0
    with (m1_dir / "specialty_genes_all.csv").open("w", encoding="utf-8",
                                                   newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SPECIALTY_COLUMNS,
                                extrasaction="ignore")
        writer.writeheader()
        for row in bundle["specialty"]:
            antibiotics = row.get("antibiotics_class") or row.get("antibiotics") or ""
            writer.writerow({
                "property": row.get("property", ""),
                "gene": row.get("gene", ""),
                "product": row.get("product", ""),
                "source": row.get("source", ""),
                "classification": row.get("classification", ""),
                "antibiotics_class": (";".join(antibiotics)
                                      if isinstance(antibiotics, list) else antibiotics),
                "identity": "" if row.get("identity") is None else row.get("identity"),
                "query_coverage": ("" if row.get("query_coverage") is None
                                   else row.get("query_coverage")),
                "patric_id": row.get("patric_id", ""),
            })
            specialty_rows += 1

    with (m1_dir / "amr_phenotypes.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=AMR_PHENOTYPE_COLUMNS,
                                extrasaction="ignore")
        writer.writeheader()
        for row in bundle.get("amr_phenotypes") or []:
            writer.writerow({k: row.get(k, "") for k in AMR_PHENOTYPE_COLUMNS})

    return {"proteins_faa": written, "feature_rows": feature_rows,
            "specialty_rows": specialty_rows}


def _write_section_csv(path: Path, columns: list[str],
                       rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({
                key: (";".join(str(v) for v in row[key])
                      if isinstance(row.get(key), list) else row.get(key, ""))
                for key in columns
            })


def write_ranked_csv(features: list[dict[str, Any]], priority: dict[str, dict[str, Any]],
                     m1_dir: str | Path) -> Path:
    """Write `proteins_ranked.csv`, ordered by rank. Both M1 routes produce this."""
    m1_dir = Path(m1_dir)
    m1_dir.mkdir(parents=True, exist_ok=True)
    by_id = {f.get("patric_id") or f.get("feature_id"): f for f in features or []}
    rows = []
    for feature_id, entry in sorted(priority.items(),
                                    key=lambda kv: kv[1].get("rank") or 10**9):
        feature = by_id.get(feature_id) or {}
        rows.append({
            "rank": entry.get("rank"), "score": entry.get("score"),
            "patric_id": feature_id, "gene": feature.get("gene"),
            "locus_tag": feature.get("refseq_locus_tag"),
            "product": feature.get("product"), "aa_length": feature.get("aa_length"),
            "categories": entry.get("categories") or [],
            "specialty_types": entry.get("specialty_types") or [],
            "selected": entry.get("selected"),
            "mechanism_hypothesis": entry.get("mechanism_hypothesis"),
        })
    path = m1_dir / "proteins_ranked.csv"
    _write_section_csv(path, SECTION_CSVS["proteins_ranked.csv"], rows)
    return path


def write_section_csvs(bundle: dict[str, Any], m1_dir: str | Path, *,
                       priority: dict[str, dict[str, Any]] | None = None) -> list[str]:
    """Write the per-section tables and return the names written."""
    m1_dir = Path(m1_dir)
    m1_dir.mkdir(parents=True, exist_ok=True)
    isolation = bundle.get("isolation") or {}
    nutrition = bundle.get("nutrition") or {}

    if priority:
        write_ranked_csv(bundle.get("features") or [], priority, m1_dir)

    content: dict[str, list[dict[str, Any]]] = {
        "amr_genes.csv": bundle.get("amr_genes") or [],
        "virulence_factors.csv": bundle.get("virulence") or [],
        "pathogenesis_host_invasion.csv": bundle.get("pathogenesis") or [],
        "close_human_pathogens.csv": (bundle.get("close_pathogens") or {}).get("species") or [],
        "growth_conditions.csv": bundle.get("growth") or [],
        "isolation_genome.csv": isolation.get("genome") or [],
        "isolation_species_distribution.csv": [
            {"field": label, "value": value, "genomes": count}
            for label, items in (isolation.get("species_distribution") or {}).items()
            for value, count in items
        ],
        "nutrition_biosynthesis.csv": nutrition.get("biosynthesis_subsystems") or [],
    }

    written: list[str] = ["proteins_ranked.csv"] if priority else []
    for name, columns in SECTION_CSVS.items():
        if name == "proteins_ranked.csv":
            continue  # already written above, for both routes
        rows = content.get(name) or []
        if not rows:
            continue
        _write_section_csv(m1_dir / name, columns, rows)
        written.append(name)
    return written


def taxonomy_tree_text(lineage: list[dict[str, Any]], genome: dict[str, Any],
                       neighbors: list[dict[str, Any]]) -> str:
    """The lineage as an indented text tree, with the genus neighbours listed under it."""
    lines = ["Taxonomic lineage (root -> leaf):"]
    for depth, node in enumerate(lineage):
        rank = f"[{node['rank']}] " if node.get("rank") else ""
        lines.append("  " * depth + "└─ " + rank + str(node.get("name")))
    lines.append("  " * len(lineage) + "└─ * " + str(genome.get("genome_name"))
                 + f"  (genome_id {genome.get('genome_id')})")
    if neighbors:
        lines += ["", "Closest reference/representative genomes in the genus:"]
        for neighbor in neighbors:
            flag = neighbor.get("reference_genome") or ""
            lines.append(f"  - {neighbor.get('genome_name')}  [{flag}]  "
                         f"({neighbor.get('species')})")
    return "\n".join(lines)

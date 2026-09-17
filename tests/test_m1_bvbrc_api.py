"""The BV-BRC Data API route for M1: client, collectors, and the shared output contract.

Every test runs offline. Responses come from `fixtures/m1/bvbrc_api/klebsiella_hs11286.json`
through a fake `requests` session, so the suite never depends on BV-BRC being up or on what
it currently holds.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from conftest import HeaderlessSession, RangeIgnoringSession
from s2f.common.http import HttpError, OfflineCacheMiss
from s2f.m1_genome import collect, parse
from s2f.m1_genome.bvbrc_api import (BvbrcApi, GenomeResolutionError, organism_from_defline,
                                     read_fasta_defline, rql_value, top_nonempty)

FIXTURE = Path(__file__).resolve().parents[1] / "fixtures/m1/bvbrc_api/klebsiella_hs11286.json"
GENOME_ID = "1125630.4"
ADHESIN = "fig|1125630.4.peg.998"        # CFA/I fimbrial minor adhesin, VFDB Adherence
SIDEROPHORE = "fig|1125630.4.peg.3382"   # yersiniabactin receptor, AMR + iron + surface
TET = "fig|1125630.4.peg.5465"           # Tet(G) efflux pump, AMR only
HYPOTHETICAL = "fig|1125630.4.peg.5521"  # hypothetical protein, no evidence


# --- the client ------------------------------------------------------------
def test_rql_quotes_multiword_values_so_they_match_exactly() -> None:
    # Unquoted, RQL treats the space as a delimiter and the match becomes loose.
    assert rql_value("Klebsiella pneumoniae") == "%22Klebsiella%20pneumoniae%22"
    assert rql_value("1125630.4") == "1125630.4"


def test_rql_is_not_re_encoded_into_the_query_string(kp_payload, make_api, fake_session) -> None:
    api, session = make_api()
    api.query("genome", f"eq(genome_id,{GENOME_ID})")
    # The parentheses and comma must survive intact or BV-BRC reads the whole thing
    # as a field name.
    assert f"eq(genome_id,{GENOME_ID})" in session.calls[0]


def test_query_all_pages_until_the_content_range_total(kp_payload, make_api, fake_session) -> None:
    api, session = make_api()
    rows = api.query_all("genome_feature", "eq(genome_id,x)", page=4)
    assert len(rows) == len(kp_payload["features"])
    assert len(session.calls) == 3  # 4 + 4 + 1 of 9


def test_query_all_stops_at_the_cap_even_if_more_exist(kp_payload, make_api, fake_session) -> None:
    api, _ = make_api()
    assert len(api.query_all("genome_feature", "eq(genome_id,x)", page=2, cap=5)) == 5


def test_a_full_page_with_no_total_is_an_error_not_a_silent_truncation(kp_payload, make_api, fake_session) -> None:
    """A truncated proteome that reports success is worse than a failed run.

    Any intermediary that strips `Content-Range` could otherwise halve a proteome while
    the run still printed "ok".
    """
    api = make_api(fake_session(cls=HeaderlessSession))[0]
    with pytest.raises(HttpError, match="Content-Range"):
        api.query_all("genome_feature", "eq(genome_id,x)", page=4)


def test_an_unknown_total_is_treated_as_unknown(kp_payload, make_api, fake_session) -> None:
    # `items 0-3/*` is legal and means the total is not known.
    api = make_api(fake_session(cls=HeaderlessSession, total_marker="items 0-3/*"))[0]
    with pytest.raises(HttpError):
        api.query_all("genome_feature", "eq(genome_id,x)", page=4)


def test_a_short_page_with_no_total_is_the_end_of_the_data(kp_payload, make_api, fake_session) -> None:
    api = make_api(fake_session(cls=HeaderlessSession))[0]
    rows = api.query_all("genome_feature", "eq(genome_id,x)", page=50)
    assert len(rows) == len(kp_payload["features"])


def test_a_server_that_ignores_range_does_not_inflate_the_result(kp_payload, make_api, fake_session) -> None:
    """Duplicated pages would inflate every count in the report.

    This is the same class of failure that pinning `eq(annotation,PATRIC)` guards
    against, arriving by a different route.
    """
    api = make_api(fake_session(cls=RangeIgnoringSession))[0]
    rows = api.query_all("genome_feature", "eq(genome_id,x)", page=4, cap=10000)
    ids = [r["patric_id"] for r in rows]
    assert len(ids) == len(set(ids))
    assert len(ids) == 4  # the one page it kept re-serving, not 10000 copies of it


def test_a_failed_facet_is_recorded_and_does_not_abort_the_run(kp_payload, make_api, fake_session) -> None:
    api, _ = make_api(fail_facets=True)
    assert api.facet("genome", "eq(species,x)", "gram_stain") == []
    assert api.failures and "gram_stain" in api.failures[0]["facet"]


def test_top_nonempty_skips_bvbrcs_placeholders_for_missing_data() -> None:
    assert top_nonempty([("", 9), ("unknown", 4), ("Negative", 2)]) == ("Negative", 2)
    assert top_nonempty([("", 1), ("not collected", 3)]) == (None, 0)
    assert top_nonempty([]) == (None, 0)


def test_facets_are_sorted_by_count_not_left_in_dict_order(kp_payload, make_api, fake_session) -> None:
    """The bug this guards against changes a scientific claim between two identical runs.

    `JsonCache` serialises with `sort_keys=True`, so a facet read back from the cache
    arrives alphabetically rather than in Solr's descending-count order. Consumers take
    the *first* value, so a live run and an `--offline` replay would report different
    "commonest host" values — a 5000-genome one live, a 4-genome one on replay.
    Sorting explicitly here makes the two identical.
    """
    session = fake_session()
    session.payload = dict(session.payload)
    session.payload["close_pathogens_facet"] = {
        "Aardvark species": 4, "Klebsiella oxytoca": 1006, "Zebra species": 900,
    }
    api, _ = make_api(session)

    live = api.facet("genome", "eq(genus,Klebsiella)", "species")
    replay = make_api(session, offline=True)[0].facet(
        "genome", "eq(genus,Klebsiella)", "species")

    assert live == replay
    assert [name for name, _ in live] == [
        "Klebsiella oxytoca", "Zebra species", "Aardvark species"]


def test_equal_counts_order_deterministically(kp_payload, make_api, fake_session) -> None:
    session = fake_session()
    session.payload = dict(session.payload)
    session.payload["close_pathogens_facet"] = {"beta": 5, "alpha": 5, "gamma": 5}
    api, _ = make_api(session)
    facets = api.facet("genome", "eq(genus,x)", "species")
    assert [name for name, _ in facets] == ["alpha", "beta", "gamma"]


def test_offline_replays_from_cache_without_touching_the_session(kp_payload, make_api, fake_session) -> None:
    api, session = make_api()
    api.query("genome", f"eq(genome_id,{GENOME_ID})")
    first = len(session.calls)

    replayed = make_api(session, offline=True)[0].query(
        "genome", f"eq(genome_id,{GENOME_ID})")
    assert replayed[0]["genome_id"] == GENOME_ID
    assert len(session.calls) == first  # nothing new went out


def test_offline_raises_on_a_query_that_was_never_cached(kp_payload, make_api, fake_session) -> None:
    api, _ = make_api(offline=True, cache_name="never-written.sqlite")
    with pytest.raises(OfflineCacheMiss):
        api.query("genome", "eq(genome_id,never.fetched)")


# --- genome resolution ------------------------------------------------------
def test_an_explicit_genome_id_is_recorded_as_such(kp_payload, make_api, fake_session) -> None:
    api, _ = make_api()
    genome, how = api.resolve_genome(genome_id=GENOME_ID)
    assert genome["genome_id"] == GENOME_ID
    assert how == f"explicit genome_id={GENOME_ID}"


def test_resolution_says_which_preference_tier_matched(kp_payload, make_api, fake_session) -> None:
    api, _ = make_api()
    _, how = api.resolve_genome(species="Klebsiella pneumoniae")
    # Which genome a species resolved to is not recoverable later, so it is recorded.
    assert "reference/representative genome" in how
    assert GENOME_ID in how


def test_an_unknown_genome_id_raises_rather_than_falling_back(kp_payload, make_api, fake_session) -> None:
    api, _ = make_api()
    with pytest.raises(GenomeResolutionError):
        api.resolve_genome(genome_id="999999.9")


def test_organism_detection_refuses_to_guess(tmp_path) -> None:
    assert organism_from_defline(
        "NZ_CP123 Klebsiella pneumoniae strain HS11286 chromosome") == "Klebsiella pneumoniae"
    # A header with no binomial must return None: a wrong organism resolves to a wrong
    # reference genome and every downstream claim inherits the error.
    assert organism_from_defline("contig_1 length=45012 cov=31.2") is None


def test_read_fasta_defline_counts_sequences_and_bases(tmp_path) -> None:
    fasta = tmp_path / "g.fna"
    fasta.write_text(">c1 Klebsiella pneumoniae\nACGT\nACGT\n>c2\nAC\n")
    assert read_fasta_defline(str(fasta)) == ("c1 Klebsiella pneumoniae", 2, 10)


# --- collectors -------------------------------------------------------------
@pytest.fixture
def bundle(kp_payload, make_api) -> dict:
    api, _ = make_api()
    return collect.collect_all(api, kp_payload["genome"], cap=500)


def test_features_are_pinned_to_one_annotation_source(kp_payload, make_api, fake_session) -> None:
    api, session = make_api()
    collect.collect_features(api, GENOME_ID, annotation="PATRIC")
    # Without eq(annotation,...) BV-BRC returns PATRIC *and* RefSeq calls, so every gene
    # appears twice and every count in the report doubles.
    assert "eq(annotation,PATRIC)" in session.calls[0]


def test_specialty_rows_are_bucketed_by_property(bundle) -> None:
    assert set(bundle["specialty_buckets"]) == {"Virulence Factor", "Antibiotic Resistance"}
    assert len(bundle["virulence"]) == 4
    assert len(bundle["amr_genes"]) == 3


def test_growth_rows_say_whether_the_value_is_this_genome_or_the_species(bundle) -> None:
    sources = {row["property"]: row["source"] for row in bundle["growth"]}
    assert sources["Oxygen requirement"] == "this genome"
    # gram_stain is null on this genome, so it may only appear as a species-typical value
    assert sources.get("Gram stain", "species-typical").startswith("species-typical")


def test_nutrition_is_labelled_an_inference_not_a_measurement(bundle) -> None:
    nutrition = bundle["nutrition"]
    assert nutrition["amino_acid_biosynthesis_count"] == 2
    assert nutrition["cofactor_vitamin_biosynthesis_count"] == 1
    assert "not a growth assay" in nutrition["basis"]


def test_close_pathogens_report_the_basis_of_the_claim(bundle) -> None:
    close = bundle["close_pathogens"]
    assert close["level"] == "genus"
    assert "host_name=Homo sapiens" in close["basis"]
    # The query species itself must not be listed as its own close relative.
    assert "Klebsiella pneumoniae" not in {r["species"] for r in close["species"]}


def test_sequences_are_fetched_once_per_md5_not_once_per_feature(kp_payload, make_api, fake_session) -> None:
    api, session = make_api()
    sequences = collect.fetch_sequences(api, kp_payload["features"], cap=500)
    assert len(sequences) == len(kp_payload["features"])
    assert sum(1 for c in session.calls if "/feature_sequence/" in c) == 1


# --- the shared output contract ---------------------------------------------
def test_api_proteins_have_exactly_the_keys_the_cga_route_produces(bundle) -> None:
    """The two M1 routes must be indistinguishable downstream.

    If this fails, M2 is about to receive a differently-shaped protein record depending
    on which route ran — the drift `docs/00a-data-contract.md` exists to prevent.
    """
    cga = parse.load_cga(Path(__file__).resolve().parents[1] / "fixtures/m1/cga_sample")
    assert set(collect.proteins_section(bundle)[0]) == set(parse.proteins_section(cga)[0])


def test_api_specialty_entries_have_the_same_shape_and_types_as_the_cga_route(
        bundle) -> None:
    """Top-level key parity is not enough: the `specialty[]` entries must match too.

    The CGA route normalizes `identity`, `coverage`, `subject_coverage` and `e_value`
    to float and carries `same_species`/`same_genus`. An API route that passed raw
    values through would hand M2 the same field with a different type — an int from one
    route, a float from the other — which is exactly the drift
    `docs/00a-data-contract.md` exists to prevent. This caught it once.
    """
    cga = parse.load_cga(Path(__file__).resolve().parents[1] / "fixtures/m1/cga_sample")
    cga_entries = [hit for p in parse.proteins_section(cga) for hit in p["specialty"]]
    api_entries = [hit for p in collect.proteins_section(bundle) for hit in p["specialty"]]
    assert cga_entries and api_entries

    assert set(api_entries[0]) == set(cga_entries[0])

    for entry in api_entries:
        for field in ("identity", "coverage", "subject_coverage", "e_value"):
            value = entry[field]
            assert value is None or isinstance(value, float), (
                f"specialty.{field} is {type(value).__name__} on the API route; "
                f"parse._as_float makes it float-or-None on the CGA route"
            )


def test_the_genome_section_records_which_route_produced_it(bundle, kp_payload) -> None:
    section = collect.genome_section(kp_payload["genome"], bundle, "explicit genome_id")
    # Everything downstream needs to know whether these features describe the submitted
    # assembly or a reference genome for the same organism.
    assert section["annotation_route"] == "api"
    assert section["cga_job_id"] is None
    assert section["resolution"] == "explicit genome_id"
    assert section["taxonomy"]["lineage_names"][0] == "Bacteria"


def test_api_subsystems_are_empty_rather_than_invented(bundle) -> None:
    # The Data API exposes subsystems per genome, not per feature. Emitting [] keeps the
    # shape; synthesising per-feature bindings would be fabricated evidence.
    assert all(p["subsystems"] == [] for p in collect.proteins_section(bundle))


def test_m1_dir_uses_the_column_names_parse_defines(bundle, tmp_path) -> None:
    counts = collect.write_m1_dir(bundle, {"fig|1125630.4.peg.998": "MKV"}, tmp_path / "m1")
    features_header = (tmp_path / "m1/genes_proteins.csv").read_text().splitlines()[0]
    specialty_header = (tmp_path / "m1/specialty_genes_all.csv").read_text().splitlines()[0]
    assert features_header == ",".join(parse.FEATURE_COLUMNS)
    assert specialty_header == ",".join(parse.SPECIALTY_COLUMNS)
    assert counts["proteins_faa"] == 1


def test_m1_dir_is_readable_by_m2s_own_loader(bundle, tmp_path) -> None:
    """The real contract test: M2 loads an API run with no knowledge that it was one."""
    from s2f.m2_triage.bvbrc_input import load_input

    sequences = {f["patric_id"]: "MKVLAAGIVRDEQ" for f in bundle["features"]}
    collect.write_m1_dir(bundle, sequences, tmp_path / "m1")
    loaded = load_input(tmp_path / "m1")
    assert len(loaded.proteins) == len(bundle["features"])
    by_id = {p.feature_id: p for p in loaded.proteins}
    assert by_id[ADHESIN].has_virulence_or_amr
    assert by_id[TET].has_virulence_or_amr
    assert not by_id[HYPOTHETICAL].has_virulence_or_amr


def test_amr_phenotypes_are_written_even_when_no_phenotype_is_resistant(bundle, tmp_path) -> None:
    collect.write_m1_dir(bundle, {}, tmp_path / "m1")
    lines = (tmp_path / "m1/amr_phenotypes.csv").read_text().splitlines()
    assert lines[0].startswith("antibiotic,resistant_phenotype")
    assert len(lines) == len(bundle["amr_phenotypes"]) + 1


# --- regressions found by the first live run against the BV-BRC API ----------------
# 2026-09-17, genome 1125630.4, on a compute node. Every bug below passed the whole
# offline suite first, for one shared reason: the fixture-backed double replies with the
# full record whatever the query selected. So a field missing from a `select(...)` list,
# or a field the live core does not define at all, is invisible offline. The assertions
# here are therefore on the *query* and on a select-honouring double, not only on the
# parsed result.

#: Fields BV-BRC's `genome` Solr core does not define. Faceting on one returns HTTP 400
#: `undefined field: "<name>"` for every genome, on every run: a guaranteed dead request
#: plus a Solr stack trace in the operator's console.
BVBRC_UNDEFINED_GENOME_FIELDS = {"cell_arrangement", "ph_range"}


def test_no_metadata_field_is_one_the_genome_core_does_not_define() -> None:
    from s2f.m1_genome.collect import GROWTH_FIELDS, ISOLATION_FIELDS

    requested = ({field for field, _ in GROWTH_FIELDS}
                 | {field for field, _ in ISOLATION_FIELDS})
    offenders = sorted(requested & BVBRC_UNDEFINED_GENOME_FIELDS)
    assert not offenders, (
        f"{offenders} are not in BV-BRC's genome Solr schema, so faceting on them returns "
        f"HTTP 400 'undefined field' for every genome. Confirmed against the live API "
        f"2026-09-17. Drop the field rather than have every run issue a dead request."
    )


def test_genetic_code_is_read_from_the_taxonomy_core_not_the_genome_record(
        api_bundle, kp_payload) -> None:
    """The genome core does not return `genetic_code`, so it must come from taxonomy.

    BV-BRC omits unset fields rather than returning null, and `genetic_code` is simply
    absent from a real `/genome/` record -- verified live against 1125630.4, which returns
    68 fields and not that one. The fixture used to supply it there, so the section was
    populated under test and `null` on every real run: a contract field present on the CGA
    route and missing on this one.
    """
    from s2f.m1_genome.collect import genome_section

    assert "genetic_code" not in kp_payload["genome"], (
        "the fixture's genome record has regained `genetic_code`. The live core does not "
        "return it, so putting it back makes this test pass while real runs write null."
    )

    section = genome_section(kp_payload["genome"], api_bundle, "explicit genome_id")
    assert section["taxonomy"]["genetic_code"] == 11, (
        "genome.taxonomy.genetic_code is not being set from the taxonomy record. The CGA "
        "route fills it from its taxon call, so losing it here breaks the claim that both "
        "M1 routes write an indistinguishable contract."
    )


def test_the_taxonomy_query_actually_selects_genetic_code(
        kp_payload, fake_session, make_api) -> None:
    """Closes the blind spot that hid the bug: the double returns unselected fields.

    `FakeSession` replies with the whole fixture row whatever `select(...)` names, so
    removing a field from a select list cannot fail a fixture-only assertion. This drives
    the collector through a double that strips anything the query did not ask for, which
    is what the real server does.
    """
    import re

    from s2f.m1_genome.collect import collect_taxonomy

    base = type(fake_session())        # FakeSession, without importing from conftest

    class SelectHonouringTaxonomy(base):
        """Returns only the fields the query selected, as the real API does."""

        def get(self, url, params=None, timeout=None, headers=None):
            response = super().get(url, params, timeout, headers)
            if "/taxonomy/" not in url:
                return response
            selected = re.search(r"select\(([^)]*)\)", url)
            if selected and isinstance(response._body, list):
                keep = {name.strip() for name in selected.group(1).split(",")}
                response._body = [{k: v for k, v in row.items() if k in keep}
                                  for row in response._body]
            return response

    api, _ = make_api(SelectHonouringTaxonomy(kp_payload))
    taxonomy = collect_taxonomy(api, kp_payload["genome"])

    assert taxonomy["genetic_code"] == 11, (
        "collect_taxonomy did not ask for `genetic_code` in its select(...) list, so a "
        "real server returns the row without it and the field silently becomes null."
    )
    assert taxonomy["lineage"], "the select-honouring double broke the lineage as well"


def test_both_m1_routes_declare_their_own_annotation_route(api_bundle, kp_payload) -> None:
    """The field that distinguishes the two routes has to be set by *both* of them.

    It was set only by the API route, so `annotation_route == "cga"` matched nothing and
    anything reading `genome` had to know to fall back to `cga_job_id`. Both report
    renderers do exactly that (`report_md.py`, `report_html.py`), which is why the
    asymmetry went unnoticed: the fallback made the output look right. That fallback is
    kept for runs written before this, but a new run on either route now states its own
    route rather than leaving it to be inferred.
    """
    from pathlib import Path

    from s2f.m1_genome import parse
    from s2f.m1_genome.collect import genome_section as api_genome_section

    api = api_genome_section(kp_payload["genome"], api_bundle, "explicit genome_id")
    assert api["annotation_route"] == "api"

    cga_sample = Path(__file__).resolve().parents[1] / "fixtures/m1/cga_sample"
    cga = parse.genome_section(parse.load_cga(cga_sample))
    assert cga.get("annotation_route") == "cga", (
        "the CGA route does not declare `annotation_route`, so nothing reading `genome` "
        "can positively identify a CGA run -- only infer one from `cga_job_id` being "
        "present. Both routes must state their own route."
    )
    assert api["annotation_route"] != cga["annotation_route"], (
        "both routes report the same annotation_route, which defeats the point of the "
        "field: a consumer cannot tell a submitted assembly from a reference genome."
    )

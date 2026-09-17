"""The gene-content phylogeny, the pathogen knowledge base, and the M1 CLI end to end."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

from s2f.m1_genome import collect
from s2f.m1_genome.pathogens import (BASE_REFERENCES, build_reference_list,
                                     resolve_disease_profile, species_disease_label)
from s2f.m1_genome.tree import build_gene_content_tree, linkage_to_newick, sanitize_newick

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
FIXTURE = FIXTURES / "m1/bvbrc_api/klebsiella_hs11286.json"

#: scipy is optional, and only the tree needs it. Scoped to those tests rather than the
#: module, so an install without scipy still exercises the CLI and the knowledge base —
#: which is exactly the configuration where the graceful-degradation path matters.
needs_scipy = pytest.mark.skipif(
    importlib.util.find_spec("scipy") is None,
    reason="scipy is an optional dependency of the gene-content tree",
)


# --- newick -----------------------------------------------------------------
def test_newick_labels_drop_characters_the_format_reserves() -> None:
    assert sanitize_newick("Klebsiella pneumoniae (strain A:1)") == \
        "Klebsiella_pneumoniae__strain_A_1_"


@needs_scipy
def test_branch_lengths_are_never_negative() -> None:
    """UPGMA float noise can make a parent marginally lower than its child.

    A negative branch length is not valid Newick and downstream tree readers reject the
    whole file, so it is clamped.
    """
    from scipy.cluster.hierarchy import linkage
    import numpy as np

    matrix = linkage(np.array([0.1, 0.1, 0.1, 0.2, 0.2, 0.2]), method="average")
    newick = linkage_to_newick(matrix, ["a", "b", "c", "d"])
    assert newick.endswith(";")
    assert "-" not in newick.replace("e-", "")


# --- the tree ---------------------------------------------------------------
@needs_scipy
def test_the_tree_is_built_from_shared_protein_families(kp_payload, make_api) -> None:
    api, _ = make_api()
    tree = build_gene_content_tree(api, kp_payload["genome"], kp_payload["neighbors"],
                                  {"species": [{"species": "Klebsiella oxytoca",
                                                "human_associated_genomes": 1006}]})
    assert tree["ok"] is True
    assert tree["n_genomes"] >= 3
    assert tree["newick"].endswith(";")
    assert "Jaccard" in tree["method"]
    # The report must not present a gene-content tree as a sequence phylogeny.
    assert "not sequence alignment" in tree["caveat"]


@needs_scipy
def test_the_query_genome_is_marked_on_the_tree(kp_payload, make_api) -> None:
    api, _ = make_api()
    tree = build_gene_content_tree(api, kp_payload["genome"], kp_payload["neighbors"], {})
    queries = [leaf for leaf in tree["leaves"] if leaf["is_query"]]
    assert len(queries) == 1
    assert queries[0]["genome_id"] == kp_payload["genome"]["genome_id"]


@needs_scipy
def test_too_few_genomes_explains_itself_rather_than_failing(kp_payload, make_api) -> None:
    api, _ = make_api()
    tree = build_gene_content_tree(api, kp_payload["genome"], [], {}, max_leaves=1)
    assert tree["ok"] is False
    assert "at least 3" in tree["reason"]


# --- disease knowledge ------------------------------------------------------
def test_an_unlisted_species_is_unknown_not_declared_harmless() -> None:
    """`None` and `False` are different claims.

    Most species are simply not in the curated table; reporting that absence as
    "not a human pathogen" would be evidence of absence from an absence of evidence.
    """
    profile = resolve_disease_profile({"species": "Bacillus subtilis", "genus": "Bacillus"})
    assert profile["human_pathogen"] is None
    assert profile["source"] == "BV-BRC genome metadata only"


def test_a_curated_species_carries_its_references() -> None:
    profile = resolve_disease_profile(
        {"species": "Chlamydia trachomatis", "genus": "Chlamydia"})
    assert profile["human_pathogen"] is True
    assert "elwell2016" in profile["refs"]
    assert "curated knowledge base" in profile["source"]


def test_every_profile_states_that_it_describes_a_species_not_an_isolate() -> None:
    for species in ["Chlamydia trachomatis", "Klebsiella pneumoniae", ""]:
        profile = resolve_disease_profile({"species": species, "genus": ""})
        assert "not this isolate" in profile["scope"]


def test_bvbrc_disease_metadata_is_used_when_nothing_is_curated() -> None:
    profile = resolve_disease_profile(
        {"species": "Klebsiella pneumoniae", "genus": "Klebsiella",
         "disease": ["Nosocomial infection"]})
    assert profile["diseases"] == ["Nosocomial infection"]
    assert profile["genus_hint"] == "pneumonia / sepsis (nosocomial)"


def test_a_leaf_label_falls_back_through_three_tiers() -> None:
    curated, _ = species_disease_label("Chlamydia trachomatis")
    assert curated.startswith("Trachoma")

    faceted, _ = species_disease_label("Klebsiella oxytoca",
                                       facet_lookup=lambda s: "wound infection")
    assert faceted == "wound infection"

    hinted, _ = species_disease_label("Klebsiella oxytoca")
    assert hinted == "pneumonia / sepsis (nosocomial)"


def test_references_are_numbered_without_gaps_and_include_the_genome_paper() -> None:
    profile = resolve_disease_profile(
        {"species": "Chlamydia trachomatis", "genus": "Chlamydia"})
    refs = build_reference_list(profile, {"publication": "22408243"})
    assert [r["n"] for r in refs] == list(range(1, len(refs) + 1))
    assert len(refs) == len(BASE_REFERENCES) + 2 + 1  # base + 2 curated + 1 PMID
    assert refs[-1]["url"].endswith("/22408243/")
    assert all(r["text"] and r["url"] for r in refs)


# --- the CLI end to end -----------------------------------------------------
def test_the_cga_route_writes_priorities_and_a_report(tmp_path) -> None:
    from s2f.common.schema import validate_report
    from s2f.m1_genome.__main__ import main

    run_dir = tmp_path / "run"
    code = main(["--run", str(run_dir),
                 "--from-cga-dir", str(FIXTURES / "m1/cga_sample"), "--html"])
    assert code == 0

    report = json.loads((run_dir / "report.json").read_text())
    validate_report(report)
    assert report["run"]["modules"]["m1"]["route"] == "cga"
    assert all(p.get("m1_priority") for p in report["proteins"])
    # M1 must never write M2's key.
    assert not any("triage" in p for p in report["proteins"])
    assert (run_dir / "m1/report.html").exists()
    for name in ("proteins.faa", "genes_proteins.csv", "specialty_genes_all.csv"):
        assert (run_dir / "m1" / name).exists()


def test_non_cds_features_do_not_consume_rank_slots(tmp_path) -> None:
    """The CGA route parses rRNA and tRNA rows too; only CDS rows become proteins.

    Ranking every feature let a non-CDS row take rank 1 and then be dropped, leaving
    `proteins[]` with no rank 1 and a hole in the report's ranked table.
    """
    from s2f.m1_genome.__main__ import _rank_for

    proteins = [{"feature_id": "fig|1.1.peg.1"}]
    features = [
        {"patric_id": "fig|1.1.rna.1", "product": "hemolysin toxin secretion system",
         "aa_length": 900, "feature_type": "rRNA"},
        {"patric_id": "fig|1.1.peg.1", "product": "hypothetical protein",
         "aa_length": 100, "feature_type": "CDS"},
    ]
    ranked = _rank_for(proteins, features, [])
    assert set(ranked) == {"fig|1.1.peg.1"}
    assert ranked["fig|1.1.peg.1"]["rank"] == 1


def test_a_failed_optional_request_still_leaves_a_complete_report(kp_payload, tmp_path, fake_session,
                                                                  monkeypatch) -> None:
    """The contract must be written before anything decorative is attempted.

    A 503 on the phylogeny's protein-family query used to abort the run with `report.json`
    holding no `genome` and no `proteins` — the opposite of what `tree.py` promises.
    """
    from s2f.common.http import HttpError
    from s2f.common.schema import validate_report
    from s2f.m1_genome.__main__ import main
    session = fake_session()
    from s2f.common.http import CachedJsonClient as RealClient

    def fake_client(**kwargs):
        kwargs["session"] = session
        kwargs["min_interval_seconds"] = 0
        return RealClient(**kwargs)

    monkeypatch.setattr("s2f.m1_genome.__main__.CachedJsonClient", fake_client)

    def exploding_pgfams(*_args, **_kwargs):
        raise HttpError("GET .../genome_feature/ failed — HTTP 503")

    monkeypatch.setattr("s2f.m1_genome.tree.fetch_pgfam_set", exploding_pgfams)

    run_dir = tmp_path / "degraded"
    assert main(["--run", str(run_dir), "--from-bvbrc-api",
                 "--genome-id", "1125630.4", "--html"]) == 0

    report = json.loads((run_dir / "report.json").read_text())
    validate_report(report)
    assert len(report["proteins"]) == len(kp_payload["features"])
    assert report["genome"]["annotation_route"] == "api"
    assert (run_dir / "m1/report.html").exists()
    assert not (run_dir / "m1/gene_content_tree.nwk").exists()


def test_the_section_tables_and_taxonomy_tree_are_written(kp_payload, tmp_path, fake_session,
                                                          monkeypatch) -> None:
    """The colleague's per-section CSVs are a deliverable, not a redundant format."""
    from s2f.m1_genome.__main__ import main
    from s2f.common.http import CachedJsonClient as RealClient

    session = fake_session()
    monkeypatch.setattr("s2f.m1_genome.__main__.CachedJsonClient",
                        lambda **kw: RealClient(**{**kw, "session": session,
                                                   "min_interval_seconds": 0}))
    run_dir = tmp_path / "tables"
    assert main(["--run", str(run_dir), "--from-bvbrc-api",
                 "--genome-id", "1125630.4", "--no-tree"]) == 0

    m1 = run_dir / "m1"
    for name in ("amr_genes.csv", "virulence_factors.csv", "close_human_pathogens.csv",
                 "growth_conditions.csv", "isolation_genome.csv",
                 "nutrition_biosynthesis.csv", "proteins_ranked.csv",
                 "pathogenesis_host_invasion.csv", "taxonomy_tree.txt"):
        assert (m1 / name).exists(), name

    # proteins_ranked.csv is ordered by rank and carries the hypothesis.
    ranked = (m1 / "proteins_ranked.csv").read_text().splitlines()
    assert ranked[0].startswith("rank,score,patric_id")
    assert [line.split(",")[0] for line in ranked[1:]] == \
        [str(i) for i in range(1, len(kp_payload["features"]) + 1)]

    tree_text = (m1 / "taxonomy_tree.txt").read_text()
    assert "Bacteria" in tree_text
    assert "Klebsiella pneumoniae subsp. pneumoniae HS11286" in tree_text


def test_sequences_are_capped_by_priority_not_by_chromosome_position(kp_payload,
                                                                     make_api) -> None:
    """A coordinate-ordered cap can drop every top candidate from M2's input."""
    from s2f.m1_genome.collect import fetch_sequences
    from s2f.m1_genome.priority import rank_proteins
    api, _ = make_api()
    ranked = rank_proteins(kp_payload["features"], kp_payload["specialty"])
    top = min(ranked, key=lambda fid: ranked[fid]["rank"])

    # The top-ranked protein is not first by coordinate, so an unprioritised cap of 1
    # would miss it.
    by_coordinate = sorted(kp_payload["features"], key=lambda f: f.get("start") or 0)
    assert by_coordinate[0]["patric_id"] != top

    got = fetch_sequences(api, kp_payload["features"], cap=1, priority=ranked)
    assert list(got) == [top]


def test_the_two_routes_are_mutually_exclusive(tmp_path) -> None:
    from s2f.m1_genome.__main__ import main

    assert main(["--run", str(tmp_path / "r"), "--from-bvbrc-api",
                 "--from-cga-dir", str(FIXTURES / "m1/cga_sample")]) == 2


def test_the_api_route_needs_something_to_resolve_from(tmp_path) -> None:
    from s2f.m1_genome.__main__ import main

    assert main(["--run", str(tmp_path / "r"), "--from-bvbrc-api"]) == 2


def test_the_api_route_writes_the_same_contract_as_the_cga_route(kp_payload, tmp_path, fake_session,
                                                                monkeypatch) -> None:
    """End to end on the API route, with the fake session standing in for BV-BRC."""
    from s2f.common.schema import validate_report
    from s2f.m1_genome import bvbrc_api
    from s2f.m1_genome.__main__ import main
    session = fake_session()
    real_client = bvbrc_api.CachedJsonClient

    def fake_client(**kwargs):
        kwargs["session"] = session
        kwargs["min_interval_seconds"] = 0
        return real_client(**kwargs)

    monkeypatch.setattr("s2f.m1_genome.__main__.CachedJsonClient", fake_client)

    run_dir = tmp_path / "api_run"
    code = main(["--run", str(run_dir), "--from-bvbrc-api", "--genome-id", "1125630.4",
                 "--html", "--no-tree"])
    assert code == 0

    report = json.loads((run_dir / "report.json").read_text())
    validate_report(report)
    assert report["genome"]["annotation_route"] == "api"
    assert report["run"]["modules"]["m1"]["route"] == "api"
    assert len(report["proteins"]) == len(kp_payload["features"])
    assert all(p.get("m1_priority") for p in report["proteins"])
    assert (run_dir / "m1/report.html").exists()
    assert (run_dir / "m1/disease_profile.json").exists()

    # And M2 reads it without knowing which route ran.
    from s2f.m2_triage.bvbrc_input import load_input
    assert len(load_input(run_dir / "m1").proteins) == len(kp_payload["features"])


# The workspace-upload filename invariants live in `test_m1_safety_invariants.py`:
# they guard a silent failure that was reverted once together with the tests beside it,
# so they are kept apart from the module they protect (docs/pitfalls.md #23).

"""Gene-content phylogeny: Jaccard distance on shared PGFam families, clustered by UPGMA.

This is **not** a phylogeny in the sequence-alignment sense and the report says so. It
clusters genomes by which protein families they encode (Snel, Bork & Huynen 1999), which
answers "which relatives have a similar gene repertoire" — a useful question when the point
is host interaction, and a fast one, since PGFam membership is already computed in BV-BRC.

The definitive tree for a run is CGA's codon tree (`genome.tree_newick` from the CGA route).
Where both exist the codon tree wins; this one is for the API route, which has no codon
tree, and as a gene-content view alongside it.

scipy and numpy are optional. Without them the function returns `{"ok": False, "reason":
...}` and the report renders a note instead of a tree — a missing figure is a far better
outcome than a failed run, and M1's actual contract (`report.json`, `<run>/m1/`) does not
depend on this at all.
"""

from __future__ import annotations

import re
from typing import Any

from .bvbrc_api import BvbrcApi, rql_value, top_nonempty
from .pathogens import CURATED_PATHOGEN_KB, species_disease_label


def fetch_pgfam_set(api: BvbrcApi, genome_id: str, *, cap: int = 8000) -> set[str]:
    """The set of PGFam ids encoded by one genome."""
    rows = api.query_all(
        "genome_feature",
        f"and(eq(genome_id,{rql_value(genome_id)}),eq(feature_type,CDS),"
        "eq(annotation,PATRIC))&select(pgfam_id)",
        cap=cap,
    )
    return {row["pgfam_id"] for row in rows if row.get("pgfam_id")}


def sanitize_newick(name: str) -> str:
    """Newick forbids parentheses, commas, colons and whitespace in a label."""
    return re.sub(r"[^A-Za-z0-9_.]", "_", str(name))


def linkage_to_newick(linkage_matrix, labels: list[str]) -> str:
    """Convert a scipy linkage matrix to a Newick string with branch lengths."""
    from scipy.cluster.hierarchy import to_tree

    tree = to_tree(linkage_matrix, rd=False)

    def walk(node, parent_distance: float) -> str:
        # A branch length is the drop in merge height from parent to child; clamped at
        # zero because floating-point noise in UPGMA can make it very slightly negative,
        # and a negative branch length is not valid Newick.
        branch = max(parent_distance - node.dist, 0.0)
        if node.is_leaf():
            return f"{sanitize_newick(labels[node.id])}:{branch:.4f}"
        return f"({walk(node.left, node.dist)},{walk(node.right, node.dist)}):{branch:.4f}"

    return f"({walk(tree.left, tree.dist)},{walk(tree.right, tree.dist)});"


def _select_genomes(api: BvbrcApi, genome: dict[str, Any], neighbors: list[dict[str, Any]],
                    close_pathogens: dict[str, Any], max_leaves: int) -> list[dict[str, Any]]:
    """Query genome first, then close human-associated species, then genus neighbours.

    The order is the point: with a leaf budget, the genomes worth spending it on are the
    ones that make the tree answer the question — which nearby species infect humans.
    """
    query_species = genome.get("species")
    chosen = [{"genome_id": genome["genome_id"], "genome_name": genome.get("genome_name"),
               "species": query_species, "is_query": True}]
    seen = {query_species}

    for row in close_pathogens.get("species") or []:
        if len(chosen) >= max_leaves:
            break
        species = row.get("species")
        if not species or species in seen:
            continue
        hits = api.query(
            "genome",
            f"and(eq(species,{rql_value(species)}),or(eq(reference_genome,Reference),"
            "eq(reference_genome,Representative)))&select(genome_id,genome_name,species)&limit(1)",
        ) or api.query(
            "genome",
            f"eq(species,{rql_value(species)})&select(genome_id,genome_name,species)&limit(1)",
        )
        if hits:
            chosen.append({**hits[0], "is_query": False})
            seen.add(species)

    for neighbor in neighbors:
        if len(chosen) >= max_leaves:
            break
        species = neighbor.get("species")
        if species in seen:
            continue
        chosen.append({"genome_id": neighbor["genome_id"],
                       "genome_name": neighbor.get("genome_name"),
                       "species": species, "is_query": False})
        seen.add(species)
    return chosen


def build_gene_content_tree(api: BvbrcApi, genome: dict[str, Any],
                            neighbors: list[dict[str, Any]],
                            close_pathogens: dict[str, Any], *,
                            max_leaves: int = 14) -> dict[str, Any]:
    """Build the tree, or explain in `reason` why it could not be built."""
    try:
        import numpy as np
        from scipy.cluster.hierarchy import dendrogram, linkage
        from scipy.spatial.distance import squareform
    except ImportError as exc:
        return {"ok": False, "reason": f"scipy/numpy not installed ({exc})"}

    query_species = genome.get("species")
    human_species = {row["species"] for row in close_pathogens.get("species") or []
                     if row.get("species")}
    if CURATED_PATHOGEN_KB.get(query_species or "", {}).get("human"):
        human_species.add(query_species)

    chosen = _select_genomes(api, genome, neighbors, close_pathogens, max_leaves)
    for entry in chosen:
        entry["pgfams"] = fetch_pgfam_set(api, entry["genome_id"])
    chosen = [entry for entry in chosen if entry["pgfams"]]
    if len(chosen) < 3:
        return {"ok": False,
                "reason": f"only {len(chosen)} genome(s) had protein-family data; "
                          "a tree needs at least 3"}

    size = len(chosen)
    distances = np.zeros((size, size))
    for i in range(size):
        for j in range(i + 1, size):
            left, right = chosen[i]["pgfams"], chosen[j]["pgfams"]
            union = len(left | right)
            distance = 1.0 - (len(left & right) / union if union else 0.0)
            distances[i, j] = distances[j, i] = distance

    labels = [entry["species"] or entry["genome_name"] for entry in chosen]
    matrix = linkage(squareform(distances, checks=False), method="average")  # UPGMA
    drawn = dendrogram(matrix, labels=labels, no_plot=True)

    def facet_lookup(species: str) -> str | None:
        value, _ = top_nonempty(
            api.facet("genome", f"eq(species,{rql_value(species)})", "disease", limit=6)
        )
        return value

    leaf_meta: dict[str, dict[str, Any]] = {}
    for entry in chosen:
        label = entry["species"] or entry["genome_name"]
        disease, is_pathogen = species_disease_label(
            entry.get("species"), facet_lookup=facet_lookup, human_species=human_species
        )
        leaf_meta[label] = {
            "genome_id": entry["genome_id"],
            "is_query": entry["is_query"],
            "is_pathogen": bool(is_pathogen) or entry.get("species") in human_species,
            "disease": disease,
            "name": label,
        }

    return {
        "ok": True,
        "method": ("Gene-content phylogeny — Jaccard distance on shared PGFam protein "
                   "families, clustered by UPGMA"),
        "caveat": ("gene content, not sequence alignment; the definitive tree for a run "
                   "is the CGA codon tree"),
        "newick": linkage_to_newick(matrix, labels),
        "icoord": drawn["icoord"],
        "dcoord": drawn["dcoord"],
        "ivl": drawn["ivl"],
        "leaf_meta": leaf_meta,
        "leaves": [leaf_meta[label] for label in drawn["ivl"]],
        "n_genomes": size,
    }

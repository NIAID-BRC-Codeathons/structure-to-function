"""Knowledge-graph assembly for M2 (issue #11).

A **small** provenance-carrying subgraph, queried from existing APIs. Nothing here ingests or
hosts a knowledge graph (pitfall #10): every source is called for the proteins we selected, the
results are capped, and what we keep is the subgraph plus the evidence that justifies it.

The useful knowledge about drugs and disease lives on characterized proteins in other organisms,
so most edges are reached through a homolog. That hop is always explicit
(``homolog_of`` carrying its identity, then ``target_of``) and never collapsed into "our protein
is a drug target" — pitfall #4 and the edge-typing rule in ``02-m2-triage.md``.

Sources, as verified against our genomes:

- **STRING** — works directly on BV-BRC locus tags, so no mapping is needed. Gives
  ``interacts_with`` edges within the genome.
- **ChEMBL** — reached through a homolog's UniProt accession. Gives ``target_of`` edges from
  compounds with measured potency.
- **BV-BRC specialty table / PDB hits** — the ``homolog_of`` edges themselves, with identity.

Open Targets is human-only (searching "gyrase" returns TOP2A), so it is reached the same way if
added later. HPIDB was unreachable when this was written.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Iterable

from ..common.http import CachedJsonClient, HttpError

STRING_API = "https://string-db.org/api/json"
CHEMBL_API = "https://www.ebi.ac.uk/chembl/api/data"

# Caps keep the subgraph small (pitfall #10). Recorded in the output so a reader knows the
# graph is bounded by choice, not by what happened to exist.
MAX_STRING_PARTNERS = 10
MAX_CHEMBL_TARGETS = 5
MAX_CHEMBL_COMPOUNDS = 10
MIN_STRING_SCORE = 700  # STRING's "high confidence" threshold, 0-1000

EDGE_INTERACTS = "interacts_with"
EDGE_HOMOLOG = "homolog_of"
EDGE_TARGET_OF = "target_of"
EDGE_CHEMBL_TARGET = "represented_by"

_SLUG = re.compile(r"[^a-z0-9]+")


def _slug(value: str) -> str:
    return _SLUG.sub("-", (value or "").lower()).strip("-")[:60]


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


@dataclass
class Node:
    id: str
    type: str
    name: str = ""
    source: str = ""
    attributes: dict[str, Any] = field(default_factory=dict)


@dataclass
class Edge:
    """One typed edge. ``evidence`` says why it exists; ``provenance`` says where it came from."""

    source_id: str
    target_id: str
    type: str
    evidence: dict[str, Any] = field(default_factory=dict)
    provenance: dict[str, Any] = field(default_factory=dict)


@dataclass
class KnowledgeGraph:
    nodes: dict[str, Node] = field(default_factory=dict)
    edges: list[Edge] = field(default_factory=list)
    source_versions: dict[str, str] = field(default_factory=dict)
    failures: list[dict[str, str]] = field(default_factory=list)

    def add_node(self, node: Node) -> str:
        existing = self.nodes.get(node.id)
        if existing is None:
            self.nodes[node.id] = node
        elif not existing.name and node.name:
            existing.name = node.name
        return node.id

    def add_edge(self, edge: Edge) -> None:
        self.edges.append(edge)

    def to_dict(self, *, caps: dict[str, Any] | None = None) -> dict[str, Any]:
        return {
            "schema_version": "0.1",
            "generated_utc": _now(),
            "note": (
                "Subgraph queried from public APIs and capped; no knowledge graph is ingested or "
                "hosted here. Edges through a homolog carry the identity that justifies them; a "
                "homolog's annotation is not our protein's function."
            ),
            "caps": caps or {},
            "source_versions": self.source_versions,
            "counts": self.summary(),
            "nodes": [
                {
                    "id": n.id,
                    "type": n.type,
                    "name": n.name,
                    "source": n.source,
                    **({"attributes": n.attributes} if n.attributes else {}),
                }
                for n in self.nodes.values()
            ],
            "edges": [
                {
                    "source_id": e.source_id,
                    "target_id": e.target_id,
                    "type": e.type,
                    "evidence": e.evidence,
                    "provenance": e.provenance,
                }
                for e in self.edges
            ],
            "failures": self.failures,
        }

    def summary(self) -> dict[str, Any]:
        node_types: dict[str, int] = {}
        for node in self.nodes.values():
            node_types[node.type] = node_types.get(node.type, 0) + 1
        edge_types: dict[str, int] = {}
        for edge in self.edges:
            edge_types[edge.type] = edge_types.get(edge.type, 0) + 1
        return {
            "nodes": len(self.nodes),
            "edges": len(self.edges),
            "node_types": dict(sorted(node_types.items())),
            "edge_types": dict(sorted(edge_types.items())),
        }

    def to_networkx(self):
        """Assemble as a NetworkX MultiDiGraph for anyone who wants to traverse it."""
        import networkx as nx

        graph = nx.MultiDiGraph()
        for node in self.nodes.values():
            graph.add_node(node.id, type=node.type, name=node.name, source=node.source, **node.attributes)
        for edge in self.edges:
            graph.add_edge(
                edge.source_id, edge.target_id, key=edge.type, type=edge.type,
                evidence=edge.evidence, provenance=edge.provenance,
            )
        return graph


def protein_node_id(feature_id: str) -> str:
    return f"brc:{feature_id}"


class KnowledgeGraphBuilder:
    """Query the sources for one set of proteins and assemble the subgraph."""

    def __init__(
        self,
        client: CachedJsonClient,
        *,
        taxon_id: int | None = None,
        max_partners: int = MAX_STRING_PARTNERS,
        max_targets: int = MAX_CHEMBL_TARGETS,
        max_compounds: int = MAX_CHEMBL_COMPOUNDS,
        min_string_score: int = MIN_STRING_SCORE,
    ) -> None:
        self.client = client
        self.taxon_id = taxon_id
        self.max_partners = max_partners
        self.max_targets = max_targets
        self.max_compounds = max_compounds
        self.min_string_score = min_string_score
        self.graph = KnowledgeGraph()

    # --- provenance --------------------------------------------------------

    def _provenance(self, source: str, method: str, params: dict[str, Any]) -> dict[str, Any]:
        return {
            "source": source,
            "source_version": self.graph.source_versions.get(source, ""),
            "retrieved_at": _now(),
            "method": method,
            "params": params,
        }

    def string_version(self) -> str:
        if "STRING" in self.graph.source_versions:
            return self.graph.source_versions["STRING"]
        try:
            payload = self.client.get_json("string_version", f"{STRING_API}/version")
        except HttpError as exc:
            self.graph.failures.append({"source": "STRING", "step": "version", "error": str(exc)})
            return ""
        version = ""
        if isinstance(payload, list) and payload:
            version = str(payload[0].get("string_version", ""))
        self.graph.source_versions["STRING"] = version
        return version

    # --- STRING ------------------------------------------------------------

    def add_string_partners(self, feature_id: str, locus_tag: str) -> int:
        """interacts_with edges. STRING accepts BV-BRC locus tags directly for our genomes."""
        if not locus_tag or not self.taxon_id:
            return 0
        self.string_version()
        params = {
            "identifiers": locus_tag,
            "species": self.taxon_id,
            "limit": self.max_partners,
            "required_score": self.min_string_score,
            "caller_identity": "s2f-m2-kg",
        }
        try:
            payload = self.client.get_json("string_partners", f"{STRING_API}/interaction_partners", params)
        except HttpError as exc:
            self.graph.failures.append(
                {"source": "STRING", "step": "interaction_partners", "feature_id": feature_id, "error": str(exc)}
            )
            return 0
        if not isinstance(payload, list):
            return 0

        source_id = self.graph.add_node(
            Node(protein_node_id(feature_id), "protein", locus_tag, "BV-BRC")
        )
        added = 0
        for row in payload[: self.max_partners]:
            partner = str(row.get("preferredName_B") or row.get("stringId_B") or "")
            if not partner:
                continue
            partner_id = f"string:{row.get('stringId_B', partner)}"
            self.graph.add_node(Node(partner_id, "protein", partner, "STRING"))
            self.graph.add_edge(
                Edge(
                    source_id,
                    partner_id,
                    EDGE_INTERACTS,
                    evidence={
                        "combined_score": row.get("score"),
                        "experimental_score": row.get("escore"),
                        "database_score": row.get("dscore"),
                        "textmining_score": row.get("tscore"),
                        "note": "STRING score is aggregated evidence, not a measured interaction",
                    },
                    provenance=self._provenance("STRING", "interaction_partners", params),
                )
            )
            added += 1
        return added

    # --- homologs ----------------------------------------------------------

    def add_homolog(
        self,
        feature_id: str,
        accession: str,
        *,
        identity: float | None,
        coverage: float | None = None,
        organism: str = "",
        via: str = "",
        name: str = "",
    ) -> str | None:
        """homolog_of edge carrying the identity that justifies everything downstream."""
        if not accession:
            return None
        source_id = self.graph.add_node(Node(protein_node_id(feature_id), "protein", "", "BV-BRC"))
        target_id = self.graph.add_node(
            Node(f"uniprot:{accession}", "protein", name, "UniProt", {"organism": organism})
        )
        self.graph.add_edge(
            Edge(
                source_id,
                target_id,
                EDGE_HOMOLOG,
                evidence={
                    "percent_identity": identity,
                    "coverage": coverage,
                    "organism": organism,
                    "via": via,
                    "note": "the homolog's annotation is not this protein's function",
                },
                provenance=self._provenance(via or "homology", "sequence comparison", {}),
            )
        )
        return target_id

    def add_named_human_homolog(
        self, feature_id: str, product: str, identity: float | None, coverage: float | None = None
    ) -> str:
        """A human homolog from BV-BRC's specialty table, which gives a name but no accession.

        Recorded because M2 scoring penalizes it (pitfall #4) and M6 must state the selectivity
        risk. The node is explicitly unaccessioned rather than guessed at.
        """
        source_id = self.graph.add_node(Node(protein_node_id(feature_id), "protein", "", "BV-BRC"))
        target_id = self.graph.add_node(
            Node(
                f"human:{_slug(product) or _slug(feature_id)}",
                "protein",
                product,
                "BV-BRC specialty genes",
                {"organism": "Homo sapiens", "accession": None},
            )
        )
        self.graph.add_edge(
            Edge(
                source_id,
                target_id,
                EDGE_HOMOLOG,
                evidence={
                    "percent_identity": identity,
                    "coverage": coverage,
                    "organism": "Homo sapiens",
                    "via": "BV-BRC specialty genes",
                    "note": "selectivity risk for an antibacterial target; scoring penalty, never a bonus",
                },
                provenance=self._provenance("BV-BRC", "specialty gene table", {"property": "Human Homolog"}),
            )
        )
        return target_id

    # --- ChEMBL ------------------------------------------------------------

    def add_chembl_for_accession(self, accession: str) -> int:
        """target_of edges for compounds measured against a homolog's accession.

        Note the filter name: ``target_components__accession`` filters; ``target_organism__*``
        is silently ignored by the API and returns everything.
        """
        if not accession:
            return 0
        params = {"target_components__accession": accession, "limit": self.max_targets, "format": "json"}
        try:
            payload = self.client.get_json("chembl_target", f"{CHEMBL_API}/target.json", params)
        except HttpError as exc:
            self.graph.failures.append(
                {"source": "ChEMBL", "step": "target", "accession": accession, "error": str(exc)}
            )
            return 0

        protein_node = f"uniprot:{accession}"
        added = 0
        for target in (payload.get("targets") or [])[: self.max_targets]:
            target_chembl_id = str(target.get("target_chembl_id") or "")
            if not target_chembl_id:
                continue
            # Guard against a filter that did not apply: the accession must really be a component.
            components = {
                str(component.get("accession"))
                for component in target.get("target_components") or []
            }
            if components and accession not in components:
                self.graph.failures.append(
                    {
                        "source": "ChEMBL",
                        "step": "target",
                        "accession": accession,
                        "error": f"{target_chembl_id} does not list {accession} as a component; filter did not apply",
                    }
                )
                continue

            target_node = self.graph.add_node(
                Node(
                    f"chembl:{target_chembl_id}",
                    "chembl_target",
                    str(target.get("pref_name") or ""),
                    "ChEMBL",
                    {"organism": target.get("organism"), "target_type": target.get("target_type")},
                )
            )
            self.graph.add_node(Node(protein_node, "protein", "", "UniProt"))
            self.graph.add_edge(
                Edge(
                    protein_node,
                    target_node,
                    EDGE_CHEMBL_TARGET,
                    evidence={"accession": accession},
                    provenance=self._provenance("ChEMBL", "target", params),
                )
            )
            added += self._add_chembl_activities(target_chembl_id, target_node)
        return added

    def _add_chembl_activities(self, target_chembl_id: str, target_node: str) -> int:
        params = {
            "target_chembl_id": target_chembl_id,
            "pchembl_value__isnull": "false",
            "limit": self.max_compounds,
            "format": "json",
        }
        try:
            payload = self.client.get_json("chembl_activity", f"{CHEMBL_API}/activity.json", params)
        except HttpError as exc:
            self.graph.failures.append(
                {"source": "ChEMBL", "step": "activity", "target": target_chembl_id, "error": str(exc)}
            )
            return 0

        added = 0
        for activity in (payload.get("activities") or [])[: self.max_compounds]:
            molecule = str(activity.get("molecule_chembl_id") or "")
            if not molecule:
                continue
            compound_node = self.graph.add_node(
                Node(f"compound:{molecule}", "compound", molecule, "ChEMBL")
            )
            self.graph.add_edge(
                Edge(
                    compound_node,
                    target_node,
                    EDGE_TARGET_OF,
                    evidence={
                        "standard_type": activity.get("standard_type"),
                        "standard_value": activity.get("standard_value"),
                        "standard_units": activity.get("standard_units"),
                        "pchembl_value": activity.get("pchembl_value"),
                        "assay_chembl_id": activity.get("assay_chembl_id"),
                        "note": "measured against this target, which is a homolog of our protein",
                    },
                    provenance=self._provenance("ChEMBL", "activity", params),
                )
            )
            added += 1
        return added

    # --- orchestration -----------------------------------------------------

    def build(self, proteins: Iterable[dict[str, Any]]) -> KnowledgeGraph:
        """Build the subgraph for the selected proteins.

        Each item: ``feature_id``, ``locus_tag``, ``product``, ``homologs`` (accession, identity,
        coverage, organism, via, name), ``human_homolog_identity``.
        """
        for protein in proteins:
            feature_id = str(protein.get("feature_id") or "")
            if not feature_id:
                continue
            self.graph.add_node(
                Node(
                    protein_node_id(feature_id),
                    "protein",
                    str(protein.get("product") or ""),
                    "BV-BRC",
                    {"locus_tag": protein.get("locus_tag"), "organism_taxon_id": self.taxon_id},
                )
            )
            self.add_string_partners(feature_id, str(protein.get("locus_tag") or ""))

            for homolog in protein.get("homologs") or []:
                accession = str(homolog.get("accession") or "")
                if not accession:
                    continue
                self.add_homolog(
                    feature_id,
                    accession,
                    identity=homolog.get("identity"),
                    coverage=homolog.get("coverage"),
                    organism=str(homolog.get("organism") or ""),
                    via=str(homolog.get("via") or ""),
                    name=str(homolog.get("name") or ""),
                )
                self.add_chembl_for_accession(accession)

            identity = protein.get("human_homolog_identity")
            if identity:
                self.add_named_human_homolog(
                    feature_id, str(protein.get("product") or ""), float(identity)
                )
        return self.graph

    def caps(self) -> dict[str, Any]:
        return {
            "max_string_partners_per_protein": self.max_partners,
            "min_string_score": self.min_string_score,
            "max_chembl_targets_per_accession": self.max_targets,
            "max_chembl_compounds_per_target": self.max_compounds,
        }

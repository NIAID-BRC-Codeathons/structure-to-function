"""Turn `report.json` into a view model the template can print (issue #66).

Every number on the page has to trace to a `report.json` field, and **nothing is computed in
the template** (`docs/06-m6-report.md`). So all arithmetic, sorting, grouping and formatting
happens here; the template only prints what this module hands it.

The report must render with sections missing — an unfinished module shows as "not run", never
as an error, and never as a zero. That distinction is the same one the pipeline makes
everywhere else: absent evidence and measured-absent evidence are different claims.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import UTC, datetime
from typing import Any

#: Sections the report knows how to show, in the order `docs/06-m6-report.md` lays out.
SECTIONS = ("run", "genome", "proteins", "structures", "kg", "ligands", "docking", "disease")

#: Fixed, per `docs/06-m6-report.md` and pitfalls #1, #2, #5, #6 and #9. Not conditional on
#: what ran: a reader must see the limits whether or not the module that earns them exists.
LIMITATIONS = [
    "Docking scores are a ranking device within one protein, not binding affinities, and are "
    "never comparable across proteins.",
    "Predicted structures carry no cofactors, metals or ligands, so a pocket that needs one is "
    "incomplete.",
    "Nothing here is experimentally validated. Candidates are computational hypotheses for "
    "wet-lab triage.",
    "No clinical or treatment implication is intended or supported.",
    "Function transferred from a homolog is the homolog's evidence, not this protein's. The "
    "identity that justified each transfer is recorded alongside it.",
]


@dataclass
class Section:
    """One section's status, so 'not run' is explicit rather than an empty table."""

    key: str
    present: bool
    count: int | None = None
    note: str = ""

    @property
    def status(self) -> str:
        return "ran" if self.present else "not run"


@dataclass
class ReportView:
    title: str = "Structure to Function"
    run_id: str = ""
    generated_at: str = ""
    organism: str = ""
    organism_source: str = ""
    taxon_id: int | None = None
    genome_id: str = ""
    sections: list[Section] = field(default_factory=list)
    counts: dict[str, Any] = field(default_factory=dict)
    proteins: list[dict[str, Any]] = field(default_factory=list)
    components: list[str] = field(default_factory=list)
    weights: dict[str, float] = field(default_factory=dict)
    structures_by_feature: dict[str, dict[str, Any]] = field(default_factory=dict)
    candidates_by_protein: list[dict[str, Any]] = field(default_factory=list)
    kg_summary: dict[str, Any] = field(default_factory=dict)
    disease_claims: list[dict[str, Any]] = field(default_factory=list)
    provenance: list[dict[str, str]] = field(default_factory=list)
    caveats: list[str] = field(default_factory=list)
    limitations: list[str] = field(default_factory=lambda: list(LIMITATIONS))


def _fmt(value: Any, digits: int = 3) -> str:
    """One place that turns a number into text, so the template never formats."""
    if value is None or value == "":
        return "—"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, (int,)) and not isinstance(value, bool):
        return str(value)
    if isinstance(value, float):
        return f"{value:.{digits}g}"
    return str(value)


def _pct(value: Any) -> str:
    if value is None or value == "":
        return "—"
    try:
        return f"{float(value) * 100:.0f}%"
    except (TypeError, ValueError):
        return "—"


def _protein_row(protein: dict[str, Any], structure: dict[str, Any] | None) -> dict[str, Any]:
    triage = protein.get("triage") or {}
    flags = protein.get("flags") or {}
    xrefs = protein.get("xrefs") or {}
    experimental = flags.get("experimental_homolog") or {}
    predicted = flags.get("predicted_model") or {}
    components = triage.get("components") or {}

    return {
        "feature_id": protein.get("feature_id", ""),
        "product": protein.get("product") or "—",
        "gene": xrefs.get("gene_name") or "—",
        "rank": triage.get("rank"),
        "score": _fmt(triage.get("score")),
        "components": {name: _fmt(value) for name, value in components.items()},
        "reason": triage.get("reason") or "",
        "uniprot": xrefs.get("uniprot") or "—",
        "same_species": bool(flags.get("same_species_hit")),
        "same_genus": bool(flags.get("same_genus_hit")),
        "homolog": {
            "entity": experimental.get("pdb_entity") or "—",
            "identity": _pct(experimental.get("identity")),
            "coverage": _pct(experimental.get("coverage")),
            "resolution": _fmt(experimental.get("resolution")),
            "holo": bool(experimental.get("holo")),
            "ligands": ", ".join(experimental.get("ligands") or []) or "—",
            "organism": experimental.get("organism") or "—",
        } if experimental else None,
        "model": {
            "entry": predicted.get("entry_id") or "—",
            "plddt": _fmt(predicted.get("mean_plddt")),
            "band": predicted.get("confidence_band") or "—",
            "coverage": _pct(predicted.get("coverage")),
            "usable": bool(predicted.get("usable_for_docking")),
        } if predicted else None,
        "structure": {
            "source": structure.get("source", "—"),
            "accession": structure.get("accession") or "—",
            "confidence": _fmt(structure.get("confidence")),
            "usable": structure.get("usable_for_docking"),
            "path": structure.get("path") or "",
            "pockets": len(structure.get("pockets") or []),
        } if structure else None,
    }


def build(report: dict[str, Any], *, top: int = 25) -> ReportView:
    """Build the whole view. Missing sections become 'not run', never empty successes."""
    view = ReportView(generated_at=datetime.now(UTC).isoformat(timespec="seconds"))

    run = report.get("run") or {}
    view.run_id = str(run.get("run_id") or "—")

    genome = report.get("genome") or {}
    # `genome.taxonomy` is a dict from M1 but the schema also allows a plain string, and the
    # committed fixture uses one. Both have to render.
    taxonomy = genome.get("taxonomy")
    if isinstance(taxonomy, str):
        taxonomy = {"scientific_name": taxonomy}
    elif not isinstance(taxonomy, dict):
        taxonomy = {}
    view.organism = (taxonomy.get("scientific_name") or "").strip() or "not determined"
    view.organism_source = str(taxonomy.get("called_by") or "")
    view.taxon_id = genome.get("taxon_id")
    view.genome_id = str(genome.get("genome_id") or "—")

    proteins = report.get("proteins") or []
    structures = report.get("structures") or []
    view.structures_by_feature = {
        str(s.get("feature_id")): s for s in structures if s.get("feature_id")
    }

    for key in SECTIONS:
        value = report.get(key)
        present = value is not None and (len(value) > 0 if isinstance(value, (list, dict)) else True)
        count = len(value) if isinstance(value, (list, dict)) else None
        view.sections.append(Section(key=key, present=bool(present), count=count))

    selected = [p for p in proteins if ((p.get("triage") or {}).get("selected"))]
    ranked = sorted(
        selected or proteins,
        key=lambda p: ((p.get("triage") or {}).get("rank") or 10**9),
    )
    view.proteins = [
        _protein_row(p, view.structures_by_feature.get(str(p.get("feature_id"))))
        for p in ranked[:top]
    ]
    if view.proteins:
        view.components = list(view.proteins[0]["components"].keys())

    view.counts = {
        "proteins": len(proteins),
        "selected": len(selected),
        "structures": len(structures),
        "usable_for_docking": sum(1 for s in structures if s.get("usable_for_docking") is True),
        "same_species_hits": sum(
            1 for p in proteins if ((p.get("flags") or {}).get("same_species_hit"))
        ),
        "shown": len(view.proteins),
    }

    # Ligands and docking, grouped by protein. Never ranked across proteins (pitfall #1).
    ligands = report.get("ligands") or []
    docking = report.get("docking") or []
    by_protein: dict[str, dict[str, Any]] = {}
    for ligand in ligands:
        feature_id = str(ligand.get("feature_id") or "")
        entry = by_protein.setdefault(feature_id, {"feature_id": feature_id, "ligands": [], "poses": []})
        entry["ligands"].append(
            {
                "ligand_id": ligand.get("ligand_id", ""),
                "score": _fmt(ligand.get("score")),
                "source": ((ligand.get("provenance") or {}).get("source")) or "—",
            }
        )
    for pose in docking:
        feature_id = str(pose.get("feature_id") or "")
        entry = by_protein.setdefault(feature_id, {"feature_id": feature_id, "ligands": [], "poses": []})
        entry["poses"].append(
            {
                "ligand_id": pose.get("ligand_id", ""),
                "score": _fmt(pose.get("score")),
                "engine": pose.get("engine") or "—",
                "control": pose.get("control") or "candidate",
            }
        )
    view.candidates_by_protein = list(by_protein.values())

    kg = report.get("kg") or {}
    if kg:
        edges = kg.get("edges") or []
        kinds: dict[str, int] = {}
        for edge in edges:
            kinds[str(edge.get("type"))] = kinds.get(str(edge.get("type")), 0) + 1
        view.kg_summary = {
            "nodes": len(kg.get("nodes") or []),
            "edges": len(edges),
            "edge_types": dict(sorted(kinds.items(), key=lambda kv: -kv[1])),
        }

    disease = report.get("disease") or {}
    view.disease_claims = [
        {
            "claim": claim.get("claim", ""),
            "evidence": ", ".join(
                list(claim.get("evidence_feature_ids") or []) + list(claim.get("evidence_pmcids") or [])
            ) or "—",
            "source_type": claim.get("source_type") or "—",
        }
        for claim in (disease.get("claims") or [])
    ]

    view.weights = dict(run.get("weights") or {})
    for name, value in (run.get("tool_versions") or {}).items():
        view.provenance.append({"what": str(name), "value": str(value)})
    if run.get("created_at"):
        view.provenance.append({"what": "run created", "value": str(run["created_at"])})
    if genome.get("cga_job_id"):
        view.provenance.append({"what": "CGA job", "value": str(genome["cga_job_id"])})

    # Caveats that depend on what this particular run contains.
    if view.counts["same_species_hits"]:
        view.caveats.append(
            f"{view.counts['same_species_hits']} proteins have a best structural hit from the "
            "same species: this organism has its own structures in the PDB, so those "
            "annotations were recovered, not discovered."
        )
    missing = [s.key for s in view.sections if not s.present]
    if missing:
        view.caveats.append(
            "Not run in this pipeline: " + ", ".join(missing) + ". Those sections are absent, "
            "which is not the same as empty."
        )
    return view

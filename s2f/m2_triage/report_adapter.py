"""Map M2's own output into the `report.json` contract (issues #2, #12).

M2 was written before the schema existed and produces TSV/JSON side files. Rather than rewrite
it, this maps those records into the sections `00-architecture.md` defines: `proteins[].xrefs`,
`annotations[]`, `flags`, `triage`, and `kg`.

M1 owns `proteins[]`. M2 only enriches records that already exist, so a feature_id M1 never wrote
is an error rather than a new row — the canonical key is the one thing every module joins on
(pitfall #11). When M1 has not run yet (M2 driven straight from BV-BRC files, as during
development), `seed_proteins` writes the minimal list first so the enrichment has somewhere to go.
"""

from __future__ import annotations

from datetime import UTC, datetime
from typing import Any, Iterable

from ..common.ids import Xrefs
from .afdb import AlphaFoldModel
from .bvbrc_input import Protein
from .score import ProteinScore


def _now() -> str:
    return datetime.now(UTC).isoformat()


def seed_proteins(proteins: Iterable[Protein]) -> list[dict[str, Any]]:
    """The minimal `proteins[]` M1 would have written, for runs where M1 has not."""
    return [
        {
            "feature_id": protein.feature_id,
            "locus_tag": protein.locus_tag or None,
            "product": protein.product or None,
            "aa_length": len(protein.sequence) or None,
            "specialty": [
                {
                    "property": hit.property_name,
                    "source": hit.source,
                    "identity": hit.identity,
                    "query_coverage": hit.query_coverage,
                }
                for hit in protein.specialty
            ],
        }
        for protein in proteins
    ]


def protein_enrichment(
    score: ProteinScore,
    *,
    xrefs: Xrefs | None = None,
    model: AlphaFoldModel | None = None,
    weights_version: str = "",
) -> dict[str, Any]:
    """The `xrefs` / `annotations` / `flags` / `triage` fields for one protein."""
    top = score.best_hit
    retrieved = _now()

    annotations: list[dict[str, Any]] = []
    if top is not None:
        annotations.append(
            {
                "source": "rcsb_sequence",
                "hit": top.entity_id,
                "description": top.description or None,
                "identity": round(top.identity, 4),
                "coverage": round(top.query_coverage, 4),
                "evalue": top.evalue,
                "retrieved_at": retrieved,
            }
        )
    if xrefs is not None and xrefs.uniprot:
        annotations.append(
            {
                "source": "uniprot",
                "hit": xrefs.uniprot,
                "description": xrefs.protein_name or xrefs.gene_name or None,
                "retrieved_at": retrieved,
                "mapping_route": xrefs.route,
            }
        )
    if model is not None and model.found:
        annotations.append(
            {
                "source": "alphafold_db",
                "hit": model.entry_id,
                "description": f"mean pLDDT {model.mean_plddt} ({model.confidence_band})",
                "coverage": model.coverage,
                "retrieved_at": retrieved,
            }
        )

    usable, reason = model.usable_for_docking() if model is not None else (None, None)

    enrichment: dict[str, Any] = {
        "xrefs": {
            "uniprot": xrefs.uniprot if xrefs else None,
            "uniparc": xrefs.uniparc if xrefs else None,
            "route": xrefs.route if xrefs else None,
            "afdb": model.entry_id if model is not None and model.found else None,
            "pdb": sorted({hit.entry_id for hit in _qualifying(score)}),
            "chembl_target": list(xrefs.chembl) if xrefs else [],
            "refseq": list(xrefs.refseq) if xrefs else [],
            "gene_name": xrefs.gene_name if xrefs else None,
            "related_uniprot": list(xrefs.related_uniprot) if xrefs else [],
        },
        "annotations": annotations,
        "flags": {
            "virulence": bool(score.components.virulence_amr),
            "amr": bool(score.components.virulence_amr),
            "essential_ortholog": bool(score.components.essential),
            "human_homolog": bool(score.flags.get("human_homolog_identity")),
            "human_homolog_identity": score.flags.get("human_homolog_identity") or None,
            # Not available until #10 lands; null means "unknown", not "false".
            "secreted": None,
            "membrane": None,
            "experimental_homolog": (
                {
                    "pdb_entity": top.entity_id,
                    "identity": round(top.identity, 4),
                    "coverage": round(top.query_coverage, 4),
                    "resolution": top.resolution,
                    "method": top.experimental_method or None,
                    "holo": bool(top.ligands),
                    "ligands": list(top.ligands),
                    "organism": top.organism or None,
                }
                if top is not None
                else None
            ),
            "predicted_model": (
                {
                    "source": "afdb",
                    "entry_id": model.entry_id,
                    "mean_plddt": model.mean_plddt,
                    "confidence_band": model.confidence_band,
                    "covered_start": model.uniprot_start,
                    "covered_end": model.uniprot_end,
                    "coverage": model.coverage,
                    "usable_for_docking": usable,
                    "reason": reason,
                    "url": model.cif_url or None,
                }
                if model is not None and model.found
                else None
            ),
        },
        "triage": {
            "score": score.score,
            "components": score.components.as_dict(),
            "weights_version": weights_version or None,
            "rank": score.rank or None,
            "selected": score.selected,
            "reason": score.reason,
            "retrieval_status": score.retrieval_status,
        },
    }
    return enrichment


def _qualifying(score: ProteinScore) -> list[Any]:
    return [score.best_hit] if score.best_hit is not None else []


def kg_section(graph_payload: dict[str, Any], *, path: str | None = None) -> dict[str, Any]:
    """The `kg` section. Edges stay inline while small; past a few thousand, point at kg.json."""
    edges = graph_payload.get("edges") or []
    if len(edges) > 3000:
        return {
            "nodes": [],
            "edges": [],
            "caps": graph_payload.get("caps") or {},
            "path": path,
            "counts": graph_payload.get("counts") or {},
            "note": "edges live in kg.json; too large to inline (02-m2-triage.md)",
        }
    return {
        "nodes": graph_payload.get("nodes") or [],
        "edges": edges,
        "caps": graph_payload.get("caps") or {},
        "path": path,
        "counts": graph_payload.get("counts") or {},
    }

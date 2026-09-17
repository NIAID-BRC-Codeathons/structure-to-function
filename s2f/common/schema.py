"""JSON Schema for each `report.json` section, plus `validate(section, obj)` (issue #2).

The contract in ``00-architecture.md``: one `report.json` per run, each module owning exactly one
top-level key and never editing another's. This module says what those keys must look like.

**Permissive on purpose, for now.** Each section requires only the fields that make a record
*identifiable and attributable* — the identity key, and where a claim came from. Everything else
is optional and unknown properties are allowed, so a module can ship fields before the shape is
settled. The shapes tighten once M3–M6 exist and we can see what they actually produce; that is a
deliberate sequence, not an oversight. Tightening bumps ``SCHEMA_VERSION`` and updates the
fixtures in the same commit (pitfall #18).

What is *not* negotiable, and is enforced from day one:

- ``feature_id`` is the canonical protein key everywhere (``common/ids.py``, pitfall #11).
- Anything retrieved from outside carries ``source`` and ``retrieved_at``.
- Anything transferred from another protein carries the identity that justified it (pitfall #4).
"""

from __future__ import annotations

from typing import Any

from jsonschema import Draft202012Validator

SCHEMA_VERSION = "0.1.0"

#: Which module owns which key. A module writes its own section and never another's.
SECTION_OWNERS = {
    "run": "M1",
    "genome": "M1",
    "proteins": "M1 (list) + M2 (enrichment)",
    "structures": "M3",
    "kg": "M2",
    "ligands": "M4",
    "docking": "M4",
    "analyses": "M3/M4/M5",
    "disease": "M5",
    "report": "M6",
}

_PROVENANCE = {
    "type": "object",
    "description": "Where a record came from. Required on anything retrieved externally.",
    "properties": {
        "source": {"type": "string", "minLength": 1},
        "source_version": {"type": ["string", "null"]},
        "retrieved_at": {"type": "string", "minLength": 1},
        "method": {"type": ["string", "null"]},
        "params": {"type": ["object", "null"]},
    },
    "required": ["source", "retrieved_at"],
    "additionalProperties": True,
}

_ANNOTATION = {
    "type": "object",
    "description": (
        "One evidence row about a protein. Several modules write here — M2 sequence hits, "
        "Foldseek (#9), eggNOG (#10) — so `source` is what keeps them apart."
    ),
    "properties": {
        "source": {"type": "string", "minLength": 1},
        "hit": {"type": ["string", "null"]},
        "description": {"type": ["string", "null"]},
        "identity": {"type": ["number", "null"], "minimum": 0, "maximum": 1},
        "coverage": {"type": ["number", "null"], "minimum": 0, "maximum": 1},
        "tm_score": {"type": ["number", "null"], "minimum": 0, "maximum": 1},
        "evalue": {"type": ["number", "null"]},
        "retrieved_at": {"type": ["string", "null"]},
    },
    "required": ["source"],
    "additionalProperties": True,
}

_XREFS = {
    "type": "object",
    "description": "Identifiers resolved by common/ids.py. Absent mapping is null, never a guess.",
    "properties": {
        "uniprot": {"type": ["string", "null"]},
        "uniparc": {"type": ["string", "null"]},
        "route": {"type": ["string", "null"], "description": "which mapping route produced it"},
        "afdb": {"type": ["string", "null"]},
        "pdb": {"type": "array", "items": {"type": "string"}},
        "chembl_target": {"type": "array", "items": {"type": "string"}},
        "refseq": {"type": "array", "items": {"type": "string"}},
        "gene_name": {"type": ["string", "null"]},
        "related_uniprot": {
            "type": "array",
            "items": {"type": "string"},
            "description": "homolog accessions; never this protein's own identity",
        },
    },
    "additionalProperties": True,
}

_FLAGS = {
    "type": "object",
    "description": (
        "Boolean-ish facts used by triage and by M3/M4. Structure evidence is kept on two "
        "separate axes: an experimental homolog (identity, coverage, ligand-bound or not) and a "
        "predicted model (pLDDT and covered range). One `structure_available` boolean would lose "
        "the distinction M3 gates on."
    ),
    "properties": {
        "virulence": {"type": ["boolean", "null"]},
        "amr": {"type": ["boolean", "null"]},
        "secreted": {"type": ["boolean", "null"]},
        "membrane": {"type": ["boolean", "null"]},
        "essential_ortholog": {"type": ["boolean", "null"]},
        "human_homolog": {"type": ["boolean", "null"]},
        "human_homolog_identity": {
            "type": ["number", "null"],
            "description": "percent identity; a penalty and a stated selectivity risk, never a bonus",
        },
        "experimental_homolog": {
            "type": ["object", "null"],
            "properties": {
                "pdb_entity": {"type": ["string", "null"]},
                "identity": {"type": ["number", "null"], "minimum": 0, "maximum": 1},
                "coverage": {"type": ["number", "null"], "minimum": 0, "maximum": 1},
                "resolution": {"type": ["number", "null"]},
                "method": {"type": ["string", "null"]},
                "holo": {"type": ["boolean", "null"], "description": "entry has a non-additive ligand"},
                "ligands": {"type": "array", "items": {"type": "string"}},
            },
            "additionalProperties": True,
        },
        "predicted_model": {
            "type": ["object", "null"],
            "properties": {
                "source": {"type": ["string", "null"]},
                "entry_id": {"type": ["string", "null"]},
                "mean_plddt": {"type": ["number", "null"]},
                "confidence_band": {"type": ["string", "null"]},
                "covered_start": {"type": ["integer", "null"]},
                "covered_end": {"type": ["integer", "null"]},
                "coverage": {"type": ["number", "null"], "minimum": 0, "maximum": 1},
                "usable_for_docking": {"type": ["boolean", "null"]},
                "reason": {"type": ["string", "null"]},
                "url": {"type": ["string", "null"]},
            },
            "additionalProperties": True,
        },
    },
    "additionalProperties": True,
}

_TRIAGE = {
    "type": "object",
    "description": "M2's ranking. `components` must be enough to recompute `score` on its own.",
    "properties": {
        "score": {"type": "number"},
        "components": {"type": "object", "additionalProperties": {"type": "number"}},
        "weights_version": {"type": ["string", "null"]},
        "rank": {"type": ["integer", "null"], "minimum": 1},
        "selected": {"type": "boolean"},
        "reason": {"type": "string"},
    },
    "required": ["score", "components", "selected", "reason"],
    "additionalProperties": True,
}

_PROTEIN = {
    "type": "object",
    "properties": {
        "feature_id": {"type": "string", "minLength": 1},
        "locus_tag": {"type": ["string", "null"]},
        "product": {"type": ["string", "null"]},
        "contig": {"type": ["string", "null"]},
        "start": {"type": ["integer", "null"]},
        "end": {"type": ["integer", "null"]},
        "strand": {"type": ["string", "null"]},
        "aa_length": {"type": ["integer", "null"]},
        "subsystems": {"type": "array"},
        "specialty": {"type": "array"},
        "xrefs": _XREFS,
        "annotations": {"type": "array", "items": _ANNOTATION},
        "flags": _FLAGS,
        "triage": _TRIAGE,
    },
    "required": ["feature_id"],
    "additionalProperties": True,
}

SECTION_SCHEMAS: dict[str, dict[str, Any]] = {
    "run": {
        "type": "object",
        "properties": {
            "run_id": {"type": "string", "minLength": 1},
            "created_at": {"type": "string", "minLength": 1},
            "schema_version": {"type": "string"},
            "input_file": {"type": ["string", "null"]},
            "tool_versions": {"type": "object"},
            "seeds": {"type": "object"},
            "caps": {"type": "object"},
            "modules": {"type": "object", "description": "per-module status and wall-clock"},
        },
        "required": ["run_id", "created_at"],
        "additionalProperties": True,
    },
    "genome": {
        "type": "object",
        "properties": {
            "taxonomy": {"type": ["object", "string", "null"]},
            "taxon_id": {"type": ["integer", "null"]},
            "genome_id": {"type": ["string", "null"], "description": "BV-BRC genome id, e.g. 243273.147"},
            "closest_genomes": {"type": "array"},
            "tree_newick": {"type": ["string", "null"]},
            "cga_job_id": {"type": ["string", "null"]},
        },
        "additionalProperties": True,
    },
    "proteins": {"type": "array", "items": _PROTEIN},
    "structures": {
        "type": "array",
        "items": {
            "type": "object",
            "properties": {
                "feature_id": {"type": "string", "minLength": 1},
                "source": {"type": "string", "enum": ["pdb", "afdb", "predicted", "esm_atlas"]},
                "accession": {"type": ["string", "null"]},
                "path": {"type": ["string", "null"]},
                "confidence": {"type": ["number", "string", "null"]},
                "holo_template": {"type": ["object", "null"]},
                "pockets": {"type": "array"},
                "usable_for_docking": {"type": ["boolean", "null"]},
                "reason": {"type": ["string", "null"]},
            },
            "required": ["feature_id", "source"],
            "additionalProperties": True,
        },
    },
    "kg": {
        "type": "object",
        "properties": {
            "nodes": {
                "type": "array",
                "items": {
                    "type": "object",
                    "properties": {
                        "id": {"type": "string", "minLength": 1},
                        "type": {"type": "string"},
                        "name": {"type": ["string", "null"]},
                        "source": {"type": ["string", "null"]},
                    },
                    "required": ["id", "type"],
                    "additionalProperties": True,
                },
            },
            "edges": {
                "type": "array",
                "items": {
                    "type": "object",
                    "properties": {
                        "source_id": {"type": "string", "minLength": 1},
                        "target_id": {"type": "string", "minLength": 1},
                        "type": {"type": "string", "minLength": 1},
                        "evidence": {"type": "object"},
                        "provenance": _PROVENANCE,
                    },
                    "required": ["source_id", "target_id", "type", "provenance"],
                    "additionalProperties": True,
                },
            },
            "caps": {"type": "object"},
            "path": {"type": ["string", "null"], "description": "set when edges live in kg.json instead"},
        },
        "additionalProperties": True,
    },
    "ligands": {
        "type": "array",
        "items": {
            "type": "object",
            "properties": {
                "ligand_id": {"type": "string", "minLength": 1},
                "feature_id": {"type": ["string", "null"], "description": "the protein it was proposed for"},
                "smiles": {"type": ["string", "null"]},
                "inchikey": {"type": ["string", "null"]},
                "score": {"type": ["number", "null"]},
                "components": {"type": "object"},
                "provenance": _PROVENANCE,
            },
            "required": ["ligand_id", "provenance"],
            "additionalProperties": True,
        },
    },
    "docking": {
        "type": "array",
        "items": {
            "type": "object",
            "properties": {
                "feature_id": {"type": "string", "minLength": 1},
                "ligand_id": {"type": "string", "minLength": 1},
                "score": {"type": ["number", "null"]},
                "engine": {"type": ["string", "null"]},
                "seed": {"type": ["integer", "null"]},
                "exhaustiveness": {"type": ["integer", "null"]},
                "pose_path": {"type": ["string", "null"]},
                "control": {
                    "type": ["string", "null"],
                    "description": "positive | decoy | null; rankings are within one protein only",
                },
            },
            "required": ["feature_id", "ligand_id"],
            "additionalProperties": True,
        },
    },
    "analyses": {
        "type": "array",
        "items": {
            "type": "object",
            "properties": {
                "feature_id": {"type": ["string", "null"]},
                "method": {"type": "string", "minLength": 1},
                "parameters": {"type": "object"},
                "result": {},
                "confidence": {"type": ["number", "string", "null"]},
            },
            "required": ["method"],
            "additionalProperties": True,
        },
    },
    "disease": {
        "type": "object",
        "properties": {
            "claims": {
                "type": "array",
                "items": {
                    "type": "object",
                    "properties": {
                        "claim": {"type": "string", "minLength": 1},
                        "evidence_feature_ids": {"type": "array", "items": {"type": "string"}},
                        "evidence_pmcids": {"type": "array", "items": {"type": "string"}},
                        "source_type": {
                            "type": ["string", "null"],
                            "description": "strain record vs genomic inference (pitfall #8)",
                        },
                        "verified": {"type": ["boolean", "null"]},
                    },
                    "required": ["claim"],
                    "additionalProperties": True,
                },
            },
            "strategies": {"type": "array"},
            "dropped_claims": {"type": ["integer", "null"]},
        },
        "additionalProperties": True,
    },
    "report": {
        "type": "object",
        "properties": {
            "rendered_path": {"type": ["string", "null"]},
            "rendered_at": {"type": ["string", "null"]},
            "counts": {"type": "object"},
        },
        "additionalProperties": True,
    },
}


class SchemaError(ValueError):
    """Raised when a section does not match its schema."""


def section_names() -> list[str]:
    return sorted(SECTION_SCHEMAS)


def validate(section: str, obj: Any) -> None:
    """Validate one section. Raises SchemaError listing every problem, not just the first."""
    schema = SECTION_SCHEMAS.get(section)
    if schema is None:
        raise SchemaError(
            f"unknown section {section!r}; known sections: {', '.join(section_names())}"
        )
    validator = Draft202012Validator(schema)
    errors = sorted(validator.iter_errors(obj), key=lambda e: list(e.absolute_path))
    if not errors:
        return
    problems = []
    for error in errors[:20]:
        where = "/".join(str(part) for part in error.absolute_path) or "(root)"
        problems.append(f"  {section}.{where}: {error.message}")
    more = "" if len(errors) <= 20 else f"\n  ... and {len(errors) - 20} more"
    raise SchemaError(f"{len(errors)} problem(s) in section {section!r}:\n" + "\n".join(problems) + more)


def is_valid(section: str, obj: Any) -> bool:
    try:
        validate(section, obj)
    except SchemaError:
        return False
    return True


def validate_report(report: dict[str, Any]) -> None:
    """Validate every known section present in a whole report. Unknown keys are reported."""
    problems: list[str] = []
    for key, value in report.items():
        if key == "schema_version":
            continue
        if key not in SECTION_SCHEMAS:
            problems.append(f"  unknown top-level key {key!r} (owners: {', '.join(SECTION_OWNERS)})")
            continue
        try:
            validate(key, value)
        except SchemaError as exc:
            problems.append(str(exc))
    if problems:
        raise SchemaError("report did not validate:\n" + "\n".join(problems))

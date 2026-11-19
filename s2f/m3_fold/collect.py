"""Collect and quality-gate existing PDB and AlphaFold DB structures selected by M2."""

from __future__ import annotations

import hashlib
import os
import re
import tempfile
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable
from urllib.parse import urlparse

from ..common.http import DownloadedFile
from .pockets import find_pockets
from .receptor import available_tool, docking_box, prepare_receptor, vina_config
from .quality import QualityThresholds, assess_structure

PDB_DOWNLOAD_URL = "https://files.rcsb.org/download/{entry_id}.cif"


def _utc_now() -> datetime:
    return datetime.now(timezone.utc)


@dataclass(frozen=True)
class StructureCandidate:
    """One existing structure chosen from M2's evidence."""

    feature_id: str
    source: str
    accession: str
    url: str
    confidence: float | None
    confidence_metric: str
    experimental_homolog: dict[str, Any] | None = None
    predicted_model: dict[str, Any] | None = None
    #: M2's best ligand- or metal-bearing qualifying homolog, which M3 reads a binding site
    #: from when it is the entry collected here (#14).
    holo_template: dict[str, Any] | None = None
    #: Set when the holo template was preferred over the best-identity hit, with what the
    #: swap cost in sequence identity.
    template_preference: dict[str, Any] | None = None

    @property
    def suffix(self) -> str:
        suffix = Path(urlparse(self.url).path).suffix.lower()
        return suffix if suffix in {".cif", ".pdb"} else ".cif"


@dataclass(frozen=True)
class CollectionSummary:
    """Structures written for one collection run."""

    structures: list[dict[str, Any]]
    selected: int
    collected: int
    prediction_required: int
    failed: int


def _pdb_entry_id(pdb_entity: str) -> str:
    """Return the four-character entry ID from values such as ``1YDX_1``."""
    return pdb_entity.split("_", 1)[0].upper()


def choose_existing(protein: dict[str, Any]) -> StructureCandidate | None:
    """Choose M2's preferred existing structure: experimental first, then AFDB.

    M2 currently records one best experimental homolog. Its ``holo`` field is preserved for
    issue 14, but there is no second experimental candidate here to compare it with.
    """
    feature_id = str(protein.get("feature_id") or "")
    if not feature_id:
        raise ValueError("selected protein is missing feature_id")

    flags = protein.get("flags") or {}
    experimental = flags.get("experimental_homolog")
    holo = flags.get("holo_template")
    holo = dict(holo) if isinstance(holo, dict) and holo.get("pdb_entity") else None

    if isinstance(experimental, dict) and experimental.get("pdb_entity"):
        chosen, preference = dict(experimental), None
        # A pocket can only be read from an entry that binds something, so a holo homolog is
        # worth more here than the last few points of sequence identity (pitfalls.md #2).
        # The swap is recorded, never silent: it trades sequence match for site evidence.
        if holo is not None and not holo.get("is_best_hit"):
            chosen = holo
            preference = {
                "policy": "prefer_holo_over_best_identity",
                "best_identity_entity": experimental.get("pdb_entity"),
                "best_identity": experimental.get("identity"),
                "chosen_entity": holo.get("pdb_entity"),
                "chosen_identity": holo.get("identity"),
                "identity_drop": holo.get("identity_drop"),
                "binds": [*(holo.get("ligands") or []), *(holo.get("metals") or [])],
            }
        entity = str(chosen["pdb_entity"])
        return StructureCandidate(
            feature_id=feature_id,
            source="pdb",
            accession=entity,
            url=PDB_DOWNLOAD_URL.format(entry_id=_pdb_entry_id(entity)),
            confidence=_number_or_none(chosen.get("resolution")),
            confidence_metric="resolution_angstrom",
            experimental_homolog=chosen,
            holo_template=holo,
            template_preference=preference,
        )

    predicted = flags.get("predicted_model")
    if isinstance(predicted, dict) and predicted.get("entry_id"):
        return StructureCandidate(
            feature_id=feature_id,
            source="afdb",
            accession=str(predicted["entry_id"]),
            url=str(predicted.get("url") or ""),
            confidence=_number_or_none(predicted.get("mean_plddt")),
            confidence_metric="mean_plddt",
            predicted_model=dict(predicted),
            holo_template=holo,
        )

    # The richer flags above are preferred, but tolerate a report containing only xrefs.
    xrefs = protein.get("xrefs") or {}
    pdb_entries = xrefs.get("pdb") or []
    if pdb_entries:
        entry_id = _pdb_entry_id(str(pdb_entries[0]))
        return StructureCandidate(
            feature_id=feature_id,
            source="pdb",
            accession=entry_id,
            url=PDB_DOWNLOAD_URL.format(entry_id=entry_id),
            confidence=None,
            confidence_metric="resolution_angstrom",
        )

    if xrefs.get("afdb"):
        return StructureCandidate(
            feature_id=feature_id,
            source="afdb",
            accession=str(xrefs["afdb"]),
            url="",
            confidence=None,
            confidence_metric="mean_plddt",
        )
    return None


def collect_existing(
    report: dict[str, Any],
    run_dir: str | Path,
    fetch: Callable[[StructureCandidate], bytes | DownloadedFile],
    *,
    limit: int = 0,
    quality_thresholds: QualityThresholds | None = None,
    receptors: bool = True,
    pdbqt_fallback: bool = False,
    find_tool: Callable[[], str | None] = available_tool,
    clock: Callable[[], float] = time.monotonic,
    now: Callable[[], datetime] = _utc_now,
) -> CollectionSummary:
    """Collect existing structures for selected proteins and return schema-ready records."""
    quality_thresholds = quality_thresholds or QualityThresholds()
    selected = [
        protein
        for protein in report.get("proteins", [])
        if (protein.get("triage") or {}).get("selected") is True
    ]
    if limit:
        selected = selected[:limit]

    run_dir = Path(run_dir)
    structures_dir = run_dir / "structures"
    records: list[dict[str, Any]] = []

    for protein in selected:
        started = clock()
        candidate = choose_existing(protein)
        if candidate is None:
            records.append(_prediction_required_record(protein, clock() - started))
            continue

        record = _candidate_record(candidate)
        try:
            if not candidate.url:
                raise ValueError(f"M2 did not provide a download URL for {candidate.accession}")
            download = _as_downloaded_file(fetch(candidate), now)
            payload = download.content
            if not payload:
                raise ValueError(f"download for {candidate.accession} was empty")
            destination = structures_dir / f"{_safe_name(candidate.feature_id)}{candidate.suffix}"
            _atomic_write(destination, payload)
        except (OSError, RuntimeError, ValueError) as exc:
            record.update(
                collection_status="failed",
                usable_for_docking=False,
                reason=f"existing structure could not be collected: {exc}",
            )
        else:
            evidence = (
                candidate.experimental_homolog
                if candidate.source == "pdb"
                else candidate.predicted_model
            )
            quality = assess_structure(candidate.source, evidence, payload, quality_thresholds)
            record["params"]["confidence_gate_applied"] = True
            pockets, pocket_summary = find_pockets(
                payload,
                source=candidate.source,
                holo_template=candidate.holo_template,
                structure_accession=candidate.accession,
                min_plddt=quality_thresholds.min_local_plddt,
            )
            record["pockets"] = pockets
            record["pocket_search"] = pocket_summary
            if receptors and quality.usable_for_docking:
                _prepare_for_docking(
                    record,
                    payload,
                    candidate=candidate,
                    run_dir=run_dir,
                    pockets=pockets,
                    allow_fallback=pdbqt_fallback,
                    find_tool=find_tool,
                )
            record.update(
                path=str(destination.relative_to(run_dir)),
                sha256=hashlib.sha256(payload).hexdigest(),
                retrieved_at=download.retrieved_at.isoformat(),
                collected_at=now().isoformat(),
                cache_hit=download.from_cache,
                collection_status="collected",
                usable_for_docking=quality.usable_for_docking,
                reason=quality.reason,
                quality_gate=quality.details,
                **_source_metadata(candidate, payload),
            )
        record["elapsed_seconds"] = round(clock() - started, 6)
        records.append(record)

    statuses = [record["collection_status"] for record in records]
    return CollectionSummary(
        structures=records,
        selected=len(selected),
        collected=statuses.count("collected"),
        prediction_required=statuses.count("prediction_required"),
        failed=statuses.count("failed"),
    )


def _candidate_record(candidate: StructureCandidate) -> dict[str, Any]:
    record: dict[str, Any] = {
        "feature_id": candidate.feature_id,
        "source": candidate.source,
        "accession": candidate.accession,
        "path": None,
        "confidence": candidate.confidence,
        "confidence_metric": candidate.confidence_metric,
        "download_url": candidate.url or None,
        "source_version": None,
        "release_date": None,
        "method": "download_existing_structure",
        "params": {
            "selection_policy": "experimental_homolog_then_afdb",
            "format": candidate.suffix.lstrip("."),
            "confidence_gate_applied": False,
        },
        "holo_template": candidate.holo_template,
        "pockets": [],
    }
    if candidate.template_preference is not None:
        record["template_preference"] = candidate.template_preference
    if candidate.experimental_homolog is not None:
        record["experimental_homolog"] = candidate.experimental_homolog
    if candidate.predicted_model is not None:
        record["predicted_model"] = candidate.predicted_model
    return record


def _prediction_required_record(protein: dict[str, Any], elapsed: float) -> dict[str, Any]:
    """Represent a structural gap until issue 13's prediction scope fills it."""
    return {
        "feature_id": str(protein.get("feature_id") or ""),
        "source": "predicted",
        "accession": None,
        "path": None,
        "confidence": None,
        "confidence_metric": None,
        "source_version": None,
        "release_date": None,
        "method": "lookup_existing_structure",
        "params": {
            "selection_policy": "experimental_homolog_then_afdb",
            "confidence_gate_applied": False,
        },
        "holo_template": None,
        "pockets": [],
        "collection_status": "prediction_required",
        "usable_for_docking": None,
        "reason": "no existing PDB or AlphaFold DB structure; prediction required",
        "elapsed_seconds": round(elapsed, 6),
    }


def _safe_name(feature_id: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", feature_id).strip("_") or "protein"


def _number_or_none(value: Any) -> float | None:
    try:
        return float(value) if value is not None else None
    except (TypeError, ValueError):
        return None


def _as_downloaded_file(
    value: bytes | DownloadedFile, now: Callable[[], datetime]
) -> DownloadedFile:
    if isinstance(value, DownloadedFile):
        return value
    return DownloadedFile(content=value, retrieved_at=now(), from_cache=False)


def _source_metadata(candidate: StructureCandidate, payload: bytes) -> dict[str, str | None]:
    if candidate.source == "pdb":
        release_date, revision_date = _pdb_revision_dates(payload)
        return {
            "source_version": f"revision:{revision_date}" if revision_date else None,
            "release_date": release_date,
        }

    predicted = candidate.predicted_model or {}
    version = predicted.get("version")
    if not version:
        match = re.search(r"model_v(\d+)", candidate.url)
        version = f"v{match.group(1)}" if match else None
    return {
        "source_version": str(version) if version else None,
        "release_date": str(predicted.get("created") or "") or None,
    }


def _pdb_revision_dates(payload: bytes) -> tuple[str | None, str | None]:
    """Read initial and latest revision dates from an RCSB mmCIF header."""
    text = payload.decode("utf-8", errors="replace")
    marker = "_pdbx_audit_revision_history.revision_date"
    if marker not in text:
        return None, None
    revision_block = text.split(marker, 1)[1].split("#", 1)[0]
    dates = re.findall(r"\b\d{4}-\d{2}-\d{2}\b", revision_block)
    return (min(dates), max(dates)) if dates else (None, None)


def _atomic_write(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile("wb", dir=path.parent, prefix=".structure-", delete=False)
    try:
        with handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(handle.name, path)
    except BaseException:
        Path(handle.name).unlink(missing_ok=True)
        raise


def _prepare_for_docking(
    record: dict[str, Any],
    payload: bytes,
    *,
    candidate: StructureCandidate,
    run_dir: Path,
    pockets: list[dict[str, Any]],
    allow_fallback: bool,
    find_tool: Callable[[], str | None],
) -> None:
    """Write the receptor and a search box per pocket (#14 part B).

    Only for a structure that passed the quality gate: preparing a receptor M3 has already
    called unusable would hand M4 something to dock that M3 said not to.
    """
    stem = run_dir / "receptors" / _safe_name(candidate.feature_id)
    entity = candidate.accession.split("_", 1)[1] if "_" in candidate.accession else ""
    prepared = prepare_receptor(
        payload,
        entity_id=entity if candidate.source == "pdb" else "",
        destination_stem=stem,
        run_dir=run_dir,
        allow_fallback=allow_fallback,
        find_tool=find_tool,
    )
    record["receptor"] = prepared.as_dict()

    docking_target = prepared.pdbqt_path or prepared.pdb_path
    for index, pocket in enumerate(pockets, 1):
        box = docking_box(payload, pocket)
        if box is None:
            continue
        pocket["box"] = box
        if not pocket.get("usable") or not docking_target:
            continue
        config = stem.parent / f"{stem.name}_pocket_{index:02d}.vina.txt"
        config.write_text(
            vina_config(box, Path(docking_target).name), encoding="utf-8"
        )
        pocket["vina_config"] = str(config.relative_to(run_dir))

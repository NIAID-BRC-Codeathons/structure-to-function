"""Provisional, auditable quality gates for structures collected by M3.

The final biological thresholds are still a team decision.  This module therefore keeps the
policy in one small dataclass and records both the thresholds and observed values in every
structure record.  The defaults reuse gates that M2 already applies instead of introducing a
second hidden policy.
"""

from __future__ import annotations

import shlex
from dataclasses import asdict, dataclass
from typing import Any


POLICY_VERSION = "provisional-2026-09-17"


@dataclass(frozen=True)
class QualityThresholds:
    """Thresholds used by the provisional M3 quality gate."""

    min_pdb_identity: float = 0.25
    min_pdb_coverage: float = 0.50
    max_pdb_resolution: float = 3.50
    min_mean_plddt: float = 70.0
    min_afdb_coverage: float = 0.80
    min_local_plddt: float = 70.0

    def as_dict(self) -> dict[str, float]:
        return asdict(self)


@dataclass(frozen=True)
class QualityDecision:
    """One pass/fail decision and the evidence that produced it."""

    usable_for_docking: bool
    reason: str
    details: dict[str, Any]


def assess_structure(
    source: str,
    evidence: dict[str, Any] | None,
    coordinate_file: bytes,
    thresholds: QualityThresholds,
) -> QualityDecision:
    """Apply the source-appropriate quality gate to one collected structure."""
    if source == "pdb":
        return _assess_pdb(evidence or {}, thresholds)
    if source == "afdb":
        return _assess_afdb(evidence or {}, coordinate_file, thresholds)
    return QualityDecision(
        False,
        f"no quality policy is defined for structure source {source!r}",
        _details(source, thresholds, observed={}, checks=[], masked_regions=[]),
    )


def _assess_pdb(evidence: dict[str, Any], thresholds: QualityThresholds) -> QualityDecision:
    identity = _number_or_none(evidence.get("identity"))
    coverage = _number_or_none(evidence.get("coverage"))
    resolution = _number_or_none(evidence.get("resolution"))
    method = str(evidence.get("method") or "").strip()

    checks = [
        _minimum_check("sequence_identity", identity, thresholds.min_pdb_identity),
        _minimum_check("sequence_coverage", coverage, thresholds.min_pdb_coverage),
    ]

    # NMR ensembles do not have an X-ray/cryo-EM-style resolution value.
    if "NMR" in method.upper():
        checks.append(
            {
                "metric": "resolution_angstrom",
                "value": resolution,
                "operator": "not_applicable",
                "threshold": thresholds.max_pdb_resolution,
                "passed": True,
                "reason": "resolution cutoff is not applicable to NMR",
            }
        )
    else:
        checks.append(_maximum_check("resolution_angstrom", resolution, thresholds.max_pdb_resolution))

    failed = [check for check in checks if not check["passed"]]
    observed = {
        "identity": identity,
        "coverage": coverage,
        "resolution_angstrom": resolution,
        "experimental_method": method or None,
    }
    details = _details("pdb", thresholds, observed, checks, masked_regions=[])
    details["pae"] = {"status": "not_applicable", "reason": "experimental PDB structure"}

    if failed:
        return QualityDecision(
            False,
            "experimental structure failed provisional quality gate: "
            + "; ".join(check["reason"] for check in failed),
            details,
        )

    resolution_text = (
        f", resolution {resolution:.2f} A" if resolution is not None else f", method {method or 'unknown'}"
    )
    return QualityDecision(
        True,
        f"experimental homolog passes provisional quality gate: identity {identity:.1%}, "
        f"coverage {coverage:.1%}{resolution_text}",
        details,
    )


def _assess_afdb(
    evidence: dict[str, Any], coordinate_file: bytes, thresholds: QualityThresholds
) -> QualityDecision:
    mean_plddt = _number_or_none(evidence.get("mean_plddt"))
    coverage = _number_or_none(evidence.get("coverage"))
    residue_confidence = _mmcif_residue_plddt(coordinate_file)
    masked_regions = _low_confidence_regions(residue_confidence, thresholds.min_local_plddt)
    confident_residues = sum(
        value >= thresholds.min_local_plddt for _, _, value in residue_confidence
    )

    checks = [
        _minimum_check("mean_plddt", mean_plddt, thresholds.min_mean_plddt),
        _minimum_check("sequence_coverage", coverage, thresholds.min_afdb_coverage),
        {
            "metric": "per_residue_plddt",
            "value": len(residue_confidence),
            "operator": ">",
            "threshold": 0,
            "passed": bool(residue_confidence),
            "reason": (
                f"per-residue pLDDT found for {len(residue_confidence)} residues"
                if residue_confidence
                else "coordinate file contains no readable per-residue pLDDT"
            ),
        },
        {
            "metric": "confident_residues",
            "value": confident_residues,
            "operator": ">",
            "threshold": 0,
            "passed": confident_residues > 0,
            "reason": (
                f"{confident_residues} residues have pLDDT >= {thresholds.min_local_plddt:g}"
                if confident_residues
                else f"no residues have pLDDT >= {thresholds.min_local_plddt:g}"
            ),
        },
    ]
    observed = {
        "mean_plddt": mean_plddt,
        "coverage": coverage,
        "assessed_residues": len(residue_confidence),
        "confident_residues": confident_residues,
        "masked_residues": len(residue_confidence) - confident_residues,
    }
    details = _details("afdb", thresholds, observed, checks, masked_regions)
    details["mask_policy"] = "exclude listed residues from pocket search; coordinates are unchanged"
    details["pae"] = {
        "status": "not_available",
        "reason": "M2 does not currently pass the AlphaFold PAE document to M3",
    }

    failed = [check for check in checks if not check["passed"]]
    if failed:
        return QualityDecision(
            False,
            "AlphaFold DB model failed provisional quality gate: "
            + "; ".join(check["reason"] for check in failed),
            details,
        )

    masked_count = observed["masked_residues"]
    mask_text = (
        f"; {masked_count} low-confidence residues masked from pocket search"
        if masked_count
        else "; no low-confidence residues require masking"
    )
    return QualityDecision(
        True,
        f"AlphaFold DB model passes provisional pLDDT/coverage gate: mean pLDDT "
        f"{mean_plddt:.1f}, coverage {coverage:.1%}{mask_text}; PAE not yet available",
        details,
    )


def _details(
    source: str,
    thresholds: QualityThresholds,
    observed: dict[str, Any],
    checks: list[dict[str, Any]],
    masked_regions: list[dict[str, Any]],
) -> dict[str, Any]:
    return {
        "policy_version": POLICY_VERSION,
        "provisional": True,
        "source": source,
        "thresholds": thresholds.as_dict(),
        "observed": observed,
        "checks": checks,
        "masked_regions": masked_regions,
    }


def _minimum_check(metric: str, value: float | None, threshold: float) -> dict[str, Any]:
    passed = value is not None and value >= threshold
    if value is None:
        reason = f"{metric} is missing"
    else:
        reason = f"{metric} {value:g} {'meets' if passed else 'is below'} minimum {threshold:g}"
    return {
        "metric": metric,
        "value": value,
        "operator": ">=",
        "threshold": threshold,
        "passed": passed,
        "reason": reason,
    }


def _maximum_check(metric: str, value: float | None, threshold: float) -> dict[str, Any]:
    passed = value is not None and value <= threshold
    if value is None:
        reason = f"{metric} is missing"
    else:
        reason = f"{metric} {value:g} {'meets' if passed else 'exceeds'} maximum {threshold:g}"
    return {
        "metric": metric,
        "value": value,
        "operator": "<=",
        "threshold": threshold,
        "passed": passed,
        "reason": reason,
    }


def _mmcif_residue_plddt(payload: bytes) -> list[tuple[str, int, float]]:
    """Read one pLDDT value per residue from an AlphaFold mmCIF atom-site loop.

    AlphaFold stores pLDDT in ``_atom_site.B_iso_or_equiv``.  CA atoms are used so each residue
    appears once.  The small parser is intentionally limited to this standard atom-site loop;
    malformed or non-AlphaFold files return no values and fail the auditable local-confidence
    check instead of silently passing.
    """
    lines = payload.decode("utf-8", errors="replace").splitlines()
    index = 0
    while index < len(lines):
        if lines[index].strip() != "loop_":
            index += 1
            continue

        index += 1
        headers: list[str] = []
        while index < len(lines) and lines[index].lstrip().startswith("_"):
            headers.append(lines[index].split()[0])
            index += 1
        if "_atom_site.B_iso_or_equiv" not in headers:
            continue

        data_tokens: list[str] = []
        while index < len(lines):
            line = lines[index].strip()
            if line == "#" or line == "loop_" or line.startswith("_"):
                break
            if line:
                try:
                    data_tokens.extend(shlex.split(line, comments=False, posix=True))
                except ValueError:
                    return []
            index += 1
        return _atom_site_plddt_rows(headers, data_tokens)
    return []


def _atom_site_plddt_rows(
    headers: list[str], tokens: list[str]
) -> list[tuple[str, int, float]]:
    width = len(headers)
    if not width or len(tokens) % width:
        return []

    positions = {name: index for index, name in enumerate(headers)}
    atom_column = positions.get("_atom_site.label_atom_id")
    sequence_column = positions.get("_atom_site.label_seq_id")
    chain_column = positions.get("_atom_site.label_asym_id")
    confidence_column = positions.get("_atom_site.B_iso_or_equiv")
    model_column = positions.get("_atom_site.pdbx_PDB_model_num")
    if None in {atom_column, sequence_column, chain_column, confidence_column}:
        return []

    residues: list[tuple[str, int, float]] = []
    for start in range(0, len(tokens), width):
        row = tokens[start : start + width]
        if model_column is not None and row[model_column] not in {"1", ".", "?"}:
            continue
        if row[atom_column].upper() != "CA":
            continue
        try:
            sequence_id = int(row[sequence_column])
            confidence = float(row[confidence_column])
        except (TypeError, ValueError):
            continue
        residues.append((row[chain_column], sequence_id, confidence))
    return residues


def _low_confidence_regions(
    residues: list[tuple[str, int, float]], threshold: float
) -> list[dict[str, Any]]:
    low = [(chain, sequence_id, value) for chain, sequence_id, value in residues if value < threshold]
    if not low:
        return []

    regions: list[dict[str, Any]] = []
    current: list[tuple[str, int, float]] = [low[0]]
    for residue in low[1:]:
        previous = current[-1]
        if residue[0] == previous[0] and residue[1] == previous[1] + 1:
            current.append(residue)
        else:
            regions.append(_masked_region(current, threshold))
            current = [residue]
    regions.append(_masked_region(current, threshold))
    return regions


def _masked_region(region: list[tuple[str, int, float]], threshold: float) -> dict[str, Any]:
    return {
        "chain": region[0][0],
        "start": region[0][1],
        "end": region[-1][1],
        "residue_count": len(region),
        "mean_plddt": round(sum(item[2] for item in region) / len(region), 2),
        "threshold": threshold,
        "action": "exclude_from_pocket_search",
    }


def _number_or_none(value: Any) -> float | None:
    try:
        return float(value) if value is not None else None
    except (TypeError, ValueError):
        return None

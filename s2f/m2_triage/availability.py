"""Which score components had no possible source on this run (issue #90).

A component with no evidence contributes 0, exactly like "checked and false" (pitfall #12). On a
genome route the BV-BRC specialty table is the source for four components; on the FASTA route
(#89) there is no BV-BRC genome, so nobody looked — and the manifest has to say so rather than
let a sequence-only ranking pass for a complete one.
"""

from __future__ import annotations

import json
from pathlib import Path

#: Components whose only default source is BV-BRC's specialty table, and the M2 provider that
#: can replace it without a genome, if any. #91 adds VFDB/CARD for the first two.
SPECIALTY_COMPONENTS: dict[str, str | None] = {
    "virulence_amr": None,
    "drug_target": None,
    "essential": "essentiality",
    "human_homolog_penalty": "human_homology",
}

NO_GENOME_ROUTES = {"fasta"}


def annotation_route(run_dir: str | Path) -> str:
    """M1's route from report.json, or "" when M1 did not run here (standalone M2)."""
    try:
        report = json.loads((Path(run_dir) / "report.json").read_text())
    except (OSError, json.JSONDecodeError):
        return ""
    return str((report.get("genome") or {}).get("annotation_route") or "")


def not_available(route: str, providers_enabled: dict[str, bool]) -> dict[str, str]:
    """`{component: reason}` for every component nothing could have measured on this route."""
    if route not in NO_GENOME_ROUTES:
        return {}
    out = {}
    for component, provider in SPECIALTY_COMPONENTS.items():
        if provider and providers_enabled.get(provider):
            continue
        hint = f"; pass --{provider.replace('_', '-')}" if provider else ""
        out[component] = f"no BV-BRC genome on the {route} route{hint}"
    return out

"""Scoring profiles: the triage policy is a choice, not a constant (issue #93).

``WEIGHTS`` in ``score.py`` encodes an antibacterial target-selection policy — human homology is
a penalty, essentiality a bonus, membrane proteins deprioritised. That is right for bacterial
drug targets and wrong for other questions (a viral protein, a single hypothetical, a human
target). The weights are not portable between genomes either: G37 and HS11286 put 32 and 10
uncharacterised proteins in the top 50 under identical weights (#61).

So the weights become a named profile, or a JSON file, and ``run.json`` records which one
produced a ranking so two runs cannot be confused. ``bacterial-target`` stays the default and is
exactly ``WEIGHTS``; changing it still needs a dated row in the change log in
``docs/02a-m2-pdb-evidence.md`` (pitfall #12). Other profiles are provisional and uncalibrated.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .score import WEIGHTS, WEIGHTS_VERSION

DEFAULT_PROFILE = "bacterial-target"


@dataclass(frozen=True)
class Profile:
    name: str
    weights: dict[str, float]
    description: str
    #: Whether a close human homolog removes a protein from selection. An antibacterial policy;
    #: meaningless when the question is not "is this a safe drug target in a pathogen".
    disqualify_human_homologs: bool = True
    version: str = ""
    source: str = "builtin"
    extra: dict[str, Any] = field(default_factory=dict)

    def as_dict(self) -> dict[str, Any]:
        return {"name": self.name, "version": self.version, "source": self.source,
                "description": self.description, "weights": dict(self.weights),
                "disqualify_human_homologs": self.disqualify_human_homologs}


def _zeroed(**weights: float) -> dict[str, float]:
    return {name: float(weights.get(name, 0.0)) for name in WEIGHTS}


PROFILES: dict[str, Profile] = {
    DEFAULT_PROFILE: Profile(
        name=DEFAULT_PROFILE, weights=dict(WEIGHTS), version=WEIGHTS_VERSION,
        description="antibacterial target selection: structural evidence, virulence/AMR, "
                    "essentiality and drug-target bonuses, membrane and human-homolog "
                    "penalties (docs/02a-m2-pdb-evidence.md)"),
    "evidence-only": Profile(
        name="evidence-only",
        weights=_zeroed(pdb_evidence=0.40, ligandable_homolog=0.10),
        version="provisional-2026-10-08", disqualify_human_homologs=False,
        description="ranks by structural evidence alone (PDB hit quality, ligand-bound "
                    "homolog); no target-selection policy, no penalties. Provisional."),
}


def load_profile(name: str = DEFAULT_PROFILE, weights_file: str | None = None) -> Profile:
    """A built-in profile, or one read from a JSON file.

    The file is ``{"name": ..., "weights": {...}, "disqualify_human_homologs": bool}``. Every
    component must be named, so a typo cannot silently drop one to 0.
    """
    if not weights_file:
        if name not in PROFILES:
            raise ValueError(f"unknown profile {name!r}; choose from {', '.join(PROFILES)}")
        return PROFILES[name]
    payload = json.loads(Path(weights_file).read_text(encoding="utf-8"))
    weights = payload.get("weights") if isinstance(payload, dict) else None
    if not isinstance(weights, dict):
        raise ValueError(f"{weights_file}: expected an object with a 'weights' object")
    missing = set(WEIGHTS) - set(weights)
    unknown = set(weights) - set(WEIGHTS)
    if missing or unknown:
        raise ValueError(f"{weights_file}: weights must name exactly the components "
                         f"{sorted(WEIGHTS)}; missing {sorted(missing)}, unknown {sorted(unknown)}")
    try:
        parsed = {key: float(weights[key]) for key in WEIGHTS}
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{weights_file}: every weight must be a number ({exc})") from None
    return Profile(
        name=str(payload.get("name") or Path(weights_file).stem), weights=parsed,
        description=str(payload.get("description") or "custom weights"),
        disqualify_human_homologs=bool(payload.get("disqualify_human_homologs", True)),
        version=str(payload.get("version") or ""), source=str(weights_file))

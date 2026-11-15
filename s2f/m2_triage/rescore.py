"""Re-rank a finished M2 run under another scoring profile, with no re-querying (issue #93).

    python -m s2f.m2_triage.rescore --run runs/kp --profile evidence-only
    python -m s2f.m2_triage.rescore --run runs/kp --weights my_weights.json --top 30

Every score component is written per protein in ``m2_pdb/proteins.tsv``, so another policy can
be applied to recorded evidence directly. The result goes to
``m2_pdb/rescored/<profile>.tsv`` with a small manifest; the run's own ranking, ``report.json``
and ``run.json`` are left untouched, because a rescore is a comparison, not a new run.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

from .profiles import DEFAULT_PROFILE, PROFILES, Profile, load_profile
from .score import DISQUALIFYING_HUMAN_HOMOLOG_IDENTITY, WEIGHTS


def _number(value: str | None) -> float:
    try:
        return float(value) if value not in (None, "") else 0.0
    except ValueError:
        return 0.0


def rescore(rows: list[dict[str, str]], profile: Profile, *, top_n: int = 50
            ) -> list[dict[str, str]]:
    """Rows re-ranked under ``profile``, with the same tie-breaks and disqualification as M2."""
    scored = []
    for row in rows:
        score = round(sum(weight * _number(row.get(name))
                          for name, weight in profile.weights.items()), 6)
        disqualified = ""
        identity = row.get("human_homolog_identity")
        if (profile.disqualify_human_homologs and identity not in (None, "")
                and _number(identity) >= DISQUALIFYING_HUMAN_HOMOLOG_IDENTITY):
            disqualified = f"human homolog {_number(identity):.1f}% identity"
        scored.append((score, disqualified, row))
    scored.sort(key=lambda item: (-item[0], -_number(item[2].get("pdb_evidence")),
                                  -_number(item[2].get("distinct_entries")),
                                  item[2].get("feature_id") or ""))
    out, selected = [], 0
    for rank, (score, disqualified, row) in enumerate(scored, start=1):
        chosen = not disqualified and selected < top_n
        selected += chosen
        out.append({"rank": str(rank), "feature_id": row.get("feature_id", ""),
                    "product": row.get("product", ""), "score": f"{score:g}",
                    "original_rank": row.get("rank", ""),
                    "original_score": row.get("triage_score", ""),
                    "selected": str(chosen), "disqualified": disqualified,
                    **{name: row.get(name, "") for name in WEIGHTS}})
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m s2f.m2_triage.rescore",
                                     description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run", required=True, help="run directory with m2_pdb/proteins.tsv")
    parser.add_argument("--profile", default=DEFAULT_PROFILE, choices=sorted(PROFILES))
    parser.add_argument("--weights", default="", metavar="JSON",
                        help="custom weights file; overrides --profile")
    parser.add_argument("--top", type=int, default=50, help="how many proteins to select")
    args = parser.parse_args(argv)

    try:
        profile = load_profile(args.profile, args.weights or None)
    except (OSError, ValueError) as exc:
        print(f"scoring profile: {exc}", file=sys.stderr)
        return 2
    source = Path(args.run) / "m2_pdb" / "proteins.tsv"
    if not source.is_file():
        print(f"no M2 output at {source}; run M2 first", file=sys.stderr)
        return 2
    with source.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    missing = [name for name in WEIGHTS if rows and name not in rows[0]]
    if missing:
        print(f"{source} lacks component columns {missing}; it predates per-component "
              "output and cannot be rescored", file=sys.stderr)
        return 2

    ranked = rescore(rows, profile, top_n=args.top)
    out_dir = source.parent / "rescored"
    out_dir.mkdir(exist_ok=True)
    target = out_dir / f"{profile.name}.tsv"
    with target.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(ranked[0]) if ranked else ["rank"],
                                delimiter="\t")
        writer.writeheader()
        writer.writerows(ranked)
    chosen = {r["feature_id"] for r in ranked if r["selected"] == "True"}
    before = {r.get("feature_id") for r in rows if r.get("selected") == "True"}
    (out_dir / f"{profile.name}.json").write_text(json.dumps({
        "profile": profile.as_dict(), "source": str(source), "top_n": args.top,
        "rescored_at": datetime.now(timezone.utc).isoformat(),
        "selected": len(chosen), "overlap_with_original_selection": len(chosen & before),
    }, indent=2) + "\n", encoding="utf-8")
    print(f"{len(rows)} proteins rescored under {profile.name!r} -> {target}")
    print(f"{len(chosen)} selected; {len(chosen & before)} of them were in the original "
          f"selection of {len(before)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Score M1's pathogenesis ranking against a curated truth set (issue #49).

    python -m s2f.evaluation --run runs/<run_id> [--truthset NAME] [--top-k 50]
    python -m s2f.evaluation --dry-run            # committed fixture, no network

Reads `proteins[]` from the run's `report.json` and writes `runs/<run_id>/eval/`:
`metrics.json`, `truth_matches.tsv`, `summary.md` and the `run.json` manifest. It never
writes a `report.json` section — see the note at the top of `metrics.py`.

Record the numbers **before** changing a weight. Pitfall #12: tuning after seeing the
results invalidates the measurement, and the weights fingerprint in `metrics.json` is what
lets the next run prove whether that happened.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from ..common.io import read_report, report_path, update_section
from ..m2_triage.score import WEIGHTS as TRIAGE_WEIGHTS
from . import agreement as agreement_mod
from .ablation import Ablation, ablate
from .metrics import BASELINES, RANDOM_SEED, Evaluation, evaluate
from .outcome import OutcomeSet, coverage_warning, outcomes_from_report
from .figures import save_figures
from .rankings import (
    SCORES, available_scores, build_ranking, population_stages, score_accessor,
    triage_components,
)
from .truthset import TruthSetError, bundled_truthsets, load_truthset, match_proteins

REPO_ROOT = Path(__file__).resolve().parents[2]
EVAL_FIXTURES = REPO_ROOT / "fixtures" / "eval"
REPORT_FIXTURE = EVAL_FIXTURES / "report.fixture.json"
DEFAULT_TRUTHSET = "kpneumoniae_hs11286"


def _seed_dry_run(run_dir: Path) -> None:
    """Populate an empty run directory from the committed evaluation fixture."""
    if report_path(run_dir).exists():
        return
    fixture = json.loads(REPORT_FIXTURE.read_text(encoding="utf-8"))
    for key, value in fixture.items():
        if key != "schema_version":
            update_section(run_dir, key, value)


def _resolve_truthset(name_or_path: str) -> Path:
    bundled = bundled_truthsets()
    if name_or_path in bundled:
        return bundled[name_or_path]
    path = Path(name_or_path)
    if path.exists():
        return path
    raise TruthSetError(
        f"no truth set {name_or_path!r}; bundled: {', '.join(sorted(bundled)) or '(none)'}"
    )


def _pct(value: float | None) -> str:
    return "n/a" if value is None else f"{value * 100:.1f}%"


def _signed(value: float | None) -> str:
    """A precision delta. `-` means the comparison was not available, not zero."""
    if value is None:
        return "&mdash;"
    return "0" if value == 0 else f"{value * 100:+.1f} pp"


def _write_matches_tsv(path: Path, result: Evaluation) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "rank", "score", "feature_id", "symbol", "label", "class",
        "match_tier", "in_top_k", "gene", "product",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t",
                                extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(result.rows)


def _write_summary(
    path: Path, result: Evaluation, *, run_id: str, dry_run: bool,
    ablation: Ablation | None = None,
    pair: "agreement_mod.Agreement | None" = None,
    outcomes: OutcomeSet | None = None,
    outcome_result: Evaluation | None = None,
    outcome_note: str | None = None,
) -> None:
    """A markdown summary short enough to paste into the issue as a comment."""
    measured, ties, leak = result.measured, result.ties, result.leakage
    lines = [
        f"# {result.score} calibration — {result.truthset}",
        "",
        f"Run `{run_id}` · {result.proteins_scored} scored proteins · top {measured.k} · "
        f"weights `{result.weights_fingerprint}`",
        "",
    ]
    if dry_run:
        lines += [
            "> **Fixture run — not a result.** These numbers come from "
            "`fixtures/eval/report.fixture.json`, a 39-protein plumbing check chosen to "
            "exercise every path in the harness. Quote figures from a full run only, the "
            "way `mgen_G37` is a shakedown genome and never a reported finding.",
            "",
        ]
    lines += [
        "## Truth-set coverage",
        "",
        f"- {result.truthset_counts['positive']} positive, "
        f"{result.truthset_counts['negative']} negative, "
        f"{result.truthset_counts['excluded']} excluded symbols curated",
        f"- {result.matched} proteins matched "
        f"({', '.join(f'{n} by {tier}' for tier, n in result.match_tiers.items()) or 'none'})",
        f"- {result.positives_in_genome} positive and {result.negatives_in_genome} negative "
        "symbols present in this genome",
        f"- {len(result.not_present)} curated symbols absent from this genome, "
        "excluded from the recall denominator",
        "",
        "## Measured",
        "",
        "| ranking | precision | recall | pos | neg | unlabelled | labelled coverage |",
        "| --- | --- | --- | --- | --- | --- | --- |",
        f"| **m1_priority** | **{_pct(measured.precision)}** | **{_pct(measured.recall)}** | "
        f"{measured.positives} | {measured.negatives} | {measured.unlabelled} | "
        f"{_pct(measured.labelled_coverage)} |",
    ]
    for name, score in result.baselines.items():
        lines.append(
            f"| baseline: {name} | {_pct(score.precision)} | {_pct(score.recall)} | "
            f"{score.positives} | {score.negatives} | {score.unlabelled} | "
            f"{_pct(score.labelled_coverage)} |"
        )
    lines += [
        "",
        "Precision is computed over the **labelled** subset of the top "
        f"{measured.k} ({measured.positives + measured.negatives} proteins), not over all "
        f"{measured.k}. Recall is over the {measured.recall_denominator} curated positives "
        "present in this genome.",
        "",
        "## Ranking behaviour",
        "",
        f"- score ceiling {ties.max_score:g}, shared by {ties.at_max} proteins",
        f"- the cut at {measured.k} falls on score {ties.cut_score:g}"
        if ties.cut_score is not None else f"- the cut at {measured.k} has no score",
        f"- {ties.at_cut_score} proteins share that score; the cut "
        + ("**splits a tie block**, so membership above the line is decided by the length "
           "tie-break, not by evidence" if ties.cut_inside_tie_block
           else "does not split a tie block"),
        "",
        "## Leakage audit",
        "",
        f"Of the {leak.labelled_in_top_k} labelled proteins in the top {measured.k}, score "
        f"points split {leak.curated_points} from curated databases against "
        f"{leak.text_points} from product text "
        f"({_pct(leak.text_share)} text-derived).",
        "",
        f"- {len(leak.curated_backed_positives)} true positives have curated-database "
        "evidence behind them",
        f"- {len(leak.text_only_positives)} true positives scored on product text alone, "
        "which is the scorer's own regex agreeing with the annotation string it read",
    ]
    if leak.text_only_positives:
        lines.append(f"  - {', '.join(sorted(leak.text_only_positives))}")
    if leak.unclassified_reasons:
        lines += [
            "",
            "> **Unrecognised score reasons** — `priority.py` emitted text this audit does "
            "not know how to attribute, so the split above is incomplete: "
            + "; ".join(leak.unclassified_reasons),
        ]
    lines += [
        "",
        "## Per class, inside the top " + str(measured.k),
        "",
        "| class | in genome | in top K |",
        "| --- | --- | --- |",
    ]
    for klass, counts in result.by_class.items():
        lines.append(f"| {klass} | {counts['in_genome']} | {counts['in_top_k']} |")
    if outcome_result is not None and outcomes is not None:
        om = outcome_result.measured
        lines += [
            "",
            "## Against M3's own outcome, no curation involved",
            "",
            f"`{result.score}` predicts \"worth folding\"; M3's quality gate reports whether a "
            "usable structure actually materialised. Same cut, objective label.",
            "",
            "| | n |",
            "| --- | --- |",
            f"| dockable (positive) | {outcomes.counts.get('dockable', 0)} |",
            f"| collected but rejected by the gate | {outcomes.counts.get('unusable', 0)} |",
            f"| no existing structure to reuse | {outcomes.counts.get('no_structure', 0)} |",
            f"| collection failed (excluded, not a negative) | {outcomes.counts.get('failed', 0)} |",
            "",
            f"**precision {_pct(om.precision)}** over the {om.positives + om.negatives} "
            f"labelled proteins in the top {om.k}: {om.positives} dockable, {om.negatives} not.",
            "",
            f"Base rate {_pct(outcomes.base_rate)} — lift "
            f"{_signed(outcomes.lift(om.precision))}. **Read the precision against the base "
            "rate, not on its own.**",
        ]
        if not outcomes.discriminates():
            lines += [
                "",
                "> **This label cannot discriminate on this run.** M3 passed "
                f"{outcomes.counts.get('dockable', 0)} of "
                f"{sum(v for k, v in outcomes.counts.items() if k != 'failed')} proteins, so "
                "almost any selection scores near-perfectly. That is expected and slightly "
                "circular: `pdb_evidence` is 0.40 of the triage score, and having a PDB hit "
                "is most of what makes a structure retrievable. The informative number here "
                "is the ablation's Δ against this label, not the headline precision.",
            ]
        if outcomes.policy_version:
            lines.append(f"Quality policy `{outcomes.policy_version}`.")
        if outcome_note:
            lines += ["", f"> {outcome_note}"]
        if outcomes.selected_but_uncollected:
            lines.append(
                f"> {len(outcomes.selected_but_uncollected)} selected proteins have no M3 "
                "record at all; M3 has not been re-run since the last M2 pass."
            )

    if ablation is not None:
        movers = ablation.load_bearing()
        lines += [
            "",
            f"## Which triage components decide the top {ablation.k}",
            "",
            "Leave-one-out: zero a component, re-apply M2's own "
            "`triage_score_from_components`, re-rank. No weight is changed and nothing is "
            "refit — this is what a different weight vector *would* have selected.",
            "",
            "| component | weight | fires on | moves out of top K | Δ precision vs truth | Δ precision vs M3 |",
            "| --- | ---: | ---: | ---: | ---: | ---: |",
        ]
        for effect in ablation.effects:
            lines.append(
                f"| `{effect.component}` | {effect.weight:+.2f} | {effect.nonzero_proteins} | "
                f"{effect.churn_at_k} | {_signed(effect.truth_delta_precision)} | "
                f"{_signed(effect.outcome_delta_precision)} |"
            )
        lines += [
            "",
            "A negative Δ means removing the component made precision *worse*, so it was "
            "earning its weight. A zero row means it changed nothing at this cut.",
            "",
            f"**{len(movers)} of {len(ablation.effects)} components move the selection at all.**",
        ]
        for note in ablation.notes:
            lines.append(f"- {note}")

    if pair is not None:
        rho = "n/a" if pair.spearman is None else f"{pair.spearman:+.3f}"
        lines += [
            "",
            "## m1_priority vs triage — the disagreement the contract predicts",
            "",
            "`00a-data-contract.md` says the two answer different questions and are expected "
            "to disagree. Measured over the "
            f"{pair.common} proteins both scored:",
            "",
            f"- Spearman rho **{rho}**",
            f"- **{pair.overlap_at_k} of {pair.k}** shared at the cut "
            f"(Jaccard {pair.jaccard_at_k:.2f})",
            f"- {len(pair.only_left)} picked by m1_priority alone, "
            f"{len(pair.only_right)} by triage alone",
        ]
        if pair.biggest_disagreements:
            lines += [
                "",
                "Where they diverge most:",
                "",
                "| m1 rank | triage rank | product |",
                "| ---: | ---: | --- |",
            ]
            for row in pair.biggest_disagreements[:8]:
                product = (row.get("product") or "")[:70]
                lines.append(
                    f"| {row['m1_priority_rank']} | {row['triage_rank']} | {product} |"
                )

    lines += [
        "",
        "---",
        "",
        f"Measured under weights `{result.weights_fingerprint}` "
        f"({json.dumps(result.weights, sort_keys=True)}). Recorded before any weight "
        "change, per pitfall #12. A later run under a different fingerprint is measuring a "
        "different scorer and the two numbers do not compare.",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def _git_state() -> dict[str, object]:
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=REPO_ROOT,
        capture_output=True, text=True, check=False,
    )
    status = subprocess.run(
        ["git", "status", "--porcelain", "--untracked-files=no"], cwd=REPO_ROOT,
        capture_output=True, text=True, check=False,
    )
    return {"commit": commit.stdout.strip() or None, "dirty": bool(status.stdout.strip())}


def _write_manifest(
    run_dir: Path, args: argparse.Namespace, result: Evaluation, *,
    started: datetime, elapsed_seconds: float, truthset_path: Path,
    input_report_sha256: str,
) -> None:
    manifest = {
        "module": "evaluation",
        "docs": "docs/08-evaluation.md",
        "issue": 49,
        "started_utc": started.isoformat(),
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": round(elapsed_seconds, 6),
        "command": getattr(args, "command", None),
        "input_report": str(report_path(run_dir)),
        "input_report_sha256": input_report_sha256,
        "truthset": {
            "name": result.truthset,
            "path": str(truthset_path),
            "sha256": hashlib.sha256(truthset_path.read_bytes()).hexdigest(),
            "counts": result.truthset_counts,
        },
        "parameters": {
            "top_k": result.measured.k,
            "baselines": list(result.baselines),
            "random_seed": args.seed,
            "include_product_matches": args.include_product_matches,
        },
        "scorer": {
            "module": "s2f.m1_genome.priority",
            "weights": result.weights,
            "weights_fingerprint": result.weights_fingerprint,
        },
        "dry_run": args.dry_run,
        "environment": {
            "python": platform.python_version(),
            "python_executable": sys.executable,
            "platform": platform.platform(),
        },
        "git": _git_state(),
        "counts": {
            "proteins_scored": result.proteins_scored,
            "matched": result.matched,
            "match_tiers": result.match_tiers,
            "positives_in_genome": result.positives_in_genome,
            "negatives_in_genome": result.negatives_in_genome,
            "not_present": len(result.not_present),
        },
    }
    out = run_dir / "eval" / "run.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def run(args: argparse.Namespace) -> int:
    started = datetime.now(timezone.utc)
    clock = time.monotonic()
    run_dir = Path(args.run)
    if args.dry_run:
        _seed_dry_run(run_dir)

    input_report = report_path(run_dir)
    if not input_report.exists():
        raise ValueError(f"{input_report} does not exist — run M1 first, or pass --dry-run")
    report = read_report(run_dir)
    proteins = report.get("proteins") or []
    if not proteins:
        raise ValueError(f"{input_report} has no proteins section")

    truthset_path = _resolve_truthset(args.truthset)
    truthset = load_truthset(truthset_path)

    # BV-BRC omits unset fields, and for some genomes `gene` is unset on every CDS — 0 of
    # 5,523 on HS11286, confirmed against the live Data API. Without this check the strong
    # tier reports "0 matched" and reads like an empty truth set rather than a missing
    # join column.
    if not any((p.get("gene") or "").strip() for p in proteins):
        print(
            f"NOTE: not one of the {len(proteins)} proteins carries a `gene` symbol, so the "
            "strong match tier has nothing to join on. BV-BRC leaves that field unset for "
            "some genomes."
            + ("" if args.include_product_matches else
               " Re-run with --include-product-matches to measure against product text "
               "instead, and read the leakage audit before quoting the number.")
        )

    matches, not_present = match_proteins(
        proteins, truthset, include_product_matches=args.include_product_matches
    )

    present = available_scores(proteins)
    if not present:
        raise ValueError(
            "this report carries neither m1_priority nor triage — run M1, and M2 for triage"
        )
    score_name = args.score or present[0]
    if score_name not in present:
        raise ValueError(
            f"{score_name!r} is not in this report; it carries: {', '.join(present)}"
        )
    read_score = score_accessor(score_name)
    result = evaluate(
        proteins, truthset, matches, not_present,
        k=args.top_k, baselines=args.baselines, seed=args.seed,
        ranker=lambda ps, n=score_name: build_ranking(ps, n).scored,
        score_of=read_score,
        score_name=score_name,
        weights=dict(TRIAGE_WEIGHTS) if score_name == "triage" else None,
    )

    # --- the three analyses that only exist once more of the pipeline has run ----------
    outcomes = outcomes_from_report(report, proteins)
    outcome_matches = outcomes.by_feature() or None
    ranking = build_ranking(proteins, score_name)
    outcome_note = (
        coverage_warning(outcomes, args.top_k, ranking.top(args.top_k))
        if outcome_matches else None
    )

    outcome_result: Evaluation | None = None
    if outcome_matches:
        outcome_result = evaluate(
            proteins, truthset, outcomes.matches, [],
            k=args.top_k, baselines=(), seed=args.seed,
            ranker=lambda ps, n=score_name: build_ranking(ps, n).scored,
            score_of=read_score, score_name=score_name,
            weights=dict(TRIAGE_WEIGHTS) if score_name == "triage" else None,
        )

    ablation: Ablation | None = None
    if any(triage_components(p) is not None for p in proteins):
        ablation = ablate(
            proteins,
            truth_matches={m.feature_id: m for m in matches} or None,
            outcome_matches=outcome_matches,
            k=args.top_k,
        )

    pair = None
    pair_common: list[str] = []
    right_ranking = None
    if len(present) > 1:
        left, right_ranking = (build_ranking(proteins, "m1_priority"),
                               build_ranking(proteins, "triage"))
        pair = agreement_mod.compare(left, right_ranking, proteins, k=args.top_k)
        pair_common = [f for f in left.order if f in right_ranking.scores]
        ranking = left

    eval_dir = run_dir / "eval"
    eval_dir.mkdir(parents=True, exist_ok=True)
    (eval_dir / "metrics.json").write_text(
        json.dumps(result.to_dict(), indent=2) + "\n", encoding="utf-8"
    )
    _write_matches_tsv(eval_dir / "truth_matches.tsv", result)
    if ablation is not None:
        (eval_dir / "ablation.json").write_text(
            json.dumps(ablation.to_dict(), indent=2) + "\n", encoding="utf-8"
        )
    if pair is not None:
        (eval_dir / "agreement.json").write_text(
            json.dumps(pair.to_dict(), indent=2) + "\n", encoding="utf-8"
        )
    if outcome_result is not None:
        (eval_dir / "outcome.json").write_text(
            json.dumps(
                {
                    "counts": outcomes.counts,
                    "base_rate": outcomes.base_rate,
                    "lift": outcomes.lift(outcome_result.measured.precision),
                    "discriminates": outcomes.discriminates(),
                    "policy_version": outcomes.policy_version,
                    "covered": outcomes.covered,
                    "selected_but_uncollected": outcomes.selected_but_uncollected,
                    "note": outcome_note,
                    "measured": outcome_result.to_dict()["measured"],
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
    figures: list = []
    if args.figures:
        outcome_payload = None
        if outcome_result is not None:
            outcome_payload = {
                "base_rate": outcomes.base_rate,
                "measured": {"precision": outcome_result.measured.precision},
            }
        figures = save_figures(
            eval_dir / "figures",
            stages=population_stages(proteins),
            ranks=(
                ([ranking.position()[f] for f in pair_common],
                 [right_ranking.position()[f] for f in pair_common])
                if pair is not None else None
            ),
            shared_at_k=pair.overlap_at_k if pair else None,
            spearman=pair.spearman if pair else None,
            k=args.top_k,
            ablation=ablation.to_dict() if ablation else None,
            outcome=outcome_payload,
        )
        if not figures:
            print("NOTE: matplotlib is not installed, so no figures were written.")

    _write_summary(eval_dir / "summary.md", result, run_id=run_dir.name,
                   dry_run=args.dry_run, ablation=ablation, pair=pair,
                   outcomes=outcomes if outcome_result else None,
                   outcome_result=outcome_result, outcome_note=outcome_note)
    _write_manifest(
        run_dir, args, result,
        started=started, elapsed_seconds=time.monotonic() - clock,
        truthset_path=truthset_path,
        input_report_sha256=hashlib.sha256(input_report.read_bytes()).hexdigest(),
    )

    measured = result.measured
    print(
        f"{result.truthset}: {result.matched} of {len(truthset)} curated symbols matched "
        f"across {result.proteins_scored} scored proteins.\n"
        f"top {measured.k}: precision {_pct(measured.precision)} "
        f"({measured.positives} positive, {measured.negatives} negative, "
        f"{measured.unlabelled} unlabelled), recall {_pct(measured.recall)} "
        f"of {measured.recall_denominator} present positives.\n"
        f"baselines: "
        + ", ".join(f"{name} {_pct(score.precision)}" for name, score in result.baselines.items())
        + f"\nceiling {result.ties.max_score:g} shared by {result.ties.at_max}; "
        f"cut {'splits' if result.ties.cut_inside_tie_block else 'does not split'} a tie block.\n"
        f"wrote {eval_dir}/summary.md"
    )
    if figures:
        print(f"figures: {len(figures)} files in {eval_dir}/figures")
    if outcome_result is not None:
        om = outcome_result.measured
        print(
            f"M3 outcome, same top {om.k}: precision {_pct(om.precision)} "
            f"({om.positives} dockable, {om.negatives} not, {om.unlabelled} uncollected); "
            f"base rate {_pct(outcomes.base_rate)}, lift {_signed(outcomes.lift(om.precision))}"
        )
        if not outcomes.discriminates():
            print(
                "  WARNING: the outcome label has almost no variance, so this precision is "
                "not evidence the ranking selects well. See summary.md."
            )
    if ablation is not None:
        movers = ablation.load_bearing()
        print(
            f"ablation: {len(movers)} of {len(ablation.effects)} components move the top "
            f"{ablation.k}"
            + (" — " + ", ".join(f"{e.component} {e.churn_at_k}" for e in movers[:4])
               if movers else "")
        )
    if pair is not None:
        rho = "n/a" if pair.spearman is None else f"{pair.spearman:+.3f}"
        print(
            f"m1_priority vs triage: spearman {rho}, "
            f"{pair.overlap_at_k} of {pair.k} shared at the cut"
        )
    if result.leakage.unclassified_reasons:
        print(
            "WARNING: priority.py emitted score reasons this audit cannot attribute; "
            "the curated/text split is incomplete. See summary.md."
        )
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="python -m s2f.evaluation", description=__doc__)
    parser.add_argument("--run", default="runs/dev", help="run directory containing report.json")
    parser.add_argument(
        "--truthset", default=DEFAULT_TRUTHSET,
        help=f"bundled truth-set name or a path to a TSV (default: {DEFAULT_TRUTHSET})",
    )
    parser.add_argument(
        "--score", default="", choices=("", *SCORES),
        help="which ranking to measure; default is m1_priority when present",
    )
    parser.add_argument("--top-k", type=int, default=50, help="rank cut to measure at (default: 50)")
    parser.add_argument(
        "--baselines", nargs="*", default=list(BASELINES), choices=list(BASELINES),
        help="trivial rankings to compare against",
    )
    parser.add_argument(
        "--include-product-matches", action="store_true",
        help="also match truth symbols as standalone tokens in product text (weaker tier)",
    )
    parser.add_argument("--seed", type=int, default=RANDOM_SEED,
                        help=f"seed for the random baseline (default: {RANDOM_SEED})")
    parser.add_argument("--figures", action="store_true",
                        help="also write PNG/SVG figures (needs matplotlib; optional)")
    parser.add_argument("--dry-run", action="store_true",
                        help="use the committed evaluation fixture; no network, no credentials")
    return parser


def main(argv: list[str] | None = None) -> int:
    raw_argv = list(sys.argv[1:] if argv is None else argv)
    args = build_parser().parse_args(raw_argv)
    args.command = [sys.executable, "-m", "s2f.evaluation", *raw_argv]
    if args.dry_run and args.run == "runs/dev":
        args.run = "runs/eval-dry-run"
    try:
        return run(args)
    except (OSError, ValueError, TruthSetError) as exc:
        print(f"Error: {exc}")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

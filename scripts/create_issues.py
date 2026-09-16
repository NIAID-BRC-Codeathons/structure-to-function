#!/usr/bin/env python3
"""Create GitHub issues from docs/issues.md.

Usage:
    python scripts/create_issues.py --dry-run          # print what would be created
    python scripts/create_issues.py                    # create them (needs gh CLI, authenticated)
    python scripts/create_issues.py --repo owner/name  # override repo

Requires the `gh` CLI (https://cli.github.com) authenticated with repo write access:
    gh auth login

Idempotent: skips any issue whose exact title already exists (open or closed).
Standard library only.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

DEFAULT_REPO = "NIAID-BRC-Codeathons/structure-to-function"
ISSUES_MD = Path(__file__).resolve().parent.parent / "docs" / "issues.md"

LABEL_COLORS = {
    "module:common": "5319e7",
    "module:m1": "1d76db",
    "module:m2": "0e8a16",
    "module:m3": "fbca04",
    "module:m4": "d93f0b",
    "module:m5": "b60205",
    "module:m6": "006b75",
    "priority:p0": "b60205",
    "priority:p1": "d93f0b",
    "priority:p2": "fef2c0",
    "type:code": "c5def5",
    "type:ops": "bfd4f2",
    "type:decision": "d4c5f9",
    "good-first-task": "7057ff",
}

CLAIM_FOOTER = (
    "\n\n---\n"
    "Unassigned on purpose: comment to claim this task, then assign yourself.\n"
)


def parse_issues(path: Path) -> list[dict]:
    text = path.read_text(encoding="utf-8")
    issues = []
    # Split on level-2 headings at line start.
    for chunk in re.split(r"^## ", text, flags=re.MULTILINE)[1:]:
        lines = chunk.splitlines()
        title = lines[0].strip()
        labels: list[str] = []
        body_start = 1
        for i, line in enumerate(lines[1:], start=1):
            if line.lower().startswith("labels:"):
                labels = [x.strip() for x in line.split(":", 1)[1].split(",") if x.strip()]
                body_start = i + 1
                break
            if line.strip():
                break
        body = "\n".join(lines[body_start:]).strip()
        issues.append({"title": title, "labels": labels, "body": body + CLAIM_FOOTER})
    return issues


def run(cmd: list[str], **kw) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, check=True, text=True, capture_output=True, **kw)


def ensure_gh() -> None:
    try:
        run(["gh", "--version"])
    except (FileNotFoundError, subprocess.CalledProcessError):
        sys.exit("gh CLI not found. Install https://cli.github.com and run `gh auth login`.")
    try:
        run(["gh", "auth", "status"])
    except subprocess.CalledProcessError:
        sys.exit("gh is not authenticated. Run `gh auth login`.")


def existing_titles(repo: str) -> set[str]:
    out = run([
        "gh", "issue", "list", "--repo", repo, "--state", "all",
        "--limit", "500", "--json", "title",
    ]).stdout
    return {row["title"] for row in json.loads(out or "[]")}


def ensure_labels(repo: str, labels: set[str]) -> None:
    for label in sorted(labels):
        color = LABEL_COLORS.get(label, "ededed")
        try:
            run(["gh", "label", "create", label, "--repo", repo, "--color", color, "--force"])
        except subprocess.CalledProcessError as exc:
            print(f"  ! label {label}: {exc.stderr.strip()}", file=sys.stderr)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", default=DEFAULT_REPO)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--file", type=Path, default=ISSUES_MD)
    args = ap.parse_args()

    issues = parse_issues(args.file)
    print(f"Parsed {len(issues)} issues from {args.file}")

    if args.dry_run:
        for issue in issues:
            print(f"\n- {issue['title']}\n  labels: {', '.join(issue['labels']) or '(none)'}")
            print("  body: " + issue["body"].splitlines()[0][:90] + " …")
        return

    ensure_gh()
    have = existing_titles(args.repo)
    ensure_labels(args.repo, {l for i in issues for l in i["labels"]})

    created = skipped = 0
    for issue in issues:
        if issue["title"] in have:
            print(f"skip (exists): {issue['title']}")
            skipped += 1
            continue
        cmd = ["gh", "issue", "create", "--repo", args.repo,
               "--title", issue["title"], "--body", issue["body"]]
        for label in issue["labels"]:
            cmd += ["--label", label]
        try:
            url = run(cmd).stdout.strip()
            print(f"created: {issue['title']}\n         {url}")
            created += 1
        except subprocess.CalledProcessError as exc:
            print(f"FAILED: {issue['title']}\n        {exc.stderr.strip()}", file=sys.stderr)

    print(f"\ncreated {created}, skipped {skipped}")


if __name__ == "__main__":
    main()

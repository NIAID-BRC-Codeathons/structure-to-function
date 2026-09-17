#!/usr/bin/env bash
# Build the environment for a run on Lambda (or any Linux node) and prove it works.
#
#   bash scripts/setup_lambda.sh                  # venv + deps + import check + tests
#   bash scripts/setup_lambda.sh --with-biolib    # also install pybiolib (DeepTMHMM)
#
# Stops at the first failure rather than reporting success over a broken install.
# Re-running reuses the venv. Prints a provider availability table at the end, so a run
# is planned against what is actually installed rather than what is hoped for.
set -euo pipefail

cd "$(dirname "$0")/.."
ROOT="$PWD"
VENV="$ROOT/.venv"
WITH_BIOLIB=0
for arg in "$@"; do
    case "$arg" in
        --with-biolib) WITH_BIOLIB=1 ;;
        *) echo "unknown argument: $arg" >&2; exit 2 ;;
    esac
done

# --- Python -----------------------------------------------------------------
# The repo targets 3.10+; jsonschema>=4.20 will not install below 3.8 at all.
PY=""
for candidate in python3.13 python3.12 python3.11 python3.10 python3; do
    if command -v "$candidate" >/dev/null 2>&1; then
        if "$candidate" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 10) else 1)'; then
            PY="$candidate"
            break
        fi
    fi
done
if [ -z "$PY" ]; then
    echo "No Python 3.10+ found. Load a module or point PATH at one first." >&2
    command -v python3 >/dev/null 2>&1 && echo "python3 is $(python3 -V 2>&1)" >&2
    exit 1
fi
echo "python: $($PY -V) ($(command -v $PY))"

# --- venv -------------------------------------------------------------------
if [ ! -d "$VENV" ]; then
    "$PY" -m venv "$VENV"
    echo "created $VENV"
else
    echo "reusing $VENV"
fi
# shellcheck disable=SC1091
source "$VENV/bin/activate"
python -m pip install --quiet --upgrade pip

python -m pip install --quiet -r requirements-common.txt -r requirements-m2.txt
python -m pip install --quiet pytest
if [ "$WITH_BIOLIB" = "1" ]; then
    python -m pip install --quiet pybiolib
fi
echo "dependencies installed"

# --- import check -----------------------------------------------------------
# Cheaper than the test suite and catches a broken install before pytest's collection does.
python - <<'PYCHECK'
import importlib
for module in (
    "s2f.common.schema", "s2f.common.io", "s2f.common.http", "s2f.common.ids",
    "s2f.m1_genome.parse", "s2f.m2_triage.bvbrc_input", "s2f.m2_triage.function",
    "s2f.m2_triage.report_adapter", "s2f.m2_triage.score", "s2f.run",
):
    importlib.import_module(module)
print("imports ok")
PYCHECK

# --- tests ------------------------------------------------------------------
python -m pytest -q

# --- what can actually annotate ---------------------------------------------
echo
echo "annotation providers on this node:"
python - <<'PYTOOLS'
import shutil
from s2f.m2_triage.function import discover_interproscan

install = discover_interproscan()
if install.available:
    print(f"  interproscan  FOUND  {install.path} ({install.version or 'version unknown'}, via {install.found_via})")
else:
    print(f"  interproscan  no     {install.error}")

for label, command, hint in (
    ("deeptmhmm", "biolib", "pip install pybiolib, then: biolib run DTU/DeepTMHMM --fasta <faa>"),
    ("signalp6", "signalp6", "free academic download from services.healthtech.dtu.dk"),
    ("psortb", "psortb", "container: brinkmanlab/psortb_commandline"),
    ("eggnog", "emapper.py", "conda install -c bioconda eggnog-mapper (~50 GB database)"),
    ("docker", "docker", "needed only for the PSORTb container"),
):
    path = shutil.which(command)
    print(f"  {label:<13} {'FOUND  ' + path if path else 'no     ' + hint}")
print()
print("Missing providers are not an error: the module falls back to its built-in")
print("sequence heuristic and records heuristic_only in run.json.")
PYTOOLS

echo
echo "Next: docs/02d-lambda-runbook.md"
echo "  source .venv/bin/activate"
echo "  python -m s2f.run --run runs/<id> --from-cga-dir data/<cga> --annotate"

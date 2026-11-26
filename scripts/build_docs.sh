#!/usr/bin/env bash
set -euo pipefail
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
if [[ -n "${BUTSON_DOCS_PYTHON:-}" ]]; then
  PYTHON_BIN="$BUTSON_DOCS_PYTHON"
elif [[ -x "$ROOT/.venv/bin/python" ]]; then
  PYTHON_BIN="$ROOT/.venv/bin/python"
else
  PYTHON_BIN=python3
fi
export PYTHONPATH="$ROOT/src:$ROOT${PYTHONPATH:+:$PYTHONPATH}"
cd "$ROOT"
BUILD_KIND="${1:-html}"
case "$BUILD_KIND" in
  html|doctest|linkcheck)
    "$PYTHON_BIN" -m sphinx -b "$BUILD_KIND" -d "docs/_build/doctrees/$BUILD_KIND" -n -W --keep-going docs/sphinx "docs/_build/$BUILD_KIND"
    ;;
  pdf)
    "$PYTHON_BIN" -m sphinx -M latexpdf docs/sphinx docs/_build -n -W --keep-going
    ;;
  all)
    "$0" html
    "$0" doctest
    "$0" pdf
    ;;
  *)
    printf 'Usage: %s [html|doctest|pdf|all|linkcheck]\n' "$0" >&2
    exit 2
    ;;
esac

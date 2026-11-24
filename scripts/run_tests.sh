#!/usr/bin/env bash
set -euo pipefail
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
export PYTHONPATH="$ROOT/src:$ROOT${PYTHONPATH:+:$PYTHONPATH}"
if [[ -n "${BUTSON_TOWERS_PYTHON:-}" ]]; then
  PYTHON_BIN="$BUTSON_TOWERS_PYTHON"
elif [[ -x "$ROOT/.venv/bin/python" ]]; then
  PYTHON_BIN="$ROOT/.venv/bin/python"
else
  PYTHON_BIN=python3
fi
cd "$ROOT"
if [[ "$#" -gt 0 ]]; then
  "$PYTHON_BIN" -m pytest -q "$@"
else
  "$PYTHON_BIN" -m pytest -q "$ROOT/tests"
fi

#!/usr/bin/env bash
set -euo pipefail
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT/formal"
python3 build.py
PYTHONPATH="$ROOT/release/lean_verification/source${PYTHONPATH:+:$PYTHONPATH}" python3 - <<'PY_AUDIT'
from pathlib import Path
from verify import parse_axioms, run
reports = parse_axioms(run(["lake", "env", "lean", "AuditAxioms.lean"], Path.cwd()))
print(f"Validated {len(reports)} theorem reports: standard axioms only.")
PY_AUDIT

PYTHONPATH="$ROOT/release/lean_sixth_root_pair/source${PYTHONPATH:+:$PYTHONPATH}" python3 - <<'PY_PAIR_AUDIT'
from pathlib import Path
from verify import parse_axioms, run
reports = parse_axioms(run(["lake", "env", "lean", "AuditSixthRootPair.lean"], Path.cwd()))
print(f"Validated {len(reports)} sixth-root matrix reports: standard axioms only.")
PY_PAIR_AUDIT

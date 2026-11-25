#!/usr/bin/env python3
"""Build local Lean modules in dependency order with bounded parallelism."""
from __future__ import annotations

import argparse
from pathlib import Path
import re
import subprocess


def module_batches(root: Path, jobs: int):
    modules = {}
    for path in sorted(root.glob("ButsonTowers/**/*.lean")) + [root / "ButsonTowers.lean"]:
        name = ".".join(path.relative_to(root).with_suffix("").parts)
        imports = set(re.findall(r"^import\s+([A-Za-z0-9_.]+)", path.read_text(), re.MULTILINE))
        modules[name] = imports
    remaining = set(modules)
    while remaining:
        ready = sorted(name for name in remaining if not (modules[name] & remaining))
        if not ready:
            raise RuntimeError("cyclic local Lean imports")
        # Check the other proofs before starting the large numerical witness.
        light = [name for name in ready if ".Data.T631" not in name]
        batch = (light or ready)[:jobs]
        yield batch
        remaining.difference_update(batch)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lake", default="lake")
    parser.add_argument("--jobs", type=int, default=3,
                        help="maximum simultaneously scheduled local modules (default:3)")
    parser.add_argument("--no-cache", action="store_true")
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    root = Path(__file__).resolve().parent
    base = [args.lake] + (["--no-cache"] if args.no_cache else []) + ["build"]
    for batch in module_batches(root, args.jobs):
        print("+ " + " ".join(base + batch), flush=True)
        subprocess.run(base + batch, cwd=root, check=True)
    subprocess.run(base, cwd=root, check=True)


if __name__ == "__main__":
    main()

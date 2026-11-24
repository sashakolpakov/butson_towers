#!/usr/bin/env python3
"""Build the source-only reproducibility ZIP from an explicit public whitelist."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import zipfile

from verify import (MATHLIB_REVISION, PAYLOAD_FILES, SOURCE_FILES, SUPPORT_FILES,
                    THEOREMS, TOOLCHAIN, digest, require)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, default=Path(__file__).resolve().parent,
                        help="repository or extracted package containing formal/ and docs/FORMALIZATION.md")
    parser.add_argument("--output", type=Path, required=True, help="output ZIP path")
    args = parser.parse_args()
    source, support = args.source_root.resolve(), Path(__file__).resolve().parent
    files = {}
    for relative in PAYLOAD_FILES:
        root = source if relative in SOURCE_FILES else support
        path = root / relative
        require(all(not part.endswith(".local") and part != ".local"
                    for part in Path(relative).parts), "local file cannot enter package")
        require(not path.is_symlink() and path.is_file(), f"missing regular source: {relative}")
        require(all(not p.is_symlink() for p in path.parents if p.is_relative_to(root)),
                f"symlinked source directory: {relative}")
        require(path.resolve().is_relative_to(root), f"source escapes package root: {relative}")
        files[relative] = path.read_bytes()
    manifest = {"schema": "butson-towers-lean-reproducibility-v1", "lean_toolchain": TOOLCHAIN,
                "mathlib_revision": MATHLIB_REVISION, "expected_theorems": list(THEOREMS),
                "scope": "partial formal construction kernel; not formal certification of every concrete matrix",
                "packaging_policy": "explicit public source whitelist; no .lake, symlinks, .local files, API credentials or submission parameters",
                "files": [{"path": path, "bytes": len(files[path]), "sha256": digest(files[path])}
                          for path in PAYLOAD_FILES]}
    files["MANIFEST.json"] = (json.dumps(manifest, indent=2) + "\n").encode()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(args.output, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for path in sorted(files):
            info = zipfile.ZipInfo("butson-towers-lean-repro/" + path, (1980,1,1,0,0,0))
            info.external_attr = 0o100644 << 16
            info.compress_type = zipfile.ZIP_DEFLATED
            archive.writestr(info, files[path])
    print(json.dumps({"archive": str(args.output), "files": len(files),
                      "sha256": digest(args.output.read_bytes()), "theorems": len(THEOREMS)}))


if __name__ == "__main__":
    main()

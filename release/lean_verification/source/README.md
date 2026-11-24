# Butson Towers: Lean verification package

This source-only package verifies **sixteen general construction theorems**:
weighted scalar Gram transport, complete framed-matrix Gram transport,
separate sign flatness and scale/order conditions, and the algebraic Spence
recurrences with their sign test. It does **not** claim that all concrete
T witnesses or eighteen recorded matrix endpoints have been formally
instantiated. See `docs/FORMALIZATION.md` for the exact boundary.

The package pins Lean **4.32.0**, mathlib commit
`81a5d257c8e410db227a6665ed08f64fea08e997`, and every transitive Git
dependency in `formal/lake-manifest.json`.

## Cold reproduction

Requirements: Python3.9+, Git, and elan/Lake on PATH. If the pinned
toolchain is not installed, install it using elan:

```sh
elan toolchain install leanprover/lean4:v4.32.0
python3 verify.py --report verification-report.json
```

The verifier runs the following in `formal/`, without changing the lock:

```sh
lake exe cache get
lake build
lake env lean AuditAxioms.lean
```

Lake automatically fetches missing dependency **source checkouts at the
manifest-locked revisions** when loading the workspace. The first command
downloads matching official precompiled mathlib artifacts to avoid a large
dependency build. **Do not run `lake update`** for this reproduction: it is
an upgrade operation and may rewrite the locked transitive revisions.

For a dependency build without precompiled artifact downloads:

```sh
python3 verify.py --from-source --report verification-report.json
```

This still needs network access to obtain absent locked source repositories
and may take substantially longer. Standard dependency caches are trusted
as imported Lean libraries; this package checks its own proof files, source
lock and axiom dependencies, rather than rebuilding Lean itself.

## Reuse an installed matching dependency cache

```sh
python3 verify.py --cache-dir /path/to/existing/.lake/packages \
  --report verification-report.json
```

The verifier checks **all nine Git HEADs** against the lock and rejects
tracked modifications before linking any cache checkout. It does not
install or upgrade packages in this mode. The supplied cache must include
the matching built mathlib dependencies. Local dependency links and proof
build outputs live in the ignored `formal/.lake/`; no machine-specific path
is part of the portable payload. An explicit Lake binary may be selected
with `--lake /path/to/lake`; its Lean version is checked.

## What constitutes success

`verify.py` checks all whitelisted file sizes and SHA256 digests, the exact
toolchain and dependency revisions, a successful Lean build, and exactly
the sixteen expected `#print axioms` reports. It accepts only the standard
axioms `propext`, `Classical.choice`, and `Quot.sound`. Missing or repeated
reports, `sorryAx`, project axioms, and unrelated theorem names cause failure.
The proof sources use no `sorry`, `admit`, project axiom or native evaluator.

The final JSON reports `success:true` and `proof_theorems_checked:16`.
The axiom-output parser has an independent small regression test:

```sh
python3 verify.py --self-test
```

## Package contents and rebuilding

The payload is an explicit whitelist of the eight formal project files,
the formalization scope document, this README, the verifier and packager,
plus `MANIFEST.json`. No `.lake` directory, compiled proof, symlink,
`.local` file, API credential or submission-parameter file is included.
The packager reads only these named public source files; it never scans a
repository for credentials or reads local configuration.

```sh
python3 build_package.py --output /tmp/butson-towers-lean-repro.zip
```

The ZIP has stable ordering, permissions and timestamps. To build from a
repository checkout while using these packaging scripts, pass
`--source-root /path/to/butson_towers`. Checksums detect accidental changes;
they are not a cryptographic signature or a substitute for Lean checking.

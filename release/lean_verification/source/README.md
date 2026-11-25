# Butson Towers: Lean verification package

This source-only package checks the **eighteen announced Hadamard orders
and all members of five Spence families**, from literal finite witnesses.
It includes kernel evaluations of both T correlations, all 27 Williamson
sources, five Paley matrices, five strong starters and the order-13 unreal
seed, connected to the general matrix construction proofs. See
`docs/FORMALIZATION.md` for the named statements and conventions.

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
python3 build.py
lake env lean AuditAxioms.lean
```

`build.py` schedules at most three independent local modules at once,
then runs a final `lake build`. This bounds memory use during the finite
kernel checks. Use `python3 build.py --jobs 1` for serial local checking.

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
the exact expected `#print axioms` reports. It accepts only the standard
axioms `propext`, `Classical.choice`, and `Quot.sound`. Missing or repeated
reports, `sorryAx`, project axioms, and unrelated theorem names cause failure.
The proof sources use no `sorry`, `admit`, project axiom or native evaluator.

The final JSON reports `success:true` and the number of checked theorem reports.
The axiom-output parser has an independent small regression test:

```sh
python3 verify.py --self-test
```

## Package contents and rebuilding

The payload is an explicit whitelist of the complete formal project,
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

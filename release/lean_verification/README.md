# Lean verification release

[Download the source-only reproducibility package](butson-towers-lean-repro.zip).
Its SHA256 is recorded in `SHA256SUMS`. The complete unpacked source is also
under [`source/`](source/README.md).

The packet pins Lean 4.32.0, mathlib commit
`81a5d257c8e410db227a6665ed08f64fea08e997` and all nine dependency revisions.
It verifies the sixteen current construction-kernel theorems, checks every
payload checksum and rejects extra theorem axioms. It contains no local
configuration, credentials, submission parameters, dependency cache or
compiled proof files.

After extracting the ZIP:

```sh
cd butson-towers-lean-repro
python3 verify.py --report verification-report.json
```

Elan/Lake, Git and Python are required. The pinned toolchain installation
and the source-only build option are documented in the packet README.
An existing matching dependency cache can be used with `--cache-dir`.
Fresh extraction and replay passed using the pinned local cache; the
recorded result is [`verification.json`](verification.json).

Rebuild the same deterministic ZIP from this repository:

```sh
python3 release/lean_verification/source/build_package.py \
  --source-root . --output /tmp/butson-towers-lean-repro.zip
```

The concrete T631/T127 witnesses and complete instantiation of all matrix
certificates remain the next formalization stage; see
[`docs/FORMALIZATION.md`](../../docs/FORMALIZATION.md).

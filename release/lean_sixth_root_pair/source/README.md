# Rational sixth-root matrix identities: Lean supplement

This small source-only package proves the ten algebraic statements in
`formal/ButsonTowers/SixthRootPair.lean`. Given PPᵀ=QQᵀ=dI,
PQᵀ+QPᵀ=dI and d≠0, it derives the identities for R=QPᵀ/d and
S=2R−I used in the detailed manuscript, including the multiplication
and conjugation formulas in coordinates a+bζ.

It supplements the separately preserved complete existence package for
the eighteen ranks and six families. It does not construct a cyclotomic
field, instantiate the literal T631-derived P,Q, or assert equality with
a saved matrix. Those distinctions are recorded in the manuscript audit.

The package pins Lean4.32.0, mathlib revision
81a5d257c8e410db227a6665ed08f64fea08e997 and all nine dependencies.
It contains only Mathlib-dependent proofs, its exact dependency lock,
checksum verifier and packager. No local configuration or compiled
artifacts enter the ZIP.

After extracting the ZIP:

```sh
python3 verify.py --report verification.json
```

For an existing matching built dependency cache:

```sh
python3 verify.py --cache-dir /path/to/.lake/packages --report verification.json
```

Use `--from-source` to build the pinned dependencies without downloading
precompiled artifacts. The default cold command obtains the official
matching cache. Missing source repositories require network access.
The verifier checks every locked Git revision and rejects tracked
modifications in an existing cache. It then builds this package's proof
sources, checks exactly ten named axiom reports and allows only
propext, Classical.choice and Quot.sound. No new axioms, placeholders
or native evaluation are permitted. Run `python3 verify.py --self-test`
for the strict report parser's regression checks.

Rebuild the deterministic source ZIP from this repository:

```sh
python3 release/lean_sixth_root_pair/source/build_package.py --source-root . \
  --output /tmp/butson-towers-sixth-root-pair.zip
```

From an extracted packet, run `python3 build_package.py --output output.zip`.
The tiny root library in this packet imports only SixthRootPair; the
repository's main library and original existence package are unchanged.

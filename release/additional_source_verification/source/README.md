# Current additional-source Lean verification

This source-only package verifies the canonical Butson Towers development
with Lean 4.32.0 and mathlib revision
`81a5d257c8e410db227a6665ed08f64fea08e997`. All nine dependency revisions
are locked. It contains no build cache, local configuration or credentials.

The additional scope is twelve literal amicable source quadruples and
sixteen finite endpoint existence theorems, six additional Spence
recurrences, and the shifted family `15131380*109^r`. Together with the
included earlier development, these prove the manuscript's 34 finite
orders and nine announced families. Additional family constructions serve
as controls. The propositions concern existence; novelty is determined by
separate literature checks.

The twelve source widths are 47, 53, 59, 67, 69, 73, 79, 93, 97, 103, 109
and 127. Width 69 uses ordinary symmetric circulants; the other sources
reverse the first circulant's columns. The six additional starter
parameters are 5, 53, 97, 109, 113 and 137. All literal signs, necessary
symmetries, complete correlation coefficients, both strong starter norms
and Paley identities are checked by kernel reduction. Generic algebraic
proofs then give the complete matrix Grams and all-parameter recurrences.

Install Python 3.10 or newer, Git and [elan](https://github.com/leanprover/elan), then:

```sh
elan toolchain install leanprover/lean4:v4.32.0
python3 verify.py --report /tmp/additional-source-verification.json
```

For matching dependency checkouts already available locally:

```sh
python3 verify.py --cache-dir /path/to/.lake/packages \
  --report /tmp/additional-source-verification.json
```

Every dependency's Git revision and tracked tree are checked before use.
Only dependency caches are reused; the package's own proof modules are
built locally with bounded parallelism. `--from-source` compiles dependencies
without downloading their compiled cache. The full proof replay can take
several minutes and needs the memory required by the pinned mathlib build.

A fast payload check, which does not claim a Lean replay:

```sh
python3 verify.py --self-test
python3 verify.py --check-only
```

Verification enforces the exact payload whitelist and every SHA256,
checks local import closure and both literal-data bindings, and runs three
strict axiom audits. Their sets contain 85 earlier reports, 10 sixth-root
matrix reports and exactly 48 additional reports. Missing, duplicate or
unexpected reports, additional output and nonstandard axioms are rejected.
Only `propext`, `Classical.choice` and `Quot.sound` are permitted.

After a successful build, rerun the additional audit directly with:

```sh
cd formal
lake env lean AuditAdditionalSources.lean
```

`formal/NEW_SOURCE_VERIFICATION.md` describes the new proof modules.
The earlier submitted release packages have not been altered. To recreate
this deterministic archive from a canonical repository checkout:

```sh
python3 build_package.py --source-root /path/to/butson_towers \
  --output /tmp/butson-towers-additional-sources.zip
```

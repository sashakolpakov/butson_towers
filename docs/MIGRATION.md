# Repository migration

This project was separated from `/Users/sasha/hadamard` on 24 November 2026
into `/Users/sasha/butson_towers`. `MIGRATION.json` records the source Git
revision and byte hashes of transferred research files, including the
uncommitted construction results present at separation.

The original layout and internal Python import names are retained so that
matrix recipes, certificate paths and mathematical source hashes continue
to replay. The installed distribution is named `butson-towers`; its small
shared arithmetic package retains the name `hadamard` for compatibility.
It is included locally and does not depend on the original checkout.

Construction engines, exact sources, certificates, manuscript, submission
tooling and dedicated tests were transferred. Shared arithmetic modules
were copied; the original finite-classification project retains its own
copies and tests. Generated search stores and historical evidence remain
ignored by Git and broad source searches.

`hexagon_key.local` and `hexagon_parameters.local` were transferred
separately with mode `0600`. They are excluded from Git and submission
archives. The public migration record contains neither their values nor
credential hashes. No remote draft or submission is part of the migration.

The local `.venv` is relocatable and contains the already verified Python
dependencies. A clean environment can instead be installed with:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -e '.[dev,search,native]'
./scripts/run_tests.sh
```

The Lean project pins its toolchain and mathlib revision separately under
`formal/`; see [formalization status](FORMALIZATION.md).

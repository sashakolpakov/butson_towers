# Construction walkthrough

The detailed exposition is now maintained in the
[Sphinx manual](sphinx/index.md), with proofs, implementation excerpts,
input/output contracts and executable examples at each step.

The main reading path is:

1. [Reconstruct a complete small example](sphinx/getting-started.md).
2. [Restrict scalars while retaining the norm](sphinx/scalar-restriction.md).
3. [Construct the disjoint design and the sixth-root/group-algebra endpoints](sphinx/construction.md).
4. [Generate exact sources and infinite output families](sphinx/families.md).
5. [Find the finite input](sphinx/exact-search.md) and [solve its signs](sphinx/sign-solvers.md).
6. [Run durable campaigns](sphinx/campaigns.md) and [verify the results](sphinx/verification.md).

Build the rendered manual from the repository root:

```sh
python -m pip install -e '.[docs,dev]' -c docs/requirements.txt
./scripts/build_docs.sh html
./scripts/build_docs.sh doctest
```

Then open `docs/_build/html/index.html`. The
[build guide](sphinx/contributing.md) also gives the printable PDF command.
The submitted manuscript and frozen verification packages retain their
separate, unchanged proof snapshots.

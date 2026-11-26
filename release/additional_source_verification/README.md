# Current additional-source verification

[Download the source-only package](butson-towers-additional-sources.zip), or
inspect the complete [unpacked source](source/README.md).
`SHA256SUMS` binds the deterministic archive.

This package includes the canonical Lean proof closure for the current
34 finite existence results and nine announced families, together with
additional source and family controls. It adds twelve literal amicable
quadruples, six strong starter/Paley pairs, sixteen finite endpoint
specializations and the shifted family `15131380*109^r`.

Lean 4.32.0, mathlib commit
`81a5d257c8e410db227a6665ed08f64fea08e997` and all dependency revisions are
pinned. The exact source whitelist, SHA256 manifest, two literal-data
checkers and three strict axiom audits are included. The audits cover
85 earlier reports, 10 sixth-root matrix reports and 48 additional reports.
Novelty labels remain separate literature claims.

After extracting the archive:

```sh
cd butson-towers-additional-sources
python3 verify.py --report /tmp/additional-source-verification.json
```

Use `--check-only` for payload/literal checks, `--self-test` for the strict
parser tests, or `--cache-dir /path/to/.lake/packages` to use validated
matching dependency checkouts. The package contains no compiled proofs,
dependency cache, Git metadata or local configuration. Existing submitted
verification packages retain their original bytes and scope.

To recreate the deterministic archive from this repository:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 \
  release/additional_source_verification/source/build_package.py \
  --source-root . --output /tmp/butson-towers-additional-sources.zip
```

The [verification record](verification.json) records a successful cold source
replay (857.71 seconds) and final archive validation, with all 143 axiom
reports. It also binds the unchanged Lean/input files across the final
documentation and formatting update.

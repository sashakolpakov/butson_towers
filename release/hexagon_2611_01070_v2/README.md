# Canceled oversized Hexagon version-2 package

This oversized package for *Butson towers and high-rank Hadamard matrices*
was **canceled by the submitter** because it included unwanted ancillary files.
The [cancellation record](cancellation-status.json) supersedes the initial
[submission receipt](submission-status.json). The replacement is the
[manuscript-only revision](../hexagon_2611_01070_manuscript/README.md), with
no ancillary files or uploaded source archive. The published
[first version](https://hexagonmath.org/2611.01070v1) remains unchanged.
The package contents and replay records below are preserved as historical data.

The revised paper proves **34 new orders and nine new infinite families**,
with the full order list and family formulas in its main theorem. A separate
table identifies positive prior constructions at all **31 recovered known
orders**. The explicit order-885924 sixth-root representative is inequivalent
to the standard Paley matrix: their ranks over the field with three elements
are 246090 and 442962. This comparison has a mathematical proof and an exact
finite certificate; the Lean existence proofs cover the 34 announced orders
and nine families, with their scope stated separately.

The release contains the exact TeX source, compiled PDF, synchronized
metadata, an explicit source archive and a member-by-member hash manifest.
Ancillary files preserve the canonical Python package paths, the complete
source data and replayable verification records. The current source-only
Lean package is pinned by immutable repository URL and SHA256 rather than
included as a nested ZIP or as files unsupported by the submission service.
No `.local` file, credential, private receipt or local configuration belongs
to this release.

The version uploader is restricted to the revision route of the original
deposit. It refuses the API operation that creates a new article. Before
committing the revision it requires an exact remote source-tree match and a
successful server PDF preview. Administrative choices are loaded from the
existing private preferences; the title, abstract and version comments come
from this revised manuscript.

## Reproduce the release

The local builder makes no network requests. It prepared this package with
Python, NumPy, SymPy and `latexmk`. After the service receipt is saved, the
builder checks the frozen bytes and refuses to overwrite them for a changed
manuscript. The independent verifier replays the submitted archive without
rewriting it:

```sh
python release/hexagon_2611_01070_v2/build.py
python release/hexagon_2611_01070_v2/verify.py --report /tmp/butson-v2-verification.json
```

The verifier extracts a fresh source tree, checks every archive member and
hash, and runs twelve exact replays. Add `--with-scalar-algebra` with
`python-flint` installed for the thirteenth check: full rational field-factor,
module-multiplicity, metric and reconstruction verification of all nineteen
additional compatible sources. Child processes use only the extracted sources
and explicitly selected third-party packages; editable-install hooks are disabled.
The saved [verification record](archive-verification.json) includes all thirteen.
It does not rerun the Lean kernel; [proof-sources.json](proof-sources.json) pins
the separately verified complete Lean source package and its 143-report audit.

For just the manuscript's finite calculations after extracting `source.zip`:

```sh
python3 anc/manuscript_data.py
```

Other canonical replays are listed with exact arguments in `anc/payload.json`.
For those commands, enter `anc/` and use `PYTHONPATH=.:src`; they require NumPy
and SymPy. The manuscript checker itself uses only the Python standard library.

[build.json](build.json) binds the warning-free 28-page PDF to the exact source.
[submission-manifest.json](submission-manifest.json) binds the archive's 81
members, metadata and reference PDF. [audit.json](audit.json) records the
publication checks and their scope. The version uploader additionally uses
`pypdf` to inspect the service's compiled preview before commit.

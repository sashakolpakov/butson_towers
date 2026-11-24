# New Hadamard orders and infinite families: Hexagon announcement

[Read the five-page announcement](announcement.pdf).
The canonical mathematical source is
[`manuscript/hadamard_7572_announcement.tex`](../../manuscript/hadamard_7572_announcement.tex);
`main.tex` is its byte-identical API snapshot. The PDF compiles with no
warnings, and all five pages have been visually reviewed.

## Source package

`source.zip` contains exactly `main.tex` and 26 ancillary files under
`anc/`. The ancillary material includes the complete sequence data, the
Python/NumPy recipe, the certificate, the dated novelty audit and the three
independent verification records. Its manifest records all member hashes.
The source archive contains no API credential, local receipt, downloaded
third-party paper, generated PDF or NPZ file.

After extracting the ZIP, reconstruct and verify the matrix with:

```sh
python3 anc/recipe.py --write h7572-reconstructed.npz
python3 anc/reuse.py
python3 anc/t127/recipe.py --write-hadamard h1524-reconstructed.npz
python3 anc/williamson_scalar_endpoints.py --output weighted-replay.json
python3 anc/williamson_scalar_endpoints.py --families --output families-replay.json
```

The replay needs Python 3.10+ and NumPy. It checks all 28,663,806 distinct row
pairs and reproduces the original packed matrix hash. The immutable
repository bundle remains in `artifacts/hadamard-7572/`; the source archive
omits its README and binary matrix, whose reference hashes remain recorded
in the original certificate.

The second command rechecks the complete P/Q/cross identities, the three
additional finite sources, and the known unreal tripling rule through
orders 39, 117 and 351. It certifies literal endpoint recipes at orders
17668, 22716, 32812, 98436, 295308 and 885924 by exact algebraic Gram
identities. It does not allocate or claim a dense audit of these six larger
matrices. Orders 17668, 22716, 32812, 98436 and 295308 are marked **new to
the best of the author's knowledge**, alongside H7572. Order 885924 is
marked **known**, by Paley's construction.
The complete portable result is retained in `reuse-verification.json`.

The third command reconstructs the periodic T127 witness from its 15 states,
checks the four complete circulant Gram identities, and writes and fully
verifies its known-order H1524 control. Its complete data, standalone script,
certificate, independent audit and dated literature audit are under
`anc/t127/`. Periodic T-matrix order 127 is marked **new to the best of the
author's knowledge**; the survey lists it as its smallest open T-matrix order.
The archive manifest records the current claims; the original dated audits
are retained as provenance.

The fourth command verifies 27 complete Williamson scalar restrictions.
Eight orders extend the earlier results: 42908, 58052, 73196, 78244, 93388,
108532, 128724 and 153964. The weighted source is over a product of
cyclotomic fields, with its exact normalized trace metric; the terminal
entry images are sign circulants. The fifth command verifies five
constructive infinite families and their first expanded sources. Four
further orders are new: 229684, 386172, 4772884 and 13001124; 18629644 is
known by Paley II. There are **18 new real orders** in total, under the
author-knowledge convention, alongside **20 recovered known orders** in
the main campaign. The source and endpoint Gram scales are retained in
each certificate; large endpoint matrices are accessed by entry and row.

The independent exact cyclotomic trace audit is also included. With SymPy
installed in addition to NumPy, replay it using:

```sh
python3 anc/williamson_cyclotomic_trace_audit.py --output trace-replay.json
```

`archive-verification.json` records an actual isolated replay from the ZIP,
with no reference matrix present. `build.json` binds the local PDF to its
exact TeX source. `submission-manifest.json` gives the source-archive
request, file size, selected entry point and hashes.

## API preparation

The current [official API](https://hexagonmath.org/api/v1/openapi.json) is
recorded by version and SHA256 in `api-contract.json`. It accepts scripts
and JSON verification inputs inside `anc/` in a source ZIP. The main entry
point after expansion is `main.tex`.

Rebuild the deterministic ZIP and check metadata readiness:

```sh
python3 release/hadamard_7572_announcement/prepare.py
```

This command makes no network request. `metadata.template.json` contains
the article's title, abstract, classifications and comments. When present,
the repository-root `hexagon_parameters.local` supplies the fixed
administrative choices authorized by the user's standing instruction:
identity, CC BY-SA 4.0, mixed text origin, human understanding `some`,
the review/provenance declaration and required confirmations. These choices
do not need to be requested again for each paper. The local file preserves
the exact reference metadata from
[2611.00913v1](https://hexagonmath.org/2611.00913v1), with that article's
geometry classifications kept separate from future article metadata.
This announcement retains `math.CO` and MSC `05B20`, `15B34`.

The API labels text origin using `no-human-author-claimed` (primarily
AI-generated text), `mixed`, or `human`; the first is a compatibility field
name, not a statement that the named human contributor is absent. The
provenance text combines the fixed review declaration with this paper's
author framing and automated mathematical work. The public work API does
not expose the confirmation booleans; their retained basis is the API's
publication requirements and the user's explicit standing authorization.

With valid metadata, `prepare.py` writes `draft-request.json`
and `draft-headers.json`, including a stable idempotency key bound to that
metadata and ZIP. The authenticated API workflow is draft creation, ZIP
registration/upload/completion/expansion, `main.tex` entry-point selection,
source check and private preview. Final commit is a separate API action.
The private preview concerns typesetting; the exact mathematical checks are
the retained local verifications.

The current preparation has fetched only public API information. It has
created no remote draft and made no submission. Readiness is a local
validation result, not permission to publish. The credential remains in
the repository-root `hexagon_key.local`; preferences are separate in
`hexagon_parameters.local`, both with mode `0600`. Existing `*.local` and
`.local/` ignore rules cover both files. Neither local file is included in
metadata or archive members, and the credential is never put in
command-line arguments or printed. Preparation reads it only in memory
for the local archive-content exclusion check.

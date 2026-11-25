# Butson towers and high-rank Hadamard matrices

Submitted to HexagonMath on 25 November 2026; currently under review.
The [submission record](submission-status.json) retains the exact package hashes;
the private receipt remains in the ignored local directory.

[Read the detailed manuscript](announcement.pdf) or download the
[TeX source and exact reconstruction recipes](source.zip).
The [canonical manuscript](../../manuscript/hadamard_7572_announcement.tex)
and `main.tex` are byte-identical. This versioned release preserves the
[previous article release](../butson_towers_v2/README.md), the
[original announcement](../hadamard_7572_announcement/README.md), and all
frozen construction bundles.

The manuscript records eighteen new real Hadamard ranks and six new
infinite output families, with individual known/new labels for the finite
ranks. Its title and concise abstract are synchronized with
`metadata.template.json` and the locally prepared API request.
This revision defines the constructions and notation explicitly and gives
detailed proofs, with the literal finite data and a standard-library
verification program included in the source archive.
`build.json` binds the reference PDF to the exact TeX source;
`submission-manifest.json` records every source-archive member and hash.

The archive includes only `main.tex` and explicit Python/JSON ancillary
files. It contains literal source data, standalone reconstruction recipes,
dated status records, exact verification records and the completed Lean
verification report. The separately downloadable
[source-only Lean reproducibility package](../lean_verification/README.md)
provides all formal proof sources. The archive's `anc/proof-sources.json`
pins that packet by immutable repository commit and SHA256; it does not
require the Hexagon service to process Lean files. The existing kernel
verification covers all eighteen announced ranks and all six infinite output
families (five using Spence recurrence and one using unreal tripling), from literal finite witnesses.

After extracting `source.zip`, the Python checks use only its files:

```sh
python3 anc/manuscript_data.py
python3 anc/recipe.py --write h7572-reconstructed.npz
python3 anc/reuse.py
python3 anc/t127/recipe.py --write-hadamard h1524-reconstructed.npz
python3 anc/williamson_scalar_endpoints.py --output weighted-replay.json
python3 anc/williamson_scalar_endpoints.py --families --output families-replay.json
python3 anc/williamson_cyclotomic_trace_audit.py --output trace-replay.json
```

The manuscript-data checker uses Python's standard library. The other
checks require Python 3.10+ and NumPy; the final independent trace audit also
uses SymPy. The H7572 recipe checks all 28,663,806 distinct row pairs of
H7572. Large endpoints are specified by literal entry/row recipes and
complete algebraic Gram proofs. The five Spence source recurrences and unreal tripling prove the
six infinite output families for every nonnegative parameter.

Rebuild this release from the repository with TeX Live/`latexmk`, Python,
NumPy and SymPy installed:

```sh
python3 release/butson_towers_v3/build.py
python3 release/butson_towers_v3/verify.py
```

The second command extracts a fresh temporary source tree, checks all
archive hashes and synchronized metadata, and executes all seven standalone
replays. Its record is `archive-verification.json`. The unchanged Lean
source packet remains separately hash checked; the editorial rebuild does
not rerun the Lean kernel. PDF page inspection is recorded separately in
`build.json` after all pages are reviewed.

Preparation uses the retained `api-contract.json` and reuses fixed
administrative defaults from the repository's private
`hexagon_parameters.local` when present. Credentials remain in the
separate `hexagon_key.local`; local files are excluded from the explicit
archive. Rebuilding makes no network requests and creates no remote draft,
submission or publication.

The [paper-to-Lean audit](../../docs/MANUSCRIPT_LEAN_AUDIT.md) checks the
actual proof chain in manuscript order. A fresh baseline kernel replay
checked 85 selected statements; a separate
[sixth-root matrix supplement](../lean_sixth_root_pair/README.md) checks
ten added identities. The general finite-field character proof and the
particular cyclotomic/Hermitian realization have the explicit scope
recorded in that audit. All eighteen ranks and all six families have the
complete existence proofs described there.

The complete review and validation record is [audit.json](audit.json).

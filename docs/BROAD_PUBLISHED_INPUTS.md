# Broad acquisition of fixed published matrices

The authors' [Aalto electronic supplement](https://wiki.aalto.fi/spaces/Butson/pages/120482304/Matrices%2Bup%2Bto%2Bmonomial%2Bequivalence)
provides complete small Butson matrices at the useful source orders 7, 9
and 13. Five raw datasets, totaling **84,288 bytes**, are now saved under
`experiments/scalar_lift/broad_catalogue_inputs/`, with pinned SHA-256
hashes. The [loader](../experiments/scalar_lift/broad_source_acquisition.py)
parses the authors' documented BH1 format, verifies every complete Gram,
and chooses a fixed monomial gauge when one passes the exact alphabet
test. No missing source is searched for or assembled from smaller ones.

| Published cell | Literal matrices checked | Gauge orbits retained |
|---|---:|---:|
| BH(7,6) | 2 | 2 |
| BH(7,12) | 4 | 4 |
| BH(9,6) | 17 | 2 |
| BH(9,12) | 65 | 15 |
| BH(13,6) | 436 | 10 |
| Total | **524** | **33** |

The [acquisition certificate](../certificates/scalar-lift/broad-source-acquisition.json)
records all 524 inputs, including the 491 excluded same-alphabet gauge
orbits. The full gauge audit visits only **11,896 row prefixes**.
These counts describe catalogue classes and selected representatives,
not 33 previously unknown mathematical objects. Some inputs may be
equivalent to ones already used in this workspace.

## Exact gauge scope

At phase six, avoiding an odd signed root relation means omitting one
antipodal residue class modulo three. At phase twelve, it means omitting
one even and one odd antipodal class modulo six. There are nine maximal
palettes of the latter kind. An even global phase shift normalizes the
missing even class to zero, leaving the three patterns `(0,1)`, `(0,3)`,
and `(0,5)`. All six row and column shifts modulo signs remain available.
The first row gauge can independently be set to zero by compensating all
column gauges. Therefore the CSP covers every root row/column gauge in
the stated sixth- or twelfth-root alphabet. Permutations preserve the
global alphabet and need no separate enumeration.

The CSP retains every successful literal representative, subtracting its
first entry as an explicit global anchor. Excluded cases may admit gauges
in a larger root alphabet; that possibility is outside this audit.
Passing this alphabet gate does not prove a normalized sign endpoint.

Tests compare the antipodal criterion with the independent exact
cyclotomic-relation engine for every nonempty phase-six and phase-twelve
palette. A separate unpruned gauge enumeration checks small examples with
both positive and negative results. Every retained input also passes a
fresh complete Gram check after its gauge is applied.

API: `sources()` lists all 524 descriptors;
`load_source('aalto-9-12-064')` returns the existing `VerifiedSource`
protocol. Every proof contains `alphabet_gate`, `monomial_gauge_audit`,
source hashes, and the explicit gauge. Root labels are the original
catalogue phases; an effective-field test should reduce a common divisor
of the actual exponents or of an explicitly chosen global anchor.

## Additional fixed Hadamard frames

Two order-764 matrices were acquired from the explicit coset seeds in
[Đoković, *Hadamard matrices of order 764 exist*](https://arxiv.org/pdf/math/0703312v1).
Both pass all 190 supplementary-difference conditions, the full integer
Gram, and an independent bit-distance check. They are saved in
`artifacts/published-h764-inputs/`; the API
`published_h764.known_h764(index=0 or1)` reconstructs them.
The [certificate](../certificates/scalar-lift/published-h764-inputs.json)
binds the full matrices to the printed seeds. They provide actual fixed
frames for the order-nine route to `6876=9*764`; compatibility with a
nonreal scalar action remains to be tested.

All five unbordered Table A1 seeds in
[the November Hexagon paper](https://hexagonmath.org/2611.00637v1)
are also acquired: orders **892, 1132, 1244, 1948, 1964**. The public TeX
source is pinned to SHA-256
`b40d2420405ce96ac56f154cb36adb071c7f67fa2f6b7bd7a23e5ae3b485e15f`.
Each source's multiplier orbits, every periodic correlation, the full
integer sign Gram, and independent bit-distance Gram are checked.
Five full bundles occupy about 1.65 MB. The API is
`published_alpoge_gs.known_gs_frame(order)` and the
[certificate is saved](../certificates/scalar-lift/published-alpoge-gs-inputs.json).
None of these particular printed quadruples has even one symmetric
circulant block, so they do not directly meet the two-symmetric-block
complex array's sufficient hypothesis. No assertion about equivalent
quadruples is made.

The apparent fifth-root route to order 23960 needs a normalized frame
of order 4792, which remains unresolved in the attributed frontier.
It cannot be supplied as an already-known frame. The registry label
`fifth-5-q20` embeds a fifth-root matrix, so its effective entry field
is Q(zeta5), not automatically the degree-eight field Q(zeta20).

Replay:

```sh
PYTHONPATH=src:. python -m experiments.scalar_lift.broad_source_acquisition --check certificates/scalar-lift/broad-source-acquisition.json
PYTHONPATH=src:. python -m experiments.scalar_lift.published_h764 --check certificates/scalar-lift/published-h764-inputs.json
PYTHONPATH=src:. python -m experiments.scalar_lift.published_alpoge_gs --check certificates/scalar-lift/published-alpoge-gs-inputs.json
```

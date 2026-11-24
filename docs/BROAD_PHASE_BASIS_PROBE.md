# Executed shared mixed-basis probe across published source classes

The [completed certificate](../certificates/scalar-lift/broad-phase-basis-probe.json)
records **512 mixed-basis geometries**, shared across **39,104 actual
source/basis checks** and **3,520 palette checks**. It covers eight
published real frames and eight target orders. Every declared case
finished; no sign Hadamard was found.

The fixed sources come from the original verified registry and the
[expanded published catalogue acquisition](BROAD_PUBLISHED_INPUTS.md).
Of the 289 necessary-test-surviving routes in the exact phase-independent
preflight, this execution covers **198 routes** compatible with the
available frames and the field `Q(zeta12)`. It retains **234 literal
embedded source/anchor matrices**, with all **238 original source/global
anchor recipes**. Only identical complete exponent matrices are merged;
no arbitrary monomial or Galois equivalence is used for that pruning.

Declared source phases 6, 12 and 24 occur. Every actual shifted entry is
embedded by the exact integer formula

`e12 = 12*((e_original+global_shift) mod declared_phase)/declared_phase`.

This includes sources declared at phase 24 whose actual shifted roots
generate phase 6 or 12. It does not misidentify a field merely from its
declared phase label, or embed fifth roots into the twelfth-root field.
The full Gram and hash of every actual embedded matrix are saved.

| Fixed published frame | Embedded source/anchor bindings | Palettes | Mixed bases | Source/basis checks |
|---|---:|---:|---:|---:|
| 668, width-four border | 98 | 10 | 64 | 6,272 |
| 668, width-twelve border | 98 | 10 | 64 | 6,272 |
| 716 | 45 | 3 | 64 | 2,880 |
| 764, first published input | 133 | 7 | 64 | 8,512 |
| 764, second published input | 133 | 7 | 64 | 8,512 |
| 892 | 48 | 6 | 64 | 3,072 |
| 1132 | 53 | 9 | 64 | 3,392 |
| 1964 | 3 | 3 | 64 | 192 |

The covered target orders are
**2676, 4676, 5892, 6876, 7924, 8684, 9308 and 11596**. They are targets
in the attributed preflight snapshot, not newly constructed Hadamard
orders. Counts in the table reuse a basis geometry across many fixed
sources; they are not claims of 39,104 separately computed field actions.

## Complete scalar maps, computed once per geometry

Each basis is a product of shuffled rational Walsh-four isometries, with
16 fixed recipes at each depth one through four. For each fixed frame F,
the exact map is

`B(z)=Q rho0(z) Q^T F`.

The output field basis is changed by `Q^T`; the input basis is that output
basis multiplied by F. Their weighted-trace metrics are respectively
`I_d` and `F^T F=d I_d`. Field parameters, faithful period-generator
definitions, the positive trace weight, every permutation, and a
nonproduct basis-column minor are included.

The primitive twelfth-root representation uses
`R=J(I-D)/2`, with `R^3=J` and `R^2-JR-I=0`. For each frame/basis,
all six representatives `B(zeta12^e)`, `e=0,...,5`, are computed exactly
once. Their six negatives supply the other roots. All source palettes
are then checked against this same complete twelve-root action map.
The code does not substitute a sixth-root-only test for a q=12 source.

Each scalar block has exact entry-category counts, its denominator and
hash, and its first failed sign entry. Every palette records a witness
for **every failed used scalar**, not only a single aggregate score.
The complete source exponent histogram gives its full implicit matrix's
exact sign-entry count. Original source positions for every root identify
the corresponding entry in the complete restricted matrix.

The original source Gram and orthogonality of the scalar action prove
the expanded ordinary Gram `n*d I` for every candidate. These candidates
fail the exact sign-entry condition. No approximate Gram or tolerance
check is used.

## Storage and duplicate scope

The driver loads one published frame at a time. For a fixed geometry it
shares `Q^T F`, streams one dense scalar block at a time, and retains only
statistics and hashes. It does not retain six order-1964 arrays or any
full expanded target matrix. The complete compact certificate is about
11.3 MB.

Literal geometry duplicate pruning found zero repeats. Source output
duplicate comparisons are conservative within each fixed literal
embedded source binding; they do not assert cross-binding or Galois
classification of every expanded output. The count 39,104 is the number
of executed source/basis checks.

## Replay and limits

```sh
PYTHONPATH=src:. ./.native-env/bin/python -m experiments.scalar_lift.broad_phase_basis_probe --check certificates/scalar-lift/broad-phase-basis-probe.json --progress
```

`SavedBroadBasisNode(path,geometry=0,node=0)` replays the fixed published
frame, faithful ambient field, all scalar blocks and their entry
statistics, every selected source's literal Gram, and every original
source/global-anchor embedding. Its `.entry(binding_id,i,j)` reconstructs
any full restricted-matrix entry while caching at most one scalar block.

Tests also check full known-order controls `BH(7,12)->H(28)` and
`BH(13,6)->H(52)`, and that a q=12 palette and a sextic palette can give
different answers in the same basis. Source alias deduplication retains
and replays every provenance recipe.

The remaining 91 preflight routes are explicitly saved as uncovered:
they require another field or a missing fixed frame width. The completed
finite basis recipes do not exhaust possible mixed unequal bases. No
new source or frame construction, SAT search, or automatic continuation
job is part of this run.

After this finite execution, the independent parity argument in
[the mixed-basis audit](https://github.com/sashakolpakov/complex-hadamard/blob/976cede671f40c810d593c8d5815e26e11a03ed6/docs/MIXED_RESTRICTION_BASIS_AUDIT.md) closes the whole
dyadic basis family at these `d=4r`, r odd >1 widths: a dyadic orthogonal
map between two order-d sign Hadamards must be a signed permutation,
which cannot have a primitive root minimal polynomial with an odd
factor. Thus deeper Walsh-only trees should be pruned without further
sampling. A useful replacement must introduce odd denominator factors
from r as well as mix coordinates globally. This later theorem is
separate from the replayed finite-run certificate.

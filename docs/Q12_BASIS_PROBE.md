# Executed mixed-basis restrictions for a genuine twelfth-root source

The [completed certificate](../certificates/scalar-lift/q12-basis-probe.json)
records **128 distinct unequal-basis scalar restrictions**, 64 for each
published real frame of orders 668 and 716. Every candidate finished;
none gives a sign Hadamard. No source-construction search or SAT search
was used.

The fixed input is the published `P7(1/12)` representative independently
verified by [the catalogue acquisition](CATALOGUE_SOURCE_ACQUISITION.md).
Its primitive twelfth-root cross ratio excludes reduction to a sixth-root
matrix by row and column phase gauges. We record an explicit global
exponent shift -1 modulo 12 before restriction, independently replay its
full Gram, and bind the resulting literal matrix to its source hash.
Its actual alphabet is

`{0,1,3,4,6,7,9,10} modulo12 = mu4 union zeta12*mu4`.

## All used scalar actions are checked

In the existing orthogonal quaternion field model let D be left
multiplication by `i+j+k`, and J right multiplication by i. The primitive
twelfth-root action is

`R = J(I-D)/2`, with `R^3=J` and `R^2-JR-I=0`.

It is not the sextic action `(I+D)/2`. The used scalar blocks are

`F, RF, JF, JRF` and their negatives,

after conjugating every scalar action by the chosen rational isometry Q.
Every candidate checks **all eight actual used roots**, using the four
opposite-root representatives to avoid redundant arithmetic. A successful
sixth-root partner alone would not certify this source.

Each Q is a product of shuffled Walsh-four layers. The executed batch has
16 fixed deterministic recipes at each depth 1, 2, 3, and 4 per published
frame. Exact layer numerators and local root numerators define

`B(z)=Q rho0(z) Q^T F`.

The output basis is the recorded weighted-trace orthonormal field basis
changed by `Q^T`; the input basis is that output basis multiplied by F.
Thus their metrics are `I_d` and `F^T F=d I_d`. The faithful auxiliary
period operator is conjugated by the same Q. The full used-root
orthogonality and the source Gram imply Gram `7d I` at every node; the
remaining test is its exact entry alphabet.

The certificate includes field parameters and proofs, the positive trace
weight, source exponents and hashes, all layer permutations, a nonproduct
basis-column minor, every used scalar's exact entry statistics and first
failure, and a hash of the complete used action map. Exact duplicate
pruning found zero repeated nodes. Entry counts for the full implicit
matrix weight each scalar-block census by its actual source frequency.

| Published frame width | Restricted source order | Nodes | Best exact sign-entry count | Total entries | All-sign endpoint |
|---:|---:|---:|---:|---:|---|
| 668 | 4676 | 64 | 8,542,658 | 21,864,976 | No |
| 716 | 5012 | 64 | 9,819,090 | 25,120,144 | No |

Order 4676 appears in the attributed unresolved-order snapshot. Order
5012 does not; it is retained as a separate arithmetic case, not a new
record target. These counts are exact diagnostics, not an assertion that
entrywise progress implies a nearby solution.

## Positive and negative controls

At width four, negate the first row of the ordinary Sylvester Hadamard
F4. Then `F4, R F4, J F4, J R F4` are all sign matrices with Gram
`4 I_4`. Applying all 49 source blocks gives an actual sign matrix of
order 28, with full integer Gram `28 I_28`. This verifies the complete
genuine-twelfth-root pipeline at a known order.

Without that row sign change, F4 and JF4 pass but RF4 and JRF4 fail.
The tests therefore check that success for only the fourth-root subgroup
does not hide failures at the other used roots. An independent dense
rational calculation at width eight verifies every local power and
two-layer conjugation, including actions needed by other q=12 palettes.

## Replay and reuse

```sh
PYTHONPATH=src:. ./.native-env/bin/python -m experiments.scalar_lift.q12_basis_probe --check certificates/scalar-lift/q12-basis-probe.json
```

`SavedQ12BasisNode(path, case=0, candidate=0)` reconstructs the source,
frame, faithful field, every used scalar block, and complete exact entry
statistics. Its `.entry(i,j)` reads any entry of the full order-4676 or
order-5012 matrix without allocating that full matrix.

The helpers `scalar_blocks(frame,layers,used_alphabet)` and
`candidate_record(source,frame_hash,frame,layers,recipe)` accept any fixed
twelfth-root alphabet. They compute additional representatives such as
exponents two and five when needed; they do not reuse a sextic-only
criterion or assume all q=12 inputs use this particular eight-root palette.
The completed batch covers its explicitly listed bases, not all unequal
input/output bases.

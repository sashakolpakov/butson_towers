# Global mixed two-basis restriction probe

The [executed certificate](../certificates/scalar-lift/alpoge-global-basis-probe.json)
contains **128 completed, distinct partner matrices**, using the verified
published order-668 and order-716 Hadamard matrices as fixed input frames.
There are no pending candidates. None is a sign partner, so no new real
Hadamard order has been constructed.

Unlike the same-basis conjugations, these coordinate pairs have the correct
expanded Euclidean Gram scale. Write F for the published frame, with
`FF^T=dI`, and R0 for the rational sextic action `(I+L_(i+j+k))/2`. For an
orthogonal rational Q, the output basis is the weighted-trace orthonormal
field basis times `Q^T`, while the input basis is the output basis times F.
The actual restricted scalar map is

`B(z)=Q rho0(z) Q^T F`, with `B(1)=F`, `B(zeta6)=G`.

Consequently `GG^T=dI` and `GF^T+FG^T=dI` hold exactly. Applied to the fixed
order-seven or order-thirteen sextic Butson sources, this map has full Gram
`ndI`. Its unresolved condition is that G have only ±1 entries. The saved
source hashes, full source proofs, faithful period-field definitions and
weighted trace basis are part of each case, rather than an inference from
the expanded order.

The two fields are the real degree-167/conductor-294923 and
degree-179/conductor-268501 period fields, each composed with Q(zeta12).
The positive trace weight and the local field basis are imported from the
independently checked [two-basis model](../experiments/scalar_lift/alpoge_mixed_basis_probe.py).
All auxiliary field operators undergo the same conjugation by Q.

The bases tested here are global isometries rather than individual low-rank
reflections. For each depth from one to four, sixteen deterministic recipes
are tested per frame. A layer partitions **all** coordinates into shuffled
quartets and applies Walsh4/2 on each quartet. Each layer is a rational
isometry. Their product Q is evaluated by quartet row sums, without cubic
matrix multiplication. Every candidate stores its actual coordinate
permutations and an explicit nonzero tensor minor in a column of the output
basis, proving that the cyclotomic and auxiliary coordinates are mixed.

Every candidate also clears the earlier necessary rank cutoff. The exact
binary rank of its integer Q numerator proves

`rank_Q(Q-I) >= rank_F2(Qnumerator)`,

because the denominator is a positive power of two. For the order-668 case,
the smallest certified lower bound is **163**, exceeding the required 84.
For order 716 it is **175**, exceeding the required 89. Thus the earlier
low-rank obstruction does not explain these failed candidates.

| Published frame | Candidates | Best sign entries in G | Total entries |
|---:|---:|---:|---:|
| 668 | 64 | 55,407 | 446,224 |
| 716 | 64 | 63,854 | 512,656 |

Deeper random layer products reduced the number of sign entries in this
batch. This is an observation about the displayed recipes, not a theorem
excluding all global isometries or unequal basis pairs. Every candidate has
an exact full entry census, a literal partner hash and a failed-entry witness.
Duplicate pruning uses literal equality of the canonical partner numerator
and denominator; no duplicates occurred.

The implementation is [alpoge_global_basis_probe.py](../experiments/scalar_lift/alpoge_global_basis_probe.py).
`SavedGlobalBasisNode(certificate_path,case=0,candidate=0)` replays a saved
partner and its fixed field/frame data.
`.entry(source_id,row,column)` accesses any entry of the complete restricted
matrix of order 4676, 8684 or 9308, without allocating that large matrix.
The six scalar blocks are `F,G,G-F,-F,-G,F-G`.

Independent small rational matrix tests check the complete product of
layers, partner conjugation, both Gram identities and the rank lower bound.
A width-four positive control confirms that a valid sign partner is accepted.

Replay:

```sh
PYTHONPATH=src:. ./.native-env/bin/python -m experiments.scalar_lift.alpoge_global_basis_probe --check certificates/scalar-lift/alpoge-global-basis-probe.json
```

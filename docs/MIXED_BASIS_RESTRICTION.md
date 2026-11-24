# Executed restrictions in mixed bases

The [completed finite batch](../certificates/scalar-lift/mixed-basis-restriction.json)
tests actual bases that mix the period and cyclotomic coordinates. There are
90 node visits, 70 distinct action/metric states and 20 exact duplicates,
with no pending recipes. It covers the order-84 control and the candidate
order 9308. Neither produces a sign matrix.

The control uses the fixed, published order-seven twelfth-root source and
the real cubic period field of conductor 31. The large case uses the
verified order-thirteen sixth-root source and the real degree-179 period
field of conductor 268501. Full source Grams and the faithful auxiliary
field definitions are replayed; no new source is constructed.

In period-first coordinates, a scalar basis has width `d=4r`. We test
`S=product(I+t E_pq)`, with one to three shears and coefficients ±1 or ±2.
Each shear changes both the period index and cyclotomic index. Every retained
nontrivial basis has an explicit nonzero 2-by-2 tensor minor in a basis
vector; it is not just a permutation of a product basis. The actual arrays
are

`A_z' = S^-1 A_z S`, and `G' = S^T G S`.

The same basis is used in the domain and codomain. Column and row updates
compute these arrays exactly, with integer overflow bounds. Inverse shears
replay every scalar and metric array. Literal equality of all used scalar
actions and the trace metric is the only duplicate test. The inverse-shear
return paths supply the twenty pruned duplicates.

Every retained node has a full matrix entry census, its source hash, basis
recipe and actual action/metric hash. The source matrix is stored implicitly
in its scalar blocks. For example, the order-9308 baseline has 86,487,788
zeros. A tested three-shear mixed basis reduces that number to 86,486,360,
and introduces entries ±2 and ±3. Thus this is an actual change of matrix
entries, while still failing the sign condition.

The implementation is [mixed_basis_restriction.py](../experiments/scalar_lift/mixed_basis_restriction.py).
`MixedRestriction(field,source,basis_shears).entry(i,j)` accesses any literal
matrix entry. `load_saved_node(certificate_path,case=1,node=-1)` reconstructs
and verifies a saved matrix directly, including its source, field, complete
entry census, scalar arrays and metric. No dense order-9308 matrix is allocated. The independent
order-84 check uses dense rational matrices: it checks every scalar
conjugation, the metric congruence, the complete weighted Gram identity and
the determinant of the complete matrix.

There is a decisive invariant beyond this finite batch. Same-basis scalar
restriction of `HH*=nI` to total order `N=nd` has

`det(A)^2=n^N`,

because its weighted Gram identity has scale n. Rational basis conjugation
does not change this determinant. A real Hadamard matrix of order N has
determinant squared `N^N`, so d>1 rules out a real Hadamard endpoint in every
same-basis restriction, including arbitrary mixed bases. For these odd-order
sources, the determinant is odd, whereas a sign matrix of order N>1 has
determinant divisible by `2^(N-1)`. That also excludes a literal all-sign
endpoint with a nonstandard metric.

This conclusion does not exclude different bases in the domain and codomain.
Such a coordinate pair transforms the restricted map as `T^-1 A S` and
changes its determinant by `det(S)/det(T)`. These unequal basis pairs are
outside this executed batch. The determinant argument must not be presented
as excluding every possible scalar-restriction formulation.

Replay:

```sh
PYTHONPATH=src:. ./.native-env/bin/python -m experiments.scalar_lift.mixed_basis_restriction --check certificates/scalar-lift/mixed-basis-restriction.json
```

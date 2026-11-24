# Completed character sources for forward scalar descent

The implementation supplies **actual known Butson matrices** of orders
5,329 and 58,081. It does not start from a missing source or an unfinished
corner. Neither source has been descended here to a real Hadamard matrix.

The primary construction is [Ong and Schmidt, Corollary 4.4, submitted
19 September 2026](https://arxiv.org/pdf/2609.22754). For a prime power
`q`, choose two nontrivial multiplicative characters `phi1, phi2`; their
orders and 6 must divide the output phase `h`. There is no character-parity
hypothesis. With any `gamma` in `mu_h`, define

\[
f(0,0)=1,\quad f(x,0)=\zeta_6\ (x\ne0),\quad
f(0,y)=\zeta_6^{-1}\ (y\ne0),\quad
f(x,y)=\gamma\phi_1(x)\phi_2(y)\ (xy\ne0).
\]

Then `H[row,column]=f(column-row)` is a `BH(q²,h)`.

The [generator](../experiments/scalar_lift/forward_character_sources.py)
implements prime fields and stores only the `q²` coefficient exponents.
`character_source(73)` uses character orders 2 and 8, giving
`BH(5329,24)`; `character_source(241)` uses orders 2 and 40, giving
`BH(58081,120)`. The `.entry(row,column)` method evaluates any entry without
materializing the large square matrix. A primitive cross-ratio certifies
the full phase in each case, even after root row/column gauges.

## Exact verification of the large inputs

Write `delta` for the identity of the additive group ring, `J` for its
all-ones element, and `A,B` for the two character arrays with value zero
at zero. All `q` periodic correlations are checked in an integral
cyclotomic basis, verifying

\[
AJ=BJ=0,\qquad AA^*=BB^*=q\delta-J.
\]

The two-dimensional coefficient array is

\[
F=\zeta_6J\otimes\delta+
\zeta_6^{-1}\delta\otimes J+\gamma A\otimes B.
\]

Its mixed terms vanish. In the tensor basis
`delta⊗delta, J⊗delta, delta⊗J, J⊗J`, the axis Gram has coefficients
`(0,q,q,-1)` and the interior Gram `(q²,-q,-q,1)`. Their sum is
`(q²,0,0,0)`. This proves **every** entry of both large Grams.

Independent controls check six literal `q²`-term correlations per large
source, and all correlations plus the dense 25×25 Gram for the order-25
control. The
[certificate](../certificates/scalar-lift/forward-character-sources.json)
records hashes and distinguishes these literal checks from the full
group-ring proof.

## Forward edges and their exact limits

The complete exponent-reduction theorem of
[Ó Catháin and Swartz](https://arxiv.org/pdf/1904.10771) supplies

\[
\mathrm{BH}(5329,24)\to\mathrm{BH}(10658,12)
\to\mathrm{BH}(21316,6),
\]
\[
\mathrm{BH}(58081,120)\to\mathrm{BH}(116162,60)
\to\mathrm{BH}(232324,30).
\]

These are theorem-backed next states; this generator does not materialize
the descended matrices. Removing the remaining odd prime factors requires
an applicable additional construction.

There is an immediate obstruction to a **uniform rational-linear final
block map on the raw alphabet**. The entries `1, zeta6, zeta6^-1` occur
for every `gamma`, and satisfy `zeta6+zeta6^-1-1=0`. If their image blocks
were all ±1, each coordinate would give a vanishing signed sum of three
signs, which is impossible. This also applies to compositions of uniform
linear scalar descents and to global root rotations of the source. It
does not exclude arbitrary source-dependent changes of bases.

The current two-cubic-coset partial frames also fail after **arbitrary
root row/column gauges** for these particular character sources. This
stronger, separately scoped gate is tiny: restrict the literal source to
coordinate grids of shapes 4×2 and 2×4. Normalize the cubic colors of the
two axes to 1 and 2. A coincident point has color 0 and an interior
difference has color `g`, the cubic color of `gamma`.

After fixing the first row gauge to zero, enumerate the `3^7=2187` row
gauges on each eight-point grid. A column admits a gauge omitting cubic
color 2 exactly when its adjusted entries do not contain all three colors.
The exact surviving row-gauge counts are:

| Coordinate grid | g=0 | g=1 | g=2 |
|---|---:|---:|---:|
| 4×2 | 0 | 0 | 2 |
| 2×4 | 0 | 2 | 0 |

Every `g` is excluded by a submatrix; the positive entries provide
nonvacuous controls. Literal submatrices of both large sources verify the
color templates. This excludes the current two-cubic-coset alphabet,
not all partial flat frames or all source-dependent descent methods.

```sh
PYTHONPATH=src:. python -m experiments.scalar_lift.forward_character_sources \
  --check certificates/scalar-lift/forward-character-sources.json
PYTHONPATH=src:. python -m pytest -q tests/test_forward_character_sources.py
```

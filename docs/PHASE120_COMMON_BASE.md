# A common-base phase-120 sign lift of width 32

There is an exact rational representation of `Q(zeta_120)` in dimension 32,
with a Hadamard frame that makes **64 of its 120 phase blocks sign matrices**.
This supplies the scalar action and frame for a conditional route

\[
BH(1499,120)\text{ in the specified alphabet}\longrightarrow H(47968).
\]

The required order-1499 source has not been constructed. The complete controls
below produce the already known real orders 480 and 1120; they establish the
composition and its literal sign entries, not a new Hadamard order.

The implementation is
[phase120_common_base.py](../experiments/scalar_lift/phase120_common_base.py),
with the exact [certificate](../certificates/scalar-lift/phase120-common-base.json)
and [tests](../tests/test_phase120_common_base.py). The construction combines the
[cubic partial frame](DYADIC_PARTIAL_LIFT.md) and
[fifth-root action](FIFTH_DYADIC_TOWER.md) over a shared field.

## The shared scalar field

Let `R` be the four-dimensional negacyclic signed shift with superdiagonal
entries 1 and bottom-left entry -1. Thus `R^4=-I`, and it represents
`zeta_8`. Write

\[
A=I+R+R^3,\qquad AA^T=3I,
\qquad
C=-\frac12\begin{pmatrix}I&A\\-A^T&I\end{pmatrix}.
\]

Then `C` is rational orthogonal, `C^2=C^T`, and `I+C+C^2=0`.
It is a rank-two action over `K=Q(zeta_8)`. Put
`R8=diag(R,R)` and `J=R8^2`, so `J^2=-I` and `J` commutes with `C`.

Let `E_b`, for `1 <= b <= 4`, be the four stored quaternary exponent
tables `FIFTH_POWER_QUATERNARY_BLOCKS`. Replacing their entries `i^e` by
`J^e` gives

\[
F^b=\frac12[J^{(E_b)_{kl}}]_{k,l=0}^3.
\]

These are the exact powers of an orthogonal action of order five. On the
32 coordinates, the cubic and eighth-root actions are four diagonal copies
of `C` and `R8`. All three actions commute. Their combination

\[
S=R_8^7 C F^4
\]

has restrictions

\[
S^{15}=R_8,\qquad S^{40}=C,\qquad S^{24}=F.
\]

Here and below these symbols denote their 32-dimensional embeddings.
The exact polynomial check gives

\[
\Phi_{120}(S)=0,\qquad
\Phi_{120}(x)=x^{32}+x^{28}-x^{20}-x^{16}-x^{12}+x^4+1.
\]

Because this polynomial is irreducible of degree 32, it is both the minimal
and characteristic polynomial. The action is faithful and has exact order
120. Its metric is the identity throughout.

This is a tensor product of ranks two and four **over `K`**; its rational
dimension is `2*4*[K:Q]=32`. It shares the eighth-root action instead of
duplicating the scalar coordinates. An ordinary tensor of the separate
width-eight cubic and fifth constructions would have width 64.

## The compatible frame is explicit

Let `P8` be the existing cubic sign frame. It satisfies

\[
P_8P_8^T=8I,\qquad P_8,\ CP_8\in\{\pm1\}^{8\times8}.
\]

Set

\[
Q=2\operatorname{diag}(P_8,P_8,P_8,P_8),\qquad P=FQ.
\]

Both frames have Gram `32I`. The raw frame `Q` is sparse, while every
eight-by-eight block of `P` is a signed permutation `J^e` applied to `P8`.
Thus `P` consists entirely of signs. For `a=0,1` and `b=0,1,2,3`,

\[
C^aF^bP=2F^{b+1}\operatorname{diag}(C^aP_8,C^aP_8,C^aP_8,C^aP_8)
\]

is also a sign matrix. The eighth-root action is a signed permutation and
preserves this property. In terms of the single primitive generator `S`,
the admitted exponents are exactly

\[
\mathcal L_P=\{e\bmod120:e\bmod3\in\{0,1\},\ e\bmod5\ne1\}.
\]

There are 64 such exponents. For the raw frame they are

\[
\mathcal L_Q=\{e\bmod120:e\bmod3\in\{0,1\},\ e\bmod5\ne0\}.
\]

The relation `P= S^24 Q` shifts a source in `L_Q` to a source in `L_P`
by subtracting 24 from every exponent. The certificate checks every phase,
including the excluded phases, so neither alphabet is merely sufficient
for this fixed action and frame: it is the complete flat alphabet.

For each frame the certificate directly checks all 120 diagonal Grams
and all 14,400 identities

\[
(S^eP)(S^fP)^T=32S^{e-f}.
\]

Consequently a restricted `BH(n,120)` becomes an actual sign matrix of
order `32n`, with Gram `32n I`, by replacing each exponent `e` with `S^eP`.
This is an exact source-to-sign morphism on the stated alphabet. Arbitrary
phase-120 sources are not covered.

## Two complete source-to-sign controls

The first control is a genuine phase-120 `BH(15,120)`. Index its rows by
`(i,a)` and columns by `(j,b)`, where `0<=i,j<3` and `0<=a,b<5`, and put

\[
e_{ia,jb}=40\,[i=j]+15aj+24\bigl(3(a+b)^2+1\bmod5\bigr)\pmod{120}.
\]

This is a Diţă construction from the cubic diagonal-omega source of order
three and a chirped fifth-root Fourier source, with an eighth-root diagonal
twist. Its entries lie in `L_Q`. The cross-ratio on rows and columns 0 and 6
has exponent 119, which is coprime to 120; hence row and column phase gauges
cannot move this control to a proper root subgroup. The complete cyclotomic
Gram and the entire integer Gram of its **480-by-480 sign lift** are checked.

The second control retains a nonzero coupled quadratic extension from
[the order-seven source](https://github.com/sashakolpakov/complex-hadamard/blob/976cede671f40c810d593c8d5815e26e11a03ed6/docs/COUPLED_SEVEN_RELATIVE_EXTENSION.md). Let `E7` be
its phase-24 branch-zero exponent matrix, and let
`f_ab=3(a+b)^2+1 mod5`. Then

\[
e_{ia,jb}=25(E_7)_{ij}+24(f_{ab}+2)\pmod{120}
\]

defines `BH(35,120)` in `L_P`. The embedding
`zeta_24 -> zeta_120^25` preserves its two admitted cubic colors. The fifth
term has residues `{0,2,4}` modulo five and avoids the excluded residue 1.
The inherited quadratic layers have `AB* != 0`. This source's full
cyclotomic Gram and its entire **1120-by-1120 integer sign Gram** are checked.

These controls prove that smaller sources, genuine cyclotomic extension,
and the common-base sign lift compose. They do not change a composite source
order into 1499. The high-order task remains a joint construction of an
order-1499 source within `L_P` or, equivalently, `L_Q`.

The [source residue engine](https://github.com/sashakolpakov/complex-hadamard/blob/976cede671f40c810d593c8d5815e26e11a03ed6/docs/PHASE120_SOURCE_RESIDUE.md) now implements
the exact quotient modulo two as two copies of `F16[pi]/(pi^4)`, with
conjugation exchanging them. It retains the same admitted alphabet and
checks the complete smaller sources in both components. Passing this
quotient check still leaves the final signs and characteristic-zero Gram
to be solved.

## Replay

```sh
PYTHONPATH=src:. python experiments/scalar_lift/phase120_common_base.py \
  --verify certificates/scalar-lift/phase120-common-base.json
PYTHONPATH=src:. python -m pytest -q tests/test_phase120_common_base.py
```

All matrix products in this certificate use exact signed 64-bit integers;
the stated matrix sizes and coefficients stay far below overflow. No
numerical near-orthogonality or prospective source completion is treated as
a matrix certificate.

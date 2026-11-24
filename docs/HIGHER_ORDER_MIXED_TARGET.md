# A strict mixed-six target at higher orders

The order-generic numerical learning experiment prescribes an entire row
comparison, not merely a six-term subsum. Its terms are a common rotation of

\[
M=(1,\zeta_{30}^{29},\zeta_{30}^{23},\zeta_{30}^{17},
   \zeta_{30}^{11},\zeta_{30}^{10}),
\]

followed by any number of antipodal pairs. The sum of `M` is zero: it is a
rotation of the familiar mixed relation

\[
\zeta_6+\zeta_6^{-1}+\zeta_5+\zeta_5^2+\zeta_5^3+\zeta_5^4=0.
\]

This strict target excludes pair/triangle decompositions at **every even
order at least eight**. It therefore retains the affine-family exclusion of
the order-eight experiment. Mere presence of a mixed-six subsum does not
have that implication at higher orders. Neither exclusion asserts that a
completion exists or that any resulting matrix is new in the literature.

## Coset obstruction

Write `mu_6` for the group of sixth roots of unity. Given a finite indexed
multiset `T` of unit complex numbers, partition its terms by their cosets in
`U(1)/mu_6`. For each coset `C`, define its subtotal

\[
\sigma_C(T)=\sum_{t\in T\cap C}t,
\]

counting multiplicities. If `T` is partitionable into rotated antipodal pairs
and equilateral triangles, then every subtotal is zero: all terms in each
such block belong to the same coset, and their sum is zero.

Adding an antipodal pair anywhere changes no subtotal. The two terms are
in the same coset because their ratio is `-1` in `mu_6`. The same is true
when adding a rotated equilateral triangle. These statements allow arbitrary
unit complex rotations, without a common root-of-unity alphabet.

For terms in `mu_30`, cosets modulo `mu_6` are indexed by exponent modulo
five. The exponents of `M` have residues

\[
(0,4,3,2,1,0)\pmod5.
\]

Four cosets consequently contain a single nonzero term. The remaining
coset has subtotal

\[
1+\zeta_{30}^{10}=1+\zeta_3=-\zeta_3^2\ne0.
\]

Thus all five nonempty coset subtotals are nonzero. Appending any number of
antipodal pairs leaves those subtotals unchanged, even if a pair overlaps
an existing coset or repeats an existing term. The resulting comparison
cannot have a pair/triangle partition. A common rotation merely translates
the cosets and multiplies their subtotals, preserving the obstruction.

The proof also applies to a tail made of any prescribed disjoint union of
rotated pairs and triangles. It applies to conjugates and Galois conjugates
of the mixed atom: multiplication of its exponents by a unit modulo 30
permutes the four singleton cosets, while the remaining two-term subtotal
is still nonzero.

## Consequence for affine Hadamard families

Consider an affine phase family

\[
H_{rj}(t)=B_{rj}\exp(i\langle v_{rj},t\rangle),
\]

orthogonal on an open real parameter domain, whose fixed seed `B` has entries
in `mu_(2^a 3^b)`. Real phase vectors are permitted; integrality is not
required. This is the same meaning of affine family as in the
[order-eight target proof](../experiments/grammar_vae/MIXED_ATOM_TARGET.md).

For a row comparison, group the terms in its identically vanishing inner
product by their frequency vectors `v_rj-v_sj`. Linear independence of
distinct exponential characters on an open domain implies that the seed
coefficients in each frequency group sum to zero. A positive vanishing sum
of roots whose orders have no prime divisors other than two and three
decomposes into antipodal pairs and equilateral triangles. This is the
two-prime case of Lam and Leung's minimal-relation structure theorem
([Corollary 3.4](https://arxiv.org/pdf/math/9511209)). One obtains an indexed
partition by repeatedly removing a minimal nonempty vanishing subsum.

Within one frequency group, every term acquires the same rotation as `t`
varies. Its seed pair/triangle partition therefore persists at every
parameter value. Repeating the argument for columns gives such a partition
for every row and column comparison throughout the family.

The property that every row and column comparison has a pair/triangle
partition is invariant under row and column phases, permutations,
transposition, and complex conjugation. Phases either cancel from a
comparison or give a common rotation; permutations reorder terms or
comparisons, while transposition exchanges row and column comparisons.

Consequently, a completed Hadamard matrix with the **strict mixed-six plus
pair-tail comparison** cannot be equivalent, under those operations, to
any point of an affine family based at a `BH(n,2^a 3^b)` seed. This is an
exclusion of a particular class of families, not a classification or a
claim of publication novelty. Nonlinear families and affine families with
other seed alphabets remain outside the exclusion.

## Why the whole comparison matters

At order twelve, let `F=F_6` and let `D` be an arbitrary diagonal unitary
matrix. Then

\[
K(D)=\begin{pmatrix}F&DF\\F&-DF\end{pmatrix}
\]

is a Hadamard matrix in an affine family based at `F_2` tensor `F_6`, a
sixth-root seed. Comparing the first columns of the two block columns
produces

\[
(a_1,\ldots,a_6,-a_1,\ldots,-a_6),
\qquad a_r=\overline{D_{rr}}.
\]

Choose the six `a_r` to be the terms of `M`. A mixed-six subsum is present,
but the other six terms are its negative, not a prescribed pair tail.
The complete comparison has six antipodal pairs. Its coset subtotals all
vanish because the second mixed block cancels the first block's subtotals.
This supplies a concrete counterexample to extending an argument based on
mixed-six *presence alone*. The strict target used by the learning problem
does not admit this counterexample.

## Numerical experiment scope

The generic problem preserves the strict row parametrization by construction.
Floating-point trigonometric evaluation approximates that exact structural
description, while the remaining orthogonality equations are only checked
numerically. A small Gram residual is therefore not an exact Hadamard
certificate. Prescribed pair-angle counts also do not certify the dimension
of a feasible family; a completed numerical matrix may lie on a lower
dimensional or singular set of those conditions.

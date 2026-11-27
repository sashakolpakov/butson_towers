# Explicit partial sign lifts at phases 24, 48, and 96

The rational actions in `JOINT_DYADIC_TOWER.md` admit explicit flat
frames. This closes the local flattening problem for the following restricted
source alphabet, while leaving the high-order source matrix unknown:

| source phase q | real block order d | admitted phases |
| --- | --- | --- |
| 24 | 8 | 16 of 24 |
| 48 | 16 | 32 of 48 |
| 96 | 32 | 64 of 96 |

The admitted exponents are exactly j congruent to 0 or 1 modulo 3. All entries
of a source must lie in this alphabet, possibly after an equivalence which
preserves its Butson property. An arbitrary BH(n,q) is not asserted to admit
such an equivalence.

## Construction and proof

Put r=4,8,16 and d=2r, q=6r. Let R be the negacyclic permutation matrix of
size r and multiplicative order 2r, so R^r=-I. Set

    J = R^(r/4) + R^(3r/4),    A = I+J,
    U = (1/2) [[I,A],[-A^T,I]],    D = diag(R,R),    S = DU.

The exact identities J^T=-J, J^2=-2I, AA^T=3I give U orthogonal,
U^2-U+I=0, and U^3=-I. D is monomial and commutes with U. The scalar-action
certificate verifies that S has minimal polynomial Phi_q in dimension d.
Consequently substitution zeta_q -> S defines a faithful rational field
representation preserving conjugation and transpose.

For r=4, enumerate sign vectors v with first coordinate +1 and Uv a sign
vector. There are exactly eight; they are mutually orthogonal. The checked
matrix P8 with those columns therefore satisfies

    P8 P8^T = 8I,       P8 and U8 P8 have entries only +1 or -1.

This is a finite exact enumeration of 128 vectors, not a floating-point
search. For general r, group coordinates by residue modulo r/4, placing the
four upper and four lower coordinates of each residue together. U becomes
I_(r/4) tensor U8. In this row order take

    P = H_(r/4) tensor P8,

where H is a Sylvester Hadamard matrix. Then P and UP are sign matrices and
PP^T=dI. Undo the row permutation to return to the original action.

Since S^jP=D^j U^jP and D is monomial, all phases j=0,1 mod 3 are flat.
For j=2 mod 3, U^jP is, up to sign, UP-P. Every entry is even, so these
phases are never flat. Thus the alphabet is verified exactly, including its
excluded third.

For a source B=(zeta_q^e_ab) in this alphabet, replace each entry by S^e_ab P.
All output entries are signs. For any source rows a,b the corresponding block
inner product is

    sum_c S^e_ac P P^T S^(-e_bc)
      = d rho(sum_c zeta_q^(e_ac-e_bc)) = dn delta_ab I.

Hence the result is a real Hadamard matrix of order dn. Both the flat-block
Gram identities and every pairwise quotient identity are checked exactly.
An order-three Butson source with diagonal a cube root and all off-diagonal
entries 1 provides an independent nontrivial substitution control; it is not
a claim of a newly constructed order.

The catalogue BH(13,6) in `bh13_6_seed.json`, with its independently verified
unreal gauge, also embeds into each source alphabet. If e is its exponent,
replace e by re+2 modulo q: the global phase shift takes the two nonzero
modulo-three classes to 0 and 1. The actual sign outputs of orders 104, 208,
and 416 have their full integer Grams checked. Their binary hashes and all
embedded source exponents are saved in the certificate. These are validation
controls at known orders.

The proof extends unchanged to r any power of two at least 4. The saved
certificates cover the three dimensions needed for the present targets.

## The remaining source equations over the dyadic field

Let eta=zeta_q^3, a primitive dth root, and choose omega to be the primitive
cube root zeta_q^(td) with td=1 mod 3. The flat alphabet is precisely
mu_d union omega mu_d. Thus each allowed source has a unique decomposition

    B = A + omega C,

where A and C have disjoint supports covering every entry, and each nonzero
entry is a dth root of unity. Neither component is required to be Hadamard.
Since Q(zeta_q) is quadratic over Q(eta), and omega^2+omega+1=0,

    BB* = AA*+CC* - AC* + omega (CA* - AC*).

Linear independence of 1 and omega over the dyadic field gives the exact
equivalence

    BB*=mI  iff  CA*=AC* and AA*+CC*-AC*=mI.

These coupled equations retain the source phases and support together. They
give a concrete norm problem along the scalar tower, with the final sign
realization already supplied by P. For each embedded order-13 control the
certificate checks both equations independently in the exact coefficient
ring Z[x]/(x^(d/2)+1).

One tempting support choice is completely excluded beyond order 3. If C is
monomial, right multiplication by C* makes C=I and A Hermitian with zero
diagonal. The norm equation becomes A^2-A=(m-1)I. Set k=sqrt(4m-3).
If k is irrational, rationality of trace forces equal multiplicities for the
two eigenvalues (1+k)/2 and (1-k)/2, giving trace m/2 instead of zero. If k
is integral, it is odd and the positive-eigenvalue multiplicity is
m(k-1)/(2k). Integrality forces k to divide m; since 4m=k^2+3, k divides 3.
Thus m is 1 or 3. The source search must permit richer support for both
cosets. This excludes only a monomial component, not the general coupled
system.

## What remains at the higher frontier

This gives explicit conditional constructions:

* a restricted BH(1499,96) produces a real matrix of order 47968;
* a restricted BH(2999,48) produces one of order 47984;
* a restricted BH(2999,96) produces one of order 95968.

The source problem retains 16, 32, or 64 choices per entry, rather than
collapsing to the old sixth-root source. In particular the determinant
obstruction for BH(2999,6) does not itself exclude a source over phase 96.
No source at these orders has been constructed here. The construction is a
proved restricted family within the broader joint source/action/frame
incidence formulation, not a classification of all its solutions.

A single circulant source is excluded at these prime orders. Do Duc and
Schmidt prove that a cyclically developed BH(Z_(p^a),q) requires p to divide
q (Theorem 3.5). Here p is 1499 or 2999 and q is 48 or 96. This obstruction
applies to that source symmetry only; it does not rule out the coupled source
equations above. [Primary paper](https://personal.ntu.edu.sg/bernhard/Publications/pub/ButsonBilinear1.pdf)

## Comparison with complementary-sequence assembly

One must track monomial representation degree as well as sequence length.
Craigen's generic iteration on s disjoint quasisymmetric real sequences uses
degree 2^(s-1). An eight-component complementary sum cannot simply be put in
an 8-variable real orthogonal design: arbitrary circulants need additional
identities involving transpose for the cross terms to vanish. The applicable
generic bound would exceed the 32-fold budget. This does not rule out a
smaller representation for a special signed subgroup.

Sources: [Craigen, 1995, Theorems 3 and 5](https://www.academia.edu/52596913/Signed_groups_sequences_and_the_asymptotic_existence_of_Hadamard_matrices);
[Egan and Ó Catháin, 2026, Theorem 2.2](https://arxiv.org/html/2610.28394v1)
for the general scalar-morphism framework. The explicit action and frame
above are established directly by the displayed identities and certificates.

Replay:

    python -m experiments.scalar_lift.dyadic_partial_lift \
      --check certificates/scalar-lift/dyadic-partial-lift.json

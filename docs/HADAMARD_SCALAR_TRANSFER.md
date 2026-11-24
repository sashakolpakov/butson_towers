# What the Hadamard condition requires of a scalar transfer

The starting point is the target identity and its entry condition. A
coefficient field or finite-field basis is useful only if it transports
both. This note gives algebraic implications and transfers; it does not
rely on a search, a completion algorithm, or an approximate Gram matrix.

The [dense sign-lift note](DENSE_SIGN_LIFTS.md) treats
`F3^times=C2` as the target multiplicative sign group and gives exact
finite-field-norm and tenth-root transfers to actual real Hadamard
matrices. The reduction-to-F3 section below concerns the additive field
equations, which are a separate issue from using its unit group.

## Start with the real target

A real Hadamard matrix of order N has entries +/-1 and satisfies

    H H^T = N I.

After normalizing its first row and column to +1, write its remaining core
as C, of order N-1. Then

    C 1 = C^T 1 = -1,
    C C^T = N I - J.

Thus an algebraic construction of a signed core with these identities
immediately supplies H by adjoining its first row and column. The rank-one
term is forced by Hadamard; it is not an error to be removed numerically.
For N=4m, the incidence matrix A=(J+C)/2 satisfies

    A 1 = (2m-1) 1,
    A A^T = m I + (m-1) J.

These are the Hadamard-design parameters 2-(4m-1,2m-1,m-1).

## What a scalar map actually transports

Let K have involution a -> bar(a), and let

    rho: K -> M_d(R),     rho(bar(a)) = rho(a)^T

be a unital algebra homomorphism. From `H H*=wI` one gets, automatically,

    rho(H) rho(H)^T = w I_(dn).

For a fixed frame Z with `Z Z^T=cI_d`, the block matrix
`W=(rho(h_ij)Z)` consequently has `W W^T=cwI`. This is the entire
orthogonality argument. A separate entrywise identity must prove that the
blocks lie in {0,+/-1}, or in {+/-1} for a dense real Hadamard target.
The Gram identity alone does not provide flat entries.

For a full cyclotomic alphabet containing an odd prime p, an additive
scalar map cannot make all its blocks dense signs. Indeed,

    1 + zeta_p + ... + zeta_p^(p-1) = 0

would imply `sum_j rho(zeta_p^j)Z=0`. An odd number of +/-1 entries cannot
sum to zero at any position. This holds for every block size, auxiliary
field extension, and basis. It excludes this full-alphabet scalar-transfer
argument, not all constructions that use finite fields or odd roots.

## A dense transfer whose hypotheses follow from quaternary flatness

For H in BH(n,4), write `H=A+iB`. The signed integer matrices A and B
have disjoint supports that partition every entry. The Hadamard equation
is equivalent to

    A A^T + B B^T = nI,
    A B^T = B A^T.

Set X=A+B, Y=A-B. Both have only +/-1 entries. Then

    R = [[X,Y],[Y,-X]]

satisfies `R R^T=2nI` identically. This is the classical quaternary-to-real
transfer. It is coefficient extension/restriction with a compatible
involution and an integral flat frame; no extra unknown matrix is needed.

For a desired real order 4m it requires an actual BH(2m,4) input. Such an
input is not supplied merely by identifying the transfer. Constructing a
family of those inputs at unresolved parameters remains a mathematical
task. The argument does not establish a new real order by itself.

## What a primitive target order forces on small scalar blocks

There is a useful restriction for a target order 4p with p>3 prime. The
following argument concerns a **uniform norm-preserving scalar transfer**:
the same unital star representation rho and fixed frame Z are used at
every entry, `Z Z^T=dI`, and the output blocks are dense signs. The source
has algebraic entries a with `a bar(a)=1`; they need not initially be
roots of unity. Both input order n and block degree d are assumed greater
than one. This is not a restriction on all Hadamard constructions or all
uses of field extensions.

Each block `rho(a)Z` is itself a real Hadamard matrix of order d. Since
`dn=4p`, the only possibilities are d=2 and d=4. A degree-one map leaves
the original-order problem; d=4p with n=1 is a direct construction of the
entire target block and is not covered by the smaller-input reduction.

Choose a used source entry a0 and absorb rho(a0) into Z. The used scalar
ratios now include 1, and

    rho(a/a0) = (rho(a)Z)(rho(a0)Z)^T/d

is rational. Below, these ratios are shown to belong to a finite matrix
group. Faithfulness of rho implies that the multiplicative group of used
scalar ratios is finite as well. A finite subgroup of a characteristic-zero
field's multiplicative group is cyclic, so the ratios generate mu_k for
some k. Consequently Q(zeta_k) embeds in M_d(Q), and phi(k) divides d.
Root-of-unity phases and their degree bound are therefore consequences of
block flatness; neither is an extra ansatz here.

For d=2 the ratios of sign Hadamard blocks are signed permutation
matrices. Their group has order eight. Combining its element orders with
the field-degree condition leaves only k=1,2,4. The real cases would
require an impossible real Hadamard input of order 2p. The remaining
possibility is a quaternary input of order 2p.

For d=4 define G to consist of signed permutation matrices and normalized
real Hadamard matrices H4/2. This is a group. Products of two normalized
sign blocks have entries in {0,+/-1/2,+/-1}. In an orthogonal row, either
one entry is +/-1 or all four entries are +/-1/2. If one row has a unit
entry, orthogonality excludes every half-valued row, so the whole matrix
is a signed permutation. This proves closure, including the mixed cases.

There are 384 signed permutations and 768 normalized sign blocks, so
`|G|=1152=2^7*3^2`. To count the second set, normalize the first row and
column of a sign block. Its last three rows are the six permutations of
the three distinct allowable sign patterns; undoing normalization has
`2^4*2^3` choices.

The degree bound gives `k in {1,2,3,4,5,6,8,10,12}`. The group order
excludes k=5,10. A zero sum of 2-power roots has even length, so an odd
input order p excludes k=1,2,4,8. A zero sum of cube roots has length
divisible by three, excluding k=3 when p>3. Only k=6 and k=12 remain.

Their flat alphabets are also forced. Write A=rho(zeta_k), and let T=Z/2.
Every A^e T is in G. For k=6, `A+A^T=I`, so A and A^2 are half-valued
sign blocks, while `A^3=-I` is a signed permutation. For k=12, the same
reasoning applied to A^2 shows that A^2 is a half-valued sign block; A
cannot be a signed permutation. Also `A^6=-I` makes A^3 skew-symmetric,
so its zero diagonal forces it to be a signed permutation.

Consequently the block type repeats modulo three. At most one of
`T,AT,A^2T` is a signed permutation: two would make A or A^2 a signed
permutation. At least one is, because otherwise the identities

    A^2 T - A T + T = 0             (k=6),
    A^4 T - A^2 T + T = 0           (k=12)

would express an odd sum of half-sign entries as zero. In the second
case, A^4 T has the same type as AT because A^3 is a signed permutation.
Exactly one residue class modulo three is therefore excluded.

Up to a common source phase, the only smaller-input cases in this framework
are:

| Block degree | Source order | Source entries |
|---|---|---|
| 2 | 2p | mu_4 |
| 4 | p | mu_6 minus mu_2 |
| 4 | p | mu_12 minus mu_4 |

This is a restriction on scalar-transfer hypotheses, not a new existence
theorem. In particular it does not say that every real Hadamard matrix
must have one of these block structures. A direct larger-degree trace
construction or a construction using a core identity remains outside it.

### Allowing a separate orthonormal basis in each block row and column

The same input restriction survives these more flexible basis choices.
Suppose normalized output blocks have the form

    T_ij = U_i rho(h_ij) V_j,

where every U_i,V_j is orthogonal and each T_ij is a sign Hadamard block
divided by sqrt(d). Dephase H so its first row and column are 1, absorbing
the removed scalar phases into U_i,V_j. The four-block product

    T_00 T_i0^T T_ij T_0j^T = U_0 rho(h_ij) U_0^T

puts all source entries in one faithfully represented field. For d=2,
the left side is a signed permutation, since it is a product of two
ratios of order-two Hadamard blocks. For d=4 it belongs to G. The finite
group and field-degree arguments above apply again, giving the same
three possible coefficient fields.

For d=4 the restricted flat alphabets follow as well. Let P be the
signed permutation subgroup of G. It has index three. The left action
`pi:G -> Sym(G/P)` has the property

    X is a flat block  <=>  pi(X)(P) != P.

Put A=rho(zeta_k), using the common conjugated representation. For k=6,
`A^3=-I` and A is not in P, so pi(A) is a three-cycle. For k=12, A^2 is
not in P, again forcing pi(A) to be a three-cycle. Label the three cosets
by 0,1,2 so pi(A) adds 1. The original blocks can be written

    T_ij = R_i A^(e_ij) C_j,
    R_i=T_i0 T_00^T,   C_j=T_0j,

with R_i,C_j in G. If `alpha_i=pi(R_i)^(-1)(0)` and
`beta_j=pi(C_j)(0)`, flatness is exactly

    e_ij + beta_j - alpha_i != 0 mod 3.

Multiplying source row i by `zeta_k^(-alpha_i)` and column j by
`zeta_k^(beta_j)` therefore gives precisely the allowed sixth- or
twelfth-root alphabet in the table. Arbitrary orthonormal bases in
individual block rows and columns do not add new possible input orders
or fields here. In this argument the three labels encode a finite group
action; the integer Gram identity still comes from the characteristic-zero
star representation.

The scope remains a scalar transfer between smaller input and output
blocks. It allows row and column basis changes, but it does not constrain
an arbitrary orthogonal change of basis on the entire target space.

## A determinant obstruction at the first reported gap

Hadamard also forces an arithmetic condition on a proposed source. For
H in BH(n,k), its determinant delta is an algebraic integer satisfying

    delta bar(delta) = n^n.

Let ell be a rational prime not dividing k. If `ell^j=-1 mod k` for some
j, then complex conjugation fixes every prime ideal P above ell in
Q(zeta_k). Since ell is unramified, taking a P-adic valuation gives

    2 v_P(delta) = n v_ell(n).

The right side must be even. This is the standard self-conjugate-prime
obstruction; see the algebraic-number-theory preliminaries in
[Armario, Egan, Kharaghani and O Cathain](https://arxiv.org/html/2412.16579v2#S2.SS3).
Passing the condition is not an existence theorem.

For a prime source order p>3, it has the following immediate consequences:

| p modulo 12 | BH(p,6) | BH(p,12) |
|---|---|---|
| 1 | Not excluded by this norm test | Not excluded by this norm test |
| 5 | Impossible | Not excluded by this norm test |
| 7 | Not excluded by this norm test | Not excluded by this norm test |
| 11 | Impossible | Impossible |

In particular, `2396=4*599` and `599=-1 mod 12`. Thus neither BH(599,6)
nor BH(599,12) exists, even before imposing the restricted entry alphabets
of the preceding section. The numerical order is the smallest gap
reported in the
[October 2026 survey](https://arxiv.org/html/2610.27175v2); current existence
claims must still be checked against later constructions.

Combining the arguments, any scalar transfer of the specified kind to
order 2396 with n,d>1 must, up to source row and column phases, start with
BH(1198,4) and use d=2. This includes separate orthonormal bases in each
block row and column. Arbitrary algebraic source phases do not enlarge
the possibilities. This does not exclude a direct field-and-trace
construction of the entire target, a construction of its signed core, or
a transfer that does not have the specified block form.

## A Gaussian-base extension for the twelfth-root case

The third case has an explicit field argument. Put

    K=Q(i),  L=K(alpha)=Q(zeta_12),  alpha^2=-3.

Use complex conjugation as involution. A positive weighted relative trace
pairing is

    h(x,y) = Tr_(L/K)((1-i alpha/3) bar(x)y)/2.

Its Gram matrix in the basis (1,alpha) is `[[1,i],[-i,3]]`. The two vectors

    1,  (alpha-i)/(i-1)

are orthonormal. Multiplication by alpha therefore has matrix

    S = [[i,1+i],[-1+i,-i]],
    S*=-S,   S^2=-3I.

Set rho(i)=iI and rho(alpha)=S. This is a star embedding of L into M_2(K).
For zeta_12=i(1-alpha)/2 its image is `i(I-S)/2`.

Every allowed scalar can be written as

    z=b(sigma+alpha)/2,    b in mu_4, sigma in {+1,-1}.

The scalar multiple `(1+i)rho(z)` is, respectively,

    b [[i,i],[-1,1]],    or    b [[-1,i],[-1,-i]].

Both are quaternary sign blocks with Gram 2I. Hence an H in BH(n,12)
whose entries avoid mu_4 gives

    G=(1+i)rho(H) in BH(2n,4),

and the ordinary quaternary transfer gives a real Hadamard matrix of
order 4n. Orthogonality is transported by the field representation and
the scalar norm `|1+i|^2=2`; flatness is the displayed two-case identity.
There is no unknown auxiliary matrix to solve for.

The same domain of eight admissible powers appears in Example 15 of
[Egan and O Cathain, *Morphisms of Butson classes*](https://arxiv.org/html/1707.08815v3#S3.SS2).
No novelty or new-order claim is made for this transfer.

Equivalently, write the input as `H=(A+alpha B)/2`. Its coefficients are
quaternary matrices with `a_ij=+/-b_ij`. Independence of 1 and alpha over
Q(i) makes the Hadamard condition exactly

    A A* + 3 B B* = 4nI,
    A B* = B A*.

Neither A nor B is assumed individually Hadamard. These two identities
and the entrywise Gaussian phase constraint are an algebraic input
construction problem. Scalar transfer alone does not establish such an
input at a previously unknown real order.

## A genuine cubic scalar extension when zeros are allowed

There is also an exact extension of the cubic cyclotomic field that
provides ordinary Euclidean orthogonality and ternary integer entries.
Let

    K=Q(alpha),   alpha^2=-3,
    L=K(beta),    beta^2=-2.

This is a degree-two extension of K, contained in a cyclotomic field of
conductor 24. Complex conjugation negates alpha and beta and fixes their
product. On L take the positive Hermitian trace pairing

    <x,y> = Tr_(L/Q)(x bar(y))/4.

In the basis (1,alpha,beta,alpha beta), its Gram matrix is
`diag(1,3,2,6)`. The exact four-square identity

    x^2+3y^2+2z^2+6t^2
      = x^2+(y+z+t)^2+(y-z+t)^2+(y-2t)^2

provides a rational orthonormal basis. One choice is

    b0 = 1,
    b1 = -alpha/3 - beta/2 + alpha beta/6,
    b2 = -alpha/3 + beta/2 + alpha beta/6,
    b3 = -alpha/3 - alpha beta/3.

Multiplication by alpha and beta in this basis is represented by

    C =  0  1  1  1       D =  0  1 -1  0
        -1  0  1 -1           -1  0  0 -1
        -1 -1  0  1            1  0  0 -1
        -1  1 -1  0            0  1  1  0.

These matrices commute, are skew-symmetric, and satisfy

    C^2=-3I,     D^2=-2I.

Thus multiplication realizes the involution as transpose. With
`zeta_6=(1+alpha)/2`, the six multiplication matrices

    rho(beta zeta_6^e) = D ((I+C)/2)^e,    e=0,...,5,

all have entries in {0,+/-1} and Gram 2I. The six cases reduce to three
matrices and their negatives using `U^2=U-I`, and can be checked directly.

Now extend H from K to L, multiply it by the **single scalar beta**, and
restrict scalars using the displayed orthonormal trace basis. If
`H H*=wI`, then

    (beta H)(beta H)* = 2wI,
    W=rho(beta H) has W W^T=2wI.

For a sixth-root source, W is an exact real W(4n,2w). This is a genuine
scalar-extension proof of the ternary weighing lift. No completion step
is part of this theorem. It guarantees zeros and signs, not dense signs;
the odd-prime obstruction explains why that distinction is unavoidable
for a complete cubic alphabet.

## What reduction to F3 implies

F3 retains three distinct symbols for -1,0,1, but reduction provides only
congruences. A real Hadamard matrix reduced modulo three satisfies

    H H^T = (N mod 3) I.

If 3 does not divide N, the rows form a basis giving an orthogonal
similitude. If 3 divides N, its row space is self-orthogonal and has
dimension at most N/2, so the rows cannot be a basis of F3^N. Neither case
promotes a modular Gram equation to an integer one.

There is a stronger exact rank statement when `v_3(N)=1`. If
`d_1 | ... | d_N` are the integer Smith factors of H, then `NH^(-1)=H^T`
implies

    d_i d_(N+1-i) = N.

Hence exactly N/2 factors are divisible by 3 in this case, so the reduced
row space is self-dual of dimension N/2. The same argument works at any
prime dividing N exactly once. This is a necessary implication of an
integer Hadamard matrix, not a sufficient modular construction.

For a concrete failure of sufficiency at an admissible real order, take

    M=J_16-2I_16,    M M^T=4I_16+12J_16.

Its entries are signs and its rows are an orthonormal basis over F3,
since `M M^T=I mod 3`. Its integer off-diagonal row inner products are
12, so M is not a real Hadamard matrix. In fact the same example works
at every N=1 mod 3 with N>4.

The coefficient field should therefore be selected from the desired norm,
involution, and flat-frame identities. A preference for F3 alone does not
supply these identities. For a dense universal scalar transfer, quaternary
flatness gives the direct working example above; cubic scalar extension
gives a weighing matrix with zeros.

## Exact algebra checks

These calculations were verified with rational and Gaussian-rational
arithmetic, without numerical optimization:

- Enumerating real order-four sign blocks gives 768; there are 384 signed
  permutations. Their union has element-order counts
  `{1:1, 2:139, 3:80, 4:228, 6:464, 8:144, 12:96}`.
- The action on the three signed-permutation cosets has image S3 and
  kernel of size 192. All 108 combinations of row action, column action
  and phase residue obey the displayed forbidden-residue condition.
- The weighted relative trace basis gives the stated matrix S exactly.
  The quaternary-flat powers of `i(I-S)/2`, after multiplication by 1+i,
  are precisely `{1,2,4,5,7,8,10,11}` modulo twelve.
- Modulo 599, `Phi_6=x^2-x+1` is irreducible and
  `Phi_12=(x^2+188x+1)(x^2-188x+1)`. The latter two factors are irreducible
  and self-reciprocal. These are explicit instances of conjugation fixing
  each prime above 599.
- The order-sixteen modular counterexample has the stated integer Gram
  matrix and is invertible over F3.

# Field and character transfers with dense sign outputs

The target alphabet is exactly the multiplicative group
`F3^times={1,-1}=C2`. Its elements are interpreted as integer signs.
Orthogonality means an integer sum of signs is zero: for any two distinct
rows, the entrywise quotients must contain each sign exactly N/2 times.
An additive equation over F3 alone does not give that integer balance.

This note considers transfers with actual +/-1 Hadamard outputs. It does
not supply a construction at a previously unknown real order.

The subsequent [inverse-lift calculation](https://github.com/sashakolpakov/complex-hadamard/blob/976cede671f40c810d593c8d5815e26e11a03ed6/docs/INVERSE_DENSE_LIFT.md) starts
with candidate target 4676, identifies four compatible stored order-seven
torsion points, and derives the precise missing order-668 block identity.
A verified order-13 point makes that lift reusable for candidate 8684,
and order-nine tripling reduces the block requirement for candidate 6876
to 764. The required order-668 rational operator is necessarily dense
with denominator 334. The larger lift pairs remain unconstructed.

The [cyclotomic tower audit](https://github.com/sashakolpakov/complex-hadamard/blob/976cede671f40c810d593c8d5815e26e11a03ed6/docs/CYCLOTOMIC_TOWERS.md) verifies iterated
transfers through Q(i), filters full cyclotomic coefficient fields using
target prime factors and trace-form obstructions, and constructs a trace
basis for a proper cyclotomic subfield of degree 167. A flat sign frame
for the larger target blocks is still required.

The [factorized trace tower](https://github.com/sashakolpakov/complex-hadamard/blob/976cede671f40c810d593c8d5815e26e11a03ed6/docs/FACTORIZED_TRACE_TOWER.md) carries metrics
without requiring an intermediate Butson alphabet. It uses a split
degree-167 field and a coupled positive trace weight; normalization is
allocated across combined stages. Its final sign frame is not constructed.

The useful established dense transfers include:

| Source | Required entry domain | Real output order |
|---|---|---|
| BH(n,4) | All fourth roots | 2n |
| BH(n,6) | Exclude +/-1 | 4n |
| BH(n,10) | Exclude +/-1 | 8n |
| BH(n,12) | Exclude the fourth roots | 4n |
| BH(n,2p), p>3 a Fermat prime | Exclude +/-1 | 2(p-1)n |

The last row is the cyclic mutually-unbiased-basis construction of Gow,
as used in Corollary 11 of
[Egan and O Cathain](https://arxiv.org/html/1707.08815v3#S3.SS1).
The tenth-root implementation below is its smallest noncubic instance.
These rows describe actual dense transfers, conditional on constructing
an input in the specified domain.

## Scalar extension through a sign character

If a group G has a nontrivial character `chi:G -> C2`, it induces a ring
homomorphism `Z[G] -> Z`, sending each group element to an integer sign
and respecting the involution `g -> g^(-1)`. Matrix norm identities can
therefore be specialized through this map.

For example, a generalized Hadamard matrix M over G has the formal Gram
identity

    M M* = nI + lambda (J-I) sum_(g in G) g,
    n=lambda |G|.

Since `sum_g chi(g)=0`, evaluating through chi gives a dense real Hadamard
matrix of order n. This is an exact scalar-specialization argument, with
the sign alphabet guaranteed by the character. Constructing the formal
input identity is essential: orthogonality under a single complex
cyclotomic character does not automatically imply it under the sign
character. No generalized-Hadamard existence assertion is assumed here.

## A finite-field extension whose norm takes values in the sign group

For E=F_(3^r), the norm

    N_(E/F3)(x)=x^((3^r-1)/2),    x != 0,

is a multiplicative character into F3^times. Interpreting it as +/-1 gives
the quadratic character chi of E. Its character identities, computed as
integer sums, are

    sum_x chi(x)=0,
    sum_x chi(x)chi(x+a)=-1       (a != 0).

The character is extended by chi(0)=0 for these identities. Set
`Q_xy=chi(y-x)`. Then `QQ^T=qI-J` for q=3^r. If r is odd, Q is skew and

    H = [[1, 1^T], [-1, I+Q]]

has only sign entries and `HH^T=(q+1)I`. If r is even, Q is symmetric and

    C = [[0,1^T],[1,Q]],     CC^T=qI,
    A=C+I, B=C-I,
    H = [[A,B],[B,-A]]

has only sign entries and `HH^T=2(q+1)I`. These are Paley constructions,
not new orders. The zeros in Q and C are removed by the displayed exact
identity before returning H; there is no weighing-matrix completion.

The replay uses F9=F3[u]/(u^2+1), where
`N(a+bu)=a^2+b^2`, and constructs a dense H20. This is a direct check of
the multiplicative interpretation of F3.

## A cyclotomic extension with dense output

There is also a dense partial transfer

    BH(n,10) with no entries +/-1
       -> BH(4n,4) -> BH(8n,2).

The first step is Example 12 of
[Egan and O Cathain, Morphisms of Butson classes](https://arxiv.org/html/1707.08815v3#S3.SS1).
The last step is Turyn's quaternary-to-real transfer. The matrices below
are an explicit, exactly checked realization; no novelty is claimed.

Let

    M = (1/2) [[-1, 1, i, i],
               [-1,-1, i,-i],
               [-1, 1,-i,-i],
               [-1,-1,-i, i]].

Direct algebra gives

    MM*=I,    I+M+M^2+M^3+M^4=0,
    every entry of 2M^k lies in {+/-1,+/-i}, k=1,2,3,4.

Thus `rho(zeta_5)=M` is a faithful star representation. A tenth root is
`zeta_10^e=(-1)^e zeta_5^(3e)`. When e is not 0 or 5 modulo 10, the block

    2(-1)^e M^(3e mod 5)

is dense and quaternary. For an admitted Hadamard input H, applying these
blocks gives `GG*=4nI` functorially. Orthogonality is a consequence of the
field representation; flatness is the four-power identity above.

### The field extension and trace basis

Put K=Q(t), where t is a primitive fifth root, and
`L=K(i)=Q(zeta_20)`. Extend the input matrix from K to L and restrict from
L to Q(i), using the positive relative Hermitian trace

    h(x,y)=Tr_(L/Q(i))(theta bar(x)y),
    theta=(3+(1+i)t^2+(1-i)t^3)/10.

In the power basis (1,t,t^2,t^3), the columns of the following matrix give
an orthonormal basis for h:

    P = [[1,-1,  -i,     i],
         [0,-2,-1-i,   1+i],
         [0,-1+i,-1-i,   2],
         [0,-1+i,   0, 1-i]].

Multiplication by t in this basis is exactly M. The weight theta is fixed
by conjugation and is totally positive; the replay verifies this by the
positive principal minors of its Hermitian multiplication matrix. The
scalar multiplier 2 has squared norm four, giving the first transfer's
Gram factor. No unknown auxiliary Hadamard matrix is assumed here.

### Restricting Gaussian scalars to signs

Use `rho(i)=J=[[0,-1],[1,0]]` and `Z=[[1,1],[1,-1]]`. All four blocks
`J^e Z`, e=0,1,2,3, have entries +/-1 and Gram 2I. Consequently every
quaternary G of order m gives a real H of order 2m with `HH^T=2mI`.
This supplies dense signs at every entry of the composed tenth-root lift.

### A fully explicit input and controls

The matrix

    H_ij=t^(3(i+j)^2+1),    i,j in F5,

is Hadamard: any distinct-row quotient is a nontrivial additive character
in j. Its exponents are in {1,3,4} modulo five, so none of its entries is
1. The tenth-root lift therefore constructs a dense real H40. Tensoring
the input with H2 supplies a genuinely signed tenth-root input and a dense
H80. These are known-order checks, not record claims.

Pure p-th-root inputs without 1 cannot be arbitrarily large. Summing all
p-1 Galois conjugates of such an input gives -J. Each conjugate has
operator norm sqrt(n), so

    n <= (p-1)sqrt(n),    hence n <= (p-1)^2.

For p=5 this excludes n>16. Large tenth-root inputs must use additional
sign structure; a large Fourier matrix and a phase gauge alone do not
supply them.

## Multiplicative projection directly to signs

For an odd prime p and a>=1, the map

    chi: mu_(2p^a) -> C2,    chi(z)=z^(p^a)

is a group retraction onto signs. It is not additive. It nevertheless has
the following useful sufficient condition:

    H in BH(N,2p^a), 4|N, N<4p
       => H^(entrywise p^a) is a real Hadamard matrix of order N.

To prove it, reduce cyclotomic integers modulo a prime above p. The p^a-th
roots reduce to 1 and the signs remain distinct. Thus the sign projection
S satisfies `SS^T=NI mod p`. The following elementary counting lemma
upgrades this to the required integer identity.

If a sign matrix of order N=4t is Hadamard modulo an odd integer m>t,
choose any three distinct rows and normalize the first to all +1. Let
a,b,c,d count the four sign patterns in the other two rows. The three
modular inner products imply

    a=b=c=d=t mod m.

Because 0<=t<m, each count is at least t. Their sum is 4t, so every count
equals t as an integer. Applying this to every triple proves real
orthogonality. The modular-Hadamard framework and pattern-counting
arguments are standard; see
[Kuperberg](https://arxiv.org/html/1409.0148v2).

The strict bound matters: `F_(2p) tensor H2` lies in BH(4p,2p), but its
sign projection has rank four, not 4p. These are exactly verified
counterexamples at p=3 and p=5.

The sufficient condition by itself does not establish a new order. In
that range, existence of the complex source is equivalent to existence
of the real target, since every real Hadamard is already a valid source.

## Relation to the order frontier

The tenth-root lift multiplies the source order by eight. In the saved
95-order historical snapshot, the only potential targets have source
orders 599,719,739,839,1031,1087. The first four and 1087 are ruled out by
the determinant norm in Q(zeta_5): a prime congruent to 2,3 or 4 modulo 5
is self-conjugate, so it cannot divide the norm to an odd exponent.

The remaining numerical target 8248=8*1031, along with 5912 and 8696, is
identified as already constructed in the
[subsequent primary frontier report](https://generalised.si/research/hadamard-matrices).
Thus this lift supplies no unsuppressed candidate in that snapshot after
the reported corrections. The snapshot is historical; this statement
does not claim a complete current novelty audit outside it.

## Replay

```sh
PYTHONPATH=src python scripts/check_dense_lifts.py --check certificates/scalar-lift/dense-lifts.json
PYTHONPATH=src python -m pytest -q tests/test_dense_lift.py
```

The [implementation](../src/hadamard/dense_lift.py) checks its source
Hadamard equations exactly and rejects real entries outside the partial
domain. The [report](../certificates/scalar-lift/dense-lifts.json) checks
the field identities, trace basis, dense entry alphabet, both integer
Gram matrices, and exact C2 quotient balance for the known-order controls.

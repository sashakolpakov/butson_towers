# Executed partial third/fifth scalar restrictions

The [actual tree](../certificates/scalar-lift/partial-fifth-descent.json)
contains **464 distinct literal Butson/sign matrices from90 fixed starting
inputs**. It evaluated463 successful scalar-edge attempts and pruned89
exact literal duplicates. There are32 partial-fifth edges,33 partial-third
edges, and384 candidate edges rejected by explicit non-root entries.
Every declared branch is finished, with zero pending nodes.

All matrices are saved under
[`artifacts/partial-fifth-descent`](../artifacts/partial-fifth-descent),
using5,371,543 bytes in the manifest's bundles. The phases actually reached
are **2,4,6,8,10,12,16,20,24,30,40,60,80,120**. The largest nonreal matrix
has order640. These are complete exponent matrices, not necessary-test
routes or missing completions.

The53 real endpoints have known orders
**28,36,40,104,312,480,640,936,1120,1280**. Every complete sign Gram was
independently checked using exact binary row distances. Newly admitted
nonreal matrices of order at most160 also received independent complete
cyclotomic Gram checks. All larger intermediate Grams follow by exact
scalar cross-Gram transport from the verified source. No new Hadamard
order is claimed.

## Fifth-root edge

For `q=5c`, with `4|c` and `5` not dividing c, the extension
`Q(zeta_q)/Q(zeta_c)` has degree4. The existing Gaussian4-by-4 matrix M
satisfies `Phi5(M)=0` and `MM*=I`; its four nonidentity powers have twice
their entries in `{1,-1,i,-i}`. All25 multiplication/adjoint pairs are
checked exactly.

Set `a=5^-1 mod c`, `b=c^-1 mod5`. The faithful field action is

`rho(zeta_q^e)=zeta_c^(a e) M^(b e mod5)`.

The cyclic vector basis identifies M with the ordinary companion matrix
of Phi5; its inverse explicitly records the output field basis. Choose
the input basis to be `2 zeta_q^s` times that output basis. The resulting
scalar blocks are

`B_e=2 rho(zeta_q^(e+s))`, `B_e B_f*=4 rho(zeta_q^(e-f))`.

They have only c-th-root entries exactly when `(e+s) mod5 !=0`. Otherwise
the block contains zeros and diagonal entries twice a root. Therefore a
source whose actual used exponents all satisfy this condition restricts
to an actual `BH(4n,c)`. The code verifies the full cyclotomic polynomial
in lower-field coefficient vectors as well as the Gaussian identities.

All five omitted-coset positions are tried as **input-basis choices**.
The starting matrix is unchanged. In particular, the raw basis for
`common120-35` fails, but shifts2 and4 succeed, yielding literal
`BH(140,24)` matrices. Its raw failure remains an explicit negative
control; it is not generalized to every shifted input basis.

## Third-root edge

For `q=3c`, `4|c` and `3` not dividing c, use

```
M3 = 1/2 [[-1+i, -1-i], [1-i, -1-i]],
P  =     [[1,1], [1,-1]].
```

The exact identities are `M3^2+M3+I=0`, `M3 M3*=I`, and `PP*=2I`.
Both P and M3P have fourth-root entries. The other power, M3^2P, has
zero entries and entries of squared absolute value2.

With `a=3^-1 mod c`, `b=c^-1 mod3`, choose input frame
`rho(zeta_q^s)P`. Then

`B_e=zeta_c^(a(e+s)) M3^(b(e+s)) P`,

and `B_e B_f*=2rho(zeta_q^(e-f))`. Its flat domain is
`b(e+s) mod3 in{0,1}`. All three input-coset positions are tested.
The cyclic basis, nine Gaussian cross Grams and full cyclotomic
polynomial are checked exactly. This gives actual edges such as
`BH(n,12)->BH(2n,4)` and `BH(n,24)->BH(2n,8)` on their verified domains.

## Actual branching and fixed inputs

At each node the breadth-first executor tries both ordinary
repeated-prime modes and every available partial-third/fifth input basis.
The phase strictly decreases on every edge. Complete arrays, phases and
orders determine the duplicate key. General rational basis equivalence
and signed permutations are not mistaken for literal equality.

The90 fixed roots comprise the seven common/fifth registry sources,
cubic phase24 levels0,1,2, the existing P7 twelfth-root representative,
all69 acquired twelfth-root catalogue representatives, and ten explicit
published-matrix gauges selected to omit a cubic residue class. Those
gauges are recorded fixed starting representatives, with verified full
Grams and raw source hashes, not operations inside the tree. The first
alphabet-admitting acquisition gauge need not fit the more specific
partial-third domain, so those additional representatives matter.

Successful towers include

`BH(15,120)->BH(30,40)->BH(120,8)->BH(240,4)->H480`,

`BH(35,120)->BH(70,40)->BH(280,8)->BH(560,4)->H1120`,

and alternate orderings through phases24 and12. Other branches stop at
complete Butson matrices of phases6,10 or30 because no further implemented
edge applies. They are retained as actual intermediate matrices.

## Replay and accessors

```sh
PYTHONPATH=src:. python -m experiments.scalar_lift.partial_fifth_descent --check certificates/scalar-lift/partial-fifth-descent.json
```

The [implementation](../experiments/scalar_lift/partial_fifth_descent.py)
exposes `scalar_certificate(q)`, `block_table(q,s)`, `restrict(H,q,s)`,
their third-root counterparts, and `SavedNode(path)`. A saved node's
`replay()` reconstructs its root and every actual ancestor edge;
`entry(i,j)` reads its literal exponent. `as_verified_source()` returns
the existing immutable `VerifiedSource` protocol for subsequent scalar
restrictions. The shared input/output basis determinant obstruction is
respected: each input/output frame and its norm factor are explicit.

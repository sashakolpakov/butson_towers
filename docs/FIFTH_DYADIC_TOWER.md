# Exact fifth-root lifts through the dyadic scalar tower

The existing fifth-root lift already contains a signed-permutation
complex structure. It therefore extends exactly to phase20,40 and80,
with real block widths8,16 and32. The admitted raw alphabet is

```text
L_q = mu_q minus mu_(q/5),     q=20,40,80,
```

equivalently all exponents not divisible by5. This alphabet is closed
under complex conjugation and under multiplication by dyadic roots.
Large compatible Butson sources remain missing; the tower itself is
not a new large Hadamard matrix.

## The commuting complex structures

The [verified dense fifth-root lift](DENSE_SIGN_LIFTS.md) uses a Gaussian
4-by-4 unitary matrix `M` with `Phi5(M)=0`, and all entries of `2M^k`
are fourth roots for `k=1,2,3,4`. Such partial fifth-root morphisms are
part of the standard Gow/sound-pair construction framework; see
[Egan–Ó Catháin, introduction](https://arxiv.org/html/2610.28394v1#S1).
No novelty of that base construction is claimed here.

Realify Gaussian scalars with

```text
J2 = [[0,-1],[1,0]],    Z=H2=[[1,1],[1,-1]].
```

Let `A=rho(M)` and `J8=I4 tensor J2`. Then `A` is rational orthogonal,
`Phi5(A)=0`, `J8^2=-I`, and `J8 A=A J8`. Its numerator `2A` is integral.

For completeness, the finite search enumerates every signed-permutation
8-by-8 matrix `J` with `J^2=-I`. Its underlying permutation must be a
fixed-point-free involution: there are105 matchings and16 orientations
per matching, hence1680 candidates. Exactly six commute with `A`.
The certificate records all six, including the canonical `J8`. This is
an exhaustive statement about this fixed action, not about all possible
fifth-root actions or frames.

## A necessary frame convention

The original frame is

```text
P0 = 2 I4 tensor Z,     P0 P0^T=8I.
```

It is sparse, with entries0 and±2. It must not be called a Hadamard
matrix. Nevertheless `A^k P0` is a sign matrix for every nonzero
`k mod5`, exactly the powers needed by the partial lift.

The frame `P=A P0` is itself sign-flat. Using it shifts all source
exponents by one fifth-root phase. Both conventions are certified, but
the **raw sparse-frame convention is preferred**: its alphabet is simply
`e mod5 != 0`, with no shifted omitted residue.

## The factorization tower

Starting with `R8=J8`, define signed-permutation roots

```text
R16 = [[0,R8],[I8,0]],
R32 = [[0,R16],[I16,0]].
```

Their squares are `diag(R8,R8)` and `diag(R16,R16)`, respectively.
They have orders8 and16 and commute with the corresponding block
diagonal fifth-root actions. This gives

```text
Q(zeta20)  subset  Q(zeta40)  subset  Q(zeta80).
```

Write `c=q/5`, `n=2c`, `R=R_n`, and let `D` be the block diagonal
extension of `A`. Choose

```text
T_q = R^(5^(-1) mod c) D^(c^(-1) mod5).
```

Explicitly these are `J8 A^4`, `R16^5 D^2`, and `R32^13 D`.
Then `T_q^5=R`, `T_q^c=D`, and `Phi_q(T_q)=0`. Since
`deg(Phi_q)=n`, this is a faithful rational representation at the minimal
field degree. Orthogonality makes it a representation preserving the
involution, with complex conjugation represented by transpose.
The field inclusions are also verified directly:
`T40^2=diag(T20,T20)` and `T80^2=diag(T40,T40)`.

For `t=c/4`, extend the raw frame to

```text
P0_q = H_t tensor P0,       P0_q P0_q^T = nI.
```

The dyadic action is signed monomial, so it preserves flatness. Thus

```text
B_e = T_q^e P0_q
```

has only±1 entries exactly when `e mod5 != 0`. Every phase has an
integral block, and each block has Gram `nI`; the excluded blocks have
zeros and entries±2. The shifted flat frame `D P0_q` instead admits
`e+c != 0 mod5`. The omitted residues for q20,40,80 are respectively
1,2,4. Translating a raw source by `e -> e-c` makes the shifted-frame
construction produce exactly the same real matrix.

The certificate verifies all phase blocks, the full cyclotomic
polynomial identities, and every one of the `q^2` cross-Gram identities

```text
B_e B_f^T = n T_q^(e-f).
```

Consequently every actual `BH(m,q)` whose entries lie in `L_q` gives an
actual real Hadamard matrix of order `nm`. Orthogonality follows by
applying the field representation to each source row inner product;
flatness follows from the exact admitted alphabet. There is no completion
or replacement of zero entries.

This also has a direct description: adjoining the dyadic scalars to the
Gaussian fifth-root construction first gives `BH(4m,c)`; the usual
power-of-two scalar lift then gives `H(2cm)`. The explicit matrices here
certify this factorization and its alphabet at each level.

## Controls and high-source requirements

The chirped Fourier source

```text
F_ij = zeta5^(3(i+j)^2+1),    i,j in F5,
```

uses only fifth-root exponents `{1,3,4}`. The verified alphabet actually
admits all four nonzero fifth classes, with arbitrary dyadic factors.
The controls apply row and column dyadic gauges to this source, verify
its cyclotomic Gram exactly, and construct sign matrices of orders40,80
and160. Both complete integer Gram matrices are checked. These are known
orders and serve only as construction controls.

The resulting high routes that pass the checked necessary source filters are:

| Missing source | Phase | Block width | Real output |
| --- | --- | --- | --- |
| 1499 | 40 | 16 | 23984 |
| 1499 | 80 | 32 | 47968 |
| 2999 | 80 | 32 | 95968 |

The tempting phase20 sources of orders1499 and2999 are impossible, as
both primes are `-1 mod20`. The phase40 source of order2999 is likewise
impossible because `2999=-1 mod40`. In each case primes above p are
self-conjugate, so the determinant norm equation
`det(B) conjugate(det(B))=p^p` would have an even valuation on the left
and an odd valuation on the right. These source exclusions do not affect
the valid scalar actions, and disappear at the indicated higher levels.

Each source must avoid the entire dyadic subgroup, not merely the scalar1.
None of these large sources is supplied by this certificate. The cubic
border obstruction does not automatically apply to this different,
conjugation-closed alphabet.

Artifacts:

- [Exact tower and validated lift function](../experiments/scalar_lift/fifth_dyadic_tower.py).
- [Matrices and complete certificate](../certificates/scalar-lift/fifth-dyadic-tower.json).
- [Finite audit, tower and output tests](../tests/test_fifth_dyadic_tower.py).

```sh
python -m experiments.scalar_lift.fifth_dyadic_tower \
  --check certificates/scalar-lift/fifth-dyadic-tower.json
```

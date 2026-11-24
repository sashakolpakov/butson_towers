# Flatness certified at the endpoint

The [executor](../experiments/scalar_lift/projective_flat_descent.py) now
constructs complete flat endpoints from verified sources. Its
[certificate](../certificates/scalar-lift/projective-flat-descent.json)
contains polynomial flatness identities, a nonflat quadratic/cubic scalar
tower with its actual metrics, and exact full Gram checks of the outputs.
The [literal matrix bundle](../artifacts/projective-flat-descent/examples.npz)
contains **five pairwise inequivalent real Hadamard matrices of order 768**,
their quaternary sources, 18 sixth-root matrices of order 42, five of order
294, a twelfth-root order-42 specialization, and the parameter-family data.

At least **four** monomial classes occur among the 18 order-42 outputs.
The five order-294 outputs are pairwise monomially inequivalent. Different
exact Haagerup multiplicity histograms prove these statements; equal
histograms are never used to assert equivalence or completeness.

The small nonmonomial action/frame pairs are the published matrices in
[Egan–Ó Catháin, Theorem 3.1](https://arxiv.org/abs/2610.28394).
The morphism and existence of these orders are known. The output families,
literal examples and certificates are new additions here; publication
priority for particular classes has not been established. The real examples
use the classical quaternary-to-real scalar transfer; the complex examples
use the published nonmonomial seventh-root action.

## A real endpoint: five inequivalent H768 matrices

Use the quadratic scalar restriction `Q(i)/Q` in the basis `(1,-i)`.
For `z=x+iy`, its raw multiplication matrix and canonical trace metric are

\[
R(z)=\begin{pmatrix}x&y\\-y&x\end{pmatrix},\qquad G=2I_2.
\]

Take `C=H2` and `D=I`. They satisfy `C^T C=G` and `DD^T=2G^-1`.
The raw matrices have zeros, while the terminal block is

\[
CR(z)D=\begin{pmatrix}x-y&x+y\\x+y&y-x\end{pmatrix}.
\]

On the projective fourth-root locus the source relations are
`g=x^2+y^2-h^2=0` and `f=xy=0`. Each terminal entry `b` satisfies
`b^2-h^2=g+2f` or `g-2f`, identically. At `h=1`, every final entry is
therefore a sign. This is the real flatness guarantee, expressed as explicit
ideal membership through the actual scalar restriction and final maps.

Construct `A=H96` from the checked `H12 tensor Sylvester H8`, and set

\[
S_{(a,b),(c,d)}=A_{ac}i^{bd+p_{bc}},\qquad
a,c=0,\ldots,95,\quad b,d=0,\ldots,3.
\]

Set `p_(0,c)=p_(b,0)=0`; the other **285 fourth-root parameters** are free.
The Fourier sum over `d` and the verified Gram of `A` prove
`S S*=384I` for every parameter choice. The terminal map consequently gives
a finite family of real `H768` matrices with flatness guaranteed. All five
saved sources and both full Grams of every real output are checked directly.

| Bundle key | Nonzero parameter exponents `(b,c,p)` | Binary code dimension |
|---|---|---:|
| `h768_example0` | none | 17 |
| `h768_example1` | `(1,1,2)` | 18 |
| `h768_example2` | `(1,1,1)` | 19 |
| `h768_example3` | `(1,1,1)`, `(2,2,2)` | 20 |
| `h768_example4` | `(1,1,1)`, `(2,2,1)` | 21 |

For a sign matrix let `b_i` encode negative entries of row `i` over `F2`.
The dimension of `span(1,b_i+b_0)` is invariant under all row and column
sign changes and permutations. Column sign changes cancel in differences;
row sign changes add the already included all-one vector; changing the
reference row preserves the span. The different dimensions therefore prove
pairwise inequivalence. Four outputs are in different classes from the
tensor control `h768_example0`. No previously unknown real order is claimed.

## Seventh-root endpoint flatness as a polynomial identity

Work over `K=Q(zeta_6)`, with `L=K(zeta_7)`. The exactly verified unitary
action `U` satisfies `Phi_7(U)=0`. A frame `Z` satisfies `ZZ*=6I`, and the
seven matrices `U^j Z` have sixth-root entries. These small identities are
checked over the rational power basis of `K`, without floating point.

Represent a source seventh root in coordinates `x_0,...,x_5` of
`1,zeta_7,...,zeta_7^5`, and homogenize with `h`. The source points at `h=1`
are the six coordinate vectors and `(-1,...,-1)`. Put

\[
\lambda=\frac{h-\sum_i x_i}{7},\qquad
g_i=x_i^2-hx_i-2h\lambda,\qquad
g_{ij}=x_ix_j-h\lambda\quad(i<j).
\]

All 21 quadrics vanish on those source points. For each of the 36 endpoint
coordinates let `c_i=(U^i Z)_(a,b)` and
`k_ij=c_i conjugate(c_j)+conjugate(c_i)c_j`. The certificate gives the exact
integer coefficients `k_ij` and verifies

\[
\left|\sum_i c_i x_i\right|^2-h^2
=\sum_i g_i+\sum_{i<j}k_{ij}g_{ij}.
\]

Indeed, `|c_i|=1` and `sum k_ij=-5`, since `-sum c_i=(U^6 Z)_(a,b)` also
has modulus one. This is an explicit ideal-membership certificate for every
pulled-back flatness equation. Independent symbolic tests expand all 72
identities for the two published pairs. Unit phases from the lower field
multiply both sides by their norm and preserve the endpoint condition.

The base-extended construction gives

\[
BH(n,q)\longrightarrow BH(6n,s),\qquad
q\mid7s,\quad 6\mid s,\quad 7\nmid s.
\]

For `q=7s`, write `zeta_(7s)^e=zeta_7^a zeta_s^b`, with
`a=e*s^-1 mod7` and `b=e*7^-1 mod s`. Its terminal block is
`zeta_s^b U^a Z`. Every used block is therefore flat by the polynomial
identity and belongs to the declared lower root alphabet.

## Nonflat intermediate matrices and the transported form

The actual field tower is `L -> E=K(eta) -> K`, of degrees two and three,
where `eta=zeta_7+zeta_7^-1`, `xi=zeta_7-zeta_7^-1`,

\[
\eta^3+\eta^2-2\eta-1=0,\qquad \xi^2=\eta^2-4.
\]

In the quadratic basis `(1,xi)`, multiplication by `zeta_7` is

\[
\frac12\begin{pmatrix}\eta&\eta^2-4\\1&\eta\end{pmatrix}.
\]

The implementation restricts this operator in the cubic basis
`(1,eta,eta^2)`. These coordinate arrays contain zeros and fractions.

The chosen incoming source form is explicitly `kappa I`, where the positive
Hermitian scalar `kappa` is derived and saved in both bases. A complete source
already satisfies `H*(kappa I)H=n(kappa I)` because it satisfies `H*H=nI` and `kappa`
is central. The executor transports this form by relative trace through
both edges. It also saves the canonical unweighted trace form separately;
the two forms are never silently identified.

Exact final maps `C,D` satisfy `C* C=G` and `D D*=6G^-1` for the transported
metric `G`. Every scalar block `C rho(zeta_7)^j D` is checked to equal
`U^j Z`. Thus ordinary orthogonality and flatness hold at the endpoint,
with no intermediate flatness assumption.

## Parameter families and exact examples

The order-42 examples descend the nine saved `BH(7,42)` representatives
through both action/frame pairs. The source Grams are rechecked. In addition,
an order-seven one-parameter source is verified as a full Laurent identity;
its order-42 image is a one-parameter flat family. Specializing its parameter
to `zeta_12` gives a literal `BH(42,12)`.

For the larger family take `a,b,c,d` in `Z/7Z` and independent unit phases
`t_(b,c)`, with `t_(0,c)=t_(b,0)=1`. The complete source is

\[
H_{(a,b),(c,d)}=\zeta_7^{ac+bd}t_{b,c}.
\]

Summing over `d` proves orthogonality when `b` differs; when `b` agrees the
phases cancel and summing over `c` proves the remaining identities. Its
endpoint has blocks `t_(b,c) U^(ac+bd) Z`, so every entry has modulus one
for all 36 independent unit phases. Both complete endpoint Laurent Grams
are also checked coefficientwise. The family is dephased throughout, and
the 36 parameters occupy disjoint nonempty entry sets. It therefore has
36 independent parameters even after removing row/column phase gauges;
permutations identify only finitely many points.

The five saved sixth-root specializations are named `zero`, `one`, `two`,
`diagonal` and `mixed`. All have full exact row and column Grams `294I`,
and all five Haagerup histograms differ.

## Replay and access

```sh
PYTHONPATH=src:. OPENBLAS_NUM_THREADS=1 .native-env/bin/python \
  -m experiments.scalar_lift.projective_flat_descent \
  --check certificates/scalar-lift/projective-flat-descent.json
```

To regenerate, replace `--check` with `--output` and the same certificate
path. Keys `h768_example*` contain literal `+1/-1` entries. The corresponding
`bh384_source*` keys contain exponents of `i`. Keys `bh42_pair*` and `bh294_*`
contain exponents of `zeta_6`; `bh42_q12` contains exponents of `zeta_12`.
The `*_base` and `*_labels` arrays describe the parameter families.

```python
import numpy as np
with np.load("artifacts/projective-flat-descent/examples.npz") as data:
    H = data["h768_example4"]  # actual +1/-1 entries
    E = data["bh294_mixed"]   # exponents of zeta_6
```

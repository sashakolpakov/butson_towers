# Foundations of Butson Towers

This project studies exact scalar restriction of complete orthogonal
designs and the additional conditions that make the final entries flat.
“Butson Towers” is our descriptive project name. The three papers below
use the language of Butson morphisms; they do not use “Butson towers” as
a named construction. The project name makes no claim to have introduced
restriction of scalars or Butson morphisms.

## The starting reference and its antecedents

The immediate starting reference for this research direction was
**Ronan Egan and Pádraig Ó Catháin,
[*On complete morphisms of Butson Hadamard matrices*](https://arxiv.org/abs/2610.28394),
arXiv:2610.28394v1, 20 October 2026**. Its
[Theorem 2.2](https://arxiv.org/html/2610.28394#S2.Thmtheorem2.2)
describes every complete morphism by a pair `(U,Z)`: `U` is unitary,
`Phi_k(U)=0`, `Z` is a Butson matrix, and all used products `U^j Z` have
the required root entries. Its
[Theorem 3.1](https://arxiv.org/html/2610.28394#S3.Thmtheorem3.1)
constructs a complete morphism `BH(n,7) -> BH(6n,6)`.

The older general framework is **Ronan Egan and Pádraig Ó Catháin,
[*Morphisms of Butson classes*](https://arxiv.org/abs/1707.08815)**,
first submitted 27 July 2017; the cited v3 is dated 25 October 2017.
It appeared in *Linear Algebra and its Applications* 577 (2019), 78–93.
It develops morphisms that reduce the order of the roots while increasing
matrix dimension, using spectral and entry conditions on sound pairs.
It treats partial as well as complete morphisms and connects earlier
quaternary-to-real and unreal-sixth-root constructions.

The direct reference for our repeated cyclotomic reduction executor is
**Pádraig Ó Catháin and Eric Swartz,
[*Homomorphisms of matrix algebras and constructions of Butson-Hadamard matrices*](https://arxiv.org/abs/1904.10771)**,
first submitted 24 April 2019; v3 is dated 16 August 2019.
It appeared in *Discrete Mathematics* 342 (2019), article 111606.
[Theorem 1.1](https://arxiv.org/html/1904.10771#S1.Thmtheorem1.1)
gives

\[
\mathrm{BH}(n,k)\longrightarrow\mathrm{BH}(pn,k/p)
\qquad\text{when }p\text{ is prime and }p^2\mid k.
\]

Section 2 explicitly embeds the cyclotomic field in a matrix algebra over
the lower field, then multiplies by a compatible Butson block.
Theorem 1.2 iterates these steps: if `k=mt` and every prime divisor of
`k` divides `t`, it gives `BH(n,k) -> BH(mn,t)`. The paper generalizes
an earlier construction of Östergård and Paavola. This is an existing
construction theorem, used here with exact source and entry certificates.

## Restriction preserves the norm with its metric

Let a complete source over `L` satisfy `H H*=w I_n`, and let `L/K` have
degree `d` and be stable under complex conjugation. For a basis `beta`,
the multiplication action and trace metric satisfy

\[
G_{ab}=\operatorname{Tr}_{L/K}(\overline{\beta_a}\beta_b),\qquad
\rho(\bar z)=G^{-1}\rho(z)^*G.
\]

The block matrix `A=rho(H)` therefore has the exact identity

\[
A(I_n\otimes G^{-1})A^*=w(I_n\otimes G^{-1}).
\]

The rational dimension has grown, while the inherited Gram scale remains
`w` in this identity. A nonidentity incoming form must itself be transferred
by trace. Trace transitivity relates a field tower to direct restriction
after the actual basis change. The
[norm-descent reference](EXACT_NORM_DESCENT_EXECUTOR.md) records that
general form of the invariant.

For example, restricting sixth roots in the basis `(1,zeta_6)` gives

\[
\rho(\zeta_6)=\begin{pmatrix}0&-1\\1&1\end{pmatrix},\qquad
G=\begin{pmatrix}2&1\\1&2\end{pmatrix}.
\]

The coordinate entries may contain zeros or have different magnitudes.
These weighted intermediate matrices retain their exact norm identity.

## Terminal flatness is an additional certificate

In coordinates where conjugation is adjoint, a frame `F` with `FF*=cI`
gives blocks `B(z)=rho(z)F`. Their full block Gram follows from the
source norm:

\[
\sum_k B(H_{ik})B(H_{jk})^*
=cw\,\delta_{ij}I_d.
\]

Every actually used block must pass its entry test. A dense sign endpoint
of order `dn` additionally needs `cw=dn`. A general source's weight and
order must not be identified to obtain that equality. Zeros in weighing
or conference inputs remain zeros under a scalar-linear reduction.

A complete Butson morphism works uniformly for every source matrix over
the entire stated root alphabet. A map certified only on an admitted
palette or particular source has that narrower domain. Our source-specific
and metric-retaining constructions carry their actual domains explicitly.
The [projective flatness executor](PROJECTIVE_FLAT_DESCENT.md) also records
polynomial endpoint conditions through nonflat intermediate stages.

## Products of fields and literal reconstruction

The Williamson constructions use commutative group algebras such as
`Q[C_w]` and `Q[G]` for explicitly given finite abelian groups. These are
products of cyclotomic number fields. Conjugation inverts group elements
and acts as complex conjugation on every field factor.

The normalized regular trace is the coefficient of the identity. It gives
metric `I` in the inverse group basis and equals the appropriately normalized
sum of the field traces, including their multiplicities. The exact change
of coordinates and reconstruction map are part of the restriction data.
Individual field components may have zeros or be nonflat; the terminal
entry test concerns the reconstructed blocks.

The [Williamson records](WILLIAMSON_ORDER_RECORDS.md) and
[Spence family proof](SPENCE_SCALAR_FAMILIES.md) retain these algebras,
field factors, dimensions, source norms and actual sign images. The
new order claims concern individually completed endpoints in the
[dated ledger](HADAMARD_HIGHER_ORDER_RECORDS.md). The five all-`r`
construction theorems and the finite members' novelty labels have
separate scopes.

## Formal verification status

The [Lean development](FORMALIZATION.md) connects literal finite witnesses
to the eighteen announced existence theorems and all six output families:
five Spence-source families and the cubic family. It proves the regular action and trace metric, the complete
orthogonal-design construction, final sign entries and the exact matrix
recurrences. Run `./scripts/check_formal.sh` for the build and fixed axiom
audit; the [portable package](../release/lean_verification/README.md)
contains the proof sources and a reproducibility verifier.

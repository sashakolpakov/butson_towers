# Forward scalar descent from a genuine Butson source

The forward construction starts with a verified Butson matrix, represents its scalars successively over smaller fields, and supplies a compatible final frame. **If the resulting blocks are signs, the output is a real Hadamard matrix.** The Gram identity is inherited exactly. The search therefore concerns admissible source matrices, scalar actions and frames that make this forward construction flat at the desired endpoint.

Intermediate matrices need not remain Butson matrices. Requiring roots at every intermediate stage is a useful sufficient subfamily, not a complete description of scalar restriction.

The current computation uses scalar restrictions only. The
[weighted divisor-field tree](WEIGHTED_SCALAR_TREE.md) now retains the
actual matrix, its tower basis, and its transported metric through odd-prime
removal and down to the rationals. Different coordinate bases are retained;
only literal action/metric duplicates are merged. The
[uniform alphabet obstruction](SCALAR_RESTRICTION_OBSTRUCTIONS.md) supplies
exact decisions independent of the chosen tower: two of the 20 fixed inputs
are excluded for uniform sign endpoints, while 18 pass that particular gate.
No source construction or frame solver is part of this continuation.

The [implemented repeated-prime executor](FORWARD_BUTSON_DESCENT.md)
replays a genuine order-1000 source through a complete real endpoint.
The [character-source registry](FORWARD_CHARACTER_SOURCES.md) also supplies
exact implicit known Butson inputs of orders 5,329 and 58,081. Their raw
alphabets fail a uniform final sign-map condition, recorded separately from
their fully proved source Grams; they are not claimed to have real endpoints.

## A sufficient forward edge, with its exact entry condition

Let \(L/K\) be a conjugation-stable extension and let

\[
 \rho:L\longrightarrow M_d(K),\qquad
 \rho(\bar z)=\rho(z)^*
\]

be a faithful unital representation. Choose a frame \(F\in M_d(K)\) with
\(FF^*=dI_d\). For a verified \(H\) with \(HH^*=nI_n\), replace each entry \(z\) by

\[
 B(z)=\rho(z)F.
\]

The block matrix \(\widehat H\) then satisfies

\[
 (\widehat H\widehat H^*)_{ij}
 =d\,\rho\left(\sum_k H_{ik}\overline{H_{jk}}\right)
 =dn\,\delta_{ij}I_d.
\]

Thus:

- If every block used by this source has entries in a declared lower root alphabet \(\mu_q\), the result is \(BH(dn,q)\).
- If every used block has entries in \(\{1,-1\}\), the result is a real Hadamard matrix of order \(dn\).

The block-entry test is necessary and sufficient for flatness of this fixed representation/frame construction. The frame itself need not be flat if the source does not use the scalar one. Its Gram, and the flatness of every phase actually used, must still be checked.

Verified edges compose: the order is multiplied by the product of their widths, and each output's actual entry alphabet becomes the next edge's input alphabet. A known source that passes every edge produces the final matrix by construction; no probabilistic assumption enters this implication.

## General scalar restrictions retain a metric

In an ordinary field basis, multiplication matrices usually satisfy

\[
 \rho(z)^*G=G\rho(\bar z)
\]

for a Hermitian trace metric \(G\), rather than for the identity metric. In source-block coordinate order, the restricted matrix then satisfies

\[
 \rho(H)(I_n\otimes G^{-1})\rho(H)^*
 =n(I_n\otimes G^{-1}).
\]

This is valid forward progress, even if its entries are sparse coefficients or sums of roots. Keep the metric, conjugation and block reconstruction data with the matrix. Do not reject a path merely because this intermediate array is not another Butson matrix.

The scalar on the right is still the original source order \(n\), even
though the restricted matrix has order \(dn\). It is not an ordinary
Hadamard Gram of scale \(dn\). The root-flat edge above supplies that
amplification through its explicitly verified norm-\(d\) frame. A general
weighted node is therefore recorded with its precise scale and is never
promoted to a Hadamard matrix from integrality or field descent alone.

One way to obtain an ordinary unitary action is to supply an actual matrix \(C\) over the lower field with \(C^*C=G\), then use \(C\rho(z)C^{-1}\). Equivalently, blocks

\[
 B(z)=C\rho(z)D,\qquad C^*C=G,\qquad DD^*=dG^{-1},
\]

have the required ordinary Gram law. Their entries must again pass the final root or sign test. A numerical complex square root of \(G\) does not establish that \(C\) lies in the required field.

For a sequence of restrictions one may retain all weighted intermediate actions and apply a single final frame to the composed representation. Requiring a separate flat frame at each stage is stronger and may miss such constructions. The [phase-120 common-base lift](PHASE120_COMMON_BASE.md) supplies a positive example with a final compatible frame.

## Full-alphabet prime-power steps

There is an unconditional root-preserving edge whenever
\(Q=pq\), where \(p\) is prime and \(p\mid q\). In this case
\([\mathbb Q(\zeta_Q):\mathbb Q(\zeta_q)]=p\). Set

\[
 R=\begin{pmatrix}
 0&0&\cdots&\zeta_q\\
 1&0&\cdots&0\\
 0&1&\cdots&0\\
 \vdots&&\ddots&\vdots
 \end{pmatrix},\qquad R^p=\zeta_q I,
 \qquad (F_p)_{rs}=\zeta_q^{(q/p)rs}.
\]

The monomial action \(R\) is unitary, \(F_pF_p^*=pI\), and
\(\Phi_Q(R)=\Phi_q(R^p)=0\). For an upper-root exponent \(e=pa+b\), \(0\leq b<p\), the block \(\zeta_q^a R^bF_p\) has exponent

\[
 a+\mathbf1_{r<b}+\frac qp\big((r-b)\bmod p\big)s
 \pmod q
\]

in row \(r\), column \(s\). Every entry is a lower root. Consequently

\[
 BH(n,pq)\longrightarrow BH(pn,q)
\]

for the entire upper root alphabet. A quadratic implementation is already independently certified in [QUADRATIC_SOURCE_INVERSE.md](https://github.com/sashakolpakov/complex-hadamard/blob/976cede671f40c810d593c8d5815e26e11a03ed6/docs/QUADRATIC_SOURCE_INVERSE.md); its forward map is `descend_source`. It uses the dyadic generator instead of the full primitive generator, so its literal blocks have different exponent formulas.

More generally, write \(Q=p^a m\), \(p\nmid m\), and choose \(\eta=\zeta_{p^a}=\zeta_Q^m\). Put \(b\equiv e m^{-1}\pmod p\), \(0\leq b<p\), and \(a_0=(e-mb)/p\). With the companion carry \(\zeta_q^m\), the block exponents are

\[
 a_0+m\mathbf1_{r<b}+\frac qp\big((r-b)\bmod p\big)s\pmod q.
\]

This variant preserves the colors from the coprime \(m\)-part, since the companion and Fourier multipliers have only \(p\)-power order. The full primitive-generator formula need not preserve those colors. Both are valid edges; their actual intermediate entry alphabets must be tracked separately.

These steps remove repeated prime factors. For a pure dyadic alphabet they finish at

\[
 BH(n,2^a)\longrightarrow BH(2^{a-1}n,2)
 =H(2^{a-1}n).
\]

This general repeated-prime construction is the established
[Ó Catháin–Swartz morphism](https://arxiv.org/abs/1904.10771), extending
Östergård–Paavola. It provides certified edges for the forward executor,
not a newly claimed construction theorem.

For phase 120, the full-alphabet quadratic chain gives

\[
 BH(n,120)\longrightarrow BH(2n,60)
 \longrightarrow BH(4n,30).
\]

It does not remove the remaining cubic and fifth phases. The fields \(\mathbb Q(\zeta_{30})\) and \(\mathbb Q(\zeta_{15})\) coincide, but the root alphabet \(\mu_{30}\) still includes its signs. Field canonicalization must not replace it with \(\mu_{15}\).

## What a final partial alphabet achieves

A scalar-linear sign map cannot be flat on an entire odd prime polygon. If it were flat on \(1,\zeta_p,\ldots,\zeta_p^{p-1}\), then

\[
 \sum_{j=0}^{p-1}\rho(\zeta_p^j)F=0
\]

would be an odd sum of signs at each matrix entry. This is impossible. The statement concerns this scalar-linear common-frame construction; it is not a nonexistence theorem for all other Butson-to-Hadamard constructions.

Partial alphabets solve this particular difficulty. The existing phase-120 representation is rational orthogonal of dimension 32 and has a frame admitting exactly 64 phases. In raw coordinates the allowed roots are

\[
 \omega^a\theta^b\eta^c,
 \quad a\in\{0,1\},\quad b\in\{1,2,3,4\},
 \quad c\in\{0,\ldots,7\},
\]

where \(\omega=\zeta_3\), \(\theta=\zeta_5\), \(\eta=\zeta_8\). The restriction factors have degrees two and four over the shared eighth-root field, followed by its rational degree four. The compatible final frame makes all used blocks literal signs. Complete source controls of orders 15 and 35 therefore give real Hadamard matrices of orders 480 and 1120, respectively.

The forward implication for the current target is exact:

\[
 \text{an admitted }BH(1499,120)
 \quad\Longrightarrow\quad H(32\cdot1499)=H(47968).
\]

The order-1499 source remains missing. Restricting a known source of another order preserves its own order factor; for example the verified order-585 input would produce order \(32\cdot585\), not 47968, under this width-32 endpoint. An alternative starting source can reach the prescribed target only if its order times the subsequent width product equals 47968.

## The source-first forward search

The descent queue begins with **actual verified sources**. An unfinished
source-construction problem belongs in a separate queue and cannot be
counted as a completed starting matrix. A descent state records the source
hash, current field, actual entry alphabet, scalar representation, adjoint
metric, accumulated dimension, reconstruction data, and exact proof status.

For each source, branch over certified restriction edges and compatible
choices of scalar bases and frames. Allowed row/column phase gauges of the
source can also be explored, with their effect on the actual entry alphabet
retained. A failed particular frame does not exclude all other frames.
Weighted intermediates continue without a premature Butson requirement.
If the final blocks are signs, emit the resulting explicit matrix or lazy
entry formula, its order, and the inherited exact Gram certificate. Compare
that real order with known constructions only after obtaining the endpoint.

If a starting source comes with a factorization, retain it as additional
data. For example, descending two actual Butson tensor factors separately
and tensoring their verified real outputs is a valid forward construction.
Such a map can depend on the two factor labels rather than only on their
complex product. It should not be conflated with one scalar-linear map on
the single cyclotomic field containing that product; the uniform odd-polygon
obstruction above does not automatically apply to the factor-dependent map.
An algebraic implementation may retain a tensor product of coefficient
algebras, together with its source Gram proof, instead of prematurely
identifying it with a single field.

The forward composition guarantees success **once** a compatible verified
source and flat endpoint have been found. Field descent alone does not prove
that every starting source or tower has such an endpoint. The fixed
order-1499 construction is one possible starting-source project; it is not
the definition of this broader forward search.

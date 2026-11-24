# A real Hadamard matrix of order 7572 by scalar reduction

An explicit sign matrix `H` of order **7572** satisfies

\[
HH^T=7572I_{7572}.
\]

The [matrix, complete sequence data, independent audits and standalone replay](../artifacts/hadamard-7572/README.md)
are retained together. Every one of the 28,663,806 distinct row pairs has
Hamming distance 3786; the diagonal entries of the Gram matrix are 7572.
The packed sign bits have SHA256
`bc3cb6fca9243a6ef01c82de2a446080c8866e0f30d76c0f23ef0508ae5e2e9b`.

The primary-literature audit on **24 November 2026** found this order still
listed as unresolved in [Egan–Kotsireas–London–Ó Catháin, Table 5](https://arxiv.org/html/2610.27175v2#S4.T5),
the [current construction frontier](https://raw.githubusercontent.com/vatskishan03/Hadamard/main/FRONTIER.md),
and the [Hadamard Matrix Explorer](https://www.math74.com/Hadamard/unknown.html).
The exact URLs, dated observations and frontier commit are recorded in
[the novelty audit](../artifacts/hadamard-7572/novelty.json).
This is an existence result at a new real Hadamard order, under the
[author-knowledge status convention](HADAMARD_HIGHER_ORDER_RECORDS.md).
Reusing this frame with further verified sources gives new orders 17668,
22716, 32812, 98436 and 295308; the completed H885924 recipe is a known-order
Paley control.

## Exact length-631 data

Let `C` fix channel 0 and cycle channels `1 -> 2 -> 3 -> 1`. The multiplier
458 has order 21 modulo 631. Enumerate the origin first and then the
nonzero multiplier orbits by their least unused representative. The states
are

```text
0 0 5 4 4 0 7 0 7 2 3 0 1 2 6 1 4 2 0 1 1 4 3 4 7 3 3 5 6 2 7
```

A state `2*j+epsilon` puts sign `(-1)^epsilon` in channel `j` at the orbit
representative. At `458^k` times that representative, rotate the channel
by `C^k`. This gives four sequences `b_i` with one signed nonzero entry in
each column. Their row sums are `(22,7,7,7)`, and exact expansion verifies
every shift:

\[
\sum_{i=0}^3\sum_{x\in\mathbb Z_{631}} b_i(x)b_i(x+s)
=631\,\delta_{s,0}.
\]

The [literal sequences](../artifacts/hadamard-7572/t631.json) and an
[independent sequence audit](../artifacts/hadamard-7572/independent-sequences.json)
accompany the recipe. Their circulants are periodic T-matrices of order 631.

Put

\[
D=\begin{pmatrix}
0&1&1&1\\-1&0&1&-1\\-1&-1&0&1\\-1&1&-1&0
\end{pmatrix},\qquad a=Db.
\]

Since `D^T=-D` and `D^T D=3I`, the summed autocorrelations of `a` are
`1893 delta_0`, and the summed mixed correlations of `a,b` vanish.
At each coordinate their supports complement each other and their nonzero
entries are signs. The real Goethals–Seidel construction therefore gives

\[
U=\operatorname{GS}(a),\quad V=\operatorname{GS}(b),\quad
UU^T=1893I,\quad VV^T=631I,\quad UV^T+VU^T=0.
\]

Define `P=U-V`, `Q=U+V` and `d=2524`. Both matrices have sign entries and
the complete exact checks give

\[
PP^T=QQ^T=dI,\qquad PQ^T+QP^T=dI.
\]

## The scalar action, trace form and endpoint

The rational operator

\[
R=QP^T/d
\]

satisfies `RR^T=I`, `R^T=I-R`, `R^2-R+I=0`, and `RP=Q`.
Consequently multiplication by `zeta_6` can act as `R` on `Q^2524`.
This is a `Q(zeta_6)`-module of rank **1262**, with rational dimension
**2524**. Complex conjugation acts by transpose.

To retain the incoming form explicitly, let `S=2R-I` and `xi=sqrt(-3)`.
The Hermitian form

\[
h(u,v)=\frac{u^Tv}{2}-\frac{u^TSv}{6}\,\xi
\]

has `Tr(h(u,v))=u^T v`. The canonical field basis `(1,zeta_6)` has trace
metric `[[2,1],[1,2]]`; the displayed `h` specifies the form in the actual
Euclidean coordinates.

Start with the complete Butson source

\[
M=\begin{pmatrix}-1&\zeta_6&\zeta_6\\
\zeta_6&-1&\zeta_6\\\zeta_6&\zeta_6&-1\end{pmatrix},
\qquad MM^*=3I_3.
\]

Its exact used alphabet is `{zeta_6,-1}`. With frame `F=P`, scalar
restriction sends these entries to `RF=Q` and `-F=-P`. Thus the completed
endpoint is the literal sign matrix

\[
H=\begin{pmatrix}-P&Q&Q\\Q&-P&Q\\Q&Q&-P\end{pmatrix}.
\]

The source scale **3** and frame scale **2524** give Gram scale **7572**,
which equals the expanded order. Only the used alphabet needs flat images:
the unused entry `zeta_6^2` would have image `Q-P=2V`, with zeros and entries
of magnitude two.

The standalone replay verifies the source, all sequence correlations,
flatness, the complete frame Grams and cross identity, equality with the
saved packed matrix, and the complete endpoint Gram. The independent audits
also read the actual saved NPZ directly and recover its scalar blocks.

## Reproduction

See [the standalone command and file descriptions](../artifacts/hadamard-7572/README.md).
Reconstruction and verification need Python and NumPy. They do not require
the stochastic search, its native executables, or the ignored research
checkpoints. The search used the auxiliary `SL(2,32768)` parity
parameterization and exact sign lifts modulo four and eight; the completed
integer witness and its proof are deterministic.

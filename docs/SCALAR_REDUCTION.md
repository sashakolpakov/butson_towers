# Butson-related designs by scalar reduction

The active higher-order direction is to reduce the scalars of **complete,
exactly verified designs**. This includes Butson matrices, complex generalized
weighing matrices and conference matrices. Here reduction means restriction
of scalars through number fields. Orthogonality is inherited; making the
endpoint have the desired root, sign or zero entries is the extra problem.

## The preserved identity

Suppose a source of order `n` over `L` satisfies `H H* = w I`, where
`w` is a positive real element of `K`, and `L/K` is stable under complex
conjugation. For a basis `beta` of degree `d`, let
`rho(z)` be multiplication by `z`, and retain the trace metric

\[
G_{ab}=\operatorname{Tr}_{L/K}(\overline{\beta_a}\beta_b),
\qquad \rho(\bar z)=G^{-1}\rho(z)^*G.
\]

The actual block-coordinate matrix `A = rho(H)`, of order `dn`, satisfies

\[
A(I_n\otimes G^{-1})A^*=w(I_n\otimes G^{-1}).
\]

This identity requires no new orthogonality search. It remains valid through
successive field restrictions, including intermediate arrays with rational
entries, zeros or sums of roots. Keep the metric and scale `w` with the
matrix. The enlarged order `dn` is a different number.

For example, sixth roots restricted to the rationals in the basis
`(1,zeta_6)` give

\[
\rho(\zeta_6)=\begin{pmatrix}0&-1\\1&1\end{pmatrix},\qquad
G=\begin{pmatrix}2&1\\1&2\end{pmatrix}.
\]

This preserves the displayed weighted identity, rather than an ordinary
Euclidean Gram. The [relative-restriction implementation](RELATIVE_SCALAR_RESTRICTION.md)
checks these relations exactly.

The [complete norm-descent reference](EXACT_NORM_DESCENT_EXECUTOR.md) also
carries an incoming nonidentity Hermitian form through successive edges.
In its equivalent column-form convention, `M* G M = w G`, the transferred
form is `G' = (I tensor G_edge) rho(G)`. Trace transitivity retains the
entire form instead of resetting it at each step.

## Obtaining the desired entries

In coordinates with a star-preserving action `rho(bar(z)) = rho(z)*`, a
fixed frame `F` with `F F* = c I_d` gives blocks `B(z) = rho(z) F` and

\[
(B(H)B(H)^*)_{ij}
=c\,\rho\!\left(\sum_k H_{ik}\overline{H_{jk}}\right)
=cw\,\delta_{ij}I_d.
\]

For a dense Butson source, `w=n`. If `c=d` and every used block contains
lower-field roots, the output is `BH(dn,q)`; signs give a real Hadamard matrix.
Thus flatness is the remaining condition once the scalar action and frame
are certified. Test only the actual used alphabet, while retaining source
gauges and any factor labels on which the map depends.

In a general trace basis the corresponding blocks are `C rho(z) D`, with
`C* C = G` and `D D* = c G^-1`. These factors must exist over the declared
field and pass the entry test. A numerical square root of `G` is insufficient.

Zeros in a weighing or conference source remain zero under a scalar-linear
map. If the resulting entries lie in `{0} union mu_q`, the output is a
weighing design of order `dn` and weight `cw`; obtaining a dense Hadamard
matrix also requires `cw=dn`. A weighted intermediate is already an exact
design even when it has no certified flat endpoint.

## Established reductions and current evidence

A new real order has now been constructed:
**[H7572](HADAMARD_7572.md)**. Exact periodic T-matrices of order 631 give
sign frames `P,Q` of width 2524, with `PP^T=QQ^T=2524I` and
`PQ^T+QP^T=2524I`. The rational action `R=QP^T/2524`, its conjugation and
trace form, and the source palette `{zeta_6,-1}` are recorded explicitly.
Restricting the verified order-3 source gives blocks `Q,-P` and a literal
sign endpoint with full integer Gram `7572I`. The matrix, standalone replay,
independent audits and dated primary-literature comparison are retained in
the [completed construction](HADAMARD_7572.md).

The [active search analysis](HADAMARD_FRAME_SEARCH.md) records the uniform
finite-field parity parameterization and complete modular sign lifts.
It also supplies a new complete [periodic T127 witness](../artifacts/t-matrices-127/README.md),
resolving the T-matrix case explicitly listed as open in the dated survey.
Its H1524 endpoint is a known-order control.
Independent prime repeats supply verified H3252 and H6564 controls.
The [repeat-target note](SCALAR_REDUCTION_REPEAT_TARGETS.md) distinguishes
these known orders, the new frame-reuse orders 17668, 22716 and 32812, and
the larger targets whose sign solver is still incomplete.
The [dated higher-order records audit](HADAMARD_HIGHER_ORDER_RECORDS.md)
states the author-knowledge status convention and records the prior
constructions located. The [completed unreal triplings](TRIPLED_UNREAL_SCALAR_ENDPOINTS.md)
also give new orders 98436 and 295308 and the known order 885924, each with
an exact endpoint recipe.

The complete T631 orthogonal design also admits weighted sources over
`Q[C_w]`, a product of cyclotomic number fields. Its entries are signed
copies of four verified Williamson group polynomials `f_j`, with
`sum_j f_j*conjugate(f_j)=4w`. The source order is 2524 and its Gram
scale is `2524w`. In the inverse group basis, the normalized regular trace
`tau(a)=[identity]a` gives metric `I_w`; multiplication conjugation becomes
transpose. Every used scalar image is a literal sign circulant. Restriction
therefore gives a sign matrix of order and Gram scale `2524w`.

The [27 literal Williamson inputs and prior-order checks](WILLIAMSON_ORDER_RECORDS.md)
add new orders **42908, 58052, 73196, 78244, 93388, 108532, 128724 and
153964**. The [five constructive Spence families](SPENCE_SCALAR_FAMILIES.md)
give literal endpoints for every integer `r>=0`, of orders
`17668*13^r`, `22716*17^r`, `78244*61^r`, `128724*101^r` and
`153964*121^r`. Their first expanded members add new orders **229684,
386172, 4772884 and 13001124**; the fifth, 18629644, is known by Paley II.
Individual finite orders carry the recorded status labels. The main
campaign now has **38 completed real orders, 18 new and 20 known recovered
orders**; the new periodic T127 ingredient is counted separately.

The [lower-order campaign](LOW_ORDER_SEARCH.md) seeks further periodic
T-matrices, starting with length 223 and Hadamard order 2676. Length 167
would supply orders 4676/8684, and length 179 would supply 7876/9308.
These are exact construction targets awaiting a completed finite witness.

The [projective flat-descent construction](PROJECTIVE_FLAT_DESCENT.md) now
supplies explicit polynomial identities guaranteeing endpoint flatness
through nonflat scalar restrictions. It constructs five pairwise inequivalent
real H768 matrices with binary-code dimensions 17 through 21, 18 `BH(42,6)` matrices
in at least four monomial classes, a one-parameter order-42 family, and a
36-parameter order-294 family with five pairwise inequivalent sixth-root
specializations. The small action/frame seeds are published, and these are
known-order controls.

The [Ó Catháin–Swartz morphism](https://arxiv.org/abs/1904.10771) gives
`BH(n,pq) -> BH(pn,q)` when `p` is prime and divides `q`. Iteration gives
`BH(n,2^a) -> H(2^(a-1)n)`. The existing executor replays the genuine
known-order control `BH(1000,8) -> BH(2000,4) -> H4000` with literal entry
access and exact Gram proofs.

Removing the last odd prime from a field can leave a weighted intermediate.
Moreover, no uniform rational-linear block map can send an alphabet with an
odd-coefficient-sum integral vanishing relation entirely to signs: modulo two,
the supposed relation would be an odd sum of all-one blocks. In particular,
a full odd-prime polygon cannot have such a sign map. This is an obstruction
to that particular class of maps, with that alphabet.

The [saved alphabet audit](SCALAR_RESTRICTION_OBSTRUCTIONS.md) excludes the
two fixed character-source alphabets and retains 18 other registry sources
by this gate. Retention is not an existence proof for a sign endpoint.
The completed [weighted inventories](WEIGHTED_SCALAR_TREE.md) have 472
coordinate matrices from those 20 inputs; all 236 rational endpoints contain
zeros. The [broader source and descent runs](BROAD_SCALAR_RESTRICTION_RUN.md)
also produced exact matrices, without a previously unknown real order.
These are finite saved runs, not exhaustion of scalar reduction or of bases.

Use [the active engines and replay commands](../experiments/scalar_lift/README.md)
to work forward from verified sources. Preserve metrics through general
restrictions, then certify any proposed endpoint frame and its actual
entries. The [retired branches](RESEARCH_HISTORY.md) are recoverable when a
specific earlier calculation is needed.

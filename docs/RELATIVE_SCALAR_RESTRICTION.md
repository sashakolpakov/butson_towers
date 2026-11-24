# Scalar restriction with the exact trace metric

`experiments/scalar_lift/relative_scalar_restriction.py` implements
\(L=\mathbb Q(\zeta_Q)\) over \(K=\mathbb Q(\zeta_q)\) for \(q\mid Q\),
including the removal of the last odd prime factors. It produces exact
coefficient matrices, not a claim that those matrices are flat.

Set \(d=\varphi(Q)/\varphi(q)\). The canonical relative basis is
\(1,\zeta_Q,\ldots,\zeta_Q^{d-1}\). Its combined rational basis has
columns \(\zeta_Q^{a+(Q/q)b}\), ordered by \((a,b)\), where
\(0\le a<d\), \(0\le b<\varphi(q)\). An exact rational basis matrix
\(E\) and its inverse recover all coefficients of all \(Q\) root powers.
The complete action table is independently compared with rational
multiplication conjugated by \(E\).

The metric is the **unnormalized** relative trace matrix

\[
 G_{ab}=\operatorname{Tr}_{L/K}(\zeta_Q^{b-a}).
\]

It is Hermitian and positive definite: under every complex embedding it
is the Gram of the relative embeddings of the basis. The executable
checks its exact inverse and the identity

\[
 \rho(z)^*G=G\rho(\bar z).
\]

Multiplication and the generator identity prove this for every scalar.
Thus a complete source \(HH^*=nI\) restricts to an actual matrix with

\[
 \rho(H)(I_n\otimes G^{-1})\rho(H)^*
 =n(I_n\otimes G^{-1}).
\]

For example, restricting sixth roots to the rationals gives

\[
 \rho(\zeta_6)=\begin{pmatrix}0&-1\\1&1\end{pmatrix},\qquad
 G=\begin{pmatrix}2&1\\1&2\end{pmatrix},\qquad
 G^{-1}=\tfrac13\begin{pmatrix}2&-1\\-1&2\end{pmatrix}.
\]

The entries are integers here, but include zeros and do not have the
ordinary Hadamard Gram. A separately verified compatible frame is still
required to obtain signs.

The tested edges 30→10, 30→6, and 30→2 have widths 2, 4, and 8 and
inverse-metric denominators 3, 5, and 15. The edges 120→40, 120→24,
and 120→8 have the same widths and denominators.

`edge(Q,q)` returns actions, metric, inverse metric, and exact basis
maps. `RationalArray` exposes integer `numerators`, a common positive
`denominator`, and `to_json()`. `trace_exponent(e)` gives the actual
relative trace of \(\zeta_Q^e\), without dividing by the relative degree.
`coordinates(e)` gives all relative-basis coefficients of that root.

`restriction_in_basis(Q,q,exponents)` also accepts a tower's actual
power basis, starting with exponent 0. Different bases retain their own
literal actions and metric. A canonical comparison
\(C=E_{\mathrm{canonical}}^{-1}E_{\mathrm{path}}\) proves the action
conjugacy and metric congruence; it does not justify discarding a
coordinate-distinct path when testing flatness.

The labels q=1 and q=2 both describe the rational field, but their
declared root labels are preserved. Identity restrictions are supported.
Exact integer bounds guard the vectorized coefficient products; no
floating-point arithmetic is used. Tests include the complete weighted
Gram of an actual order-3 source and an independent path-basis
congruence check.

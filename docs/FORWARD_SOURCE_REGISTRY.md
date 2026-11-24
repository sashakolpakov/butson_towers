# Verified inputs for forward restriction

`experiments/scalar_lift/forward_source_registry.py` supplies 20 complete
starting matrices for the [forward engine](FORWARD_SCALAR_DESCENT.md).
It contains no unfinished order-1499 or order-2999 source.

| Family | Source orders | Root phases | Proof |
|---|---|---|---|
| Conference | 1000, 2000 | 8, 16 | Actual skew-conference identities, exact block Gram, existing source hash |
| Shifted sixth-root sources | 3, 13 | 6 | Complete exact Gram; global shift −1 gives exponents in {0,1,3,4} |
| Mixed cubic/fifth | 15, 35 | 120 | Complete exact cyclotomic Gram |
| Restricted tripling from 13 | 13, 39, 117, 351, 1053 | 24 | Exact seed Gram and coefficient equations, verified tensor recurrence |
| Fifth chirp | 5 | 20 | Complete exact cyclotomic Gram |
| Fifth tensor | 40 | 40, 80 | Complete exact cyclotomic Gram |
| Forward fifth descendants | 80, 160 | 40, 20 | Explicit repeated-prime descent, independent complete Gram |
| Character sources | 5329, 58081 | 24, 120 | Exact character correlations and full group-ring Gram identity |

`sources()` returns small JSON descriptors without loading matrices.
`load_source(descriptor_or_id)` checks the canonical descriptor and
returns a cached `VerifiedSource` with order, phase, source hash, exact
used alphabet, proof, and `entry(i,j)` using ordinary linear indices.
Dense exponent arrays are read-only and are limited to order 2000.

The IDs `unreal3-q6` and `unreal13-q6` use a global exponent shift of
−1 from the saved unreal sixth-root inputs. Add 1 modulo 6 to recover
the original inputs for constructions expecting that convention. The
registry hashes bind the shifted arrays, and both shifts are recorded
explicitly in their proofs.

The character matrices remain implicit. Their hashes bind the complete
small convolution coefficient table, the phase, dimensions, and the
linear-index convention. They are explicitly **not** presented as
hashes of unallocated dense matrices. The exact group-ring proof covers
every entry of each full Gram.

For the tripling family, write an unreal sixth-root source as
\((A+i\sqrt3 B)/2\), with \(A,B\) sign matrices. Its complete Gram is
equivalent to \(AA^T+3BB^T=4nI\) and \(AB^T=BA^T\). With the fixed
skew conference \(C_3\), the implemented recurrence is

\[
 A'=B\otimes J_3,\qquad B'=A\otimes I_3+B\otimes C_3.
\]

The identities \(C_3C_3^T=3I-J\), \(C_3J=JC_3=0\), and
\(J^2=3J\) prove the next full Gram. Every literal tensor expansion and
root reconstruction is checked. This avoids substituting an expensive
large direct Gram for an exact construction proof.

These entries are starting substrates, not novelty certificates.
Conference inputs give known real orders 4000, 8000, and 16000 under
dyadic descent. The other available terminal lifts also produce known
control orders. Character inputs offer genuinely large complex sources,
but a valid final sign frame is still required. No uniform sign lift is
assumed for their entire odd-factor alphabet.

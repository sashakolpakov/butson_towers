# Exact exclusion of eight entire odd-denominator basis families

The [executed certificate](../certificates/scalar-lift/odd-denominator-basis-audit.json)
closes a complete, non-dyadic family of mixed scalar bases for each of
the eight supplied published frames. This is an exact finite-field
argument over **all `2^d` sign choices per frame**, not a sampled search.
No new Hadamard matrix was found or excluded outside these families.

Fix a published sign Hadamard F of order `d=4r`, with r an odd prime.
For every sign diagonal E, define the rational orthogonal matrix

`Q_E = F E F^T/d`.

Unlike Walsh-only bases, this family can introduce the necessary odd
denominator r and can have large rank of `Q_E-I`. It therefore tests a
different family from the previously excluded dyadic and low-rank
ansatzes. It uses only already verified F as basis data.

Let `T_e=2rho0(zeta12^e)`, an integer local quaternion action tensored
with the period identity. The used two-basis scalar block is

`G_e=Q_E rho0(zeta12^e) Q_E^T F = F E K_e E/(2d)`,

where `K_e=F^T T_e F`. A necessary condition for integer, hence sign,
entries is

`F E K_e=0 mod r`.

All four non-two-power root representatives `e=1,2,4,5` are covered;
their negatives give exponents 7,8,10,11. No sixth-root-only test is
substituted for the genuine twelfth-root cases.

## The kernel-support argument closes six frames

Modulo r, the exact computation gives `rank(F)=d/2` and `F F^T=0`.
Always `image(K_e) subset image(F^T)=kernel(F)`. If
`rank(K_e)=d/2`, the image is the entire kernel. The necessary condition
then says E must preserve `kernel(F)`.

After putting F in reduced row echelon form `[I A]` and permuting its
columns accordingly, that condition is exactly

`epsilon_pivot_i = epsilon_free_j whenever A[i,j] != 0`.

The connected components of this bipartite support graph give all
possible diagonal sign stabilizers. For each frame in the following
table the graph is connected. Thus E is globally `+I` or `-I`, giving
`Q_E=+I` or `-I`; conjugation returns the already failing baseline.

| Published frame | r | rank(F) | rank(K1) | rank(K2) | Support components |
|---|---:|---:|---:|---:|---:|
| 668, width-twelve border | 167 | 334 | 334 | 334 | 1 |
| 764, first input | 191 | 382 | 382 | 382 | 1 |
| 764, second input | 191 | 382 | 382 | 382 | 1 |
| 892 | 223 | 446 | 446 | 446 | 1 |
| 1132 | 283 | 566 | 566 | 566 | 1 |
| 1964 | 491 | 982 | 982 | 982 | 1 |

For these full-image cases, `K4=K2-2dI` has the same reduction as K2.
Also `K5=-K1^T` has the same maximal rank and its image is again the
entire kernel. This transfers the support-graph exclusion to both
remaining representatives. Every baseline is dyadic; the independent
[Hadamard parity argument](https://github.com/sashakolpakov/complex-hadamard/blob/976cede671f40c810d593c8d5815e26e11a03ed6/docs/MIXED_RESTRICTION_BASIS_AUDIT.md) excludes its
odd-part primitive root actions at these widths.

## Small linear systems close both rank-deficient frames

The remaining frames genuinely have deficient K images:

| Frame | rank(F) | rank(K1) | rank(K2) |
|---|---:|---:|---:|
| 668, width-four border | 334 | 330 | 332 |
| 716 | 358 | 350 | 356 |

The full-image argument is not applied to these cases. Instead,
`F E K_e=0 mod r` is linear in the signs epsilon. Its equations from
F rows 0,1,2 are the stacked matrix

`[K_e^T diag(F_0); K_e^T diag(F_1); K_e^T diag(F_2)] epsilon=0`.

For each of `e=1,2,4,5`, and both frames, this system has nullity three.
A sign solution is determined by its three free-coordinate signs, so
there are only **eight assignments to check exhaustively**. The
certificate records every assignment that actually gives a full sign
vector in the kernel.

For e=1 only the two globally constant vectors remain. For e=2,4,5
there are six sign vectors: the constants and four vectors with only
two minority signs. The transpose-root case e=5 was computed explicitly;
its small kernel was not inferred from K1 merely having the same rank.

To exclude the nonconstant survivors, put `A0=rho0(zeta12^e)F`, an
integer matrix. If G_e were signs, then

`rank_F2(A0-J) <= rank_Q(G_e-A0) <= 2 rank_Q(Q_E-sigma I)`.

Here sigma can be chosen as the majority sign in E. Since Q_E is
similar to E, that last rank is exactly the number of minority signs,
zero or two. The actual baseline binary ranks are:

| Frame | e=1 | e=2 | e=4 | e=5 |
|---|---:|---:|---:|---:|
| 668, width-four border | 167 | 166 | 167 | 166 |
| 716 | 179 | 178 | 179 | 178 |

Every value exceeds four, so all surviving sign vectors are excluded.
This completes the entire sign-diagonal family for both frames without
enumerating `2^668` or `2^716` cases.

## Native computation and replay

The audit uses FLINT exact modular matrices. It verifies the complete
column factorization `F=V C` from the RREF, where C has full row rank.
It then uses

`rank(K_e)=rank(V^T T_e V)`

to reduce the large K rank computation to a matrix of order `d/2`.
It saves the source frame proof/hash, RREF hash and pivots, connected
component labels, a nonzero spanning tree, compressed matrix hashes,
exact ranks, small-system hashes, every surviving sign point, and the
final binary-rank exclusions.

A separate four-coordinate finite-field control compares the graph,
kernel-image, and direct stabilizer conditions on all 16 sign choices.
This validates the linear-algebra criterion; it is not claimed as a
new Hadamard construction. Tests also replay a complete deficient-frame
audit and check that no rank deficiency is silently treated as full rank.

```sh
PYTHONPATH=src:. ./.native-env/bin/python -m experiments.scalar_lift.odd_denominator_basis_audit --check certificates/scalar-lift/odd-denominator-basis-audit.json --progress
```

All eight families are closed, with no unfinished audit. Rational
mixed bases outside `F diag(epsilon) F^T/d` remain outside this result.

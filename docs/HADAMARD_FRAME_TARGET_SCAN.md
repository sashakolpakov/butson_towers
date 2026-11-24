# Compact scalar-frame target scan

The [finite parameter certificate](../certificates/scalar-lift/hadamard-frame-target-scan.json)
lists every positive DC configuration in the following class: prime
`7 <= p <= 20000`, an odd multiplier order `m` divisible by three, the
same fixed-head/three-cycle channel action, and at most **50 nonorigin
orbits**. It contains **28 charts**, including small controls. The largest
surviving prime is **1951**. A surviving chart is a necessary-condition
candidate, not an endpoint or a Hadamard nonexistence result.

## Why the finite arithmetic scan is complete

Write the common row sums as

\[
a=1+mh,\qquad b=(m/3)k,\qquad a^2+3b^2=p,
\qquad N=(p-1)/m.
\]

Here `k != 0`, since a prime cannot be a square. The integers `p,m` are
odd, so `N` is even and `h,k` have the same parity. Expanding the norm gives

\[
N=2h+m(h^2+k^2/3).
\]

For these integers and `m>=3`,

\[
3mh^2+mk^2+6h+6-4m\geq0.
\]

For `h=-1,0,1`, this follows from `k^2>=1,4,1` respectively. For
`|h|>=2`, use `m>=3` in the quadratic lower bound. Consequently

\[
m\leq\frac{3(N+2)}4,\qquad
p\leq1+\frac{3N(N+2)}4.
\]

Thus **N<=50 implies p<=1951**, even without the original prime cutoff.
Equality is realized arithmetically by `p1951,m39,h=-1,k=1`, with row sums
`(-38,13,13,13)`. Scanning larger primes cannot produce a smaller chart in
this particular family. This bound does not concern other twists, groups,
or unrestricted Hadamard matrices.

## Knownness gates and remaining research targets

Turyn's theorem supplies Williamson matrices of order `p` when `2p-1` is
a prime power. The prime-power parameter is automatically `1 mod4` for
odd `p`. Cooper–Wallis then combines these with known T-matrices of orders
`3,7,9,13`, giving every corresponding order `4*p*w` already. The
certificate records the exact factorization and this gate separately for
every chart. The hypotheses are stated in the
[2026 survey, §2.1](https://arxiv.org/html/2610.27175v2#S2.SS1)
and its [Cooper–Wallis construction](https://arxiv.org/html/2610.27175v2#S4).

This eliminates the small-source new-order rationale for the compact
271, 307, 331 and 547 charts, though they remain useful construction
controls. Further explicit Paley and Yamada checks give the current scope:

| Chart | Nonorigin orbits | Parity field | Remaining small-source order with unverified prior status |
|---|---:|---|---:|
| p919, m27 | 34 | GF(2^17) | 25732, source7 |
| p397, m9 | 44 | GF(2^22) | 14292, source9 |
| p1453, m33 | 44 | GF(2^22) | 75556, source13 |
| p1951, m39 | 50 | GF(2^25) arithmetic candidate | Higher-order audit still required |

**Unverified prior status is not a novelty claim.** The currently checked
unresolved-order inventories stop below these three orders. The concrete
primary record for each must be assessed after a completed endpoint is
available; absence from a lower-order table proves nothing about it.

The explicitly unresolved small-source targets at p223 and p283 need
74 and 94 nonorigin orbits in this twist, respectively. They exceed the
present single-word sign representation and practical finite sign filter.
Even multiplier orders acting trivially on colors at inversion have a
modulo-two obstruction; a transposition at inversion forces the special
DC form `p=1+2s^2`. In particular, using order6 to shrink the p223 chart
does not evade that obstruction.

## Ingredient orders are a separate opportunity

Ranking only new real matrix orders misses open reusable ingredients.
The [2026 survey, Table3](https://arxiv.org/html/2610.27175v2#S4.T3)
explicitly lists **127** as its smallest open T-matrix order, even though
the immediate Hadamard orders from its small-source products are known.
The p127/m9 chart has only 14 nonorigin orbits and field GF(128).

A bounded pilot produced exact T127 in 0.002314 seconds of native search,
with states, all four literal sequences, a direct circulant Gram audit and
the full H1524 scalar endpoint retained in
[its certificate](../certificates/scalar-lift/hadamard-frame-p127.json).
The later-record/priority audit is separate from this exact existence proof.

The [active engine guide](HADAMARD_FRAME_SEARCH.md) describes the verified
degree17/22 ports and the conditional sampler for the headless p1453
stratum. No finite timeout establishes nonexistence, and no parity map
alone establishes a completed frame.

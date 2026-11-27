# Compact scalar-frame target scan

All construction searches are paused. The next
[selection](LOW_ORDER_SEARCH.md) matches unresolved endpoint orders below
50000 with admissible multiplier symmetries, then compares their reduced
sign systems before allocating search time. Earlier runtime results are
retained privately.

The [finite parameter certificate](../certificates/scalar-lift/hadamard-frame-target-scan.json)
lists every positive DC configuration in the following class: prime
`7 <= p <= 20000`, an odd multiplier order `m` divisible by three, the
same fixed-head/three-cycle channel action, and at most **50 nonorigin
orbits**. It contains **28 charts**, including small controls. The largest
surviving prime is **1951**. A surviving chart is a necessary-condition
candidate, not an endpoint or a Hadamard nonexistence result.

**26 November 2026 update.** The
[expanded 87-chart multiplier scan](../certificates/scalar-lift/multiplier-target-scan-20261126.json)
extends to 160 nonorigin orbits, with the proved prime bound 19441. The
50-orbit bound below describes the original scan. The p1951 parity map is
now verified; p397 has been removed
from the missing-order shortlist by a positive classical construction.

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

## Searching at the scale of the successful case

The completed p631 search used 30 nonzero orbits and 13 free sign bits
in each normalized sign system. Its recorded winning run took 2347.74
seconds. The p271 and p547 controls took 0.037877 and 25.2677 seconds of
native search; these timings exclude subsequent endpoint verification.
The [run record](HADAMARD_FRAME_SEARCH.md) gives the seeds and counters.
A successful run time also depends on how frequently the sampled coordinate
patterns admit exact signs, which varies between charts.

The arithmetic scanner exhausts this particular three-cycle action at
N<=30: the prime bound is 721 and its two independent enumerations agree
on exactly 16 charts, one each at

```
7, 13, 19, 31, 37, 43, 61, 67, 73, 79, 109, 127, 271, 331, 547, 631.
```

For every listed prime other than 631, Williamson-type matrices of that
prime order are already known. Turyn's family covers
7,13,19,31,37,61,79,271,331,547. The remaining five are covered by
[Kharaghani–Mohammadian–Tayfeh-Rezaie–Zaitsev, Definition 1.1 and §2.3](https://arxiv.org/html/2605.08661#S2.SS3):
near-Williamson matrices give Williamson-type matrices after reversing the
columns of the first matrix. Classical T-sequences exist at each of our
27 retained odd widths, all at most 63, by the base-sequence construction
in [Đoković, (2.3) and §3](https://arxiv.org/html/0912.5091#S3).
Cooper–Wallis with these two ingredients therefore already supplies every
corresponding order 4*p*w. This comparison concerns existence at the
output orders; it does not classify the T-inputs or Hadamard equivalence
classes.

Consequently finding another new order with as few orbit variables as 631
requires broadening this coordinate action, the underlying group, or the
source inputs. A larger multiplier within the same action can still leave
more independent orbits and a larger sign system.

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
controls. The p397/m9 chart is also a control: Paley gives a symmetric
conference matrix of order 198, and Seberry–Yamada's
[Corollary 29(ii), p.120](https://documents.uow.edu.au/~jennie/WEBPDF/139_1990.pdf#page=24)
gives symmetric Williamson-type matrices of order 397. Cooper–Wallis,
[Theorem 3 and Corollary 4, p.274](https://doi.org/10.1017/S0004972700045019),
then gives `4*397*w` for every odd `w<=19`. In particular, **H14292 is
known**. The valid GF(2^22) parity map remains an engine control.

The remaining compact targets in the 26 November check are:

| Chart | Nonorigin orbits | Parity field | Candidate orders `4*p*w`, with Williamson width `w` |
|---|---:|---|---:|
| p919, m27 | 34 | GF(2^17) | 25732, w7; 40436, w11 |
| p1453, m33 | 44 | GF(2^22) | 63932, w11; 75556, w13 |
| p1951, m39 | 50 | GF(2^25) | 23412, w3; 54628, w7; 101452, w13 |

No prior construction was located for these targets. Each requires a
completed T-input and the recorded Williamson source. Width 11 uses the
four-variable orthogonal design and weighted group-algebra restriction;
it is not an additional admitted sixth-root source for a P/Q frame alone.

The [p1951 certificate](../certificates/scalar-lift/hadamard-frame-p1951.json)
verifies the coordinate-parity map for `mu=1127` and row sums
`(-38,13,13,13)`, including every expanded nonzero shift. The generic
normalized affine sign-lift estimate is 23 free bits; each sampled axis
word still needs exact elimination, counts, mod-8 and integer tests.
This supplies a search route, not a completed T-matrix or endpoint.

The larger charts now have verified parity maps:

| Chart | Coordinate sums | Coefficient probes | Expanded parity comparisons | Free sign bits in four checked patterns |
|---|---|---:|---:|---:|
| [p3571, m51, mu292](../certificates/scalar-lift/hadamard-frame-p3571.json) | 52,17,17,17 | 9,871 | 35,239,470 | 33 |
| [p4219, m57, mu24](../certificates/scalar-lift/hadamard-frame-p4219.json) | -56,19,19,19 | 11,027 | 46,511,886 | 35 |

Each expanded comparison uses an actual length-p coordinate assignment
and one nonzero cyclic shift. Every constant, linear and quadratic probe
was checked at every nonzero shift. Both binary maps have full rank and
both field moduli are irreducible. These certify the parity reductions;
neither chart has a signed T-input. The modern arbitrary-size sign masks
are required: these charts exceed the legacy engine's 63-orbit limit.
Their measured sign dimensions make them substantially more expensive
than the completed p631 search, despite the larger multipliers. The
[higher-order audit](HADAMARD_HIGHER_ORDER_RECORDS.md) records the positive
known-order exclusions separately from these incomplete searches.

The explicitly unresolved small-source targets at p223 and p283 need
74 and 94 nonorigin orbits in this twist, respectively. They require the
modern arbitrary-size orbit masks and substantially larger sign systems
than the successful p631 run.
Even-order multiplier subgroups with signed-permutation coordinate actions
are excluded at every prime length p>3. Pairing a nonzero correlation
under x -> -x-s first forces exactly two coordinate sequences away from
zero, exchanged by inversion. Their circulant matrix would have weight
(p-1)/2=s^2 and zero diagonal. A prime divisor of s is self-conjugate
modulo p, contradicting that circulant Gram identity. The
[complete proof](EVEN_MULTIPLIER_OBSTRUCTION.md) includes the genuine
p=3 exception. In particular, order 6 cannot shrink the p223 chart.

For odd multiplier order, a signed permutation fixing e0 acts on the
remaining three coordinates. Their signed-permutation group has order 48,
so its odd cyclic subgroups have order 1 or 3. An order-three action becomes
the ordinary three-cycle after changing coordinate signs. Thus identity
and three-cycle coordinate actions exhaust this particular symmetry model
at prime lengths greater than 3; composite lengths and more general
coordinate actions remain outside the obstruction.

## A different multiplier action at p463

The identity coordinate action requires b(mu*x)=b(x): each of the four
coordinate sequences is constant, including its sign, on multiplication
orbits. This allows odd multiplier orders that are not divisible by three.
Normalize b(0)=e0 and let a_i be the signed number of nonzero orbits
assigned to coordinate i. The expanded row sums and their necessary norm
equation are

\[
(1+ma_0,ma_1,ma_2,ma_3),\qquad
N=2a_0+m\sum_{i=0}^3 a_i^2,\qquad N=(p-1)/m.
\]

For p463, multiplication by mu158 has order 11. The choice
(a0,a1,a2,a3)=(-1,1,1,1) gives 42 nonzero orbits and row sums
(-10,11,11,11), whose squared norm is 463. A completed T463 would give
**H9260=4*463*5** using the retained Williamson input of width five.
Order 9260 is explicitly unresolved in the
[author-maintained frontier](https://github.com/vatskishan03/Hadamard/blob/main/FRONTIER.md),
checked 26 November 2026.

This case needs its own orbit coefficients and four separate sign-count
conditions. The three-cycle map and its two-count sign system do not apply.
It remains a construction target awaiting a complete signed witness.

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

# Repeating the H7572 scalar reduction

The width-2524 scalar frame gives three further **new real Hadamard orders:
17668, 22716 and 32812**. Unreal-source triplings add **98436 and 295308**;
the completed order 885924 is known. Status labels follow the
[author-knowledge convention and dated literature ledger](HADAMARD_HIGHER_ORDER_RECORDS.md).

The new [executor](../experiments/scalar_lift/reuse_hadamard_7572_frame.py)
and [exact certificate](../certificates/scalar-lift/hadamard-7572-frame-reuse.json)
give complete sources, gauges, entry access, and algebraic endpoint Grams.
The literal entry recipes specify the matrices without large dense allocations.

## Completed reuse of the same frame

Let P,Q be the certified sign matrices of width d=2524 from
[H7572](HADAMARD_7572.md). They obey

    PP^T=QQ^T=dI,       PQ^T+QP^T=dI.

Thus R=QP^T/d is orthogonal, R+R^T=I, and R^2-R+I=0. The field action is
Q(zeta6)/Q, with degree two and module rank 1262. Conjugation is transpose.
The canonical trace basis is (1,zeta6), with metric [[2,1],[1,2]]. In the
actual rational coordinates, the Hermitian form is

    h(u,v)=u^T v/2-(u^T S v)*sqrt(-3)/6,       S=2R-I,
    Tr(h(u,v))=u^T v.

The four flat entry images are

    1 -> P,       zeta6 -> Q,       -1 -> -P,       -zeta6 -> -Q.

The unused roots zeta6^2 and zeta6^5 map to Q-P and P-Q, which contain
zeros and entries of magnitude two. Sources with the full six-root palette
are therefore not silently admitted.

| Complete active source | Source scale | Actual used exponents modulo 6 | Frame scale | Completed endpoint order and Gram scale |
|---|---:|---|---:|---:|
| `catalogue-p7-q6` | 7 | 0,1,3,4 | 2524 | 17668 |
| `aalto-9-6-000` | 9 | 0,1,3,4 | 2524 | 22716 |
| `unreal13-q6` | 13 | 0,1,3,4 | 2524 | 32812 |

For each source, the executor independently checks its complete exact row
and column Grams in Z[zeta6], reconstructs its actual row/column/global
gauge from the original literal data, and records both raw and final hashes.
The order-13 source uses its recorded unreal gauge followed by global
shift -1. The order-9 source uses its complete root-gauge witness and the
recorded global anchor. Its catalogue class 004 is also admitted, but a
second complex input does not establish a different real equivalence class.
Exactly classes 000 and 004 among the 17 active BH(9,6) representatives pass
the original sixth-root gauge admission test.

For B(z)=rho(z)P, the full endpoint Gram is proved by the exact identity

    sum_k B(H_ik) B(H_jk)^T
       =d*rho(sum_k H_ik conjugate(H_jk))=d*n*delta_ij*I.

The source scale n and frame scale d are separately recorded; dn equals
the expanded order only after the sign-entry test. The executor replays all
full P/Q/cross integer identities once and checks literal endpoint row
controls. Its public `endpoint(source_id).entry(i,j)` and `.row(i)` specify
every matrix entry without allocating the full endpoint.

```sh
PYTHONPATH=src:. .native-env/bin/python -m experiments.scalar_lift.reuse_hadamard_7572_frame --check certificates/scalar-lift/hadamard-7572-frame-reuse.json
PYTHONPATH=src:. .native-env/bin/python -m unittest tests.test_hadamard_7572_frame_reuse
```

## What the primary records establish

The [2026 survey, Table 5](https://arxiv.org/html/2610.27175v2#S4.T5) and
[current construction frontier](https://raw.githubusercontent.com/vatskishan03/Hadamard/main/FRONTIER.md)
explicitly listed 7572 as unresolved before this construction. They cover
orders only to 10000 and 11996 respectively. The
[research repository's higher-order constructions](https://github.com/vatskishan03/Hadamard/blob/d885022144ba6bf4efa8315e14871a3736947527/hadamard/README.md#orders-beyond-the-tables)
and the [higher-order records audit](HADAMARD_HIGHER_ORDER_RECORDS.md)
record the subsequent construction checks and catalog coverage. No prior
construction of 17668, 22716 or 32812 was located in those checks.

A focused check on 24 November 2026 rules out the direct Paley prime-power
tests and the prime-power hypothesis in Yamada's order `4(q+2)` construction
for these three numbers:

| Endpoint order N | N-1 | N/2-1 | N/4-2 |
|---:|---|---|---|
| 17668 | 3^2 * 13 * 151 | 11^2 * 73 | 5 * 883 |
| 22716 | 5 * 7 * 11 * 59 | 41 * 277 | 7 * 811 |
| 32812 | 3 * 10937 | 5 * 17 * 193 | 59 * 139 |

The hypotheses are stated in the
[survey's Theorem 2 and Table 4](https://arxiv.org/html/2610.27175v2#S4.T4)
and [Đoković's Theorems 2.2(b) and 2.3](https://www.m-hikari.com/imf-2011/61-64-2011/djokovicIMF61-64-2011.pdf#page=3).
None of these parameters is a prime power.
The [1987 primary table](https://documents.uow.edu.au/~jennie/WEBPDF/110_1987.pdf)
mentions `9*631=5679` among cases where `16*5679` was known and `8*5679`
was then unknown.

H15144 and all its successive doublings were already known: Cati and
Pasechnik's [Table 4](https://arxiv.org/pdf/2411.18897v2#page=18) records
1893(3), meaning H(2^3*1893). Doubling H7572 therefore produces no new order.
More generally, ordinary Kronecker products with real Hadamards are stated
as closure consequences, not independent new-order discoveries.

The [completed order-13 tripling reductions](TRIPLED_UNREAL_SCALAR_ENDPOINTS.md)
give new orders 98436 and 295308 and the known order 885924. Their phase-24
exponents are 4e+2; the global shift -6 followed by division by four gives
the admitted sixth-root exponents e-1. The executor replays the complete
inductive source proof, actual palettes and exact algebraic endpoint Grams
through source order 351. The next finite execution at source order 1053,
which would give H2657772, remains a target.

## Generalization of the successful parity coordinates

For an odd prime p, suppose s=ord_p(2) is odd and m is an odd multiple of
three dividing p-1. Use the same four-color twist c(mu*x)=omega*c(x), with
mu of order m. Put

    g=gcd(s,m),       e=s/g,       h=lcm(s,m),
    r=(p-1)/(2h),     N=(p-1)/m=2re.

The group generated by 4 and mu has order h and does not contain -1, so its
frequency orbits occur in r opposite pairs. The compatibility on each
orbit is a Frobenius eigenspace t*F_(2^(2e)). Conjugation is the 2^s
Frobenius, whose restriction to that subfield is the 2^e Frobenius.
The scalar t can be normalized to t*bar(t)=1: its norm lies in F_(2^e),
where the needed inverse square root exists uniquely.

Writing each Fourier pair as t*(a0+omega*a1), t*(b0+omega*b1), the parity
condition is independently

    a0*b1+a1*b0=1       over F_(2^e).

Thus the exact parity locus is a product of r copies of SL(2,2^e).
This is a parity parametrization, not a theorem that the subsequent sign
and integer equations have a solution. The DC condition must also pass:

    head^2+3*tail^2=p,
    head=1 (mod m),       tail=0 (mod m/3).

| p | m, mu | Row sums | Nonorigin orbits | Correlation classes | Parity factors | Role |
|---:|---|---|---:|---:|---|---|
| 223 | 3,183 | -14;3,3,3 | 74 | 37 | SL(2,2^37) | H2676 and H11596 remain explicitly unresolved; the active search uses arbitrary-size orbit masks |
| 271 | 15,90 | -14;5,5,5 | 18 | 9 | SL(2,2^9) | Independent repeat completed; H3252 already has Paley's construction with prime 3251 |
| 631 | 21,458 | 22;7,7,7 | 30 | 15 | SL(2,2^15) | Completed H7572 witness |
| 919 | 27,207 | -26;9,9,9 | 34 | 17 | SL(2,2^17) | Targets H25732/H40436 with Williamson widths 7/11; H11028 is already Paley with prime 11027 |
| 1951 | 39,1127 | -38;13,13,13 | 50 | 25 | SL(2,2^25) | Verified parity map; targets H23412/H54628/H101452 from sources 3/7/13 |

The width-11 target uses the complete T-input's four-variable orthogonal
design and the verified Williamson source over its group algebra. It does
not follow from admitting a new order-11 sixth-root source to P/Q alone.

The earlier arithmetic scans are retained in the ignored
`artifacts/scalar-frontier/repeat-targets` store; the current public scan
and p1951 map are linked below.
The two known-order labels use [Paley's primary construction as restated in
the 2026 survey, Theorem 2](https://arxiv.org/html/2610.27175v2#S2.Thm2),
not an absence argument. Under the odd-s, three-cycle hypotheses, p=223 is
the remaining direct small unresolved-order target; its 74-axis size makes
it a further solver task, not a quick rerun of the 631 engine.

**26 November 2026 update.** The
[expanded 87-chart scan](../certificates/scalar-lift/multiplier-target-scan-20261126.json)
extends the arithmetic shortlist. The
[p1951 certificate](../certificates/scalar-lift/hadamard-frame-p1951.json)
verifies the exact coordinate-parity map, including all expanded shifts.
Its generic normalized affine sign-lift estimate is 23 free bits; each
sampled axis word needs exact elimination, counts, mod-8 and integer tests.
No completed T1951 or endpoint is claimed.
The charts `p3571,m51` and `p4219,m57` now have verified
[degree-35](../certificates/scalar-lift/hadamard-frame-p3571.json) and
[degree-37](../certificates/scalar-lift/hadamard-frame-p4219.json) parity
maps. All coefficient probes pass at every expanded nonzero shift:
35,239,470 and 46,511,886 comparisons, respectively. Four sign preparations
per chart leave 33 and 35 free bits. No signed input has been found; these
larger charts are lower priorities than the smaller validated searches.

H14292, previously a p397/source-9 target, is **known**. Paley's symmetric
conference matrix of order 198 gives symmetric Williamson-type matrices
of order 397 by Seberry–Yamada,
[Corollary 29(ii), p.120](https://documents.uow.edu.au/~jennie/WEBPDF/139_1990.pdf#page=24).
Cooper–Wallis, [Theorem 3 and Corollary 4, p.274](https://doi.org/10.1017/S0004972700045019),
then settles `4*397*w` for every odd `w<=19`. The p397 map remains a valid
engine control. This positive construction supersedes the earlier lookup
gap in the [dated order audit](HADAMARD_HIGHER_ORDER_RECORDS.md).

The [generalized engine](HADAMARD_FRAME_SEARCH.md) has now repeated the
construction at p=271, with the [completed exact certificate](../certificates/scalar-lift/hadamard-frame-p271.json).
A bounded run found the states

```text
0 1 7 6 3 1 2 6 2 0 0 3 7 7 6 4 5 6 1
```

in about 0.038 seconds; full frame/cross and H3252 integer Grams passed.
The endpoint packed-sign SHA256 is
`599b87c00f75a0b29118443f8d5fcaab9073056452c55bc1071c173910aef4fd`.
This is a successful known-order method control at a different prime.
Its 36-column map was checked on every constant, singleton and
pair basis input (667 quadratic-coefficient checks), both binary inverse
compositions, and 24 full expanded periodic-parity controls.

## The even-order-of-two port at p=547

The even-s case has a different Fourier involution. Its concrete p=547
port is now exactly validated in the
[parity-map certificate](../certificates/scalar-lift/hadamard-frame-p547.json).
Here m=21, mu=169, ord_547(2)=546 and ord_547(4)=273. The chart has 26
nonorigin orbits and 13 correlation classes, with target row sums
(-20,7,7,7), head sign net -1 and tail sign net +1.

The ambient field is F_(2^546). The coefficient/root involution is the
2^273 Frobenius, which now sends zeta547 to its inverse. The eigenspace
condition

    z_(+/-1)^(2^26)=omega^-1*z_(+/-1)

gives t*F_(2^26); choose t with t*bar(t)=1. The parity equation is

    N(a)+N(b)=1,       N(x+y*omega)=x^2+x*y+y^2 over F_(2^13).

It is a sum of norms, not the alternating pairing used for p631. Nevertheless
an exact linear coordinate change turns it into the same determinant form.
For determinant coordinates (P,Q,U,V), set

    a0=P+U+Q,       a1=V+P,
    b0=U+Q,         b1=V+P+Q.

Then N(a)+N(b)=P*V+Q*U. The exported 52-column inverse color map already
contains this change, so the same SL(2,8192) sampler applies. All 1379
constant/singleton/pair coefficient checks, 104 inverse-basis compositions,
and 24 full 547-shift parity controls passed. This certifies the entire
parity identity; signs and exact integer correlations remain separate.

The p547 repeat is now completed. Seed 20261124548 found the states

```text
0 1 1 0 6 7 5 7 4 3 6 7 1 6 7 1 7 2 4 4 0 3 0 5 2 4 2
```

after 25.2677 seconds of native search. The complete frame, cross identity,
and actual H6564 integer Gram passed in a further 6.0158 seconds. Its packed
sign SHA256 is
`ff886685971143ef09c602b29e73463c1d023b43c4fef340747a5894ec8e0997`.
The [public p547 certificate](../certificates/scalar-lift/hadamard-frame-p547.json)
records the completed known-order repeat and actual states.

This gives a width-2188 frame, but all three immediate source-order targets
are already known. H6564 is Paley I (prime 6563), and H28444 is Paley II
(prime 14221). **H15316 is also known**, using the prime-power version of
the Paley-block product: 4*547-1=2187=3^7 is a prime power congruent to 3
modulo 8. Paley's H2188 can be arranged as sixteen circulant blocks of
size 547; the standard T-matrices of order seven then give

    4*547*7=15316.

The prime-power hypothesis is explicit in
[the 2026 survey, Table 4 and its following discussion](https://arxiv.org/html/2610.27175v2#S4.T4).
Its product statement cites Seberry–Yamada (1990), Theorem 4, and its
Paley-block statement cites Yamada (1991) and Schmidt (1999). The
[research repository's construction](https://github.com/vatskishan03/Hadamard/blob/d885022144ba6bf4efa8315e14871a3736947527/hadamard/PALEY_BLOCKS.md#construction)
gives the same sixteen-block/T-matrix argument. Thus this use is not
restricted to the case where 4*547-1 is prime.

The p547 witness is a successful independent prime repeat at a known real
Hadamard order.

# How the H7572 search succeeded, and what can be repeated

The completed [H7572 matrix](HADAMARD_7572.md) came from an exact periodic
T631 sequence, found by generating parity-compatible axes and then exhausting
a small affine space of signs. The active implementation is
[`hadamard_frame_search.py`](../experiments/scalar_lift/hadamard_frame_search.py)
and its [C++ kernel](../experiments/scalar_lift/hadamard_frame_search.cpp).
They regenerate their coefficient tables and need no ignored research files.
The original successful sources and the immutable matrix bundle are retained
unchanged.

## Why these variables were useful

The construction starts with the complete `BH(3,6)` having diagonal `-1`
and off-diagonal `zeta_6`. Orthogonality transfers through its scalar action;
the remaining task is a compatible flat frame. Four periodic sequences
`b_i` with exactly one signed nonzero per coordinate supply such a frame
when their summed autocorrelation is `631 delta_0`.

With the skew conference matrix `D4`, put `a=D4*b`. Disjoint supports ensure
that every entry of `a-b` and `a+b` is a sign. The exact correlation identities
then give `P=GS(a-b)`, `Q=GS(a+b)` with

\[
PP^T=QQ^T=2524I,\qquad PQ^T+QP^T=2524I.
\]

Thus `R=QP^T/2524` is the sixth-root scalar action and `F=P` is the frame.
The source scale 3 and frame scale 2524 give the final Gram scale/order 7572.
The field, conjugation, trace form, actual alphabets and full integer proof
are recorded in the [standalone certificate](../artifacts/hadamard-7572/certificate.json).

The multiplier `mu=458` has order 21 modulo 631. It fixes channel 0 and
cycles the other three. This replaces 631 four-channel entries by 31 signed
axes: origin plus 30 nonzero orbits. Origin and first head are positive;
the required row sums are `(22,7,7,7)`. There are 15 nonzero correlation
classes. The underlying twisted length-631 chart was proposed previously;
the completed witness and its exact endpoint are the result here.

## The search reductions

For fixed axes, write the signs as `s_i=1-2x_i`. Each correlation is a
weighted signed graph sum

\[
C_\ell(x)=C^{\rm abs}_\ell
       -2\sum_{i<j}w_{ij,\ell}[x_i\ne x_j].
\]

Diagonal terms are sign independent. In particular, odd `Cabs` cannot be
repaired by changing signs. Generating axes with all 15 unsigned correlations
even avoids spending time on that obstruction.

The 60 binary axis coordinates have an invertible binary change of basis to
`(p,q,b0,b1)` in `GF(2^15)^4`, using modulus `0xe7d5`. In these coordinates
the parity equations are

\[
p\,b1+q\,b0=1.
\]

This change of basis was obtained by exact finite-field Fourier evaluation
of the orbit sums at a primitive631st root. The multiplier covariance gives
linear eigenspace conditions on those evaluations. Pairing opposite
frequencies turns the correlation parity equations into a bilinear trace
condition; expressing its eigenspace in a subfield basis gives the displayed
determinant equation. This is the step that exposes the finite group,
instead of treating the parity equations as unrelated Boolean constraints.
The retained maps specify the resulting binary basis explicitly.

The public tool checks this identity completely: equality of all **1831**
constant, linear and quadratic coefficients, including invertibility of
both the 60-bit coordinate map and the 15-bit equation map. Channel equality
over four colors is quadratic in their two binary coordinates, so these
coefficient checks establish the identity on all `2^60` inputs.

Uniform rejection-free sampling of this parity locus is simple. Choose
`(p,q)` uniformly with at least one nonzero, and choose free `w`. If `p!=0`,
set `b0=w`, `b1=(1+q*w)/p`; otherwise set `b0=1/q`, `b1=w`.
Every determinant-one matrix has exactly one representation of this form.
The first-head gauge and optional head-count filter are applied afterward.

Modulo four, the correlation equations become linear:

\[
\sum_i \deg_{i,\ell}x_i=C^{\rm abs}_\ell/2\pmod2.
\]

The fixed signs and parities of the required head/tail negative counts are
included in the same binary solve. Typical rank is 18 among 31 variables,
leaving **13 free signs**, or 8192 affine sign assignments per axis word.
Exact negative counts are checked later; the linear count parities alone
are insufficient.

On this affine kernel, every `C_l` is divisible by four and `C_l/4 mod2`
is quadratic in the free signs. The implementation obtains its constant and
linear coefficients together using two bitplanes of weighted-cut counts
modulo four. An edge's odd weight contributes the upper outer product of
its free-sign flip mask to the quadratic coefficients. Thirteen variables
need only 78 quadratic bits. Six low variables fit in a 64-bit truth word;
seven high variables select 128 words. This tests all 8192 signs using
exact bit operations before exact DC counts and integer correlations.

For the length 547 port, the Fourier calculation initially gives a sum of two
norms `N(a)+N(b)=1`, where `N(x,y)=x^2+xy+y^2`. It reduces to the same sampler:
set `P=a 0+b0`, `Q=a 1+b1`, `U=b0+Q`, `V=b1+P+Q`. Direct expansion gives
`N(a)+N(b)=P*V+Q*U`. The inverse is
`(a 0,a 1,b0,b1)=(P+U+Q,V+P,U+Q,V+P+Q)`.
Composing the binary basis with this invertible linear map avoids a separate
norm-sphere search; its full coefficient identity is verified again.

This avoided dependence on a numerical objective or a nearby incomplete
matrix. It also made every positive candidate independently checkable.
It does not prove that every larger parameter choice has a small sign
kernel or a comparable finite-field parameterization.

## The successful finite run and an independent repeat

The winning run used seed **2026112416**, head counts **5,7,9**, and elapsed
time **2347.74 seconds**. Its saved counters were:

| Stage | Count |
|---|---:|
| Field parameters sampled | 273,324,655 |
| Axis words sent to sign lifting | 52,029,943 |
| Sum of affine sign-kernel sizes | 426,230,308,864 |
| Modulo-eight survivors | 13,018,215 |
| Exact-DC survivors | 5,202,822 |
| Completed exact endpoint | H7572 |

The kernel-size counter sums the sizes of the affine spaces presented to
the packed filters; the final successful kernel stops at its first exact
zero. It does not count that many full integer correlation checks.
These are finite-run counters, not independent random trials or an exhaustion
of all axes. The successful native kernel's SHA256 was
`28d3ce519e1ac9066b070b835f6a82badbd1d300e298c15f4e1e9d9f736b1d87`;
the sampler source SHA256 was
`d5b1b14202e52e39fddb5cf81d91cdd3a271d4bdf4f3230beed62afb93de09aa`.
The active extraction retains those formulas and checks its generated631
coefficient table against the successful table's deterministic hash.
Later ports use 128-bit color coordinates. Retained input maps reach field
degree 22; the implementation accepts degrees up to25. Binary-polynomial
[Rabin witnesses](../certificates/scalar-lift/hadamard-frame-fieldhood.json)
certify every retained modulus. Large fields use polynomial Euclid for
inverses, with the multiplication identity checked at each use.

An independent bounded repeat with seed **2026112421** used 45.0008 seconds,
lifted **1,562,644** axis words and covered **12,801,220,608** affine signs.
It retained 392,241 modulo-eight and 156,214 exact-DC survivors, with no
second endpoint. This run neither excludes further T631 sequences nor
establishes a new equivalence class.

The method then **succeeded at a different prime**: length 271, multiplier 90
of order 15, row sums `(-14,5,5,5)`, and parity field `GF(512)` with modulus
`0x295`. The independently derived 36-column map passes all 667 quadratic
coefficient checks. With seed **20261124271**, the extracted active engine
found an exact T271 after **2956** axis lifts and **0.037877 seconds** of
native search; complete coefficient, frame, cross and assembled **H3252**
integer verification took another **1.008 seconds**. Its packed sign-bit
hash is `599b87c00f75a0b29118443f8d5fcaab9073056452c55bc1071c173910aef4fd`.
Order 3252 was already known (Paley I applies to the prime 3251); no new
equivalence-class claim is made.
The [map, states and exact repeat certificate](../certificates/scalar-lift/hadamard-frame-p271.json)
make this a reproducible independent-prime positive control.

A second port **also completed**: length 547, multiplier 169 of order 21,
row sums `(-20,7,7,7)`, and parity field `GF(8192)` with modulus `0x34bd`.
Seed **20261124548** found T547 after **25.2677 seconds**:

| Stage | Count |
|---|---:|
| Field parameters sampled | 4,219,348 |
| Axis words sent to sign lifting | 1,056,115 |
| Sum of affine sign-kernel sizes | 2,166,233,088 |
| Modulo-eight survivors | 265,614 |
| Exact-DC survivors | 120,076 |

The same source of order 3, with frame dimension/scale **2188**, gives a
fully verified real **H6564**. Complete coefficient, frame, cross and
assembled matrix integer verification took **6.015846 seconds**. Its
packed sign-bit hash is
`ff886685971143ef09c602b29e73463c1d023b43c4fef340747a5894ec8e0997`.
The [completed p547 certificate](../certificates/scalar-lift/hadamard-frame-p547.json)
retains the map, all 27 states, source, scalar action, incoming trace form,
frame data, actual alphabets, search counters and source hashes.
Order 6564 was already known by Paley I at the prime 6563. Novelty of T547
and equivalence classes has not been assessed. This is a successful repeat
of the construction method, with no additional order-novelty claim.

Compact run summaries and controls are in the
[engine certificate](../certificates/scalar-lift/hadamard-frame-search.json).
Its original source hashes retain historical scope; the current active
source binding and the degree 17/22 control certificates are distinct.

The p127/m9 chart also produced completed periodic **T127** in 0.002314
seconds of native search. The [p127 certificate](../certificates/scalar-lift/hadamard-frame-p127.json)
retains its 15 states, four literal sequences and a direct disjoint-circulant
Gram audit. Full scalar-frame and H1524 integer checks also pass. The
2026 survey explicitly lists T127 as its smallest open periodic T-matrix
order; ingredient priority is audited separately from the already known
Hadamard order 1524. This illustrates the value of checking open ingredients
alongside unresolved real matrix orders.

## Commands and current limits

Run the controls, including all 150 solutions in the complete fixed-`k=3`
T19 chart with origin/first head positive and row sums `(4,1,1,1)`,
a full H228 integer Gram, all 8192 masks of the planted T631 sign kernel,
and the independently found H3252 with its complete integer Gram:

```sh
PYTHONPATH=src:. python -m experiments.scalar_lift.hadamard_frame_search \
  --self-test --workdir /tmp/hadamard-frame-controls
```

Run a bounded independent search:

```sh
PYTHONPATH=src:. python -m experiments.scalar_lift.hadamard_frame_search \
  --seconds 60 --seed 2026112422 --heads 5,7,9 --workdir /tmp/hadamard-frame-run
```

`--heads all` removes the optional head-count selection. Each native zero
undergoes complete coefficient, frame, cross and assembled endpoint integer
verification before the Python runner saves a completed endpoint certificate.
The only build dependency beyond Python/NumPy is a C++17 compiler.

`--lift-pool pool.json` accepts finite axis pools with fields `p,m,mu`,
`head_net,tail_net,fixed_first_head_negative,patterns`; each pattern is an
integer axis list. It saves every exact solution and every modulo-eight
survivor. This part applies to prime cyclic charts with a multiplier of
order divisible by three, at most 63 orbit variables, and sign kernels no
larger than the stated native enumeration limit. The current stochastic
SL2 sampler ships with the 60-bit map for length 631. `--parity-map map.json`
can supply another determinant-coordinate map with fields `p,m,mu`,
`field_degree,field_modulus,inverse_columns`, `head_net,tail_net` and
`fixed_first_head_negative`. The tool proves the full quadratic identity
before running it. The current sampler requires one SL2 factor of field
degree at most 25; the largest validated supplied maps have degree 22. It does
not infer new field maps automatically. The exact
[compact target scan and bound](HADAMARD_FRAME_TARGET_SCAN.md) describe this
parameter scope.

Length 223 has 74 nonzero multiplication orbits, a parity field
`GF(2^37)`, and a much larger sign kernel. It exceeds this original engine's
word and finite-enumeration limits. The current
[length-223 implementation](LOW_ORDER_SEARCH.md) supplies the larger exact
map, native determinant sampler, sign filter and resumable campaign while
retaining the same reduction chain. The compact
length 919 chart is now ported: its 68-bit color map and all 32768 signs of
a native kernel pass independent exact checks. The
length 547 map has now been derived and passes all 1379 quadratic coefficient
checks; its [52-column map and completed endpoint](../certificates/scalar-lift/hadamard-frame-p547.json)
can be replayed with `--parity-map`. The deterministic planted-state test
reconstructs H6564 and verifies its full integer Gram without repeating the
stochastic search.

## The even-degree port and conditional headless sampler

At p397 the parity field has degree 22. A basis `(1,omega)` is dependent
over that field, so the port uses frequency orbits represented by 1 and 1/2,
and a basis `(1,gamma)` with `gamma+bar(gamma)=1` in its quadratic extension.
After a linear normalization and a Frobenius square, parity is again
`P*V+Q*U=1`. Its 88-bit map passes all 3917 quadratic coefficient checks.
The native filter matches all 1,048,576 sign masks of a 20-variable kernel,
including both surviving masks, using 211 coefficients recomputed from
expanded integer correlations. Recursive truth words handle larger blocks.

The p1453 row sums `(1,22,22,22)` permit **zero nonorigin head entries**.
The ordinary first-head gauge requires at least one nonorigin head.
`--headless` covers the zero-head stratum using a different verified gauge:
choose a positive orbit at `r` and set `b'(x)=C^k*b(r*x)`, rotating its
channel to 1. Positive tail net guarantees such an orbit. Multiplication
by `r` commutes with the multiplier, fixes the origin and preserves all
correlations; the common positive tail sums are unchanged by the channel
cycle. No global sign flip is used.

For each uniformly sampled nonzero SL2 first row `(p,q)`, the color word is
affine in the free field parameter `w`. Two linear equations impose the
first-tail gauge. Bit masks then exclude every zero color and emit **all**
admissible `w`. Each admissible axis word has the same visit intensity;
the emitted stream is correlated, rather than independent uniform draws.

The complete p13 control gives 12 conditional axis words and 18 sign
solutions, equal to independent enumeration, and a full H156 Gram check.
For an actual p1453 field row, all 4,194,304 values of `w` were checked:
the same nine survivors were returned by the conditional filter in 0.000640
seconds versus 0.1354 seconds naively. A zero-p row also matched completely,
with eleven survivors. A ten-second calibration lifted 10044 headless axis
words. The [headless control certificate](../certificates/scalar-lift/hadamard-frame-headless.json) retains these checks and scopes.

Four bounded long streams, two at 919, one at 397 and one headless 1453,
finished without an exact zero. Their
[terminal campaign certificate](../certificates/scalar-lift/hadamard-frame-campaign-20261124.json)
binds the frozen sources and actual result files:

| Combined counter | Value |
|---|---:|
| Axis words sent to sign lifting | 54,695,329 |
| Sum of affine sign-kernel sizes | 5,955,345,580,032 |
| Modulo-eight survivors | 14,235,675 |
| Exact-DC survivors | 4,770,525 |

Summed native wall time was 5400.001 seconds, with concurrent streams; this
is not a CPU-time measurement. The kernel-size sum is not a count of full
integer correlation evaluations. These finite negative runs are separate
from the positive T127 pilot and imply no nonexistence result.

```sh
PYTHONPATH=src:. python -m experiments.scalar_lift.hadamard_frame_search \
  --headless --parity-map certificates/scalar-lift/hadamard-frame-p1453.json \
  --seconds 60 --seed 202611261455 --workdir /tmp/hadamard-frame-headless
```

Source hashes are captured before the native build and search. Every exact
T candidate still receives full coefficient/frame/cross/endpoint integer
verification before a completed scalar endpoint is reported.

The completed frame can also be reused with other admitted sources; see
[`reuse_hadamard_7572_frame`](../experiments/scalar_lift/reuse_hadamard_7572_frame.py).
Every resulting order still needs its own exact source proof, alphabet test
and order-novelty audit.

# Lower-order construction search

The first target is **2676**, followed by the remaining targets in increasing
order. The submitted manuscript and its frozen construction/proof packages
are preserved separately from this search.

| Missing design | Direct target orders |
|---|---|
| Periodic T-matrices of length 223 | 2676, 4460, 11596 |
| Periodic T-matrices of length 167 | 4676, 8684 |
| Periodic T-matrices of length 179 | 7876, 9308 |

A length 223 witness gives the four-variable orthogonal design proved in the
manuscript. The already verified Williamson inputs of orders 3, 5, 13 then
produce 2676, 4460, 11596. Miyamoto's construction would also give 10708 from
2676, with prime parameter 2677. These orders occur in the
[author-maintained unresolved frontier](https://github.com/vatskishan03/Hadamard/blob/main/FRONTIER.md),
checked 25 November 2026. Every completed endpoint will be compared again
with primary records before adding a new-order claim.

## Exact target

Let $p=223$ and let $e_0,e_1,e_2,e_3$ be the standard basis of
$\mathbb Z^4$. Multiplication by $\mu=183$ on $\mathbb Z/p\mathbb Z$
has order three. Let $\Pi e_0=e_0$, $\Pi e_1=e_2$, $\Pi e_2=e_3$,
and $\Pi e_3=e_1$. We seek a function

\[
b:\mathbb Z/p\mathbb Z\longrightarrow\{\pm e_0,\pm e_1,\pm e_2,\pm e_3\}
\]

such that $b(0)=e_0$, $b(\mu x)=\Pi b(x)$, and

\[
C(s):=\sum_{x\in\mathbb Z/p\mathbb Z}b(x)^{\mathsf T}b(x+s)
=\begin{cases}p,&s=0,\\0,&s\ne0.\end{cases}
\]

The four coordinate sequences then have disjoint supports and entries
$0,1,-1$. Their circulant matrices $T_0,T_1,T_2,T_3$ satisfy
$\sum_i T_iT_i^{\mathsf T}=pI_p$, since its entry at displacement $s$
is exactly $C(s)$. These are the required periodic T-matrices.

There are 74 nonzero multiplication orbits, each of size three. Since
$C(s)=C(-s)=C(\mu s)$, 37 nonzero shift representatives suffice for
searching; the final verifier checks all 223 shifts.

The sum $\sum_x b(x)$ has the form $(a,t,t,t)$ with $a\equiv1\pmod3$.
Summing $C(s)$ over all $s$ gives $a^2+3t^2=223$. The integer solutions
are $(a,t)=(\pm14,\pm3)$. The congruence forces $a=-14$; changing the
signs of the last three coordinates permits $t=3$. Thus the normalized
coordinate sums are $(-14,3,3,3)$.

Let $h$ be the number of nonzero multiplication orbits whose values use
coordinate zero. Their signs must sum to $(-14-1)/3=-5$. The remaining
$74-h$ orbit signs must sum to $3$. Consequently $h$ is odd and
$5\leq h\leq71$. The numbers of positive and negative signs are

\[
(n_{0,+},n_{0,-},n_{\ne0,+},n_{\ne0,-})
=\left(\frac{h-5}{2},\frac{h+5}{2},\frac{77-h}{2},\frac{71-h}{2}\right).
\]

There is a negative coordinate-zero orbit. Rescaling the domain places
one at the residue $1$, which is the gauge used by the sign solver.
At $h=5$ all five coordinate-zero signs are negative. After fixing the
negative sign at $1$, the choices of positive coordinate-zero signs at
$h=7$ and $h=9$ number $\binom61=6$ and $\binom82=28$.

## Verified parametrization and sign solver

[`large_parity_map.py`](../experiments/scalar_lift/large_parity_map.py)
constructs the 148-bit coordinate map by recovering the field structure of
the correlation quadratic forms and splitting two hyperbolic planes.
The [retained map](../certificates/scalar-lift/t223-parity-map.json) uses
GF(2^37). All 11,027 constant, linear and quadratic coefficients are checked.
A second check expands all coefficient probes and compares all 222 nonzero
shifts: 2,447,994 parity correlations. These finite checks establish the
map used by the sampler; signs are solved separately.

More precisely, let $F:\mathbb F_2^{148}\to\mathbb F_2^{37}$ give the
37 unsigned correlation parities, encoding each coordinate label by two
bits. The certificate records invertible binary linear maps $T$ and $L$
for which

\[
F\bigl(T(P,Q,U,V)\bigr)=L(PV+QU+1),
\qquad P,Q,U,V\in\mathbb F_{2^{37}}.
\]

Each side is quadratic over $\mathbb F_2$. Checking its constant term,
148 linear coefficients and $\binom{148}{2}$ quadratic coefficients
proves the identity on every input. The rank checks therefore identify
*all* admissible coordinate assignments with determinant-one two-by-two
matrices over $\mathbb F_{2^{37}}$. The sampling step chooses from this
set; it does not exhaust it.

[`quadratic_sign_lift.py`](../experiments/scalar_lift/quadratic_sign_lift.py)
eliminates the modulo-four equations with arbitrary-size integer bitsets,
forces signs implied by exact counts, and derives every quadratic
modulo-eight equation. The new
[native filter](../experiments/scalar_lift/quadratic_sign_lift_bits.cpp)
checks 64 assignments per word with bounded memory. Returned modular roots
must satisfy the exact sign counts and every integer correlation, followed
by an independent check at all 223 shifts.

The CP-SAT and optional CryptoMiniSat implementations provide independent
solver controls. Complete small controls include all 90 T19 coordinate
patterns and 150 full sign solutions, the planted T127/T631 witnesses,
and a negative-coordinate-sum T31 witness. Native Boolean controls compare
with complete enumeration and cover all 2^32 assignments of a planted
32-variable test.

## Run a bounded search

Use the repository Python environment with NumPy, SymPy and OR-Tools, and a
Clang/GCC compiler supporting C++17. Generated runs belong in the ignored
artifact store; each worker needs its own output directory.

```sh
PYTHONPATH=src:. .venv/bin/python -m experiments.scalar_lift.large_parity_map \
  --output /tmp/t223-parity-map.json
PYTHONPATH=src:. .venv/bin/python -m experiments.scalar_lift.low_order_campaign \
  --config certificates/scalar-lift/t223-parity-map.json \
  --output artifacts/scalar-frontier/low-order-search/h5 \
  --seconds 900 --heads 5 --seed 20261125223
```

Here `--heads` is the number $h$ defined above. Use `--heads 7` or
`--heads 9` for those other counts. The campaign enumerates
the exact coordinate-zero sign choices, saves a checkpoint after each coordinate
pattern, and writes a verified witness only after all expanded
correlations vanish. Its reported affine-assignment count refers to the
modular filter, with exact integer checks on survivors.

An independent
[local search](../experiments/scalar_lift/t_sequence_local_search.py)
changes the signed coordinate states while preserving exact counts. Its
move set includes full-state swaps and balanced sign exchanges. Weighted
objectives and exact
[neighborhood repair](../experiments/scalar_lift/t_sequence_neighborhood.py)
provide further searches around small correlation residuals. Every
checkpoint records its actual residuals. Only a zero is a T-matrix witness.

The local search can also change $h$: replacing an oppositely signed pair
using coordinates $1,2,3$ by $(+e_0,-e_0)$ preserves all four coordinate
sums, as does the reverse replacement. The two original coordinates may
be different, because each nonzero-coordinate orbit contributes its sign
once to each of the last three sums.

The [Fourier projection search](../experiments/scalar_lift/t_sequence_projection.py)
alternates a discrete assignment step with a spectral norm step. With the
unnormalized Fourier transform
$\widehat b(k)=\sum_x b(x)e^{-2\pi ikx/p}$, the desired correlation
identity is equivalent to $\|\widehat b(k)\|^2=p$ at every frequency.
The zero-frequency vector is fixed at $(-14,3,3,3)$. The discrete step
maximizes the inner product with a signed-basis function having the exact
orbit counts; it uses an assignment algorithm for fixed $h$, or an exact
dynamic program when $h$ may vary. Floating point is used to choose search
iterates. Every saved best candidate is rechecked using integer sums at
all shifts, and only an exact zero is a completed design.

Run a longer checkpointed projection campaign with:

```sh
PYTHONPATH=src:. .venv/bin/python -m experiments.scalar_lift.low_order_campaign \
  --engine projection --config certificates/scalar-lift/t223-parity-map.json \
  --output artifacts/scalar-frontier/low-order-search/projection \
  --seconds 3600 --trial-seconds 60 --seed 22320261126
```

It varies the projection order, relaxation parameter and restart interval.
Most trials start afresh; every fifth trial may use the retained best
integer candidate. For a previous checkpoint, add `--initial PATH`.
`campaign.json` records progress, `best.json` retains the best exact
residuals, and `verified-witness.json` is written only after every
correlation is zero. To run local search with changing coordinate counts,
use `t_sequence_local_search --free-heads` and a fresh `--report PATH`.

Focused verification:

```sh
PYTHONPATH=src:. .venv/bin/python -m pytest -q \
  tests/test_large_parity_map.py tests/test_quadratic_sign_lift.py \
  tests/test_boolean_quadratic_sat.py tests/test_t_sequence_neighborhood.py \
  tests/test_t_sequence_projection.py tests/test_t_sequence_local_search.py \
  tests/test_low_order_campaign.py
```

A timeout applies only to its run, and a fully checked coordinate pattern
applies only to that pattern and its stated gauges. Neither excludes a
T-design or a Hadamard matrix at the target order.

## Bounded results, 25 November 2026

Three 900-second sign-search batches at $h=5,7,9$ completed 360, 282,
and 224 coordinate patterns respectively. Including their final partial
patterns, the native filter covered 5,052,591,767,552 affine assignments.
There were 34 modulo-eight roots and seven survivors of the exact count
conditions; none satisfied every integer correlation. These are searches
inside the multiplier-covariant chart defined above.

The first 32 fixed-count Fourier trials also found no zero. The best local
candidate has $E=\sum_{\ell=1}^{37}C(s_\ell)^2=42$, with largest absolute
correlation two. Allowing coordinate-count changes visited $5\leq h\leq31$
but did not improve it in a 180-second run. It is an incomplete candidate,
and **no new order is claimed by these runs**.

The full repository suite passed **359 tests and four subtests** after
adding these engines. The submitted paper's 18 completed new orders and
six families are unchanged.

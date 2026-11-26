# Lower-order construction search

The first target is **2676**, followed by the remaining targets in increasing
order. The submitted manuscript and its frozen construction/proof packages
are preserved separately from this search.

The construction is the manuscript's restriction of scalars with a certified
flat endpoint. The exact search described here supplies a new compatible
frame, using the same determinant-coordinate and sign-lifting method that
found T631 and hence H7572.
The [end-to-end walkthrough](CONSTRUCTION_WALKTHROUGH.md) explains both the
finite search and the subsequent scalar-restriction proof.

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

The CP-SAT implementation provides an independent solver control. Complete
small controls include all 90 T19 coordinate patterns and 150 full sign
solutions, the planted T127/T631 witnesses,
and a negative-coordinate-sum T31 witness. Native Boolean controls compare
with complete enumeration and cover all 2^32 assignments of a planted
32-variable test.

## Run the exact campaign

Use the repository Python environment with NumPy, SymPy and OR-Tools, and a
Clang/GCC compiler supporting C++17. Each worker needs its own ignored output
directory. The runner freezes its source files and parity map before detaching.
The native sampler checks its field arithmetic and binary coordinate map
against the Python reference. Its tagged 256-bit generator state and draw
count are saved in the same ledger as the search cursor.

```sh
PYTHONPATH=src:. .venv/bin/python -m experiments.scalar_lift.exact_t_campaign \
  --config certificates/scalar-lift/t223-parity-map.json \
  --output artifacts/scalar-frontier/low-order-search/exact223-h5 \
  --seconds 259200 --heads 5,7,9 --branch-seconds 120 --nice 10 \
  --sampler native --branch-policy all \
  --seed 20261126001 --detach
```

Here `--heads` is the number $h$ defined above. The `all` branch policy
checks every coordinate-zero sign assignment with the exact required count.
The `affine` policy instead leaves these signs in the complete affine system.
Use separate directories and seeds for independent workers.
The counts listed by `--heads` are visited in turn. A complementary worker
uses `--heads 15,17,19,21 --branch-policy affine`: it checks the full sign
kernel once per coordinate pattern instead of separately preparing tens of
thousands of exact coordinate-zero sign choices. `--sampler python` retains
the original reference generator, with its own explicitly distinct stream.

The SQLite ledger stores sampled coordinate patterns, completed branch
checks and sampler state. Restarting with the same configuration and frozen
sources resumes the ledger. An interrupted branch remains unfinished;
`status.json` reports the actual completed work. Creating a `STOP` file in
the worker directory requests a stop after the current branch returns.
The run log and source manifest are kept beside the ledger.

A candidate becomes a witness only after every expanded integer correlation
vanishes and the exact coordinate sums hold. The worker then independently
checks the orthogonal-design identities and the verified Williamson inputs,
and writes literal endpoint recipes. Literature novelty and new Lean
literals are assessed separately after a completed witness exists.

Focused verification:

```sh
PYTHONPATH=src:. .venv/bin/python -m pytest -q \
  tests/test_large_parity_map.py tests/test_quadratic_sign_lift.py \
  tests/test_exact_t_campaign.py tests/test_exact_parity_sampler.py
```

A timeout applies only to its run, and a fully checked coordinate pattern
applies only to that pattern and its stated gauges. Neither excludes a
T-design or a Hadamard matrix at the target order.

## Unrestricted sign solvers for lengths 167 and 179

The [unrestricted parity module](../experiments/scalar_lift/unrestricted_parity_map.py)
uses the same determinant parametrization without a multiplier symmetry.
The [167 certificate](../certificates/scalar-lift/unrestricted-parity-p167.json)
has field degree 83, with 55,279 coefficients and 9,176,314 expanded parity
comparisons checked. The [179 certificate](../certificates/scalar-lift/unrestricted-parity-p179.json)
has degree 89, with 63,547 coefficients and 11,311,366 expanded comparisons.

Fix the value at residue zero to $+e_0$. Up to coordinate signs,
permutations and cyclic translation, the complete row-sum types are:

| Length | Row-sum types |
|---:|---|
| 167 | (9,6,5,5), (9,7,6,1), (9,9,2,1), (10,7,3,3), (11,6,3,1) |
| 179 | (9,7,7,0), (9,8,5,3), (9,9,4,1), (11,7,3,0), (12,5,3,1), (13,3,1,0) |

The campaign computes these lists by complete integer four-square
enumeration. A positive first coordinate guarantees a positive entry in
that coordinate, which can be translated to residue zero.

For each sampled coordinate pattern, the
[shared sign model](../experiments/scalar_lift/unrestricted_t_sign_model.py)
eliminates the exact modulo-four equations and imposes all four negative-sign
counts. Typical affine dimensions are 80 and 86. The two backends are:

- [SAT](../experiments/scalar_lift/unrestricted_t_sign_sat.py): native XORs,
  modulo-eight quadratic constraints, and exact binary-adder encodings of
  all four counts and every integer correlation; requires `cryptominisat5`.
- [CP-SAT](../experiments/scalar_lift/unrestricted_t_sign_cp.py): the same
  constraints in OR-Tools, using one solver worker and arbitrary-size
  binary affine masks.

Both reconstruct every returned sign word independently and check all $p$
integer shifts. `UNSAT` excludes only the given coordinate pattern, row
sums and recorded gauges. `UNKNOWN` is an unfinished search. Solver
implementation and successful model construction alone supply no endpoint.

Run the full exact pipeline with source snapshots and a SQLite ledger:

```sh
PYTHONPATH=src:. .venv/bin/python -m experiments.scalar_lift.unrestricted_t_campaign \
  --config certificates/scalar-lift/unrestricted-parity-p167.json \
  --output artifacts/scalar-frontier/unrestricted-p167 \
  --seconds 259200 --pattern-seconds 900 --backends sat,cp \
  --seed 20261126167 --nice 10 --detach
PYTHONPATH=src:. .venv/bin/python -m experiments.scalar_lift.unrestricted_t_campaign \
  --config certificates/scalar-lift/unrestricted-parity-p179.json \
  --output artifacts/scalar-frontier/unrestricted-p179 \
  --seconds 259200 --pattern-seconds 900 --backends cp,sat \
  --seed 20261126179 --nice 10 --detach
```

Each worker cycles through every row-sum type and alternates the backends
between complete passes. The default gives each fixed-pattern solve up to
15 minutes within a three-day campaign. A heartbeat updates `status.json`
while a solver is running. The ledger distinguishes exact exclusions from
timeouts. Interrupted work can be resumed with its frozen launch command;
a committed positive witness is recovered before any further search.

An exact witness triggers the independent orthogonal-design/source checks
and endpoint recipes already used by the length-223 campaign. The finite
controls include the complete T7 census (180 sign solutions), independent
fixed-pattern T19 enumeration, planted T127/T631 recovery, and real 80/86-bit
model-building checks. Replay them with:

```sh
PYTHONPATH=src:. .venv/bin/python -m pytest -q \
  tests/test_unrestricted_parity_map.py tests/test_unrestricted_t_sign_model.py \
  tests/test_unrestricted_t_sign_sat.py tests/test_unrestricted_t_sign_cp.py \
  tests/test_unrestricted_t_campaign.py
```

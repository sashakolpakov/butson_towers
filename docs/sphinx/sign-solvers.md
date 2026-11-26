(sign-solvers)=
# Solving the exact sign equations

Every backend solves the same kind of fixed-coordinate problem described in
{ref}`exact-search`: choose the remaining signs with the prescribed coordinate
sums and zero nontrivial integer correlations. The parity parametrization
chooses coordinates; the sign solver cannot compensate for a different
coordinate pattern. A verified positive result is a literal finite input to
the construction. A negative or unfinished solve has only its recorded
fixed-pattern scope.

The retained length-631 witness was found with the original C++ search.
The larger native filter extends the same modular-lifting method to the
length-223 sign kernels. SAT and CP-SAT provide exact constraint models for
the unrestricted 167 and 179 kernels, where representative affine dimensions
are 80 and 86. They avoid explicitly generating all $2^{80}$ or $2^{86}$
assignments; this does not guarantee that a bounded solve will finish.

| Implementation | Prepared variables | Method after modulo four | Final gate |
|---|---|---|---|
| Original C++ | Up to 63 multiplication-orbit indices; configured nullity limit | Native modulo-eight filtering and finite sign enumeration | Exact counts, integer correlations and expanded shifts |
| Native filter | Arbitrary-size preparation; native enumeration limited to 60 free bits | Packed, bit-sliced exhaustive modulo-eight filter | Recheck every root, exact counts, integer correlations and expanded shifts |
| SAT | Actual sites; Python masks have no 64-bit bound | Native-XOR SAT, exact count and correlation adders | Independent integer reconstruction and expanded shifts |
| CP-SAT | Actual sites; Python masks have no 64-bit bound | OR-Tools Boolean, XOR and integer constraints | Independent integer reconstruction and expanded shifts |

The last two rows are implemented solvers, not new witnesses. No new
length-167, length-179 or length-223 input follows from their successful model
construction, validation tests or timed `UNKNOWN` runs.

(sign-solvers-shared-model)=
## The shared model

The unrestricted public preparation API is

```python
system = unrestricted_t_sign_model.prepare(
    axes, data_or_p, target_row, fixed_signs=None, force_extremes=True
)
```

The last two arguments are keyword-only. `axes[i]` is the coordinate label
at actual residue $i$, in $\{0,1,2,3\}$, and `axes[0]` must be zero.
`data_or_p` is an odd prime or a chart from `periodic_chart(p)`.
`target_row` contains the four integer sums. It must have squared norm $p$.
`fixed_signs` maps site indices to negative-sign bits, with zero meaning
positive. The origin bit is fixed to zero unconditionally.

The result retains the affine particular solution and basis, each site's
affine form, the four site groups and exact negative counts, the symmetric
pair weights and diagonal constants, the target sums and the fixed signs.
`InfeasibleAxes` reports inconsistent parity, support, fixed-sign or
modulo-four arithmetic. Malformed input is a `ValueError`, rather than a
mathematical exclusion.

For a reconstructed negative-sign mask $x$, `correlations(system, x)` uses
the exact integers

```{math}
C_\ell=D_\ell+\sum_{i<j}w_{ij,\ell}(1-2x_i)(1-2x_j).
```

`quadratic_mod8(system)` derives $C_\ell/4\pmod2$ on the affine space.
The derivation and its divisibility conditions appear in
{ref}`exact-search-modular-lifting`.

(sign-solvers-native)=
## Original C++ search and the larger native filter

The original
{repo-file}`C++ engine <experiments/scalar_lift/hadamard_frame_search.cpp>`
stores orbit masks in machine words. It checks the table dimensions,
eliminates the modulo-four equations, and records a nullity-limit skip
explicitly. Its usual limit is 24 free bits; the successful 631 chart's
retained input has 13. Reaching that configured limit is not an exclusion.
Within a completed allowed kernel, the acceptance code uses every required
exact condition:

```{literalinclude} ../../experiments/scalar_lift/hadamard_frame_search.cpp
:language: cpp
:start-at:             auto accept=[&](uint64_t x){
:end-at:                 return cfg.stopFirst;
```

The newer Python preparation removes the orbit-mask size restriction. Its
`solve_bits` wrapper sends only the reduced Boolean quadratic system to
{repo-file}`quadratic_sign_lift_bits.cpp <experiments/scalar_lift/quadratic_sign_lift_bits.cpp>`.
The native enumerator still has a declared limit of 60 free bits. This is an
implementation boundary, not a claim that such an enumeration is practical
near that boundary.

### Bit slicing and exact Gray updates

Split the free variables into low variables $z$ and high variables $h$.
For fixed $h$, an equation is a fixed quadratic in $z$, plus a linear form
in $z$ whose coefficients depend linearly on $h$, plus a Boolean constant
quadratic in $h$. Bit slicing represents the values at 64 choices of $z$
in one machine word. The default uses 11 low variables, hence 2,048
assignments in 32 words per block. Block widths from 6 to 12 are supported.

The high assignments are traversed in reflected Gray order
$h(t)=t\oplus(t\mathbin{\gg}1)$. At step $t>0$, the flipped bit is
$u=\operatorname{ctz}(t)$. For a quadratic polynomial, its derivative in
$h_u$ is independent of $h_u$ and linear in the other high bits. Between
consecutive flips of $u$, all lower Gray bits repeat and precisely one
higher bit $v$ changes. Thus that derivative changes by the single Hessian
coefficient indexed by $u,v$. Packing coefficients across equations makes
this an exact word-XOR update:

```{literalinclude} ../../experiments/scalar_lift/quadratic_sign_lift_bits.cpp
:language: cpp
:start-at:             unsigned u=__builtin_ctzll(step);
:end-at:                 derivative[u*constantWords+w]^=hessian[(u*high+v)*constantWords+w];
```

Every high assignment and every low lane is still tested unless a declared
time or root limit stops the run. There is no random pruning. The block
width changes performance and enumeration order, while preserving the full
set being enumerated. The Python wrapper reevaluates every returned root
against the Boolean equations before its integer checks.

The symmetric entrypoint is

```python
result = quadratic_sign_lift.solve_bits(
    axes, data, head_net=1, tail_net=1, first_negative=0,
    seconds=10, max_roots=100000, workdir=None, fixed_signs=None,
    block_variables=11,
)
```

Its count parameters refer to signed *orbit* counts. They are not the four
site sums used by the unrestricted API. For the 223 chart they must be
`head_net=-5`, `tail_net=3`, `first_negative=1`; see
{ref}`exact-search-counts`. `iter_head_sign_branches` generates complete
head-sign branches when explicit branching is desired.

Here `seconds` bounds the native enumeration phase. Python preparation,
compilation when needed, native input parsing and final exact checks add to
the wrapper's reported `elapsed_seconds`.

### Replay a retained small sign kernel

This command reads the existing length-127 input, discards its signs before
solving, recovers them with the native filter, and verifies all 127 shifts.
It compiles a small C++17 binary in a temporary directory. It does not
materialize a large Hadamard endpoint or run a stochastic search.

```sh
PYTHONPATH=src:. .venv/bin/python - <<'PY'
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from experiments.scalar_lift import hadamard_frame_search as frame
from experiments.scalar_lift import quadratic_sign_lift as lift

states = json.loads(Path("artifacts/t-matrices-127/t127.json").read_text())["states"]
data = frame.chart(127, 9, 22)
with TemporaryDirectory() as workdir:
    result = lift.solve_bits([s // 2 for s in states], data,
                             seconds=5, workdir=workdir)
assert result["complete_for_fixed_axes"]
assert [item["states"] for item in result["solutions"]] == [states]
print("T127 signs recovered; every expanded shift verified")
PY
```

(sign-solvers-sat)=
## Native-XOR SAT

The SAT backend requires the `cryptominisat5` executable. Its campaign API is

```python
result = unrestricted_t_sign_sat.solve_rows(
    axes, data, target_row, seconds=120, seed=0,
    max_solutions=1, workdir=None, **options
)
```

`data` may also be the prime length. `solve(system, ...)` consumes an already
prepared model. The default has `exact_correlations=True`, `mod8=True` and
`all_solutions=False`. All four exact counts and every integer correlation
are encoded before solving. In `solve(system, ...)`, the budget starts before
SAT encoding; `solve_rows` first prepares the affine model, then calls that
budgeted routine. Its returned `elapsed_seconds` measures the latter call.
Native attempt times and process-limit outcomes are retained in `attempts`.
Process startup, serialization and shutdown can add overhead to the requested
limit; the recorded elapsed wall time reports what actually occurred. Seeds
are reduced modulo 2,147,483,647 for the native integer argument.

### Affine links and Boolean products

Each site sign bit is linked by a native XOR to its affine free-variable
form. For the modulo-eight equation
$c+\sum_u a_uy_u+\sum_{u<v}b_{uv}y_uy_v=0$, one cached auxiliary literal
represents each needed product $y_uy_v$. If $z=a\land b$, the clauses are

```{math}
(\neg z\lor a)\land(\neg z\lor b)
\land(z\lor\neg a\lor\neg b).
```

They force both implications of $z\leftrightarrow(a\land b)$. A native XOR
then imposes the whole quadratic equation. The caches share identical
products between equations without changing their values.

`XorCNF.xor` cancels duplicate variables and transfers each negated literal
to the right-hand parity. A one-variable XOR becomes a unit clause; an empty
XOR is either a tautology or an empty contradiction clause. The emitted
CryptoMiniSat `x` rows use odd parity, so an even right-hand side is written
by negating one literal. The signed-literal tests check this convention
through the actual executable.

### Pair cuts turn correlations into exact integer sums

Introduce the cut bit $t_{ij}=x_i\oplus x_j$. With
$A_\ell=D_\ell+\sum_{i<j}w_{ij,\ell}$, the correlation equation becomes

```{math}
C_\ell=0\quad\Longleftrightarrow\quad
\sum_{i<j}w_{ij,\ell}t_{ij}=A_\ell/2.
```

The parity stage ensures the right side is integral. Each pair cut is linked
to the already defined site literals by a ternary XOR, avoiding a repeated
long affine XOR for every pair. This is the exact encoding used by the
backend:

```{literalinclude} ../../experiments/scalar_lift/unrestricted_t_sign_sat.py
:language: python
:pyobject: add_correlation
```

`XorCNF.weighted_exact` encodes a nonnegative weighted sum using Boolean
binary adders. Complementary literals first contribute their common weight
to the constant target. A common divisor of all remaining weights is
cancelled only if it also divides the target. Binary weight bits are then
placed into columns. Three literals $a,b,c$ in a column produce a same-column
sum bit $a\oplus b\oplus c$ and a next-column carry
$\operatorname{maj}(a,b,c)$, because

```{math}
a+b+c=(a\oplus b\oplus c)+2\operatorname{maj}(a,b,c).
```

Two literals use XOR and AND in the same way. Finally every output bit is
fixed to the target bit, retaining all carries. The resulting predicate is
the exact integer sum, rather than just its low bits. The identical primitive
encodes the four negative-sign counts with unit weights.

The complete model builder interleaves those pieces as follows:

```{literalinclude} ../../experiments/scalar_lift/unrestricted_t_sign_sat.py
:language: python
:pyobject: encode_sign_system
```

### Optional lazy correlations and complete blocking

With `exact_correlations=False`, the initial model retains the exact counts
and modulo-eight conditions. Every returned assignment is checked against
all integer correlations. For a nonzero residual, the backend adds the full
exact equation for every violated shift and solves again within the remaining
budget. Each addition is necessary for every true solution of that fixed
model, so a subsequent `UNSAT` is sound. A modular survivor is never returned
as a successful T-sequence.

For `all_solutions=True`, each verified solution is blocked by its complete
free-variable assignment, rather than by an auxiliary adder assignment. Set
`max_solutions=0` to remove the solution-count limit. A final unsatisfiable
residual model proves that enumeration is complete. A timeout or a requested
solution limit leaves it incomplete.

(sign-solvers-cp)=
## OR-Tools CP-SAT

The CP-SAT backend requires OR-Tools in the Python environment. Its API is

```python
result = unrestricted_t_sign_cp.solve(
    axes, data_or_p, target_row, seconds=30, seed=0,
    max_solutions=1, workdir=None,
    all_solutions=False, add_mod8=True, affine_links=True,
    fixed_signs=None, fixed_free=None, log_search=False,
)
```

Arguments after `workdir` are keyword-only. The implementation creates site
Boolean variables, affine free variables, cached pair cuts and cached
quadratic products. It imposes the same four exact counts and all exact
weighted cut equations. Native CP-SAT XOR constraints impose the affine
links and modulo-eight equations; multiplication equalities impose the
quadratic products. It validates the model before calling the solver.

```{literalinclude} ../../experiments/scalar_lift/unrestricted_t_sign_cp.py
:language: python
:pyobject: build_model
```

The solver uses one search worker and subtracts model-construction time from
the requested budget. `fixed_free` can additionally fix selected affine bits;
a completeness claim then applies to that smaller constrained model.
`add_mod8=False` omits redundant modular quadratics while retaining all
integer equations. Omitting affine links is supported only when the
modulo-eight equations are also omitted; the site-level exact counts,
correlations and recorded fixed signs remain imposed.

Every solution callback reconstructs its negative mask, checks all four
counts and compressed integer correlations, and calls the independent
all-shift verifier. A callback that stops after the requested number of
solutions does not certify exhaustion, even if the returned word is valid.

(sign-solvers-status)=
## What the statuses prove

Always inspect `solutions`, their `expanded_verification`, the completeness
flag and the recorded scope together.

| Result | Meaning |
|---|---|
| Verified nonempty `solutions` | These literal words satisfy the requested counts and every integer shift. Other signs or coordinates may remain unsearched. |
| `INFEASIBLE_ARITHMETIC` | Necessary exact count, fixed-sign or modulo-four arithmetic is inconsistent for this input. |
| SAT `status="UNSAT"`, with completeness true | No sign word satisfies this fixed coordinate pattern, row sums and supplied gauges. |
| SAT `status="UNKNOWN"` | The bounded search has not established a solution or exclusion. |
| CP-SAT `INFEASIBLE`, with completeness true | The full recorded CP model is infeasible. Any `fixed_free` restrictions are part of that scope. |
| Native `COMPLETE`, with completeness true | Every assignment of the prepared affine kernel passed through the modular filter, and all returned roots received the exact gates. |
| Native `UNFINISHED`, root limit, time limit or stopped solution enumeration | The stated enumeration remains incomplete. |

`complete_for_fixed_axes` is the common completeness field. The CP backend
also exposes `complete_for_fixed_constraints`. The SAT result's `status`
summarizes whether any verified solutions were found; its `solver_status`
records the last native attempt. In particular, complete all-solution
enumeration may have `status="SAT"` and final `solver_status="UNSAT"`:
solutions were found, then the blocked residual problem was exhausted.
CP-SAT preserves its native names such as `FEASIBLE`, `OPTIMAL`, `INFEASIBLE`
and `UNKNOWN`; its explicit completeness flag accounts for callback limits.

None of these fixed-pattern outcomes is a nonexistence theorem for a
T-matrix length or a Hadamard order. A campaign samples many coordinate
patterns and records which solves completed. Sampling the full proven
parity locus is not the same as exhausting it. The submitted manuscript's
frozen finite witnesses and all-parameter construction proofs retain their
own scope; ongoing runs do not revise them merely by timing out.

(sign-solvers-replay)=
## Replay both unrestricted backends on one small pattern

The following command is deterministic and normally finishes in well under
a second once dependencies are installed. It compares the complete solution
sets for a fixed T7 coordinate pattern. It runs no frontier search.

```sh
PYTHONPATH=src:. .venv/bin/python - <<'PY'
from experiments.scalar_lift import unrestricted_t_sign_sat as sat
from experiments.scalar_lift import unrestricted_t_sign_cp as cp

axes = [0, 3, 0, 2, 2, 2, 1]
target = (2, 1, 1, 1)
a = sat.solve_rows(axes, 7, target, seconds=5,
                   all_solutions=True, max_solutions=0)
b = cp.solve(axes, 7, target, seconds=5,
             all_solutions=True, max_solutions=0)
assert a["complete_for_fixed_axes"] and b["complete_for_fixed_axes"]
left = {tuple(item["states"]) for item in a["solutions"]}
right = {tuple(item["states"]) for item in b["solutions"]}
assert left == right and len(left) == 2
assert all(item["expanded_verification"][
    "every_expanded_periodic_correlation_verified"
] for result in (a, b) for item in result["solutions"])
print("Both exact backends recover the same two complete T7 solutions")
PY
```

Focused controls are kept beside the implementations:

- {repo-file}`Modular and native controls <tests/test_quadratic_sign_lift.py>`:
  complete symmetric T19 sign sets, retained T127/T631 recovery,
  arbitrary-width elimination, every supported block width, packed-equation
  boundaries, and a planted complete $2^{32}$ assignment check.
- {repo-file}`Unrestricted model controls <tests/test_unrestricted_t_sign_model.py>`:
  the complete length-7 signed-word census, actual 80/86-bit affine models,
  exact correlation comparisons and corruption rejection by the verifier.
- {repo-file}`SAT controls <tests/test_unrestricted_t_sign_sat.py>`:
  exhaustive small integer-adder assignments, native XOR semantics, complete
  T7 sets, planted 127/631 signs, complete small UNSAT, lazy rejection of a
  modulo-eight survivor with nonzero integer correlation, and the real CLI
  with a large production seed.
- {repo-file}`CP-SAT controls <tests/test_unrestricted_t_sign_cp.py>`:
  complete T7 and fixed-pattern T19 sign sets, optional modular encodings,
  actual 80/86-bit models and incomplete-budget/solution-limit outcomes.

For a short local verification, the two fixed-pattern backend suites can be
run without starting a campaign:

```sh
PYTHONPATH=src:. .venv/bin/python -m pytest -q \
  tests/test_unrestricted_t_sign_sat.py tests/test_unrestricted_t_sign_cp.py
```

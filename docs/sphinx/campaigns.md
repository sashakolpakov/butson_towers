(campaigns)=
# Exact search campaigns and recoverable computation

The current selection matches eight unresolved orders below 10000 with seven
strongly symmetric Williamson inputs. The
{repo-file}`source model and selection <docs/WILLIAMSON_SYMMETRY_SEARCH.md>`
derive every sign from the occupied coordinate on each orbit. CP-SAT and
an independent native axis sampler and annealer check the remaining integer
conditions. Runs for six orders are active; 3004 and 5692 remain queued.
The earlier eleven one-hot T-input workers are paused.

## Earlier one-hot T-input engines

Their interfaces below remain exact replay references.

The earlier campaigns searched for the finite T-input required by the construction.
It fixes occupied coordinates, solves their signs, and accepts a design
only after an independent integer check at every cyclic shift. The later
orthogonal-design and source proofs then give actual Hadamard endpoints.
A finite search result always records its coordinate pattern, requested
sums and any additional gauges.

Three runners use the same exact modular strategy:

| Runner | Coordinate pattern | Sign solver |
|---|---|---|
| {repo-file}`223 worker <experiments/scalar_lift/exact_t_campaign.py>` | The multiplier-covariant length-223 chart | Native packed Boolean mod8 filtering, followed by exact counts and correlations |
| {repo-file}`Identity multiplier worker <experiments/scalar_lift/identity_t_campaign.py>` | Constant coordinate vectors on multiplier orbits at lengths 223 or 491 | Native mod8 filtering with four separate exact sign counts and every expanded integer shift |
| {repo-file}`167/179 worker <experiments/scalar_lift/unrestricted_t_campaign.py>` | One coordinate label at every site of lengths 167 or 179 | CryptoMiniSat or OR-Tools CP-SAT with all exact counts and correlations |

The submitted manuscript's literal witnesses and frozen packages are
separate from these campaigns. Their existence proofs do not depend on a
campaign finishing.

## Why the filtering stages preserve the problem

For fixed coordinate labels $a_i\in\{0,1,2,3\}$, write a sign as
$(-1)^{x_i}$ with $x_i\in\mathbb F_2$. In an orbit chart $i$ labels an
orbit representative; in an unrestricted chart it labels an actual site.
Each correlation has the form

```{math}
:label: eq-campaign-cuts
C_\ell=C_\ell^{\mathrm{abs}}
-2\sum_{i<j}w_{ij,\ell}(x_i\mathbin{\oplus}x_j).
```

The diagonal contributions are sign independent and included in
$C_\ell^{\mathrm{abs}}$. Coordinate labels first pass a complete parity
parametrization. In the retained determinant coordinates this is

```{math}
:label: eq-campaign-determinant
PV+QU=1\quad\text{over }\mathbb F_{2^k}.
```

The fields have degrees $k=37$ for both 223 actions, $35$ for the identity
491 chart, $83$ for
unrestricted length 167, and $89$ for unrestricted length 179. The stored
binary maps are checked through every constant, linear and quadratic
coefficient. Independent geometry checks compare those probes with the
expanded cyclic correlations. A parity word is a coordinate assignment,
not a signed T-witness.

Modulo four, {eq}`eq-campaign-cuts` becomes a binary linear system. Its
solution is $x=x^{(0)}+By$. On this affine space, each $C_\ell/4$ modulo
two is a Boolean quadratic polynomial in $y$. The mod8 stage solves those
quadratics. Exact negative-sign counts remain necessary: for coordinate
support size $n_j$ and requested sum $r_j$, the count is
$(n_j-r_j)/2$. Summing all correlations also requires
$\sum_jr_j^2=p$.

The unrestricted backends encode the full integer equation

```{math}
:label: eq-campaign-integer-cut
\sum_{i<j}w_{ij,\ell}(x_i\mathbin{\oplus}x_j)
=C_\ell^{\mathrm{abs}}/2
```

and all four exact counts. The affine and mod8 constraints are additional
exact cuts. Their Python integer masks support the 80- and 86-dimensional
free-sign systems without a machine-word truncation. The native 223
filter instead evaluates sign assignments in packed words; its survivors
still undergo the complete integer tests.

The shared implementations are
{repo-file}`quadratic_sign_lift.py <experiments/scalar_lift/quadratic_sign_lift.py>` and
{repo-file}`unrestricted_t_sign_model.py <experiments/scalar_lift/unrestricted_t_sign_model.py>`.
The unrestricted positive verifier reconstructs actual signed vectors and
checks all $p$ shifts, independently of the compressed cut equations.

## Run the multiplier-covariant length-223 search

Install the Python environment from the repository root:

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -e '.[dev,native,search]'
```

The exact 223 runner also needs a Clang/GCC C++17 compiler. This example
requests a three-day search session in a fresh ignored directory:

```sh
PYTHONPATH=src:. .venv/bin/python -m experiments.scalar_lift.exact_t_campaign \
  --config certificates/scalar-lift/t223-parity-map.json \
  --output artifacts/scalar-frontier/example-223 \
  --seconds 259200 --branch-seconds 120 \
  --heads 5,7,9 --branch-policy all --sampler native \
  --seed 20261126001 --nice 10 --detach
```

Here $h$ is the number of nonzero multiplier orbits using coordinate zero.
The required row sums are $(-14,3,3,3)$, so $h$ is odd with
$5\leq h\leq71$. `--heads` selects the strata, scheduled in round-robin
order. The `all` branch policy enumerates exact coordinate-zero sign
choices and solves the remaining signs for each choice. The `affine`
policy leaves those signs in a single larger affine system. A completed
branch has only its stated fixed-sign scope.

The native sampler uses exact binary-field arithmetic and a packed color
map. Its tagged xoshiro256** state is seeded through SplitMix64 and is
stored for exact stream replay. `--sampler python` selects the distinct
Python `Random` stream. The two states are not interchangeable; the
sampler choice is included in the run fingerprint. Phase-correct domain
normalization is applied after sampling.

See {repo-file}`exact_parity_sampler.py <experiments/scalar_lift/exact_parity_sampler.py>`
and its native companion for that interface. Their arithmetic,
coordinate conversion, filtering and serialized stream are checked
against independent Python controls.

## Identity coordinate action

The {repo-file}`identity geometry <experiments/scalar_lift/identity_multiplier.py>`
requires $b(\mu x)=b(x)$, fixing every coordinate label. It differs from
the three-cycle action above. If the four signed orbit counts are $a_j$,
the row sums are $(1+ma_0,ma_1,ma_2,ma_3)$ and must have squared norm $p$.
This necessary condition is followed by four exact sign counts and the
full integer correlation checks.

The {repo-file}`223 certificate <certificates/scalar-lift/identity-parity-p223.json>`
validates the full parity map. Its saved minimum-head pattern for
$(-14,3,3,3)$ has 30 free sign bits, compared with 32 in the earlier
three-cycle control. This comparison concerns the recorded affine systems;
neither is a signed witness. A focused three-day campaign uses

```sh
PYTHONPATH=src:. .venv/bin/python -B -m experiments.scalar_lift.identity_t_campaign \
  --config certificates/scalar-lift/identity-parity-p223.json \
  --output artifacts/scalar-frontier/example-identity223 \
  --target-row=-14,3,3,3 --head-offsets=0,2,4 \
  --seconds 259200 --branch-seconds 120 --branch-policy all \
  --sampler native --seed 20261127223001 --nice 10 --detach
```

Offsets are added to $|a_0|$, so this selects $h=5,7,9$.
`--target-row` is repeatable; its omission selects all compatible row types.
The {repo-file}`491 certificate <certificates/scalar-lift/identity-parity-p491.json>`
uses $m=7$, $\mu=138$ and compatible row sums $(1,21,7,0)$, checking every
parity coefficient and expanded-shift comparison. A signed witness would
yield order $4\cdot491\cdot3=5892$. No such witness is claimed.

## Run unrestricted length 167 or 179

The unrestricted runner does not impose multiplier covariance. It cycles
through the four-square row-sum types and the chosen backends. Each job
is one fixed coordinate word, target row and backend attempt.

```sh
PYTHONPATH=src:. .venv/bin/python -m experiments.scalar_lift.unrestricted_t_campaign \
  --config certificates/scalar-lift/unrestricted-parity-p167.json \
  --output artifacts/scalar-frontier/example-167 \
  --seconds 259200 --pattern-seconds 900 \
  --backends sat,cp --seed 20261126010 --nice 10 --detach
```

For length 179, use
`certificates/scalar-lift/unrestricted-parity-p179.json` and a different
output directory. CryptoMiniSat requires the external `cryptominisat5`
executable on `PATH`; OR-Tools is installed by the `search` extra. Use
`--backends cp` when only CP-SAT is available, or `--backends sat` for
CryptoMiniSat alone. Backends run sequentially within one worker, not as
multiple solver threads.

{repo-file}`unrestricted_t_sign_sat.py <experiments/scalar_lift/unrestricted_t_sign_sat.py>`
uses native XOR constraints and exact Boolean adders.
{repo-file}`unrestricted_t_sign_cp.py <experiments/scalar_lift/unrestricted_t_sign_cp.py>`
uses Boolean XOR/product constraints and exact integer equalities. Both
expose the campaign's fixed-axis solve interface and independently verify
positive states before returning them. Complete small controls compare
their solution sets with independent finite enumeration; true-width
controls also exercise the 80/86-bit models.

## What is committed, and what resumes

Each output directory has a single-writer lock and a SQLite ledger using
WAL mode and `synchronous=FULL`. JSON status records are replaced
atomically. The ledger is the durable record; the status file is a readable
summary.

For the three-cycle 223 runner, the ledger stores coordinate patterns, a head-sign
branch cursor, skipped-branch counts, attempts, sampler state and totals.
A completely checked branch advances the cursor. If a per-branch limit
interrupts a branch while session time remains, that branch is recorded
unfinished and may be skipped to explore further branches. Its pattern
is marked `partial`, not complete. A branch interrupted at the session
deadline stays active for replay on resume. The native solver does not
serialize its internal Gray-code computation at every word.

The identity runner commits accepted patterns and sampler state together.
Each unfinished branch remains pending and is retried on resume; it is
never converted into a completed exclusion. Only a pattern with every
required branch complete advances the completed-pattern total.

For the unrestricted runner, insertion of a coordinate/target/backend job
and the updated sampler state are committed together. Completed jobs are
recorded as `excluded`, `unknown` or `witness`. A job still `active` after
a crash is retried; a bounded `unknown` job remains an unfinished search
result. Neither SAT nor CP-SAT search trees are serialized. Retrying an
interrupted solve rebuilds its exact model and restarts that solve.

A successful result is additionally checked by the runner. Its
`pending_witness` record commits in the same transaction as the job or
branch result. On restart, the worker verifies that record again before
publishing `verified-witness.json` and rebuilding `endpoints.json`.
Consequently a crash between the database commit and file publication
does not lose a found design.

```{list-table} Worker files and their roles
:header-rows: 1
:widths: 35 65

* - File
  - Meaning
* - `frozen/`
  - The explicit source set, parity configuration and source manifest used by the detached worker
* - `launch.json`
  - The command and process identifier for a particular launch
* - `work.sqlite3`
  - Committed work and sampler state; keep its WAL companions with a live database
* - `worker.log`
  - Worker/solver output
* - `status.json`
  - Current phase and measured progress
* - `verified-witness.json`
  - A literal signed-coordinate word that passed the independent integer check
* - `endpoints.json`
  - Independently checked source/design identities and literal real-matrix recipes derived from the witness
```

## Frozen sources and controlled resource use

`--detach` copies the runner's explicit `SOURCE_NAMES` and parity map into
`frozen/`, records their digests, and launches from that copy. Changes in
the working checkout therefore do not change a running worker. The ledger
fingerprint binds the configuration, source hashes, seed and search policy;
unrestricted fingerprints also bind the backend sequence and row types.
A different source/configuration identity requires a fresh run directory.
An existing frozen directory rejects a conflicting supplied parity map.

Resume an output directory by rerunning its saved launch command, or the
same runner options against that frozen identity. The three-cycle and
unrestricted time budgets apply to each new session; the identity runner
uses a cumulative budget across resumed sessions. Committed totals carry
across sessions. Remove a
previously requested `STOP` file only when intending to resume the worker.

Both native search and each SAT/CP-SAT invocation use one solver thread.
The runners have separate campaign and branch/pattern limits. Exact
startup validation, model construction and final verification can add
overhead; the recorded elapsed times state what actually ran. A positive
witness can end the search early. `--max-patterns` counts completed coordinate patterns in the 223
runner and checked jobs in the unrestricted runner; partial 223 patterns
do not advance its completed-pattern count. `--nice 10` lowers scheduling
priority; it is not a memory or CPU quota. Use a separate directory and seed for each independent worker and
choose the number of workers to fit the machine.

Creating `STOP` in the output directory requests a stop after the current
bounded branch or solver call returns. Detached execution keeps a process
independent of the terminal; it does not prevent the operating system from
sleeping. The runners make no submission or publication requests.

Generated ledgers, logs and solver caches belong under an ignored output
path such as the examples above, away from frozen releases and literal
source fixtures. Both `.gitignore` and `.ignore` exclude `*.local` and
`.local/`. API credentials and saved administrative parameters remain in
those private files and are not part of a campaign's source whitelist.

## Interpret outcomes at their actual scope

| Recorded outcome | Valid conclusion |
|---|---|
| Exact verified solution | This literal word is a periodic T-design with the recorded row sums; the separate endpoint reconstruction proves the derived real matrices |
| Complete UNSAT/infeasible result | The fixed coordinate word has no signs satisfying the recorded sums, gauges and exact constraints |
| Complete native enumeration without a zero | The fully enumerated affine branch has no solution; a whole coordinate pattern is complete only when all its required branches are complete |
| UNKNOWN, interrupted or partial | That computation did not finish its stated search; it excludes no remaining sign assignment |

Assignment counters describe the native filter's work, including any replay
of interrupted branches. They are not counts of full integer endpoint
checks or necessarily of distinct assignments. A parity-map certificate
proves a coordinate parametrization; it proves no signed T-design by itself.

A completed T-input of length $p$ combines with a verified source quadruple
of width $w$ to give a real endpoint of order $4pw$. A new T-input needs its
literal integer certificate and endpoint proof; a new-order label then
requires the dated literature check. Its new Lean literal bridge must also
be added and successfully checked before claiming formal verification of
that new result. These stages do not alter the manuscript's already
completed thirty-four new orders and nine stated families.

The {repo-file}`low-gap scanner <experiments/scalar_lift/low_gap_target_scan.py>`
also screens both odd coordinate actions at prime lengths, with at most
100 nonzero orbits and the retained literal source widths. Its
{repo-file}`certificate <certificates/scalar-lift/low-gap-target-scan.json>`
recorded nine necessary-arithmetic survivors among 119 pinned frontier
entries below 12000. Subsequent records identify 5164 as known, so the
dated certificate is retained for replay rather than current gap selection.
It constructs no parity map or signed design.
Replay with `python -m experiments.scalar_lift.low_gap_target_scan --check`.

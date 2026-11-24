# Running the forward search

The search starts with completed, verified Butson matrices. It follows
certified scalar restrictions, checks final sign frames against each actual
entry alphabet, and records explicit endpoint recipes. The source Gram,
scalar representation, adjoint metric, and frame proof are separate records.

The implementation is
[`forward_descent_search.py`](../experiments/scalar_lift/forward_descent_search.py).
Its source registry is [documented here](FORWARD_SOURCE_REGISTRY.md), and the
[frame/gauge search](FORWARD_FRAME_SEARCH.md) verifies every actual sign block.

## Registered maps and completed first run

The initial 18-source run evaluated 238 states, with 108 general weighted
frame problems retained as open. It produced 20 distinct real orders:

```
40, 80, 104, 208, 312, 480, 624, 640, 936, 1120,
1280, 1872, 2560, 2808, 4000, 5616, 8000, 8424, 16000, 16848
```

None intersects the attributed unresolved comparison set. These outputs
reproduce known constructions; no priority is claimed. The larger
character sources remain complete Butson inputs, but fail the registered
two-cubic-color frame domains even after row/column root gauges. Other
source-dependent constructions remain outside that scoped exclusion.

The expanded registry also contains the completed sixth-root sources of
orders 3 and 13, so it can connect a search for a new final frame to a
verified starting matrix directly. The current report is
[`forward-descent-search.json`](../certificates/scalar-lift/forward-descent-search.json).
Its counts supersede the historical first-run counts above.
That expanded run verified 20 sources, evaluated 244 states, retained 110
open general-frame records, and saved 59 endpoint recipes with 22 distinct
known real orders. It completed in about four seconds and uses 576 KiB on
disk. The two added real orders are 12 and 52.

The finite comparison cache records 156 attributed unresolved orders,
including the previously considered higher-order candidates. None of the
first 18 source orders reaches these comparison orders through the current
power-of-two total widths. This is a limitation of the registered inputs
and maps; repeating them cannot settle those candidates. A new terminal
frame can have a different width, which is why the remaining-frame search
must be an actual computational branch.

## Persistent state and explicit outputs

SQLite stores verified source hashes, lazy descent recipes, exact gate
results, edges, and terminal recipes. It does not store full descendant
matrices in its frontier. Each evaluation, its children, and its result
commit together. Interruption or a disk-cap failure leaves that evaluation
pending for retry. The database allows one writer and has an explicit page
cap; rollback-journal headroom is reserved. Reports are replaced atomically.

The registered graph includes both distinct generator conventions at
repeated-prime steps, and all registered partial frames under global
phase/Galois transformations. At initial sources it additionally checks
row/column gauges with a bounded dense solver or a certified implicit
obstruction. Timed-out or omitted gauge searches remain explicitly scoped.

General weighted restrictions are recorded with their integral trace
metric and verified companion-action identity. Finishing the registered
maps does not exhaust the choice of final frames or bases. A signed odd
scalar relation, when present, excludes a uniform linear sign frame only
for that actual alphabet; it does not silently exclude changed row/column
gauges or factor-dependent constructions.

`endpoint_entry(database, endpoint_id, row, column)` returns a literal sign.
It checks the current source hash, scalar edge hashes, and frame proof
against the saved recipe before evaluating the nested blocks. Tests
materialize an emitted H40 and independently check its complete integer
Gram.

## Commands

```sh
PYTHONPATH=src:. ./.native-env/bin/python -m experiments.scalar_lift.forward_descent_search \
  --database experiments/scalar_lift/search-state/forward-descent-v2.sqlite3 \
  --seconds 60 --max-db-mib 256 \
  --report certificates/scalar-lift/forward-descent-search.json
```

The same command resumes. A source or frame catalogue change requires a new
checkpoint identity, so a completed old cursor cannot silently omit added
branches. For detached execution, use the existing supervised launcher
with `--engine forward`. `status`, `stop`, and `reap` use its job directory.

The source-first final-frame search uses a separate `forward-frame` engine
and checkpoint. It fixes a verified sixth-root source and searches exact
sign-frame equations. See its report for each disjoint partition and solver
status. `UNKNOWN` is an open case, never a proof of infeasibility.

The active bounded batch was launched with:

```sh
PYTHONPATH=src:. ./.native-env/bin/python -m experiments.scalar_lift.vanishing_search_job start \
  --engine forward-frame --job-dir runs/forward-frame \
  --database experiments/scalar_lift/search-state/forward-frame-t223.sqlite3 \
  --seconds 3600 --max-db-mib 256 \
  --report certificates/scalar-lift/forward-design-frame-search.json

PYTHONPATH=src:. ./.native-env/bin/python -m experiments.scalar_lift.vanishing_search_job status \
  --job-dir runs/forward-frame
```

It covers 192 disjoint cases, at most 15 native seconds per case, with one
solver worker. Model construction and verification have their own overhead;
the supervisor enforces the one-hour overall bound. The process may finish
earlier when all scheduled attempts finish or an exact witness is verified.
See [the frame equations and partition](FORWARD_DESIGN_FRAME_SEARCH.md).
The SQLite checkpoint preserves case results and candidate witnesses; it
does not serialize the CP solver's internal search tree. The batch runs on
the local machine and does not prevent sleep.

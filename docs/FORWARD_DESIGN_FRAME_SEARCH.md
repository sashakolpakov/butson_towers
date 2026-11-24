# Forward search for a new final scalar frame

This worker begins with the **complete verified** unreal sixth-root Butson
matrix of order 13. It searches for a new order-892 scalar action and frame.
A successful branch gives an explicit Hadamard matrix of order
\(13\cdot892=11596\), listed as unresolved in the dated primary frontier.
The separate complete order-3 source would also give order 2676. No such
large frame or new order is currently claimed.

The registry stores globally shifted sources. The worker explicitly adds
one to each exponent modulo six, checks the complete source Grams again,
and records both the registry and recovered-source hashes. The recovered
alphabet is contained in \(\{\zeta_6,\zeta_6^2,-\zeta_6,-\zeta_6^2\}\).

For a periodic one-hot quadruple \(b\) of length \(p\), let \(a=D_4b\),
where \(D_4\) is the real skew conference matrix of order four. Zero summed
autocorrelations give the orthogonal design \(OD(4p;3p,p)\). The exact
Goethals–Seidel array supplies

\[
P=GS(a-b),\qquad Q=GS(a+b),\qquad
PP^T=QQ^T=dI,\quad PQ^T+QP^T=dI,\quad d=4p.
\]

Thus \(R=QP^T/d\) is orthogonal and satisfies \(R^2-R+I=0\).
The frame \(Z=P-Q\) has \(ZZ^T=dI\), and the used scalar blocks are
\(RZ=P,R^2Z=Q,R^4Z=-P,R^5Z=-Q\). This is a forward scalar restriction
of an already completed source. The action and frame are recovered jointly
from a completed design; an unfinished Butson matrix is not a starting input.

The default declared chart is the existing multiplier-order-3 T223 model.
Its three normalized row-sum types are
\((-14,3,3,3),(-5,6,9,9),(13,3,3,6)\). Fixing the eight-state values at
orbits 2 and 3 partitions this normal form into exactly \(3\cdot8^2=192\)
disjoint cells. Each initially receives one deterministic seed-0 attempt
of at most 15 seconds, with one CP-SAT worker. This is a restricted design
chart, not an exhaustive search over all possible order-892 frames.
Optional T491 contributes 64 cells and outputs 25532 and 5892.

SQLite records the exact chart, source identities, implementation hashes,
OR-Tools version, compressed base model, cell-model hashes, solver settings,
statuses and witnesses. A solver `INFEASIBLE` result closes only that exact
cell; it is labeled as a solver conclusion without an exported standalone
UNSAT proof. `UNKNOWN` remains open. Resume continues unattempted cells;
it does not repeatedly launch the same unresolved cell and seed.
`--retry-seed` explicitly schedules a new seed for open cells. A killed
attempt is recorded as interrupted because the native search tree was not
saved. Raw positive candidates are committed before independent verification,
so an interrupted verification can resume without solving the cell again.

Every positive candidate must pass all expanded integer autocorrelations
and the complete sign-pair Gram identities. The large output is stored as
a lazy entry formula, avoiding a dense order-11596 allocation. Its full
Gram follows from the independently verified source and sign-pair laws.
The worker stops after its first verified frame so that the actual output
and prior construction status can be reviewed.

Run through the existing bounded launcher with `--engine forward-frame`.
The engine also accepts `--database`, `--seconds`, `--report`,
`--max-db-mib`, `--cell-seconds`, `--include-491`, `--retry-seed` and
`--attempts`. The solver budget is checked per cell; the detached supervisor
supplies the overall wall-clock watchdog. No sustained worker is started by
the module import or its tests.

Tests use a complete T7 control and independently verify the actual output
matrices of orders 364 and 84. They also check the finite partition,
checkpoint locking, unknown-result retention, explicit retry accounting,
and recovery of interrupted attempts and saved positive candidates.

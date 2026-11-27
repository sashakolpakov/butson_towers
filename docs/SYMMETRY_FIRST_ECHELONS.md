# Symmetry-first production in three order ranges

All search campaigns are **paused for strategy review** as of 27 November
2026. The ten most recent workers produced no verified source before they
were stopped. Their complete inputs, seeds, checkpoints and frozen engines
are preserved in a losslessly verified private archive. The current approach
took too long; these cases remain available for another attempt.

## What made the completed constructions productive

The successful T631 run took 2347.74 seconds. Its order-21 multiplier left
30 nonzero residue orbits. A finite-field parameterization generated
coordinate patterns satisfying the binary correlation conditions, after
which a normalized affine system left 13 free sign bits per pattern.
The exact modular filters and final integer test completed the source.
These are search times after the engine and reduction had been developed;
they are not a timer for the whole research project.

All 34 new real Hadamard orders in the paper reuse that one T631 input.
The later productivity came from combining it with more complete scalar
sources, including explicit recurrences, rather than finding a fresh
T-input for each order. The next campaign should therefore rank a missing
ingredient by both its reduced search space and the number of gaps it
would fill with already verified partners.

## The ranking error to correct

A large multiplier and a small sign system do not by themselves make a
search small. In the paused dense Williamson route, four choices remain
for the coordinate index on each nonzero orbit. With $r$ such orbits,
the unfiltered space has $4^r=2^{2r}$ coordinate patterns, even though the
signs are then forced. Thus the length-463 chart with $r=21$ starts with
$2^{42}$ patterns before further filtering. Eliminating its signs did not
supply the finite-field parameterization used for T631.

For every proposed source, record the multiplier and its coordinate action,
the number of orbits, the number of binary parameters after every exact
reduction, and the dimension of any remaining sign system. Distinguish a
proved parameterization from an estimate based on generic rank. Record
all surviving row-sum and orbit-count conditions. Benchmark the same
formulas on a complete positive control before estimating cost at a new
length. An integer norm test remains necessary at the end.

## Three target ranges

The reproducible [frontier audit](https://github.com/sashakolpakov/butson_towers/blob/main/docs/ECHELON_FRONTIER.md) starts with primary
construction records and removes positive constructions before considering
a search. Its compact snapshot retains hashes rather than a long generated
list of larger candidates.

| Range | Source-listed gaps | Additional orders surviving the current positive construction gates |
|---|---:|---:|
| Below 10,000 | 81 | 0 |
| 10,000 through 50,000 | 37 | 4,102 |
| Above 50,000 through 100,000 | 0 | 5,981 |

The 81 entries come from the current Hadamard Matrix Explorer table. The
37 further entries come from the author-maintained frontier through
11,996. The other survivors need individual literature checks; surviving
this finite set of construction gates alone does not establish unknown
status. Keep that check with the proposed endpoint, rather than labelling
the entire larger pool unknown.

The order ranges are reporting groups, not separate searches. One small
source can have useful endpoints in several groups. Start with the smallest
gaps when their reduced costs are comparable, and retain a cheaper source
in a higher group when it offers a substantially smaller search or more
outputs.

The bounded T-input audit screens both odd multiplier actions through
74 nonzero orbits. It cross-checks two independent enumerations of the
integer row-norm conditions, considers all 46 retained finite Williamson
widths and the admitted finite sixth-root sources, and includes powers of
two in the final orders. Its selected outputs include the source-listed
gaps 2676, 4460, 5892, 9260 and 11596. This is a bounded selection rule;
the other frontier entries remain available to other source classes or
better symmetry reductions.

Replay the compact inventories without starting a search:

```sh
python -m experiments.scalar_lift.echelon_frontier --check
python -m experiments.scalar_lift.echelon_multiplier_scan --check
```

## Shortlist before any restart

The following examples pass the exact row-norm and signed orbit-count
conditions of the stated dense symmetric Williamson charts. Here $p$ is
the component length, $m$ the multiplier order, and the multiplier cycles
coordinates 1, 2 and 3 while fixing coordinate 0 unless stated otherwise.
The current complete T-inputs supply the endpoint factors.

| Proposed Williamson source | $m$ | Nonzero orbits $r$ | Unfiltered coordinate bits $2r$ | Conditional endpoints |
|---|---:|---:|---:|---|
| $p=463$, identity action | 22 | 21 | 42 | 9260; also 16668 and 20372 |
| $p=919$ | 54 | 17 | 34 | 25732 and 40436 |
| $p=1801$ | 72 | 25 | 50 | 36020, 50428, 64836, 79244 and 93652 |
| $p=1951$ | 78 | 25 | 50 | 23412 and 54628 |
| $p=2269$ | 84 | 27 | 54 | 45380, 63532 and 99836 |

Only 9260 in this table belongs to the current below-10,000 frontier; the
higher entries currently carry the larger-pool evidence class. No source
or endpoint in this table has been constructed by this audit. In particular,
the smaller raw coordinate count at 919 is a reason to examine its exact
algebraic reduction, not a runtime prediction or an existence proof.

The parallel T-input audit retains the successful fixed-coordinate and
three-cycle actions. For example, T919 with multiplier order 27 has
34 nonzero orbits and a verified binary parameterization; T1453 with
multiplier order 33 has 44. Their dense Williamson alternatives and scalar
partners must be compared on the same record before choosing an engine.
The [earlier chart analysis](https://github.com/sashakolpakov/butson_towers/blob/main/docs/HADAMARD_FRAME_TARGET_SCAN.md) proves that simply
increasing the prime cannot produce another new-order case with at most
30 orbits in that same T-input three-cycle class and the retained partners.
An effective improvement must change the source, the symmetry action or
the algebraic reduction.

Before restarting, complete the source-status check, prove or replay the
chosen symmetry and parity reductions, and list every compatible endpoint
through 100,000. Prioritize sources with a small actual parameter space
and several remaining endpoints. When a source is found, replay its full
integer norm, then construct and check every listed endpoint from that
source before allocating another long search.

## Exact endpoint obligations and current stopping point

For a Williamson source of length $p$ and a complete T-input of length $t$,
the orthogonal design has order $4t$ over $\mathbb Q[C_p]$ and Gram scale
$4tp$. Under the full regular scalar action, the normalized trace metric
in the group basis is $I_p$ and conjugation becomes transpose. The actual
blocks are signed copies of the four symmetric sign circulants. Thus the
endpoint is flat, with order and Gram scale both $4tp$. Retain both factors
$\mathbb Q\times\mathbb Q(\zeta_p)$; do not replace the regular action by
only its cyclotomic factor. Other scalar routes must retain their own
source, action, metric, frame factor and actual entry alphabet.

The mathematical engines, positive controls, certificates and manuscript
remain available. Generated unsuccessful-run artifacts have been removed
from the active tree, with complete retry states privately archived.
There are no active workers and no automatic restart queue. No new matrix,
family or manuscript claim results from this strategy audit.

(results)=
# Completed results and construction targets

The current manuscript proves **34 new real Hadamard orders** and **nine
infinite families**, alongside the separate new periodic T127 ingredient.
The source inventory reconstructs **31 known orders**, giving **65 distinct
completed real orders**. The submitted 25 November package is preserved
with its earlier 18-order scope; the current PDF includes the later results.
Version 1 is [published on HexagonMath](https://hexagonmath.org/2611.01070v1);
version 2 was submitted on 26 November 2026 (Toronto) as a revision of that
article. The {repo-file}`version 2 receipt <release/hexagon_2611_01070_v2/submission-status.json>`
records service status **under review**. Version 1 remains the published record.
The finite claims are existence constructions; novelty labels refer to
the dated literature ledger.

All 31 positive prior constructions and the order-885924 matrix comparison
are given separately in {ref}`known-orders`.

The complete finite list is in the manuscript's main theorem and the
{repo-file}`source inventory <docs/SOURCE_INVENTORY.md>`.

The source recipes and individual status evidence are retained in the
{repo-file}`dated construction ledger <docs/HADAMARD_HIGHER_ORDER_RECORDS.md>`
and its {repo-file}`machine-readable union <certificates/scalar-lift/completed-order-inventory.json>`.
A repeated construction of the same order contributes once to these counts.

## Which objects are new

The new finite T631 input supplies the reusable disjoint design. The small
Butson and Williamson inputs and the classical source recurrences retain
their original attribution. The output orders and the finite T127 ingredient
have their own novelty claims. In particular, a classical source family and
a new output family are different statements; see {ref}`families`.

The T127 input gives a real order-1524 matrix by width-three substitution,
but that real order is known. Its separate T-matrix result is not counted as
a new real Hadamard order. The project does not classify all matrices at any
of the high orders, nor claim inequivalence from a changed search seed.

The group-algebra restrictions can use known sources even when their complex
field components are nonflat. What matters is their exact norm together with
the actual reconstructed sign blocks; see {ref}`construction-group-algebra`.

## The current lower-order ingredients

A completed periodic T-input of length $p$ combines with a verified source
width $w$ to give order $4pw$. The following ingredients remain construction
targets until their actual finite witnesses and endpoint checks exist:

| Target order | T-input length $p$ | Source width $w$ |
|---:|---:|---:|
| 2676 | 223 | 3 |
| 4460 | 223 | 5 |
| 4676 | 167 | 7 |
| 7876 | 179 | 11 |
| 8684 | 167 | 13 |
| 9308 | 179 | 13 |
| 9260 | 463, identity multiplier | 5 |
| 11596 | 223 | 13 |

Several rows share a missing ingredient. There is one length-223 input
problem for the three rows using 223. Increasing a run's elapsed time,
checking a parity parametrization, or constructing a SAT model does not
complete any of these rows. {ref}`campaigns` explains how the runtime reports
record actual witnesses, exclusions and unfinished work.

The {repo-file}`construction scope <docs/RESULTS_SCOPE.md>` defines the
project boundary. Finite phase/torsion classifications in the earlier
`complex-hadamard` repository have separate objects and claims.

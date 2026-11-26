(results)=
# Completed results and construction targets

The submitted manuscript dated 25 November 2026 contains eighteen real
Hadamard orders labelled new to the best of the author's knowledge, six
infinite output families, and a separate new periodic T127 ingredient. Its
main campaign also reconstructs twenty known orders, giving 38 distinct
completed real orders. The finite claims are existence constructions;
novelty labels refer to the dated literature ledger.

The eighteen orders are

**7572,17668,22716,32812,42908,58052,73196,78244,93388,98436,
108532,128724,153964,229684,295308,386172,4772884,13001124.**

The source recipes and individual status evidence are retained in the
{repo-file}`dated construction ledger <docs/HADAMARD_HIGHER_ORDER_RECORDS.md>`
and its {repo-file}`machine-readable record <certificates/scalar-lift/higher-order-records.json>`.
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
| 11596 | 223 | 13 |

Several rows share a missing ingredient. There is one length-223 input
problem for the three rows using 223. Increasing a run's elapsed time,
checking a parity parametrization, or constructing a SAT model does not
complete any of these rows. {ref}`campaigns` explains how the runtime reports
record actual witnesses, exclusions and unfinished work.

The {repo-file}`construction scope <docs/RESULTS_SCOPE.md>` defines the
project boundary. Finite phase/torsion classifications in the earlier
`complex-hadamard` repository have separate objects and claims.

# Provenance of targets through order 100,000

The target inventory has three disjoint ranges: orders below 10,000,
orders from 10,000 through 50,000, and orders above 50,000 through 100,000.
The public records checked on 27 November 2026 have different coverage.

The [Hadamard Matrix Explorer, version 61](https://www.math74.com/Hadamard/index.html)
covers orders through 10,000. Its index and
[unknown-order table](https://www.math74.com/Hadamard/unknown.html)
independently list the same 81 orders without a construction in the program.
The source is an implementation inventory, so its word “unknown” does not
establish nonexistence or rule out an unrecorded published construction.

[Kumar's current author frontier](https://github.com/vatskishan03/Hadamard/blob/0a73016a2a6e80bfd4a36b530ad1e5ed02387947/FRONTIER.md)
contains 120 orders through 11,996. Its 37 orders above 10,000 are retained
as source-listed unresolved targets. The file retrieved on 27 November has
the same SHA-256 digest as the repository's earlier pinned copy. Below
10,000, the newer Explorer inventory replaces this older list.

Neither source supplies a complete current inventory through 100,000.
The larger ranges are therefore screened by positive construction tests:
all positive Explorer entries, the completed Butson Towers ledger, direct
Paley I and II constructions for prime powers, the recorded Williamson-type
397 construction with each odd T-length at most 19, and all Kronecker products
of these orders within the bound. These are sufficient construction tests.
An order that passes none of them remains unassessed by these tests.
Further construction theorems or individual literature records can remove it.

| Order range | Current Explorer gaps | Additional author-listed gaps | Other candidates after these positive tests |
|---|---:|---:|---:|
| Below 10,000 | 81 | 0 | 0 |
| 10,000 through 50,000 | 0 | 37 | 4,102 |
| Above 50,000 through 100,000 | 0 | 0 | 5,981 |

These counts describe the starting pool before symmetry selection and
before additional positive construction tests. They are not counts of new
matrices or a claim that 10,201 Hadamard orders are unknown.

For Paley I the tested parameter is an odd prime power `q = 3 mod 4`,
giving order `q+1`; for Paley II it is an odd prime power `q = 1 mod 4`,
giving order `2(q+1)`. The
[2026 survey, Theorem 2 and the preceding product construction](https://arxiv.org/html/2610.27175v2#S2)
state these results. The Williamson-type 397 case follows from the Paley symmetric conference
matrix of order 198, [Seberry–Yamada, Corollary 29(ii), p.120](https://documents.uow.edu.au/~jennie/WEBPDF/139_1990.pdf#page=24),
and [Cooper–Wallis, Theorem 3 and Corollary 4, p.274](https://doi.org/10.1017/S0004972700045019).
It removes the already known order 14,292 as well as the other corresponding
products. The repository previously recorded this case in the
[compact target scan](HADAMARD_FRAME_TARGET_SCAN.md).

Each positive test retains the parameter and its prime
and exponent, or the two smaller factors for a Kronecker product. One
increasing pass over orders computes the whole bounded product closure,
because both nontrivial factors of an order are smaller than that order.

The compact [snapshot](../experiments/scalar_lift/echelon_frontier_snapshot.json)
stores input digests, source provenance, the 81 and 37 source-listed orders,
and the counts and hashes of the larger generated pools. It omits the large
pool lists. The code reconstructs those lists without network access or a
matrix search:

```python
from experiments.scalar_lift.echelon_frontier import build_report, candidate_status

report = build_report(include_candidate_pool=True)
for order in report['candidate_pool']:
    evidence = candidate_status(order, report)
    # A symmetry selector can now attach its exact orbit and sign dimensions.
```

Replay the saved counts and hashes with:

```sh
python -m experiments.scalar_lift.echelon_frontier --check
```

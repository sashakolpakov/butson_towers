# Butson Towers

Exact scalar restrictions of verified Butson-related designs, with transported
trace metrics and certified flat endpoints.

Start with a complete source and its exact Gram identity. Restrict scalars
through number fields or products of number fields, carrying conjugation,
the metric, the actual entry alphabet and the Gram scale. Intermediate
matrices may be weighted and nonflat. A final compatible frame and entry
test certify a real Hadamard or lower-field Butson endpoint.

Read the [foundations and attribution](docs/FOUNDATIONS.md),
[scalar-reduction identities](docs/SCALAR_REDUCTION.md), and
[active engine guide](experiments/scalar_lift/README.md).

## Completed results

The main construction campaign has **38 distinct completed real Hadamard
orders: 18 new to the best of the author's knowledge and 20 known recovered
orders**, following the [dated literature ledger](docs/HADAMARD_HIGHER_ORDER_RECORDS.md).
The new orders are

**7572, 17668, 22716, 32812, 42908, 58052, 73196, 78244, 93388, 98436,
108532, 128724, 153964, 229684, 295308, 386172, 4772884 and 13001124.**

The complete [H7572 construction](docs/HADAMARD_7572.md) supplies periodic
T-matrices of order 631, a compatible rational sixth-root action and a sign
frame of width 2524. Further complete sources give
[frame-reuse endpoints](docs/SCALAR_REDUCTION_REPEAT_TARGETS.md),
[unreal-source tripling endpoints](docs/TRIPLED_UNREAL_SCALAR_ENDPOINTS.md), and
[Williamson scalar endpoints](docs/WILLIAMSON_ORDER_RECORDS.md).
Large outputs have literal entry recipes and complete algebraic Gram proofs;
they do not require large dense matrix files.

The [five constructive infinite families](docs/SPENCE_SCALAR_FAMILIES.md)
provide literal Hadamard matrices for every integer `r>=0` at orders

| Family parameter | Real endpoint order |
|---:|---|
| 13 | `17668*13^r` |
| 17 | `22716*17^r` |
| 61 | `78244*61^r` |
| 101 | `128724*101^r` |
| 121 | `153964*121^r` |

Existence holds throughout each family. The finite members audited in the
ledger have individual known/new labels; the infinite families do not
claim every member's order as new.

There is also a separate new [periodic T127 ingredient](artifacts/t-matrices-127/README.md).
Its real endpoint H1524 is a known-order control. The
[search analysis](docs/HADAMARD_FRAME_SEARCH.md) explains the finite-field
parity coordinates and exact modular sign lifts. These are existence
results, with their scope recorded in [the results index](docs/RESULTS_SCOPE.md).

## Foundations

The immediate starting reference for this project was Egan–Ó Catháin,
[*On complete morphisms of Butson Hadamard matrices*](https://arxiv.org/abs/2610.28394),
20 October 2026. The older foundations include their
[*Morphisms of Butson classes*](https://arxiv.org/abs/1707.08815), 27 July 2017,
and Ó Catháin–Swartz,
[*Homomorphisms of matrix algebras and constructions of Butson-Hadamard matrices*](https://arxiv.org/abs/1904.10771),
24 April 2019. The [foundations note](docs/FOUNDATIONS.md) distinguishes
these contributions from the present constructions.

“Butson Towers” is this project's descriptive name. Intermediate stages
need not have Butson entries, and a tower may include a product of
cyclotomic fields with its exact reconstruction map.

## Verification and formalization

The completed constructions have executable exact proofs, source hashes,
literal matrix recipes and certificates. The [active engine guide](experiments/scalar_lift/README.md)
gives the replay commands. Run the repository suite with

```sh
./scripts/run_tests.sh
```

The [Lean project](docs/FORMALIZATION.md) proves the **eighteen announced
orders and all five infinite families**, from literal finite witnesses.
It verifies T127/T631 correlations in the kernel, the regular scalar
action, orthogonal-design substitution, the actual Spence recurrence,
and the order-13 cubic recurrence. Build and audit all final statements
with `./scripts/check_formal.sh`.
A separate [Lean reproducibility package](release/lean_verification/README.md)
contains the locked sources, checksum verifier and successful extracted replay.

## Manuscript and HexagonMath tooling

The [announcement and source package](release/hadamard_7572_announcement/README.md)
are maintained alongside the construction proofs. Local administrative
defaults live in `hexagon_parameters.local`; credentials live separately
in `hexagon_key.local`. Both files remain private, outside Git and release
archives. Local preparation and saved defaults do not authorize a remote
draft, submission or publication.

Audit the Git index, reachable history and release ZIPs without printing
private values:

```sh
.venv/bin/python scripts/check_private_files.py --root .
```

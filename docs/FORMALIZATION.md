# Formal verification

The Lean project proves existence of the thirty-four Hadamard orders in the
manuscript and every member of its nine stated families: eight from Spence sources
and one from the cubic recurrence. The statements
construct rational sign matrices with the exact Gram identity and the
specified positive order; no dense high-order Gram calculation is needed.

## Final statements

`ButsonTowers.HadamardExists n` means that a matrix indexed by `Fin n`
has entries in `{1,-1}` and satisfies `H * H.transpose = n • 1`, with
`0 < n`. The named `hadamard_<order>` theorems and all nine stated all-iteration
family theorems are collected in `Endpoints.lean`, `NearWilliamsonEndpoints.lean`,
`AdditionalSpenceFamilies.lean` and `AdditionalResults.lean`.
Novelty labels belong to the dated literature ledger; the formal
statements prove existence.

The nine stated families are

| Character order | Hadamard order, for every `r : ℕ` |
|---|---|
| 5 | `7572 * 5^r` |
| 13 | `17668 * 13^r` |
| 17 | `22716 * 17^r` |
| Cubic recurrence | `32812 * 3^r` |
| 61 | `78244 * 61^r` |
| 101 | `128724 * 101^r` |
| 121 | `153964 * 121^r` |
| 137 | `174156 * 137^r` |
| 109, stage r+1 | `15131380 * 109^r` |

The finite data are literal Lean objects. Sound finite checkers connect
integer computations to rational periodic correlations, matrix products
and sign tests. Their successful evaluations are proved by kernel
reduction, without an external evaluator. Balanced lookup trees make the
T631 all-shift calculation practical. A source-file hash is never used as
a mathematical premise.

The [manuscript-to-Lean audit](MANUSCRIPT_LEAN_AUDIT.md) follows the
rewritten proof in order and records the exact scope and conventions.

## Proof chain

* `SixthRootPair.lean`: the rational sixth-root matrix action, its
  polynomial and adjoint identities, and multiplication/conjugation in
  scalar coordinates, from the explicit paired-matrix assumptions.
* `ScalarFrame.lean`: inverse-metric adjoints, weighted Gram transport,
  full framed-matrix Grams, preservation of zero, separate flatness and
  the final equality between Gram scale and expanded order.
* `GroupAlgebra.lean` and `FiniteGroupAlgebra.lean`: the faithful rational
  finite-abelian-group regular action, multiplication, inversion as
  transpose, normalized trace and the identity trace metric in the
  inverse-group basis.
* `GoethalsSeidel.lean` and `OrthogonalDesign.lean`: the actual block-array
  Gram identity, all four coefficient norms and six mixed identities,
  and the disjoint complete signed supports furnished by a T tuple.
  `ReflectionConvention.lean` identifies the reversal used in Python
  by an explicit coordinate translation for periods 127 and 631.
* `OrthogonalSubstitution.lean`: substituting an amicable sign quadruple
  gives a sign matrix with the required Gram. The scale/order equality
  remains explicit. Pairwise amicability suffices for the substituted
  matrices.
* `SpenceMatrix.lean` and `SpenceHadamard.lean`: the literal tensor
  recurrence, sign entries, symmetry, commutativity, both strong norms,
  and the growing index type and its cardinality at every iteration.
* `UnrealTripling.lean`: the cubic recurrence for the order-13 unreal
  pair, including signs, amicability and norm. Two iterations give the
  order-117 substitution used for order 295308.
* `Data/`: complete T127/T631 witnesses, all 27 original Williamson inputs,
  twelve additional compatible quadruples, eleven strong starters and their
  Paley matrices, and the order-13 unreal pair.
  For character order 121 the index group is `ZMod 11 × ZMod 11`.
* `Existence.lean`, `Construction.lean` and `Endpoints.lean`: reindexing
  the constructed matrices to `Fin n` and instantiating the concrete
  sources in the final existence and family statements.

The trace and basis-change identities in `TraceTowers.lean` provide the
metric formulation of successive scalar restrictions. The endpoint
proofs use the explicit regular representations, so they do not require
a separate cyclotomic decomposition theorem or an unproved assertion
that a group algebra is a field.

## Reproduction and trust record

The project pins Lean **4.32.0** and mathlib commit
`81a5d257c8e410db227a6665ed08f64fea08e997`, with all nine dependency
revisions locked. From the repository root:

```sh
./scripts/check_formal.sh
```

This builds the project and checks three exact axiom-report lists:
85 original statements, ten rational sixth-root matrix identities, and
48 additional source, endpoint and family statements. The 143 selected
reports include the finite T certificates, construction theorems, every
announced order and all announced families. The only permitted dependencies are
Lean's standard `propext`, `Classical.choice` and `Quot.sound` axioms.
There are no project-specific axioms, proof placeholders or native
computation certificates. Every numeric T block certificate is axiom-free;
the coverage and semantic bridges use only the standard axioms.

The [current source-only package](../release/additional_source_verification/README.md)
contains all current proof dependencies, the literal input bindings and all
three strict audit sets. The earlier
[source-only reproducibility package](../release/lean_verification/README.md)
contains the complete proof sources, dependency lock, exact audit list
and a checksum verifier. Its recorded replay builds a fresh extraction
using a matching pinned dependency cache. A cold download route and a
full dependency source-build route are also provided.

The [sixth-root supplement](../release/lean_sixth_root_pair/README.md)
provides a separate source-only reproducer for those ten identities.
The original complete endpoint archive remains byte-for-byte unchanged.

## Conventions

The regular representation has entries `rho(f)[i,j] = f(j-i)`; its
involution inverts group elements. The normalized identity-coefficient
trace gives metric `I`. An unnormalized regular trace would give `|G|I`.
The group algebra may have zero divisors, and no field cancellation is
used.

For a metric identity `rho(a)^T G = G rho(star a)`, row-Gram transport
uses `G^-1`. Nondegeneracy suffices for that identity; positivity and a
flat endpoint require their own hypotheses. Weighted intermediates
need not have root entries. The final sign proof and equality with the
expanded order are separate from inherited orthogonality.

The formal scope is the mathematical construction and its concrete
witnesses. Literature priority, the correctness of search software and
the general number-field factor descriptions are independent of these existence
proofs.

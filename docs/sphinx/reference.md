(reference)=
# Objects, interfaces and data contracts

The interfaces below implement the objects used in the preceding arguments.
Importing them does not start a search. A constructor that stores a matrix
recipe does not by itself certify arbitrary supplied arrays: use the named
verifiers, retain their literal inputs, and keep the returned evidence with
an output. The {ref}`getting-started` example shows this sequence in full.

(reference-restriction)=
## Exact coefficients and a declared scalar edge

`RationalArray` stores a common positive denominator and an integer numerator
array. When representing a lower-field scalar, its last axis contains
coefficients in that field's declared power basis. `regular_matrix` expands
these entries into rational multiplication matrices, which is useful for
independently checking identities over the lower field.

A `RelativeRestriction(Q, q)` describes the declared root labels and basis for
$\mathbb Q(\zeta_Q)/\mathbb Q(\zeta_q)$. The constructor verifies the basis,
all scalar powers, conjugation, trace metric and its inverse. Its `width`
is the relative degree. This implementation requires positive integer labels
with $q\mid Q$. A supplied `basis_exponents` has exactly `width` entries,
starts with exponent zero and must give a nonsingular combined rational
basis. Root labels are retained even when their cyclotomic
fields have a smaller conductor. For example the labels6 and3 generate the
same field, but their exponent alphabets are different data.

```{eval-rst}
.. automodule:: experiments.scalar_lift.relative_scalar_restriction
   :members: RationalArray, regular_matrix, conjugate_transpose, RelativeRestriction, edge, restriction_in_basis
   :member-order: bysource
   :undoc-members:
```

An edge's `metric` is the unnormalized relative trace form. The ordinary row
Gram uses its inverse; the transported column form uses the metric itself.
The derivation and a checked example appear in {ref}`restriction-row-column`.
The JSON metadata retains both matrices, the combined rational basis, scalar
actions, trace automorphisms and a definition hash.

## A whole operator with its incoming form

The `ExactNormState` interface transports an existing nonidentity form.
`order` is the current operator's side length over the declared field;
`source_order` is the fixed norm scale admitted from the original Butson
source. The latter does not grow merely because an edge expands dimensions.
`relative_degree` records the accumulated degree. `digest` binds the complete
state, including its incoming form and provenance.

```{eval-rst}
.. automodule:: experiments.scalar_lift.exact_norm_descent
   :members: ExactNormState, from_butson, positive_control
   :member-order: bysource
   :undoc-members:
```

`from_butson` checks the actual complete source before constructing its root
state. `restrict` verifies the child with the transferred form. `endpoint`
audits supplied rational input/output bases and checks normalization and signs
separately; it does not solve for such bases. These are the contracts used
in {ref}`restriction-towers` and {ref}`restriction-flatness`.

(reference-endpoints)=
## T-inputs, coefficient matrices and literal endpoints

A generating-row array has shape $(4,p)$, exact integer dtype and entries in
$\{0,1,-1\}$. `verify_t_rows` checks every support position and every cyclic
shift. `verify_orthogonal_design` checks all diagonal and paired polynomial
Gram coefficients and the entrywise support partition. Its coefficient array
has shape $(4,4p,4p)$.

A source-row array has shape $(4,w)$ and sign entries. `verify_williamson`
checks symmetry, amicability, the summed matrix Gram and the
cyclic-algebra norm. The source array then has shape $(4,w,w)$.
Commutation follows from the circulant construction; pairwise amicability
is also checked explicitly.

```{eval-rst}
.. automodule:: experiments.scalar_lift.williamson_scalar_endpoints
   :members: verify_t_rows, coefficient_sequences, verify_orthogonal_design, verified_t_design, verify_williamson, verified_williamson, ScalarEndpoint, ImplicitScalarEndpoint, row_controls
   :member-order: bysource
   :undoc-members:
```

`ScalarEndpoint(X,A)` stores these arrays and computes an exact entry or row
of $\sum_jX_j\otimes A_j$. It assumes that callers obtained them from the
verifiers. `ImplicitScalarEndpoint` replaces the dense source array by an
entry/row recipe. Neither class allocates the full endpoint at construction.
`row_controls` checks selected rows and pairs; the full algebraic Gram proof
comes from the verified input identities. A dense full audit is an additional
operation with its own recorded scope.

## Symmetry-selected Williamson sources

The {repo-file}`current source model <docs/WILLIAMSON_SYMMETRY_SEARCH.md>`
uses even multiplier orbits. `chart` checks every expanded coefficient
table; `coefficients` derives the unique signs for a supplied axis word.
`verify_source` requires literal integer sign rows, symmetry and every
integer correlation. Its result is a complete source certificate, whereas
`select_charts` checks necessary arithmetic before a construction search.

```{eval-rst}
.. automodule:: experiments.scalar_lift.williamson_symmetry
   :members: chart, coefficients, expand, verify_source, coefficient_record, select_charts
   :member-order: bysource
```

The CP model retains the exact pair contributions and the forced-sign XOR
identities. A positive result is independently expanded and verified.
The native launcher copies its coefficient input, source and binary into
an immutable snapshot; its checkpoint retains the generator state and
current annealing word. Random and annealing evaluation counts can include
repeated words. Both interfaces accept explicit time and seed arguments.

```{eval-rst}
.. automodule:: experiments.scalar_lift.williamson_axis_cp
   :members: integer_nets, build, solve
   :member-order: bysource
```

```{eval-rst}
.. automodule:: experiments.scalar_lift.williamson_axis_native
   :members: row_targets, serialize, prepare
   :member-order: bysource
```

## Recursive source entries

The additional compatible finite sources use a separate verifier with an
explicit first-matrix reversal flag. It checks all full matrix norms and
pairwise identities before constructing a `ScalarEndpoint`.

```{eval-rst}
.. automodule:: experiments.scalar_lift.near_williamson_scalar_endpoints
   :members: literal_matrices, verify_quadruple, verified_source, endpoint, certificate
   :member-order: bysource
```

`SpenceSource(q,r)` validates the supported parameter and its starter, then
evaluates recursive sign entries without materializing all source rows.
`order` is $((q+1)/2)q^r$ and `group_dimensions` specifies the actual finite
abelian group. For $q=121$ each additive coordinate is represented over
$\mathbb F_{11}$; see {ref}`families` for the finite-field arithmetic.

```{eval-rst}
.. automodule:: experiments.scalar_lift.williamson_spence_family
   :members: field_subtract, character, starter, SpenceSource, cyclotomic_factors, finite_member
   :member-order: bysource
   :undoc-members:
```

`entries(i,j)` returns all four source entries. `entry(j,i,k)` selects one
source matrix. `row` and `coefficient_rows` materialize what their names
request and can become large; entry access is the appropriate interface for
large parameters.

```{eval-rst}
.. automodule:: experiments.scalar_lift.additional_spence_families
   :members: starter, paley, AdditionalSpenceSource, source, endpoint, finite_member
   :member-order: bysource
```

`source(q,r)` dispatches between the original five parameters and the six
additional prime parameters. `finite_member` checks all group coefficients
of the three source norms with an integer NTT; it does not allocate a dense
Hadamard endpoint. Its certificate distinguishes these complete finite
source checks from the all-parameter recurrence proof.

(reference-parity)=
## Coordinate charts and determinant maps

The symmetric chart records multiplication orbits and their coordinate-cycle
phases. The unrestricted chart records actual sites. The distinction controls
both sign variables and which shift representatives suffice during solving.
Final verifiers reconstruct all residues in either case.

```{eval-rst}
.. automodule:: experiments.scalar_lift.large_parity_map
   :members: BinaryField, recover_map, validate_map, validate_expanded_geometry, DeterminantSampler, normalize_head
   :member-order: bysource
   :undoc-members:
```

`inverse_columns` stores an invertible binary coordinate map; the field's
polynomial modulus is encoded as an integer bitmask. `validate_map` checks
irreducibility, both ranks and every constant/linear/quadratic coefficient.
`validate_expanded_geometry` independently checks those coefficients against
actual expanded shifts. Sampling is admitted only after these checks.

```{eval-rst}
.. automodule:: experiments.scalar_lift.exact_parity_sampler
   :members: NativeDeterminantSampler
   :member-order: bysource
   :undoc-members:
```

The native sampler uses fields of degree at most 63 and four machine words
for the coordinate labels. Its tagged 256-bit generator state is distinct
from Python `Random` state. A saved draw count and state are both required to
resume the exact stream. The unrestricted83/89-degree maps use arbitrary-size
Python field arithmetic:

```{eval-rst}
.. automodule:: experiments.scalar_lift.unrestricted_parity_map
   :members: periodic_chart, recover_map, validate_map, validate_expanded_geometry, DeterminantSampler
   :member-order: bysource
   :undoc-members:
```

The identity-action chart instead keeps every coordinate fixed under the
multiplier. Its sign variables index whole orbits, including the singleton
at residue zero. `prepare` converts the four requested row sums into four
exact orbit-sign counts, accounting for that singleton separately. Its
`groups` therefore contain orbit indices; the unrestricted model below
uses actual residue indices. `validate_map` checks both the complete
quadratic identity and its expansion at every nonzero shift.

```{eval-rst}
.. automodule:: experiments.scalar_lift.identity_multiplier
   :members: chart, target_rows, normalize_head, validate_map, prepare, solve_bits, verify_states
   :member-order: bysource
```

The bounded target scanner checks row-sum norms and signed orbit counts
before constructing a parity map. It binds its frontier, completed-order
ledger and source-width tables by hashes. Its output contains necessary
arithmetic candidates, with no signed design or endpoint assertion.

```{eval-rst}
.. automodule:: experiments.scalar_lift.low_gap_target_scan
   :members: identity_dc_rows, c3_dc_rows, scan_targets, build_report
   :member-order: bysource
```

(reference-sign-model)=
## Prepared sign systems

All sign variables are negative-sign bits: $x_i=0$ means $+1$ and $x_i=1$
means $-1$. The unrestricted model's dictionary retains the following exact
data:

| Field | Meaning |
|---|---|
| `axes` | Coordinate label 0–3 at each actual residue |
| `data` | Complete chart and shift coefficients |
| `weights[ell,i,j]` | Symmetric nonnegative pair weights; sum only over $i<j$ |
| `diagonal[ell]` | Sign-independent contribution from loops |
| `particular` | Packed particular solution of the modulo-four linear equations |
| `basis[k]` | Packed original-variable change for free bit $y_k$ |
| `forms[i]` | Mask giving $x_i=x_i^{(0)}\mathbin{\oplus}\operatorname{parity}(\mathrm{forms}[i]\mathbin{\&}y)$ |
| `groups` | Four lists of actual sites, one per coordinate |
| `negative_counts` | Four exact negative counts determined by supports and target sums |
| `target_row` | Four integer sums, with squared norm $p$ |
| `forced_signs` | Explicit gauges and signs forced by exact counts |

```{eval-rst}
.. automodule:: experiments.scalar_lift.unrestricted_t_sign_model
   :members: prepare, exact_count_mask, states_from_assignment, verify_states
   :member-order: bysource
   :undoc-members:
```

`InfeasibleAxes` is a specific arithmetic exclusion. Other `ValueError`
exceptions identify malformed inputs; a campaign must not silently interpret
those as nonexistence proofs. `verify_states` reconstructs every integer
entry independently of the solver's cut variables or cached Gram weights.

(reference-solvers)=
## Exact sign backends

The original symmetric preparation has two count groups, dictated by its
coordinate-cycle symmetry. Its full affine masks use arbitrary-size Python
integers; complete native enumeration is limited to 60 free bits.

```{eval-rst}
.. automodule:: experiments.scalar_lift.quadratic_sign_lift
   :members: prepare, quadratic_mod8, sign_mask, correlations, verify_states, iter_head_sign_branches, solve_bits
   :member-order: bysource
   :undoc-members:
```

The two unrestricted backends use the four-count model. Both enforce full
integer correlations by default. Their positive outputs contain literal
`states` and `expanded_verification`; their completeness flags concern only
the given pattern and constraints. SAT may optionally add exact correlations
lazily, with every rejected candidate checked and sound exact cuts added.
CP-SAT reports its own solver status strings. Interpret statuses together
with `solutions` and the explicit completeness flags, as in
{ref}`sign-solvers-status`.

```{eval-rst}
.. automodule:: experiments.scalar_lift.unrestricted_t_sign_sat
   :members: solve, solve_axes, solve_rows
   :member-order: bysource
   :undoc-members:
```

```{eval-rst}
.. automodule:: experiments.scalar_lift.unrestricted_t_sign_cp
   :members: build_model, solve
   :member-order: bysource
   :undoc-members:
```

## Durable workers

These functions coordinate the preceding exact operations. Their SQLite
ledgers distinguish completed exclusions, unfinished work and positive
witnesses. `--detach` captures source files before starting an independent
process; a saved launch command resumes those sources. Full command examples
and recovery semantics are in {ref}`campaigns`.

```{eval-rst}
.. automodule:: experiments.scalar_lift.exact_t_campaign
   :members: run, verified_endpoints
   :member-order: bysource
   :undoc-members:
```

```{eval-rst}
.. automodule:: experiments.scalar_lift.unrestricted_t_campaign
   :members: row_types, run
   :member-order: bysource
   :undoc-members:
```

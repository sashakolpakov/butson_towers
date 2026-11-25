# Manuscript and Lean correspondence audit

The detailed manuscript and the Lean source were reviewed in the order of
its mathematical proof. The review checked the actual endpoint dependencies,
not only the presence of similarly named lemmas. No incorrect construction,
unproved endpoint assumption or finite-input mismatch was found.

The paper states eighteen existence results and six all-parameter families.
Those are the formal claims. Literature novelty is recorded separately.
The trace interpretation and the alternative sixth-root representatives have
the precise formal scope described below.

## Ordered proof correspondence

| Manuscript step | Lean definitions and proved statements | Audit result |
|---|---|---|
| Two signed-coordinate functions and all their correlations | `Data/TreeLiterals.lean`, `Data/FiniteChecks.lean`: `dotCode_eq_sum`, `unitCode_onehot`, `code_corr_sound`; `Data/T127.lean`; `Data/T631Blocked.lean` and its four chunk modules | Forward permutation is e1→e2→e3→e1, with e0 fixed. All expanded literals agree with the printed representative tables. Lean checks all 127 and 631 shifts, independently of the paper's7 and15 shift-class reductions. |
| Circulants, convolution, transpose | `GroupAlgebra.lean`: `regular_mul`, `regular_mul_comm`, `regular_involution`, `regular_gram_apply` | Entries are exactly f(j−i). The opposite sign of an autocorrelation argument is harmless, since the sum is invariant under negating the shift. |
| GS Gram calculation | `GSAlgebra.blocks_gram`; `GoethalsSeidel.blocks_gram`, `gs_gram` | All block signs, transpose positions and required commutation hypotheses match the printed six cancellations. |
| Quaternion coefficient construction | `OrthogonalDesignCore.combine`, `combine_dot`, `combine_correlations`; `OrthogonalDesign.combined_gram`, `coefficient_gram`, `coefficient_cross` | The four components are exactly the four rows of the paper's M(x)b. Norm and mixed identities use the same specialization and polarization. |
| Final signs | `combine_basis_apply`, `combine_basis_onehot`, `coefficients_signed_partition`; `tensor_substitute_flat` | Uniqueness of the occupied coefficient is proved separately from orthogonality. |
| Amicable substitution and order | `OrthogonalSubstitution.lean`: `tensor_substitute_gram`, `tensor_substitute_hadamard`; `Construction.exists_from_t` | Requires signs, pairwise A_i A_jᵀ=A_j A_iᵀ and the summed Gram. No extra symmetry/commutation assumption is hidden. The scale/cardinality equality is explicit. |
| Group-algebra scalar interpretation | `FiniteGroupAlgebra.lean`: `regularAlgHom`, `regularAlgHom_injective`, `regular_star`, `regularAlgHom_trace`, `regularAlgHom_normalized_trace`, `group_algebra_identity_trace_metric` | Faithful unital rational action; inverse-group basis; inversion becomes transpose; normalized trace gives metric I. No field cancellation or cyclotomic decomposition is assumed. |
| Trace metric and row Gram | `TraceTowers.trace_metric_adjoint`, `inverse_trace_metric_adjoint`, `weighted_trace_gram_transport`; `ScalarFrame.weighted_gram_block` | The row Gram uses the inverse metric on the correct side. Invertibility is supplied explicitly by the two inverse identities. |
| Successive restrictions | `finite_free_trace_transitivity`, `tower_action_coordinates`, `tower_trace_metric`, `basis_gram_congruence` | Actual product-basis coordinates and traces are proved. Pair-index ordering is stated explicitly; swapping the pair coordinates is an indexing permutation. |
| Rational sixth-root pair algebra | `SixthRootPair.lean`: `transpose_mul_of_mul_transpose`, `action_mul_transpose`, `action_add_transpose`, `action_transpose`, `action_polynomial`, `action_mul_left`, `imaginaryAction_transpose`, `imaginaryAction_square`, `coordinateAction_mul`, `coordinateAction_conjugate` | Derives the ten matrix identities from precisely the paired-matrix hypotheses and nonzero scale. |
| Final scalar entry test | `ScalarFrame.frame_gram_block`, `framed_matrix_flat`, `framed_matrix_gram`, `framed_matrix_hadamard` | The scalar representation's adjoint property, frame norm, actual sign-entry test and final scale equality are separate hypotheses. |
| Thirteen finite Williamson inputs | `Data/RowChecks.lean`; `Data/Williamson/Chunk0.lean` through `Chunk2.lean`; `Endpoints.williamson*_quad` | All thirteen printed quadruples agree entry for entry with Lean. Sign, symmetry, every cyclic correlation and amicability are proved from the literals. The larger archive retains 27 inputs. |
| Cubic seed and recurrence | `Data/Unreal13.lean`; `UnrealTripling.nextA`, `nextB`, `nextA_sign`, `nextB_sign`, `next_amicable`, `next_norm`, `tower`, `quad_amicable` | All338 seed signs and both169-entry pair identities match. Iteration proves signs, norm, amicability and the index cardinality for every r. |
| Five strong Spence starters | `Data/StrongChecks.quad_matrices_strong`; `Data/StrongStarters.lean` | The literals give both M0=4I and N0=4qI. A weaker Williamson norm is never substituted for these two hypotheses. |
| Five character matrices | `Data/Paley.lean`, `PaleyChecks.lean`, `PaleyPairChecks.lean` | All printed rows agree, including all 121 entries at q=121. The latter group is `ZMod 11 × ZMod 11`, encoded a+11b, with componentwise addition. |
| Spence recurrence | `SpenceMatrix.stepPlus`, `stepMinus`, `pairStep_sign`, `pairStep_symmetric`, `pairStep_commute`, `step_minus_norm`, `step_plus_norm`, `strong_step`, `family_strong`, `family_norm`, `family_card`; `SpenceHadamard.family_amicable` | The same two strong norms are carried through induction. The condition 2w0=q+1 is explicit. Signs, symmetry, commutativity and growing order are proved. |
| Eighteen ranks and six families | `Endpoints.lean`: `hadamard_*`, `spence_family_{13,17,61,101,121}`, `unreal_cubic_family`, `announced_orders_exist` | Thirteen finite endpoints use Williamson substitution; four use Spence at r=1;295308 uses the cubic family at r=2. All depend on the proved construction chain. |

## Exact conventions and proof scope

The paper and saved Python recipes use reflection u+v=−1. The core Lean
GS construction uses u+v=0. `ReflectionConvention.lean` proves

```
gsAt(-1,a)((r,u),(c,v)) = gs(a)((r,u+q),(c,v+q)),  2q=1.
```

The audited concrete translations are q=316 modulo 631 and q=64 modulo 127:
`python_reflection_631` and `python_reflection_127`. Thus the conventions
are related by simultaneous row and column permutations. The formal
existence claims do not assert equality with the saved packed matrix.

The paper proves character identities by counting solutions over a finite
field. Lean proves exactly the required identities of each of the five
finite integer matrices from its literal entries, using kernel reduction.
It does not need, or claim, a general formal quadratic-character theorem.
The orbit encoding in the paper is likewise expanded into literal inputs;
Lean verifies every shift rather than assuming the encoding's covariance.

The formal rank7572 endpoint uses the printed order3 Williamson quadruple.
The cubic endpoint uses the printed real quadruple (A_r,B_r,B_r,B_r).
The paper proves these same routes, then separately explains their scalar
interpretations. The group-algebra interpretation is formally justified,
but endpoint proofs use its equivalent concrete tensor formulas directly.
They do not require an unproved decomposition of a group algebra into fields.

The paper's particular sixth-root field embedding, faithful module-dimension
calculation and explicit Hermitian form are mathematical proofs in the text.
The original formal development contains the general scalar-frame theorem;
it does not construct that particular field or Hermitian form. The added
`SixthRootPair.lean` proves the rational matrix identities underlying that
illustration. The paragraph transferring an arbitrary incoming Hermitian matrix is proved
in the manuscript; the existing formal trace module proves scalar trace
metrics, their coordinate changes and tower identities, rather than that
full matrix-valued statement as a named theorem. This generality is not a
premise of any endpoint proof.

This scope distinction does not weaken any of the eighteen
formal existence statements or six all-parameter theorems.

## Finite input checks

An independent review parsed the actual TeX tables and compared them with
the Lean lookup literals: both T functions, all thirteen Williamson
quadruples, all five strong starters, every character row and the order-13
pair agree, including ordering and signs. Every printed correlation and
both strong starter norms were recomputed with integer arithmetic.

Reproduce the executable provenance and finite arithmetic checks with

```sh
python3 scripts/check_lean_literal_data.py
python3 scripts/verify_manuscript_data.py
```

These checks are additional source comparisons. Lean never assumes their
outputs or file hashes as mathematical premises. Its finite proofs use
`decide +kernel`, and not native evaluation or new axioms.

A fresh extraction of the original package rebuilt successfully and checked
all 85 selected reports in 180.4 seconds. A separate fresh extraction of the
supplement built and checked its ten reports in 6.67 seconds. Only the standard
axioms propext, Classical.choice and Quot.sound occur. The exact JSON records
are retained as `release/butson_towers_v3/lean-rerun.json` and
`release/lean_sixth_root_pair/verification.json`.
The original package's narrative scope string mentions five Spence families;
its exact audited theorem list also contains `unreal_cubic_family`. The
actual certified total is six families, as recorded in the current scope
document, without rewriting the frozen older archive.
 The current `./scripts/check_formal.sh` also passed both exact report lists
(85 + 10 statements). The repository suite passed 326 tests and 4 subtests.

The original proof archive, original matrix
bundles and earlier article releases remain unchanged.

import ButsonTowers.GroupAlgebra
import ButsonTowers.ScalarFrame
import Mathlib.Algebra.MonoidAlgebra.Basic
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.NormNum

namespace ButsonTowers.GroupAlgebra

open scoped BigOperators
variable {G : Type*} [AddCommGroup G] [Fintype G] [DecidableEq G]

theorem regular_injective : Function.Injective (regular : (G → ℚ) → Matrix G G ℚ) := by
  intro f g h
  funext x
  simpa [regular] using congrArg (fun M : Matrix G G ℚ => M 0 x) h

theorem algebra_coeff_mul (f g : AddMonoidAlgebra ℚ G) :
    (fun x => (f*g).coeff x) = convolution f.coeff g.coeff := by
  funext x
  rw [AddMonoidAlgebra.coeff_mul_apply_left,
    Finsupp.sum_fintype _ _ (by intro; simp)]
  simp [convolution, sub_eq_add_neg, add_comm]

@[simp] theorem algebra_one_coeff (x : G) :
    (1 : AddMonoidAlgebra ℚ G).coeff x = if x=0 then 1 else 0 := by
  rw [AddMonoidAlgebra.one_def, AddMonoidAlgebra.coeff_single]
  simp [Finsupp.single_apply, eq_comm]

/-- The faithful rational algebra homomorphism for the actual finite group algebra. -/
noncomputable def regularAlgHom : AddMonoidAlgebra ℚ G →ₐ[ℚ] Matrix G G ℚ where
  toFun f := regular f.coeff
  map_zero' := by ext i j; simp [regular]
  map_one' := by
    ext i j
    simp [regular, algebra_one_coeff,
      Matrix.one_apply, sub_eq_zero, eq_comm]
  map_add' f g := by ext i j; simp [regular]
  map_mul' f g := by
    change regular (fun x => (f*g).coeff x) = _
    rw [algebra_coeff_mul, regular_mul]
  commutes' q := by
    ext i j
    simp [regular, algebra_one_coeff,
      Algebra.algebraMap_eq_smul_one, Matrix.smul_apply, Matrix.one_apply,
      sub_eq_zero, eq_comm]

theorem regularAlgHom_injective :
    Function.Injective (regularAlgHom : AddMonoidAlgebra ℚ G →ₐ[ℚ] Matrix G G ℚ) := by
  intro f g h
  have hc := regular_injective h
  exact AddMonoidAlgebra.ext (Finsupp.ext (fun x => congrFun hc x))

noncomputable def ofFunction (f : G → ℚ) : AddMonoidAlgebra ℚ G :=
  .ofCoeff (Finsupp.equivFunOnFinite.symm f)

@[simp] theorem ofFunction_coeff (f : G → ℚ) (x : G) : (ofFunction f).coeff x = f x := by
  simp [ofFunction]

@[simp] theorem regularAlgHom_ofFunction (f : G → ℚ) : regularAlgHom (ofFunction f) = regular f := by
  ext i j
  simp [regularAlgHom, regular]

noncomputable def algebraInvolution (f : AddMonoidAlgebra ℚ G) : AddMonoidAlgebra ℚ G :=
  .ofCoeff (Finsupp.equivFunOnFinite.symm (involution f.coeff))

@[simp] theorem algebraInvolution_coeff (f : AddMonoidAlgebra ℚ G) (x : G) :
    (algebraInvolution f).coeff x = f.coeff (-x) := by
  simp [algebraInvolution, involution]

theorem regular_algebraInvolution (f : AddMonoidAlgebra ℚ G) :
    regularAlgHom (algebraInvolution f) = (regularAlgHom f).transpose := by
  ext i j
  simp [regularAlgHom, regular, algebraInvolution_coeff]

noncomputable instance groupAlgebraStarRing : StarRing (AddMonoidAlgebra ℚ G) where
  star := algebraInvolution
  star_involutive f := by
    apply regularAlgHom_injective
    change regularAlgHom (algebraInvolution (algebraInvolution f)) = regularAlgHom f
    rw [regular_algebraInvolution, regular_algebraInvolution, Matrix.transpose_transpose]
  star_mul f g := by
    apply regularAlgHom_injective
    change regularAlgHom (algebraInvolution (f*g)) =
      regularAlgHom (algebraInvolution g * algebraInvolution f)
    rw [regular_algebraInvolution, map_mul, Matrix.transpose_mul, map_mul,
      regular_algebraInvolution, regular_algebraInvolution]
  star_add f g := by
    apply regularAlgHom_injective
    change regularAlgHom (algebraInvolution (f+g)) =
      regularAlgHom (algebraInvolution f + algebraInvolution g)
    rw [regular_algebraInvolution, map_add, Matrix.transpose_add, map_add,
      regular_algebraInvolution, regular_algebraInvolution]

theorem regular_star (f : AddMonoidAlgebra ℚ G) :
    regularAlgHom (star f) = (regularAlgHom f).transpose := regular_algebraInvolution f

def coefficientTrace (f : AddMonoidAlgebra ℚ G) : ℚ := f.coeff 0

theorem regularAlgHom_trace (f : AddMonoidAlgebra ℚ G) :
    Matrix.trace (regularAlgHom f) = (Fintype.card G : ℚ) * coefficientTrace f :=
  regular_trace f.coeff

theorem regularAlgHom_normalized_trace (f : AddMonoidAlgebra ℚ G) :
    Matrix.trace (regularAlgHom f) / (Fintype.card G : ℚ) = coefficientTrace f := by
  rw [regularAlgHom_trace]
  have hn : (Fintype.card G : ℚ) ≠ 0 := by exact_mod_cast Fintype.card_ne_zero
  field_simp

noncomputable def basisElement (j : G) : AddMonoidAlgebra ℚ G := ofFunction (inverseBasis j)

theorem group_algebra_identity_trace_metric (i j : G) :
    coefficientTrace (star (basisElement i) * basisElement j) = if i=j then 1 else 0 := by
  unfold coefficientTrace
  rw [congrFun (algebra_coeff_mul (star (basisElement i)) (basisElement j)) 0]
  change convolution (fun x => (algebraInvolution (basisElement i)).coeff x)
    (fun x => (basisElement j).coeff x) 0 = _
  simp only [algebraInvolution_coeff, basisElement, ofFunction_coeff]
  change convolution (involution (inverseBasis i)) (inverseBasis j) 0 = _
  exact inverse_basis_trace_metric (G := G) i j

end ButsonTowers.GroupAlgebra

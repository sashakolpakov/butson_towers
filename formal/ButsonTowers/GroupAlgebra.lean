import Mathlib.LinearAlgebra.Matrix.Circulant
import Mathlib.LinearAlgebra.Matrix.Trace
import Mathlib.LinearAlgebra.Matrix.Notation
import Mathlib.Data.ZMod.Basic
import Mathlib.Tactic.Abel
import Mathlib.Tactic.FinCases
import Mathlib.Tactic.Ring
import Mathlib.Tactic.NoncommRing
import Mathlib.LinearAlgebra.Matrix.ConjTranspose

/-! The finite additive-group regular action in the inverse group basis.
The coefficient functions use convolution, not pointwise multiplication. -/

namespace ButsonTowers.GroupAlgebra

open scoped BigOperators

variable {G : Type*} [AddCommGroup G] [Fintype G] [DecidableEq G]

def regular (f : G → ℚ) : Matrix G G ℚ := fun i j => f (j-i)

def convolution (f g : G → ℚ) : G → ℚ := fun s => ∑ x, f x * g (s-x)

def involution (f : G → ℚ) : G → ℚ := fun x => f (-x)

@[simp] theorem involution_involution (f : G → ℚ) : involution (involution f) = f := by
  ext x
  simp [involution]

def coefficientOne : G → ℚ := fun x => if x = 0 then 1 else 0

@[simp] theorem regular_apply (f : G → ℚ) (i j : G) :
    regular f i j = f (j-i) := rfl

theorem regular_eq_circulant (f : G → ℚ) :
    regular f = Matrix.circulant (fun x => f (-x)) := by
  ext i j
  simp [regular]

@[simp] theorem regular_zero : regular (0 : G → ℚ) = 0 := by ext; rfl

@[simp] theorem regular_add (f g : G → ℚ) : regular (f+g) = regular f + regular g := by
  ext; rfl

@[simp] theorem regular_neg (f : G → ℚ) : regular (-f) = -regular f := by ext; rfl

@[simp] theorem regular_smul (c : ℚ) (f : G → ℚ) :
    regular (c • f) = c • regular f := by ext; rfl

theorem regular_mul (f g : G → ℚ) :
    regular (convolution f g) = regular f * regular g := by
  ext i j
  simp only [regular, convolution, Matrix.mul_apply]
  refine Fintype.sum_equiv (Equiv.addRight i) _ _ ?_
  intro x
  simp only [Equiv.coe_addRight, add_sub_cancel_right]
  congr 2
  abel

theorem regular_mul_comm (f g : G → ℚ) :
    regular f * regular g = regular g * regular f := by
  simpa only [regular_eq_circulant] using
    Matrix.circulant_mul_comm (fun x => f (-x)) (fun x => g (-x))

@[simp] theorem regular_involution (f : G → ℚ) :
    regular (involution f) = (regular f).transpose := by
  ext i j
  simp [regular, involution]

@[simp] theorem regular_one : regular (coefficientOne : G → ℚ) = 1 := by
  ext i j
  simp [regular, coefficientOne, Matrix.one_apply, sub_eq_zero, eq_comm]

theorem regular_trace (f : G → ℚ) :
    Matrix.trace (regular f) = (Fintype.card G : ℚ) * f 0 := by
  simp [Matrix.trace, regular]

/-- The inverse group basis `β_j = c^(-j)` has the identity trace metric. -/
def inverseBasis (j : G) : G → ℚ := fun x => if x = -j then 1 else 0

theorem inverse_basis_trace_metric (i j : G) :
    convolution (involution (inverseBasis i)) (inverseBasis j) 0 =
      if i = j then 1 else 0 := by
  classical
  simp [convolution, involution, inverseBasis, neg_eq_iff_eq_neg,
    neg_inj, eq_comm]

theorem regular_gram_apply (f g : G → ℚ) (i j : G) :
    (regular f * (regular g).transpose) i j =
      ∑ x, f x * g (x+(i-j)) := by
  simp only [Matrix.mul_apply, Matrix.transpose_apply, regular_apply]
  refine Fintype.sum_equiv (Equiv.subRight i) _ _ ?_
  intro x
  simp only [Equiv.subRight_apply]
  congr 2
  abel

theorem sum_regular_gram {ι : Type*} [Fintype ι]
    (b : ι → G → ℚ) (w : ℚ)
    (h : ∀ s, ∑ c, ∑ x, b c x * b c (x+s) = if s = 0 then w else 0) :
    (∑ c, regular (b c) * (regular (b c)).transpose) =
      w • (1 : Matrix G G ℚ) := by
  ext i j
  simp only [Matrix.sum_apply, regular_gram_apply]
  rw [h (i-j)]
  by_cases hij : i = j <;> simp [Matrix.smul_apply, Matrix.one_apply, hij, sub_eq_zero]

end ButsonTowers.GroupAlgebra

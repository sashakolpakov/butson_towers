import ButsonTowers.GroupAlgebra

namespace ButsonTowers.GoethalsSeidel

open scoped BigOperators
open GroupAlgebra

variable {G : Type*} [AddCommGroup G] [Fintype G] [DecidableEq G]

def reflection : Matrix G G ℚ := fun i j => if j = -i then 1 else 0

theorem eq_neg_swap (i j : G) : (i = -j) ↔ (j = -i) := by
  constructor <;> intro h <;> simpa using congrArg Neg.neg h.symm

@[simp] theorem reflection_transpose : (reflection : Matrix G G ℚ).transpose = reflection := by
  ext i j
  simp only [reflection, Matrix.transpose_apply, eq_neg_swap i j]

@[simp] theorem reflection_sq : (reflection : Matrix G G ℚ) * reflection = 1 := by
  ext i j
  simp only [reflection, Matrix.mul_apply, ite_mul, one_mul, zero_mul]
  simp [Matrix.one_apply, eq_comm]

theorem reflection_regular (f : G → ℚ) :
    (reflection : Matrix G G ℚ) * regular f = regular (involution f) * reflection := by
  ext i j
  simp only [reflection, regular, involution, Matrix.mul_apply, ite_mul,
    mul_ite, one_mul, mul_one, zero_mul, mul_zero]
  simp_rw [eq_neg_swap j]
  simp [sub_eq_add_neg, add_comm]

theorem reflection_sandwich (f : G → ℚ) :
    (reflection : Matrix G G ℚ) * (regular f * reflection) = regular (involution f) := by
  rw [← mul_assoc, reflection_regular, mul_assoc, reflection_sq, mul_one]

theorem reflection_square_right (C : Matrix G G ℚ) :
    (reflection : Matrix G G ℚ) * (reflection * C) = C := by
  rw [← mul_assoc, reflection_sq, one_mul]

theorem regular_mul_right_comm (f g : G → ℚ) (C : Matrix G G ℚ) :
    regular f * (regular g * C) = regular g * (regular f * C) := by
  rw [← mul_assoc, regular_mul_comm f g, mul_assoc]

def blockAdjoint (B : Matrix (Fin 4) (Fin 4) (Matrix G G ℚ)) :
    Matrix (Fin 4) (Fin 4) (Matrix G G ℚ) := fun i j => (B j i).transpose

def flatten (B : Matrix (Fin 4) (Fin 4) (Matrix G G ℚ)) :
    Matrix (Fin 4 × G) (Fin 4 × G) ℚ := fun i j => B i.1 j.1 i.2 j.2

theorem flatten_mul (B C : Matrix (Fin 4) (Fin 4) (Matrix G G ℚ)) :
    flatten (B*C) = flatten B * flatten C := by
  ext ⟨i,a⟩ ⟨j,b⟩
  simp [flatten, Matrix.mul_apply, Fintype.sum_prod_type, Matrix.sum_apply]

@[simp] theorem flatten_adjoint (B : Matrix (Fin 4) (Fin 4) (Matrix G G ℚ)) :
    flatten (blockAdjoint B) = (flatten B).transpose := rfl

def blocks (a : Fin 4 → G → ℚ) : Matrix (Fin 4) (Fin 4) (Matrix G G ℚ) :=
  let A := regular (a 0)
  let B := regular (a 1)
  let C := regular (a 2)
  let D := regular (a 3)
  let J : Matrix G G ℚ := reflection
  !![A, B*J, C*J, D*J;
     -(B*J), A, -(D.transpose*J), C.transpose*J;
     -(C*J), D.transpose*J, A, -(B.transpose*J);
     -(D*J), -(C.transpose*J), B.transpose*J, A]

def gs (a : Fin 4 → G → ℚ) : Matrix (Fin 4 × G) (Fin 4 × G) ℚ := flatten (blocks a)


end ButsonTowers.GoethalsSeidel

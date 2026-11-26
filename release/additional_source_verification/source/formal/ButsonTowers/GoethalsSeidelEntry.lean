import ButsonTowers.GoethalsSeidelCore

namespace ButsonTowers.GoethalsSeidel
open scoped BigOperators
open GroupAlgebra
variable {G : Type*} [AddCommGroup G] [Fintype G] [DecidableEq G]

def channel : Fin 4 → Fin 4 → Fin 4 :=
  !![0,1,2,3; 1,0,3,2; 2,3,0,1; 3,2,1,0]

def entrySign (r c : Fin 4) : ℚ :=
  if (r=1 ∧ (c=0 ∨ c=2)) ∨ (r=2 ∧ (c=0 ∨ c=3)) ∨
      (r=3 ∧ (c=0 ∨ c=1)) then -1 else 1

def entryArgument (r c : Fin 4) (u v : G) : G :=
  if r=c then v-u else if r=0 ∨ c=0 then -v-u else u+v

theorem channel_involutive (r : Fin 4) : Function.Involutive (channel r) := by
  intro c
  fin_cases r <;> fin_cases c <;> rfl

theorem entrySign_sign (r c : Fin 4) : entrySign r c = 1 ∨ entrySign r c = -1 := by
  unfold entrySign
  split_ifs
  · exact Or.inr rfl
  · exact Or.inl rfl


theorem regular_reflection_apply (f : G → ℚ) (i j : G) :
    (regular f * (reflection : Matrix G G ℚ)) i j = f (-j-i) := by
  simp only [Matrix.mul_apply, reflection, mul_ite, mul_one, mul_zero,
    eq_neg_swap j]
  simp [regular]

theorem transpose_regular_reflection_apply (f : G → ℚ) (i j : G) :
    ((regular f).transpose * (reflection : Matrix G G ℚ)) i j = f (i+j) := by
  rw [← regular_involution, regular_reflection_apply]
  simp [involution, sub_eq_add_neg, add_comm]


theorem gs_apply (a : Fin 4 → G → ℚ) (r c : Fin 4) (u v : G) :
    gs a (r,u) (c,v) = entrySign r c * a (channel r c) (entryArgument r c u v) := by
  fin_cases r <;> fin_cases c
  · change a 0 (v-u) = 1 * a 0 (v-u)
    exact (one_mul _).symm
  · change ((regular (a 1)) * (reflection : Matrix G G ℚ)) u v = 1 * a 1 (-v-u)
    rw [regular_reflection_apply, one_mul]
  · change ((regular (a 2)) * (reflection : Matrix G G ℚ)) u v = 1 * a 2 (-v-u)
    rw [regular_reflection_apply, one_mul]
  · change ((regular (a 3)) * (reflection : Matrix G G ℚ)) u v = 1 * a 3 (-v-u)
    rw [regular_reflection_apply, one_mul]
  · change -(((regular (a 1)) * (reflection : Matrix G G ℚ)) u v) = -1 * a 1 (-v-u)
    rw [regular_reflection_apply, neg_one_mul]
  · change a 0 (v-u) = 1 * a 0 (v-u)
    exact (one_mul _).symm
  · change -(((regular (a 3)).transpose * (reflection : Matrix G G ℚ)) u v) = -1 * a 3 (u+v)
    rw [transpose_regular_reflection_apply, neg_one_mul]
  · change ((regular (a 2)).transpose * (reflection : Matrix G G ℚ)) u v = 1 * a 2 (u+v)
    rw [transpose_regular_reflection_apply, one_mul]
  · change -(((regular (a 2)) * (reflection : Matrix G G ℚ)) u v) = -1 * a 2 (-v-u)
    rw [regular_reflection_apply, neg_one_mul]
  · change ((regular (a 3)).transpose * (reflection : Matrix G G ℚ)) u v = 1 * a 3 (u+v)
    rw [transpose_regular_reflection_apply, one_mul]
  · change a 0 (v-u) = 1 * a 0 (v-u)
    exact (one_mul _).symm
  · change -(((regular (a 1)).transpose * (reflection : Matrix G G ℚ)) u v) = -1 * a 1 (u+v)
    rw [transpose_regular_reflection_apply, neg_one_mul]
  · change -(((regular (a 3)) * (reflection : Matrix G G ℚ)) u v) = -1 * a 3 (-v-u)
    rw [regular_reflection_apply, neg_one_mul]
  · change -(((regular (a 2)).transpose * (reflection : Matrix G G ℚ)) u v) = -1 * a 2 (u+v)
    rw [transpose_regular_reflection_apply, neg_one_mul]
  · change ((regular (a 1)).transpose * (reflection : Matrix G G ℚ)) u v = 1 * a 1 (u+v)
    rw [transpose_regular_reflection_apply, one_mul]
  · change a 0 (v-u) = 1 * a 0 (v-u)
    exact (one_mul _).symm

theorem gs_add (a b : Fin 4 → G → ℚ) : gs (a+b) = gs a + gs b := by
  ext ⟨r,u⟩ ⟨c,v⟩
  simp only [gs_apply, Pi.add_apply, Matrix.add_apply, mul_add]

end ButsonTowers.GoethalsSeidel

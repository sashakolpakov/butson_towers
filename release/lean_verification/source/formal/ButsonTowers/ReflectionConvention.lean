import ButsonTowers.GoethalsSeidelEntry

/-! Exact relation between reflection about zero and the literal Python
convention `J i j = 1` iff `i+j=-1`.  Both row and column indices are
translated by the same group element. -/

namespace ButsonTowers.GoethalsSeidel

open GroupAlgebra
variable {G : Type*} [AddCommGroup G] [Fintype G] [DecidableEq G]

def reflectionAt (t : G) : Matrix G G ℚ := fun i j => if j=t-i then 1 else 0

theorem eq_reflection_swap (t j i : G) : (j=t-i) ↔ (i=t-j) := by
  constructor <;> intro h <;> rw [h] <;> abel

theorem regular_reflectionAt_apply (f : G → ℚ) (t i j : G) :
    (regular f * (reflectionAt t : Matrix G G ℚ)) i j = f (t-j-i) := by
  simp only [Matrix.mul_apply, reflectionAt, mul_ite, mul_one, mul_zero,
    eq_reflection_swap t j]
  simp [regular]

theorem transpose_regular_reflectionAt_apply (f : G → ℚ) (t i j : G) :
    ((regular f).transpose * (reflectionAt t : Matrix G G ℚ)) i j = f (i+j-t) := by
  rw [← regular_involution, regular_reflectionAt_apply]
  change f (-(t-j-i)) = f (i+j-t)
  congr 1
  abel

def blocksAt (t : G) (a : Fin 4 → G → ℚ) :
    Matrix (Fin 4) (Fin 4) (Matrix G G ℚ) :=
  let A := regular (a 0)
  let B := regular (a 1)
  let C := regular (a 2)
  let D := regular (a 3)
  let J := reflectionAt t
  !![A, B*J, C*J, D*J;
     -(B*J), A, -(D.transpose*J), C.transpose*J;
     -(C*J), D.transpose*J, A, -(B.transpose*J);
     -(D*J), -(C.transpose*J), B.transpose*J, A]

def gsAt (t : G) (a : Fin 4 → G → ℚ) :
    Matrix (Fin 4 × G) (Fin 4 × G) ℚ := flatten (blocksAt t a)

def entryArgumentAt (t : G) (r c : Fin 4) (u v : G) : G :=
  if r=c then v-u else if r=0 ∨ c=0 then t-v-u else u+v-t

theorem gsAt_apply (t : G) (a : Fin 4 → G → ℚ) (r c : Fin 4) (u v : G) :
    gsAt t a (r,u) (c,v) = entrySign r c * a (channel r c) (entryArgumentAt t r c u v) := by
  fin_cases r <;> fin_cases c
  · change a 0 (v-u) = 1 * a 0 (v-u)
    exact (one_mul _).symm
  · change ((regular (a 1)) * (reflectionAt t : Matrix G G ℚ)) u v = 1 * a 1 (t-v-u)
    rw [regular_reflectionAt_apply, one_mul]
  · change ((regular (a 2)) * (reflectionAt t : Matrix G G ℚ)) u v = 1 * a 2 (t-v-u)
    rw [regular_reflectionAt_apply, one_mul]
  · change ((regular (a 3)) * (reflectionAt t : Matrix G G ℚ)) u v = 1 * a 3 (t-v-u)
    rw [regular_reflectionAt_apply, one_mul]
  · change -(((regular (a 1)) * (reflectionAt t : Matrix G G ℚ)) u v) = -1 * a 1 (t-v-u)
    rw [regular_reflectionAt_apply, neg_one_mul]
  · change a 0 (v-u) = 1 * a 0 (v-u)
    exact (one_mul _).symm
  · change -(((regular (a 3)).transpose * (reflectionAt t : Matrix G G ℚ)) u v) = -1 * a 3 (u+v-t)
    rw [transpose_regular_reflectionAt_apply, neg_one_mul]
  · change ((regular (a 2)).transpose * (reflectionAt t : Matrix G G ℚ)) u v = 1 * a 2 (u+v-t)
    rw [transpose_regular_reflectionAt_apply, one_mul]
  · change -(((regular (a 2)) * (reflectionAt t : Matrix G G ℚ)) u v) = -1 * a 2 (t-v-u)
    rw [regular_reflectionAt_apply, neg_one_mul]
  · change ((regular (a 3)).transpose * (reflectionAt t : Matrix G G ℚ)) u v = 1 * a 3 (u+v-t)
    rw [transpose_regular_reflectionAt_apply, one_mul]
  · change a 0 (v-u) = 1 * a 0 (v-u)
    exact (one_mul _).symm
  · change -(((regular (a 1)).transpose * (reflectionAt t : Matrix G G ℚ)) u v) = -1 * a 1 (u+v-t)
    rw [transpose_regular_reflectionAt_apply, neg_one_mul]
  · change -(((regular (a 3)) * (reflectionAt t : Matrix G G ℚ)) u v) = -1 * a 3 (t-v-u)
    rw [regular_reflectionAt_apply, neg_one_mul]
  · change -(((regular (a 2)).transpose * (reflectionAt t : Matrix G G ℚ)) u v) = -1 * a 2 (u+v-t)
    rw [transpose_regular_reflectionAt_apply, neg_one_mul]
  · change ((regular (a 1)).transpose * (reflectionAt t : Matrix G G ℚ)) u v = 1 * a 1 (u+v-t)
    rw [transpose_regular_reflectionAt_apply, one_mul]
  · change a 0 (v-u) = 1 * a 0 (v-u)
    exact (one_mul _).symm

def translatedGS (q : G) (a : Fin 4 → G → ℚ) :
    Matrix (Fin 4 × G) (Fin 4 × G) ℚ :=
  fun i j => gs a (i.1,i.2+q) (j.1,j.2+q)

theorem translated_gs_reflection (q : G) (a : Fin 4 → G → ℚ) :
    translatedGS q a = gsAt (-(q+q)) a := by
  ext ⟨r,u⟩ ⟨c,v⟩
  simp only [translatedGS, gs_apply, gsAt_apply]
  congr 2
  unfold entryArgument entryArgumentAt
  split_ifs <;> abel

def translationIndex (q : G) : (Fin 4 × G) ≃ (Fin 4 × G) :=
  Equiv.prodCongr (Equiv.refl _) (Equiv.addRight q)

theorem gsAt_eq_reindex (q : G) (a : Fin 4 → G → ℚ) :
    gsAt (-(q+q)) a =
      Matrix.reindex (translationIndex q).symm (translationIndex q).symm (gs a) := by
  rw [← translated_gs_reflection]
  rfl

theorem python_reflection_eq_reindex [One G] (q : G) (hq : q+q=1)
    (a : Fin 4 → G → ℚ) :
    gsAt (-1) a =
      Matrix.reindex (translationIndex q).symm (translationIndex q).symm (gs a) := by
  simpa [hq] using gsAt_eq_reindex q a

theorem python_reflection_127 (a : Fin 4 → ZMod 127 → ℚ) :
    gsAt (-1) a =
      Matrix.reindex (translationIndex (64 : ZMod 127)).symm
        (translationIndex (64 : ZMod 127)).symm (gs a) :=
  python_reflection_eq_reindex 64 (by decide) a

theorem python_reflection_631 (a : Fin 4 → ZMod 631 → ℚ) :
    gsAt (-1) a =
      Matrix.reindex (translationIndex (316 : ZMod 631)).symm
        (translationIndex (316 : ZMod 631)).symm (gs a) :=
  python_reflection_eq_reindex 316 (by decide) a

end ButsonTowers.GoethalsSeidel

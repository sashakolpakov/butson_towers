import ButsonTowers.Data.LiteralTree
import ButsonTowers.ScalarFrame
namespace ButsonTowers.Data
set_option maxRecDepth 1000000
set_option maxHeartbeats 0
def unreal13ATree : LiteralTree ℤ := (.node 84 (.node 42 (.node 21 (.node 10 (.node 5 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1)))))) (.node 5 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1))))))) (.node 10 (.node 5 (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1)))))) (.node 5 (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1)))))))) (.node 21 (.node 10 (.node 5 (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1)))))) (.node 5 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1))))))) (.node 10 (.node 5 (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1)))))) (.node 5 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1))))))))) (.node 42 (.node 21 (.node 10 (.node 5 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1)))))) (.node 5 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1))))))) (.node 10 (.node 5 (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))))) (.node 5 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))))))) (.node 21 (.node 10 (.node 5 (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1)))))) (.node 5 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1))))))) (.node 11 (.node 5 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))))) (.node 5 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1))))))))))
def unreal13BTree : LiteralTree ℤ := (.node 84 (.node 42 (.node 21 (.node 10 (.node 5 (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1)))))) (.node 5 (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1))))))) (.node 10 (.node 5 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1)))))) (.node 5 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1)))))))) (.node 21 (.node 10 (.node 5 (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1)))))) (.node 5 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1))))))) (.node 10 (.node 5 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1)))))) (.node 5 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1))))))))) (.node 42 (.node 21 (.node 10 (.node 5 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))))) (.node 5 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1))))))) (.node 10 (.node 5 (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1)))))) (.node 5 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1)))))))) (.node 21 (.node 10 (.node 5 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1)))))) (.node 5 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1))))))) (.node 11 (.node 5 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1)))))) (.node 5 (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1))))))))))
def unreal13AInt : Matrix (Fin 13) (Fin 13) ℤ := fun i j => treeAt unreal13ATree (13*i.val+j.val)
def unreal13BInt : Matrix (Fin 13) (Fin 13) ℤ := fun i j => treeAt unreal13BTree (13*i.val+j.val)
def unreal13A : Matrix (Fin 13) (Fin 13) ℚ := fun i j => (unreal13AInt i j : ℚ)
def unreal13B : Matrix (Fin 13) (Fin 13) ℚ := fun i j => (unreal13BInt i j : ℚ)
theorem unreal13A_int_sign : ∀ i j, unreal13AInt i j=1 ∨ unreal13AInt i j= -1 := by decide +kernel
theorem unreal13B_int_sign : ∀ i j, unreal13BInt i j=1 ∨ unreal13BInt i j= -1 := by decide +kernel
theorem unreal13_int_amicable : unreal13AInt * unreal13BInt.transpose = unreal13BInt * unreal13AInt.transpose := by decide +kernel
theorem unreal13_int_norm : unreal13AInt * unreal13AInt.transpose + (3 : ℤ) • (unreal13BInt * unreal13BInt.transpose) =
    (52 : ℤ) • (1 : Matrix (Fin 13) (Fin 13) ℤ) := by decide +kernel

theorem unreal13A_sign : IsSignMatrix unreal13A := by
  intro i j
  change (unreal13AInt i j : ℚ)=1 ∨ (unreal13AInt i j : ℚ)= -1
  obtain h | h := unreal13A_int_sign i j
  · left; exact_mod_cast h
  · right; exact_mod_cast h

theorem unreal13B_sign : IsSignMatrix unreal13B := by
  intro i j
  change (unreal13BInt i j : ℚ)=1 ∨ (unreal13BInt i j : ℚ)= -1
  obtain h | h := unreal13B_int_sign i j
  · left; exact_mod_cast h
  · right; exact_mod_cast h

theorem unreal13_amicable : unreal13A * unreal13B.transpose = unreal13B * unreal13A.transpose := by
  ext i j
  have h := congrArg (fun M : Matrix (Fin 13) (Fin 13) ℤ => M i j) unreal13_int_amicable
  simp only [Matrix.mul_apply, Matrix.transpose_apply] at h ⊢
  simp only [unreal13A, unreal13B]
  exact_mod_cast h

theorem unreal13_norm : unreal13A * unreal13A.transpose + (3 : ℚ) • (unreal13B * unreal13B.transpose) =
    (4*(13 : ℚ)) • (1 : Matrix (Fin 13) (Fin 13) ℚ) := by
  ext i j
  have h := congrArg (fun M : Matrix (Fin 13) (Fin 13) ℤ => M i j) unreal13_int_norm
  simp only [Matrix.add_apply, Matrix.smul_apply, smul_eq_mul, Matrix.mul_apply,
    Matrix.transpose_apply, Matrix.one_apply] at h ⊢
  simp only [unreal13A, unreal13B]
  norm_num at h ⊢
  exact_mod_cast h
#print axioms unreal13_int_norm
#print axioms unreal13_amicable
#print axioms unreal13_norm
end ButsonTowers.Data

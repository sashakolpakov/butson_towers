import ButsonTowers.ScalarFrame

namespace ButsonTowers

variable {ι κ : Type*} [Fintype ι] [DecidableEq ι]
  [Fintype κ] [DecidableEq κ]

def HadamardExists (n : ℕ) : Prop :=
  ∃ H : Matrix (Fin n) (Fin n) ℚ, IsHadamard H

theorem hadamard_reindex (H : Matrix ι ι ℚ) (e : κ ≃ ι)
    (h : IsHadamard H) : IsHadamard (H.submatrix e e) := by
  refine ⟨?_,?_,?_⟩
  · rw [Fintype.card_congr e]
    exact h.1
  · intro i j
    exact h.2.1 (e i) (e j)
  · have ht : (H.submatrix e e).transpose = H.transpose.submatrix e e := rfl
    rw [ht,Matrix.submatrix_mul_equiv,h.2.2]
    ext i j
    simp [Matrix.submatrix_apply,Matrix.smul_apply,Matrix.one_apply,
      Fintype.card_congr e]

theorem hadamard_exists_of_card {n : ℕ} (H : Matrix ι ι ℚ)
    (h : IsHadamard H) (hcard : Fintype.card ι = n) : HadamardExists n := by
  let e := (Fintype.equivFinOfCardEq hcard).symm
  exact ⟨H.submatrix e e,hadamard_reindex H e h⟩

end ButsonTowers

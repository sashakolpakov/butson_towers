import ButsonTowers.Data.RowChecks

namespace ButsonTowers.Data
open scoped BigOperators

-- Avoid a Lean 4.32 simplifier-procedure diagnostic; ordinary rewrite lemmas
-- and the kernel suffice for this generic matrix bridge.
set_option simprocs false

def pairRows (a : Fin 4 → ℕ → ℤ) (plus : Bool) : Fin 4 → ℕ → ℤ :=
  if plus then ![fun x => a 0 x + a 1 x, fun x => a 2 x + a 3 x, fun _ => 0, fun _ => 0]
  else ![fun x => a 0 x - a 1 x, fun x => a 2 x - a 3 x, fun _ => 0, fun _ => 0]

theorem pairRows_matrices {p : ℕ} (a : Fin 4 → ℕ → ℤ) (plus : Bool) :
    rowMatrices p (pairRows a plus) =
      if plus then ![rowMatrices p a 0 + rowMatrices p a 1,
        rowMatrices p a 2 + rowMatrices p a 3, 0, 0]
      else ![rowMatrices p a 0 - rowMatrices p a 1,
        rowMatrices p a 2 - rowMatrices p a 3, 0, 0] := by
  funext i
  ext x y
  fin_cases i <;> cases plus <;>
    simp [rowMatrices, pairRows, GroupAlgebra.regular, rationalRows]

theorem quad_matrices_strong {p q : ℕ} [NeZero p] {a : Fin 4 → ℕ → ℤ}
    (hsign : checkQuadSign p a = true)
    (hsymm : checkQuadSymmetric p a = true)
    (hminus : checkQuadCorr p (pairRows a false) 4 = true)
    (hplus : checkQuadCorr p (pairRows a true) (4*(q : ℤ)) = true) :
    SpenceMatrix.Strong q 0 (rowMatrices p a) := by
  refine ⟨quad_matrices_sign hsign, quad_matrices_symmetric hsymm,
    quad_matrices_commute a, ?_, ?_⟩
  · have h := quad_matrices_norm hminus
    rw [pairRows_matrices] at h
    simp only [Bool.false_eq_true, ↓reduceIte] at h
    simpa [Fin.sum_univ_succ, Matrix.transpose_sub, quad_matrices_symmetric hsymm,
      pow_two] using h
  · have h := quad_matrices_norm hplus
    rw [pairRows_matrices] at h
    simp only [↓reduceIte] at h
    simpa [Fin.sum_univ_succ, Matrix.transpose_add, quad_matrices_symmetric hsymm,
      pow_two] using h

#print axioms quad_matrices_strong
end ButsonTowers.Data

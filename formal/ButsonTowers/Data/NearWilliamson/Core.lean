import ButsonTowers.Data.RowChecks
import ButsonTowers.Data.LiteralTree
import ButsonTowers.ReflectionConvention
import ButsonTowers.OrthogonalSubstitution

/-! Near-Williamson rows become an amicable quadruple by reversing the
columns of the first circulant. The reversal is exactly x+y=-1, matching
the literal finite source recipe. Only rows 1, 2 and 3 require symmetry. -/

namespace ButsonTowers.Data
open scoped BigOperators
open GroupAlgebra GoethalsSeidel

variable {p : ℕ} [NeZero p]

def tailRows (a : Fin 4 → ℕ → ℤ) : Fin 4 → ℕ → ℤ :=
  fun i x => if i = 0 then 0 else a i x

def nearMatrices (p : ℕ) [NeZero p] (a : Fin 4 → ℕ → ℤ) :
    Fin 4 → Matrix (ZMod p) (ZMod p) ℚ :=
  fun i => if i = 0 then rowMatrices p a i * reflectionAt (-1 : ZMod p) else rowMatrices p a i

theorem reflectionAt_transpose (t : ZMod p) :
    (reflectionAt t).transpose = reflectionAt t := by
  ext i j
  simp only [reflectionAt, Matrix.transpose_apply, eq_reflection_swap t i j]

theorem reflectionAt_sq (t : ZMod p) :
    reflectionAt t * reflectionAt t = (1 : Matrix (ZMod p) (ZMod p) ℚ) := by
  ext i j
  simp only [reflectionAt, Matrix.mul_apply, ite_mul, one_mul, zero_mul]
  simp [Matrix.one_apply, eq_reflection_swap]

theorem reflectionAt_regular_transpose (t : ZMod p) (f : ZMod p → ℚ) :
    reflectionAt t * (regular f).transpose = regular f * reflectionAt t := by
  ext i j
  simp only [reflectionAt, regular, Matrix.mul_apply, Matrix.transpose_apply,
    ite_mul, mul_ite, one_mul, zero_mul, mul_one, mul_zero]
  simp_rw [eq_reflection_swap t j]
  simp only [Finset.sum_ite_eq', Finset.mem_univ, if_true]
  congr 1
  abel

theorem near_matrices_sign {a : Fin 4 → ℕ → ℤ}
    (h : checkQuadSign p a = true) : ∀ i, IsSignMatrix (nearMatrices p a i) := by
  intro i x y
  by_cases hi : i = 0
  · subst i
    simp only [nearMatrices, if_true, rowMatrices]
    rw [regular_reflectionAt_apply]
    simp only [rationalRows]
    obtain hs | hs := quad_sign_sound h 0 (-1-y-x)
    · left; exact_mod_cast hs
    · right; exact_mod_cast hs
  · simpa only [nearMatrices, hi, if_false] using quad_matrices_sign h i x y

theorem near_tail_symmetric {a : Fin 4 → ℕ → ℤ}
    (h : checkQuadSymmetric p (tailRows a) = true) (i : Fin 4) (hi : i ≠ 0) :
    (rowMatrices p a i).transpose = rowMatrices p a i := by
  have hh := quad_matrices_symmetric h i
  have he : rowMatrices p (tailRows a) i = rowMatrices p a i := by
    ext x y
    simp only [rowMatrices, regular, rationalRows, tailRows, hi, if_false]
  rw [he] at hh
  exact hh

theorem near_matrices_amicable {a : Fin 4 → ℕ → ℤ}
    (h : checkQuadSymmetric p (tailRows a) = true) :
    ∀ i j, nearMatrices p a i * (nearMatrices p a j).transpose =
      nearMatrices p a j * (nearMatrices p a i).transpose := by
  intro i j
  have hr (k : Fin 4) :
      reflectionAt (-1 : ZMod p) * (rowMatrices p a k).transpose =
        rowMatrices p a k * reflectionAt (-1) :=
    reflectionAt_regular_transpose (-1) (rationalRows a k)
  by_cases hi : i = 0 <;> by_cases hj : j = 0
  · subst i; subst j; rfl
  · simp only [nearMatrices, hi, hj, if_true, if_false, Matrix.transpose_mul,
      reflectionAt_transpose]
    rw [mul_assoc, hr j, ← mul_assoc, quad_matrices_commute a 0 j,
      mul_assoc, ← hr 0]
  · symm
    simp only [nearMatrices, hi, hj, if_true, if_false, Matrix.transpose_mul,
      reflectionAt_transpose]
    rw [mul_assoc, hr i, ← mul_assoc, quad_matrices_commute a 0 i,
      mul_assoc, ← hr 0]
  · simp only [nearMatrices, hi, hj, if_false, near_tail_symmetric h i hi,
      near_tail_symmetric h j hj]
    exact quad_matrices_commute a i j

theorem near_matrices_norm {a : Fin 4 → ℕ → ℤ} {w : ℤ}
    (h : checkQuadCorr p a w = true) :
    (∑ i : Fin 4, nearMatrices p a i * (nearMatrices p a i).transpose) =
      (w : ℚ) • (1 : Matrix (ZMod p) (ZMod p) ℚ) := by
  rw [← quad_matrices_norm h]
  apply Finset.sum_congr rfl
  intro i _
  by_cases hi : i = 0
  · simp only [nearMatrices, hi, if_true, Matrix.transpose_mul,
      reflectionAt_transpose]
    rw [mul_assoc, ← mul_assoc (reflectionAt (-1 : ZMod p)), reflectionAt_sq,
      one_mul]
  · simp only [nearMatrices, hi, if_false]

theorem near_matrices_quad {a : Fin 4 → ℕ → ℤ}
    (hs : checkQuadSign p a = true)
    (ht : checkQuadSymmetric p (tailRows a) = true)
    (hn : checkQuadCorr p a (4*p) = true) :
    AmicableQuad (nearMatrices p a) (4*(p : ℚ)) := by
  refine ⟨near_matrices_sign hs, near_matrices_amicable ht, ?_⟩
  simpa only [Int.cast_mul, Int.cast_ofNat, Int.cast_natCast] using near_matrices_norm hn

theorem ordinary_matrices_quad {a : Fin 4 → ℕ → ℤ}
    (hs : checkQuadSign p a = true)
    (ht : checkQuadSymmetric p a = true)
    (hn : checkQuadCorr p a (4*p) = true) :
    AmicableQuad (rowMatrices p a) (4*(p : ℚ)) := by
  refine ⟨quad_matrices_sign hs, ?_, ?_⟩
  · intro i j
    rw [quad_matrices_symmetric ht i, quad_matrices_symmetric ht j]
    exact quad_matrices_commute a i j
  · simpa only [Int.cast_mul, Int.cast_ofNat, Int.cast_natCast] using quad_matrices_norm hn

end ButsonTowers.Data

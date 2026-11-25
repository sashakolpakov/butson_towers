import ButsonTowers.OrthogonalDesignCore
import ButsonTowers.GoethalsSeidel

namespace ButsonTowers.OrthogonalDesign
open scoped BigOperators
open GroupAlgebra GoethalsSeidel
variable {G : Type*} [AddCommGroup G] [Fintype G] [DecidableEq G]

theorem combined_gram (x : Fin 4 → ℚ) (b : Fin 4 → G → ℚ) (p : ℚ)
    (h : ∀ s, ∑ c, ∑ u, b c u * b c (u+s) = if s=0 then p else 0) :
    gs (combine x b) * (gs (combine x b)).transpose =
      (p*(∑ r, (x r)^2)) • (1 : Matrix (Fin 4 × G) (Fin 4 × G) ℚ) := by
  apply gs_gram
  exact sum_regular_gram _ _ (combine_correlations x b p h)

theorem coefficient_gram (b : Fin 4 → G → ℚ) (p : ℚ)
    (h : ∀ s, ∑ c, ∑ u, b c u * b c (u+s) = if s=0 then p else 0)
    (i : Fin 4) :
    coefficients b i * (coefficients b i).transpose =
      p • (1 : Matrix (Fin 4 × G) (Fin 4 × G) ℚ) := by
  simpa [coefficients, basis, Pi.single_apply] using combined_gram (basis i) b p h

theorem coefficient_cross (b : Fin 4 → G → ℚ) (p : ℚ)
    (h : ∀ s, ∑ c, ∑ u, b c u * b c (u+s) = if s=0 then p else 0)
    (i j : Fin 4) (hij : i ≠ j) :
    coefficients b i * (coefficients b j).transpose +
      coefficients b j * (coefficients b i).transpose = 0 := by
  have hnorm : (∑ k, ((basis i + basis j) k)^2) = 2 := by
    have he : ∀ k, ((basis i + basis j) k)^2 =
        (basis i k)^2 + (basis j k)^2 + 2*(basis i k*basis j k) := by
      intro k
      simp only [Pi.add_apply]
      ring
    simp_rw [he, Finset.sum_add_distrib, ← Finset.mul_sum]
    simp [basis, hij, eq_comm]
    norm_num
  have hh := combined_gram (basis i + basis j) b p h
  rw [hnorm, combine_add, gs_add] at hh
  change (coefficients b i + coefficients b j) *
    (coefficients b i + coefficients b j).transpose = _ at hh
  have hexpand : (coefficients b i + coefficients b j) *
      (coefficients b i + coefficients b j).transpose =
      (coefficients b i * (coefficients b i).transpose +
        coefficients b j * (coefficients b j).transpose) +
      (coefficients b i * (coefficients b j).transpose +
        coefficients b j * (coefficients b i).transpose) := by
    rw [Matrix.transpose_add]
    noncomm_ring
  rw [hexpand, coefficient_gram b p h, coefficient_gram b p h] at hh
  have hsum : p • (1 : Matrix (Fin 4 × G) (Fin 4 × G) ℚ) + p • 1 = (p*2) • 1 := by
    rw [← add_smul]
    congr 1
    ring
  rw [hsum] at hh
  exact add_left_cancel (hh.trans (add_zero _).symm)

end ButsonTowers.OrthogonalDesign

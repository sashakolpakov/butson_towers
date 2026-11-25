import ButsonTowers.Existence
import ButsonTowers.OrthogonalDesign
import ButsonTowers.SpenceHadamard

namespace ButsonTowers
open scoped BigOperators

variable {p : ℕ} [NeZero p]
variable {ι κ : Type} [Fintype ι] [DecidableEq ι] [Nonempty ι]
  [Fintype κ] [DecidableEq κ] [Nonempty κ]

theorem exists_from_t
    (b : Fin 4 → ZMod p → ℚ)
    (hcorr : ∀ s, (∑ i, ∑ x, b i x * b i (x+s)) = if s=0 then (p : ℚ) else 0)
    (honehot : ∀ x, ∃ c, (b c x=1 ∨ b c x= -1) ∧ ∀ d, d≠c → b d x=0)
    (A : Fin 4 → Matrix ι ι ℚ) (hA : AmicableQuad A (4*(Fintype.card ι : ℚ))) :
    HadamardExists (4*p*Fintype.card ι) := by
  apply hadamard_exists_of_card (tensorSubstitute (OrthogonalDesign.coefficients b) A)
  · apply tensor_substitute_hadamard _ _ (p : ℚ) (4*(Fintype.card ι : ℚ))
      (OrthogonalDesign.coefficient_gram b p hcorr)
      (OrthogonalDesign.coefficient_cross b p hcorr)
      (OrthogonalDesign.coefficients_signed_partition b honehot) hA
    simp only [Fintype.card_prod, Fintype.card_fin, ZMod.card]
    push_cast
    ring
  · simp [Nat.mul_assoc]

theorem spence_exists_from_t
    (b : Fin 4 → ZMod p → ℚ)
    (hcorr : ∀ s, (∑ i, ∑ x, b i x * b i (x+s)) = if s=0 then (p : ℚ) else 0)
    (honehot : ∀ x, ∃ c, (b c x=1 ∨ b c x= -1) ∧ ∀ d, d≠c → b d x=0)
    (A : Fin 4 → Matrix ι ι ℚ) (Q : Matrix κ κ ℚ)
    (hA : SpenceMatrix.Strong (Fintype.card κ) 0 A) (hQ : SpenceMatrix.IsPaley Q)
    (hcard : 2*Fintype.card ι=Fintype.card κ+1) (r : ℕ) :
    HadamardExists (4*p*Fintype.card ι*(Fintype.card κ)^r) := by
  have h := exists_from_t b hcorr honehot (SpenceMatrix.family A Q r)
    (SpenceMatrix.family_amicable A Q hA hQ hcard r)
  simpa [SpenceMatrix.family_card, Nat.mul_assoc] using h

end ButsonTowers

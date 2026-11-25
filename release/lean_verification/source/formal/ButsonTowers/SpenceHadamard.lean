import ButsonTowers.SpenceMatrix
import ButsonTowers.OrthogonalSubstitution

namespace ButsonTowers.SpenceMatrix

open scoped BigOperators

variable {ν ι κ : Type} [Fintype ν] [DecidableEq ν]
  [Fintype ι] [DecidableEq ι] [Fintype κ] [DecidableEq κ]

theorem family_amicable
    (A : Fin 4 → Matrix ι ι ℚ) (Q : Matrix κ κ ℚ)
    (hA : Strong (Fintype.card κ) 0 A) (hQ : IsPaley Q)
    (hcard : 2 * Fintype.card ι = Fintype.card κ + 1) (r : ℕ) :
    AmicableQuad (family A Q r) (4 * (Fintype.card (FamilyIndex ι κ r) : ℚ)) := by
  have h := family_strong A Q hA hQ r
  refine ⟨h.sign, ?_, family_norm A Q hA hQ hcard r⟩
  intro k l
  rw [h.symmetric k, h.symmetric l]
  exact h.commute k l

theorem family_hadamard [Nonempty ν] [Nonempty ι] [Nonempty κ]
    (X : Fin 4 → Matrix ν ν ℚ) (p : ℕ)
    (hX : ∀ i, X i * (X i).transpose = (p : ℚ) • (1 : Matrix ν ν ℚ))
    (hcross : ∀ i j, i ≠ j → X i * (X j).transpose + X j * (X i).transpose = 0)
    (hpart : SignedPartition X) (hwidth : Fintype.card ν = 4*p)
    (A : Fin 4 → Matrix ι ι ℚ) (Q : Matrix κ κ ℚ)
    (hA : Strong (Fintype.card κ) 0 A) (hQ : IsPaley Q)
    (hcard : 2 * Fintype.card ι = Fintype.card κ + 1) (r : ℕ) :
    IsHadamard (tensorSubstitute X (family A Q r)) := by
  apply tensor_substitute_hadamard X (family A Q r) (p : ℚ)
    (4 * (Fintype.card (FamilyIndex ι κ r) : ℚ)) hX hcross hpart
    (family_amicable A Q hA hQ hcard r)
  rw [Fintype.card_prod, hwidth]
  push_cast
  ring

theorem family_endpoint_card (r : ℕ) :
    Fintype.card (ν × FamilyIndex ι κ r) =
      Fintype.card ν * Fintype.card ι * (Fintype.card κ)^r := by
  rw [Fintype.card_prod, family_card, Nat.mul_assoc]

end ButsonTowers.SpenceMatrix

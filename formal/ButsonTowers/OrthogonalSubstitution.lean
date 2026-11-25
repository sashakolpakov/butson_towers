import ButsonTowers.ScalarFrame
import Mathlib.Tactic.Linarith

namespace ButsonTowers
open scoped BigOperators Kronecker

variable {ι κ η : Type*} [Fintype ι] [DecidableEq ι]
  [Fintype κ] [DecidableEq κ] [Fintype η] [DecidableEq η]

def tensorSubstitute (X : η → Matrix ι ι ℚ) (A : η → Matrix κ κ ℚ) :
    Matrix (ι × κ) (ι × κ) ℚ := ∑ k, X k ⊗ₖ A k

def SignedPartition (X : η → Matrix ι ι ℚ) : Prop :=
  ∀ i j, ∃ k, (X k i j = 1 ∨ X k i j = -1) ∧ ∀ l, l ≠ k → X l i j = 0

def AmicableQuad (A : η → Matrix κ κ ℚ) (scale : ℚ) : Prop :=
  (∀ k, IsSignMatrix (A k)) ∧
  (∀ k l, A k * (A l).transpose = A l * (A k).transpose) ∧
  (∑ k, A k * (A k).transpose) = scale • (1 : Matrix κ κ ℚ)

omit [Fintype ι] [DecidableEq ι] [Fintype κ] [DecidableEq κ] in
lemma transpose_tensor (M : Matrix ι ι ℚ) (N : Matrix κ κ ℚ) :
    (M ⊗ₖ N).transpose = M.transpose ⊗ₖ N.transpose := by
  rfl

omit [Fintype ι] [DecidableEq ι] in
lemma skew_double_sum (B : η → η → Matrix ι ι ℚ)
    (hB : ∀ i j, i ≠ j → B i j + B j i = 0) :
    (∑ i, ∑ j, B i j) = ∑ i, B i i := by
  have hp : ∀ i j, B i j + B j i = if i = j then B i i + B i i else 0 := by
    intro i j
    by_cases h : i = j
    · subst j; simp
    · simp [h, hB i j h]
  have htwo : (∑ i, ∑ j, B i j) + (∑ i, ∑ j, B i j) =
      (∑ i, B i i) + (∑ i, B i i) := by
    calc
      _ = ∑ i, ∑ j, (B i j + B j i) := by
        simp only [Finset.sum_add_distrib]
        rw [Finset.sum_comm (f := fun i j => B j i)]
      _ = ∑ i, ∑ j, (if i = j then B i i + B i i else 0) := by
        apply Finset.sum_congr rfl
        intro i _
        apply Finset.sum_congr rfl
        intro j _
        exact hp i j
      _ = _ := by simp [Finset.sum_add_distrib]
  ext a b
  have he := congrArg (fun M : Matrix ι ι ℚ => M a b) htwo
  simp only [Matrix.add_apply, Matrix.sum_apply] at he ⊢
  linarith

theorem tensor_substitute_gram
    (X : η → Matrix ι ι ℚ) (A : η → Matrix κ κ ℚ) (p scale : ℚ)
    (hX : ∀ i, X i * (X i).transpose = p • (1 : Matrix ι ι ℚ))
    (hcross : ∀ i j, i ≠ j → X i * (X j).transpose + X j * (X i).transpose = 0)
    (hamicable : ∀ i j, A i * (A j).transpose = A j * (A i).transpose)
    (hnorm : (∑ i, A i * (A i).transpose) = scale • (1 : Matrix κ κ ℚ)) :
    tensorSubstitute X A * (tensorSubstitute X A).transpose =
      (p*scale) • (1 : Matrix (ι × κ) (ι × κ) ℚ) := by
  let B := fun i j => (X i * (X j).transpose) ⊗ₖ (A i * (A j).transpose)
  have hb : ∀ i j, i ≠ j → B i j + B j i = 0 := by
    intro i j hij
    dsimp [B]
    rw [hamicable j i, ← Matrix.add_kronecker, hcross i j hij]
    simp
  calc
    _ = ∑ i, ∑ j, B i j := by
      unfold tensorSubstitute
      simp only [Matrix.transpose_sum, Matrix.sum_mul, Matrix.mul_sum, transpose_tensor]
      simp_rw [← Matrix.mul_kronecker_mul]
      exact Finset.sum_comm
    _ = ∑ i, B i i := skew_double_sum B hb
    _ = ∑ i, (p • (1 : Matrix ι ι ℚ)) ⊗ₖ (A i * (A i).transpose) := by
      simp only [B, hX]
    _ = p • (∑ i, (1 : Matrix ι ι ℚ) ⊗ₖ (A i * (A i).transpose)) := by
      rw [Finset.smul_sum]
      simp only [Matrix.smul_kronecker]
    _ = p • ((1 : Matrix ι ι ℚ) ⊗ₖ (∑ i, A i * (A i).transpose)) := by
      congr 1
      ext ⟨i,a⟩ ⟨j,b⟩
      simp only [Matrix.sum_apply, Matrix.kroneckerMap_apply, Finset.mul_sum]
    _ = (p*scale) • (1 : Matrix (ι × κ) (ι × κ) ℚ) := by
      rw [hnorm, Matrix.kronecker_smul, Matrix.one_kronecker_one, smul_smul]

omit [Fintype ι] [DecidableEq ι] [Fintype κ] [DecidableEq κ] [DecidableEq η] in
theorem tensor_substitute_flat
    (X : η → Matrix ι ι ℚ) (A : η → Matrix κ κ ℚ)
    (hpart : SignedPartition X) (hsign : ∀ k, IsSignMatrix (A k)) :
    IsSignMatrix (tensorSubstitute X A) := by
  intro ⟨i,a⟩ ⟨j,b⟩
  obtain ⟨k,hk,hzero⟩ := hpart i j
  have hentry : tensorSubstitute X A (i,a) (j,b) = X k i j * A k a b := by
    unfold tensorSubstitute
    simp only [Matrix.sum_apply, Matrix.kronecker_apply]
    apply Finset.sum_eq_single k
    · intro l _ hl
      simp [hzero l hl]
    · simp
  rw [hentry]
  rcases hk with hk | hk <;> rcases hsign k a b with ha | ha <;>
    simp [hk, ha]

theorem tensor_substitute_hadamard [Nonempty ι] [Nonempty κ]
    (X : η → Matrix ι ι ℚ) (A : η → Matrix κ κ ℚ) (p scale : ℚ)
    (hX : ∀ i, X i * (X i).transpose = p • (1 : Matrix ι ι ℚ))
    (hcross : ∀ i j, i ≠ j → X i * (X j).transpose + X j * (X i).transpose = 0)
    (hpart : SignedPartition X)
    (hA : AmicableQuad A scale)
    (hscale : p*scale = (Fintype.card (ι × κ) : ℚ)) :
    IsHadamard (tensorSubstitute X A) := by
  refine ⟨Fintype.card_pos,tensor_substitute_flat X A hpart hA.1,?_⟩
  rw [tensor_substitute_gram X A p scale hX hcross hA.2.1 hA.2.2,hscale]

end ButsonTowers

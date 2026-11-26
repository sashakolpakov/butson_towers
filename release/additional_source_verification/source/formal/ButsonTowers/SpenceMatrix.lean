import ButsonTowers.ScalarFrame
import ButsonTowers.Spence
import Mathlib.Tactic.FinCases
import Mathlib.Tactic.Module

namespace ButsonTowers.SpenceMatrix

open Kronecker

variable {ι κ : Type*} [Fintype ι] [DecidableEq ι]
  [Fintype κ] [DecidableEq κ]

def ones : Matrix κ κ ℚ := fun _ _ => 1

structure IsPaley (Q : Matrix κ κ ℚ) : Prop where
  symmetric : Q.transpose = Q
  diagonal : ∀ i, Q i i = 0
  off_sign : ∀ i j, i ≠ j → Q i j = 1 ∨ Q i j = -1
  square : Q * Q = (Fintype.card κ : ℚ) • 1 - ones
  left_zero : ones * Q = 0
  right_zero : Q * ones = 0

structure Strong (q r : ℕ) (A : Fin 4 → Matrix ι ι ℚ) : Prop where
  sign : ∀ k, IsSignMatrix (A k)
  symmetric : ∀ k, (A k).transpose = A k
  commute : ∀ k l, A k * A l = A l * A k
  minus : (A 0 - A 1)^2 + (A 2 - A 3)^2 = (4 * (q : ℚ)^r) • 1
  plus : (A 0 + A 1)^2 + (A 2 + A 3)^2 = (4 * (q : ℚ)^(r+1)) • 1

theorem ones_square : (ones : Matrix κ κ ℚ) * ones = (Fintype.card κ : ℚ) • ones := by
  ext i j
  simp [ones, Matrix.mul_apply]

def stepPlus (A B : Matrix ι ι ℚ) (Q : Matrix κ κ ℚ) : Matrix (ι × κ) (ι × κ) ℚ :=
  (1/2 : ℚ) • ((A-B) ⊗ₖ ones + (A+B) ⊗ₖ (Q+1))

def stepMinus (A B : Matrix ι ι ℚ) (Q : Matrix κ κ ℚ) : Matrix (ι × κ) (ι × κ) ℚ :=
  (1/2 : ℚ) • ((A-B) ⊗ₖ ones + (A+B) ⊗ₖ (Q-1))

def step (A : Fin 4 → Matrix ι ι ℚ) (Q : Matrix κ κ ℚ) :
    Fin 4 → Matrix (ι × κ) (ι × κ) ℚ :=
  ![stepPlus (A 0) (A 1) Q, stepMinus (A 0) (A 1) Q,
    stepPlus (A 2) (A 3) Q, stepMinus (A 2) (A 3) Q]

def pairStep (A B : Matrix ι ι ℚ) (Q : Matrix κ κ ℚ) (plus : Bool) :
    Matrix (ι × κ) (ι × κ) ℚ :=
  (1/2 : ℚ) • ((A-B) ⊗ₖ ones + (A+B) ⊗ₖ (Q + if plus then 1 else -1))

theorem step_is_pair (A : Fin 4 → Matrix ι ι ℚ) (Q : Matrix κ κ ℚ) (k : Fin 4) :
    ∃ a b : Fin 4, ∃ plus : Bool, step A Q k = pairStep (A a) (A b) Q plus := by
  fin_cases k
  · exact ⟨0,1,true,rfl⟩
  · exact ⟨0,1,false,by simp [step, pairStep, stepMinus, sub_eq_add_neg]⟩
  · exact ⟨2,3,true,rfl⟩
  · exact ⟨2,3,false,by simp [step, pairStep, stepMinus, sub_eq_add_neg]⟩

theorem pairStep_sign (A B : Matrix ι ι ℚ) (Q : Matrix κ κ ℚ) (plus : Bool)
    (hA : IsSignMatrix A) (hB : IsSignMatrix B) (hQ : IsPaley Q) :
    IsSignMatrix (pairStep A B Q plus) := by
  intro ⟨i,a⟩ ⟨j,b⟩
  have entry : pairStep A B Q plus (i,a) (j,b) =
      Spence.entryStep (A i j) (B i j) (Q a b) (decide (a=b)) plus := by
    by_cases hab : a = b <;> cases plus <;>
      simp [pairStep, Spence.entryStep, Matrix.kroneckerMap_apply,
        Matrix.smul_apply, ones, Matrix.one_apply, hab] <;> ring
  rw [entry]
  apply Spence.entry_step_sign _ _ _ _ _ (hA i j) (hB i j)
  by_cases hab : a = b
  · subst b
    simpa using hQ.diagonal a
  · simpa [hab] using hQ.off_sign a b hab

theorem pairStep_symmetric (A B : Matrix ι ι ℚ) (Q : Matrix κ κ ℚ) (plus : Bool)
    (hA : A.transpose = A) (hB : B.transpose = B) (hQ : Q.transpose = Q) :
    (pairStep A B Q plus).transpose = pairStep A B Q plus := by
  ext ⟨i,a⟩ ⟨j,b⟩
  have hA' := congrArg (fun M : Matrix ι ι ℚ => M i j) hA
  have hB' := congrArg (fun M : Matrix ι ι ℚ => M i j) hB
  have hQ' := congrArg (fun M : Matrix κ κ ℚ => M a b) hQ
  simp only [Matrix.transpose_apply] at hA' hB' hQ'
  cases plus <;>
    simp [pairStep, Matrix.kroneckerMap_apply, Matrix.transpose_apply,
      Matrix.smul_apply, ones, Matrix.one_apply, hA', hB', hQ', eq_comm] <;> ring

theorem kronecker_commute {A B : Matrix ι ι ℚ} {C D : Matrix κ κ ℚ}
    (hab : Commute A B) (hcd : Commute C D) : Commute (A ⊗ₖ C) (B ⊗ₖ D) := by
  change (A ⊗ₖ C) * (B ⊗ₖ D) = (B ⊗ₖ D) * (A ⊗ₖ C)
  rw [← Matrix.mul_kronecker_mul, ← Matrix.mul_kronecker_mul, hab.eq, hcd.eq]

theorem smul_commute {μ : Type*} [Fintype μ] [DecidableEq μ]
    {A B : Matrix μ μ ℚ} (a b : ℚ) (h : Commute A B) : Commute (a • A) (b • B) := by
  change (a • A) * (b • B) = (b • B) * (a • A)
  simp only [Matrix.smul_mul, Matrix.mul_smul, smul_smul]
  rw [h.eq, mul_comm a b]

theorem pairStep_commute (A B C D : Matrix ι ι ℚ) (Q : Matrix κ κ ℚ)
    (u v : Bool) (hac : Commute A C) (had : Commute A D)
    (hbc : Commute B C) (hbd : Commute B D) (hQ : IsPaley Q) :
    Commute (pairStep A B Q u) (pairStep C D Q v) := by
  have hjq : Commute (ones : Matrix κ κ ℚ) Q := by
    change ones * Q = Q * ones
    rw [hQ.left_zero, hQ.right_zero]
  have hcoef (w : Bool) : Commute (ones : Matrix κ κ ℚ) (Q + if w then 1 else -1) := by
    cases w
    · exact hjq.add_right (Commute.neg_one_right _)
    · exact hjq.add_right (Commute.one_right _)
  have hcoefs : Commute (Q + if u then 1 else -1) (Q + if v then 1 else -1) := by
    have hconst (w : Bool) : Commute Q (if w then 1 else -1) := by
      cases w
      · exact Commute.neg_one_right _
      · exact Commute.one_right _
    have hconsts : Commute (if u then (1 : Matrix κ κ ℚ) else -1)
        (if v then 1 else -1) := by
      cases u
      · exact Commute.neg_one_left _
      · exact Commute.one_left _
    exact ((Commute.refl Q).add_right (hconst v)).add_left
      ((hconst u).symm.add_right hconsts)
  unfold pairStep
  apply smul_commute (1/2) (1/2)
  apply Commute.add_left
  · apply Commute.add_right
    · exact kronecker_commute ((hac.sub_right had).sub_left (hbc.sub_right hbd))
        (Commute.refl _)
    · exact kronecker_commute ((hac.add_right had).sub_left (hbc.add_right hbd)) (hcoef v)
  · apply Commute.add_right
    · exact kronecker_commute ((hac.sub_right had).add_left (hbc.sub_right hbd)) (hcoef u).symm
    · exact kronecker_commute ((hac.add_right had).add_left (hbc.add_right hbd)) hcoefs

theorem step_difference (A B : Matrix ι ι ℚ) (Q : Matrix κ κ ℚ) :
    stepPlus A B Q - stepMinus A B Q = (A+B) ⊗ₖ (1 : Matrix κ κ ℚ) := by
  ext ⟨i,a⟩ ⟨j,b⟩
  simp [stepPlus, stepMinus, Matrix.kroneckerMap_apply, Matrix.smul_apply]
  ring

theorem step_sum (A B : Matrix ι ι ℚ) (Q : Matrix κ κ ℚ) :
    stepPlus A B Q + stepMinus A B Q = (A-B) ⊗ₖ ones + (A+B) ⊗ₖ Q := by
  ext ⟨i,a⟩ ⟨j,b⟩
  simp [stepPlus, stepMinus, Matrix.kroneckerMap_apply, Matrix.smul_apply]
  ring

theorem squared_tensor_sum (U V : Matrix ι ι ℚ) (Q : Matrix κ κ ℚ)
    (hQ : IsPaley Q) :
    (U ⊗ₖ ones + V ⊗ₖ Q)^2 =
      (Fintype.card κ : ℚ) • ((V*V) ⊗ₖ (1 : Matrix κ κ ℚ)) +
      ((Fintype.card κ : ℚ) • (U*U) - V*V) ⊗ₖ ones := by
  rw [pow_two, add_mul, mul_add, mul_add]
  simp only [← Matrix.mul_kronecker_mul, ones_square, hQ.square,
    hQ.left_zero, hQ.right_zero, Matrix.kronecker_zero, add_zero, zero_add]
  ext ⟨i,a⟩ ⟨j,b⟩
  simp [Matrix.kroneckerMap_apply, Matrix.smul_apply]
  ring

theorem step_minus_norm (A : Fin 4 → Matrix ι ι ℚ) (Q : Matrix κ κ ℚ) :
    ((step A Q) 0 - (step A Q) 1)^2 + ((step A Q) 2 - (step A Q) 3)^2 =
      ((A 0 + A 1)^2 + (A 2 + A 3)^2) ⊗ₖ (1 : Matrix κ κ ℚ) := by
  simp only [step, Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val,
    Fin.reduceFinMk, step_difference, pow_two, ← Matrix.mul_kronecker_mul, one_mul]
  rw [Matrix.add_kronecker]

theorem step_plus_norm (A : Fin 4 → Matrix ι ι ℚ) (Q : Matrix κ κ ℚ)
    (hQ : IsPaley Q) :
    ((step A Q) 0 + (step A Q) 1)^2 + ((step A Q) 2 + (step A Q) 3)^2 =
      (Fintype.card κ : ℚ) •
        (((A 0 + A 1)^2 + (A 2 + A 3)^2) ⊗ₖ (1 : Matrix κ κ ℚ)) +
      ((Fintype.card κ : ℚ) • ((A 0 - A 1)^2 + (A 2 - A 3)^2) -
        ((A 0 + A 1)^2 + (A 2 + A 3)^2)) ⊗ₖ ones := by
  simp only [step, Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val,
    Fin.reduceFinMk, step_sum, squared_tensor_sum _ _ Q hQ, pow_two]
  ext ⟨i,a⟩ ⟨j,b⟩
  simp [Matrix.kroneckerMap_apply, Matrix.smul_apply]
  ring

theorem strong_step {r : ℕ} (A : Fin 4 → Matrix ι ι ℚ) (Q : Matrix κ κ ℚ)
    (hA : Strong (Fintype.card κ) r A) (hQ : IsPaley Q) :
    Strong (Fintype.card κ) (r+1) (step A Q) := by
  constructor
  · intro k
    obtain ⟨a,b,u,hu⟩ := step_is_pair A Q k
    rw [hu]
    exact pairStep_sign _ _ _ _ (hA.sign a) (hA.sign b) hQ
  · intro k
    obtain ⟨a,b,u,hu⟩ := step_is_pair A Q k
    rw [hu]
    exact pairStep_symmetric _ _ _ _ (hA.symmetric a) (hA.symmetric b) hQ.symmetric
  · intro k l
    obtain ⟨a,b,u,hu⟩ := step_is_pair A Q k
    obtain ⟨c,d,v,hv⟩ := step_is_pair A Q l
    rw [hu,hv]
    exact (pairStep_commute _ _ _ _ _ _ _ (hA.commute a c) (hA.commute a d)
      (hA.commute b c) (hA.commute b d) hQ).eq
  · rw [step_minus_norm, hA.plus, Matrix.smul_kronecker, Matrix.one_kronecker_one]
  · rw [step_plus_norm A Q hQ, hA.minus, hA.plus]
    simp only [Matrix.smul_kronecker, Matrix.one_kronecker_one, smul_smul]
    have hpow : (Fintype.card κ : ℚ) * (4*(Fintype.card κ : ℚ)^r) =
        4*(Fintype.card κ : ℚ)^(r+1) := by rw [pow_succ]; ring
    rw [hpow, sub_self, Matrix.zero_kronecker, add_zero]
    congr 1
    rw [pow_succ]
    ring

theorem Strong.sum_squares {q r : ℕ} {A : Fin 4 → Matrix ι ι ℚ}
    (hA : Strong q r A) :
    (A 0)^2 + (A 1)^2 + (A 2)^2 + (A 3)^2 =
      (2*(q:ℚ)^r*((q:ℚ)+1)) • (1 : Matrix ι ι ℚ) := by
  calc
    _ = (1/2 : ℚ) • (((A 0 - A 1)^2 + (A 2 - A 3)^2) +
        ((A 0 + A 1)^2 + (A 2 + A 3)^2)) := by
      simp only [pow_two, sub_mul, mul_sub, add_mul, mul_add]
      module
    _ = (1/2 : ℚ) • (((4*(q:ℚ)^r) + (4*(q:ℚ)^(r+1))) • (1 : Matrix ι ι ℚ)) := by
      rw [hA.minus,hA.plus,add_smul]
    _ = _ := by
      rw [smul_smul]
      congr 1
      rw [pow_succ]
      ring

def FamilyIndex (ι κ : Type) : ℕ → Type
  | 0 => ι
  | r+1 => FamilyIndex ι κ r × κ

instance familyFintype (ι κ : Type) [Fintype ι] [Fintype κ] :
    (r : ℕ) → Fintype (FamilyIndex ι κ r)
  | 0 => inferInstanceAs (Fintype ι)
  | r+1 => by
    letI := familyFintype ι κ r
    exact inferInstanceAs (Fintype (FamilyIndex ι κ r × κ))

instance familyDecidableEq (ι κ : Type) [DecidableEq ι] [DecidableEq κ] :
    (r : ℕ) → DecidableEq (FamilyIndex ι κ r)
  | 0 => inferInstanceAs (DecidableEq ι)
  | r+1 => by
    letI := familyDecidableEq ι κ r
    exact inferInstanceAs (DecidableEq (FamilyIndex ι κ r × κ))

instance familyNonempty (ι κ : Type) [Nonempty ι] [Nonempty κ] :
    (r : ℕ) → Nonempty (FamilyIndex ι κ r)
  | 0 => inferInstanceAs (Nonempty ι)
  | r+1 => by
    letI := familyNonempty ι κ r
    exact inferInstanceAs (Nonempty (FamilyIndex ι κ r × κ))

theorem family_card (ι κ : Type) [Fintype ι] [Fintype κ] (r : ℕ) :
    Fintype.card (FamilyIndex ι κ r) = Fintype.card ι * (Fintype.card κ)^r := by
  induction r with
  | zero =>
    rw [pow_zero, Nat.mul_one]
    exact @Fintype.card_congr (FamilyIndex ι κ 0) ι
      (familyFintype ι κ 0) ‹Fintype ι› (Equiv.refl ι)
  | succ r ih =>
    change Fintype.card (FamilyIndex ι κ r × κ) = _
    rw [Fintype.card_prod, ih, pow_succ, Nat.mul_assoc]

def family {ι κ : Type} [Fintype ι] [DecidableEq ι] [Fintype κ] [DecidableEq κ]
    (A : Fin 4 → Matrix ι ι ℚ) (Q : Matrix κ κ ℚ) :
    (r : ℕ) → Fin 4 → Matrix (FamilyIndex ι κ r) (FamilyIndex ι κ r) ℚ
  | 0 => A
  | r+1 => step (family A Q r) Q

theorem family_strong {ι κ : Type} [Fintype ι] [DecidableEq ι]
    [Fintype κ] [DecidableEq κ]
    (A : Fin 4 → Matrix ι ι ℚ) (Q : Matrix κ κ ℚ)
    (hA : Strong (Fintype.card κ) 0 A) (hQ : IsPaley Q) :
    ∀ r, Strong (Fintype.card κ) r (family A Q r) := by
  intro r
  induction r with
  | zero => exact hA
  | succ r ih => exact strong_step _ _ ih hQ

theorem family_norm {ι κ : Type} [Fintype ι] [DecidableEq ι]
    [Fintype κ] [DecidableEq κ]
    (A : Fin 4 → Matrix ι ι ℚ) (Q : Matrix κ κ ℚ)
    (hA : Strong (Fintype.card κ) 0 A) (hQ : IsPaley Q)
    (hcard : 2 * Fintype.card ι = Fintype.card κ + 1) (r : ℕ) :
    ∑ k : Fin 4, family A Q r k * (family A Q r k).transpose =
      (4 * (Fintype.card (FamilyIndex ι κ r) : ℚ)) • 1 := by
  have h := family_strong A Q hA hQ r
  simp only [h.symmetric]
  have hc : 2 * (Fintype.card ι : ℚ) = (Fintype.card κ : ℚ) + 1 := by exact_mod_cast hcard
  have hn : 2*(Fintype.card κ : ℚ)^r*((Fintype.card κ : ℚ)+1) =
      4*(Fintype.card (FamilyIndex ι κ r) : ℚ) := by
    rw [family_card]
    push_cast
    rw [hc.symm]
    ring
  have hs := h.sum_squares
  rw [hn] at hs
  simpa [Fin.sum_univ_succ, pow_two, add_assoc] using hs

end ButsonTowers.SpenceMatrix

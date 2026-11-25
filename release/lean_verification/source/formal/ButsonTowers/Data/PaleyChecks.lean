import ButsonTowers.Data.FiniteChecks
import ButsonTowers.GroupAlgebra
import ButsonTowers.SpenceMatrix

namespace ButsonTowers.Data

open scoped BigOperators

theorem paley_of_integer_rows {G : Type} [AddCommGroup G] [Fintype G] [DecidableEq G]
    (r : G → ℤ) (hz : r 0 = 0) (hsym : ∀ x, r (-x) = r x)
    (hsign : ∀ x, x ≠ 0 → r x = 1 ∨ r x = -1)
    (hsum : ∑ x, r x = 0)
    (hcorr : ∀ s, ∑ x, r x * r (x+s) =
      if s = 0 then (Fintype.card G : ℤ)-1 else -1) :
    SpenceMatrix.IsPaley (GroupAlgebra.regular (fun x => (r x : ℚ))) := by
  let Q := GroupAlgebra.regular (fun x => (r x : ℚ))
  have hsymQ : Q.transpose = Q := by
    ext i j
    change (r (i-j) : ℚ) = (r (j-i) : ℚ)
    have h := hsym (j-i)
    have hh : -(j-i) = i-j := by abel
    rw [hh] at h
    exact_mod_cast h
  have hright : Q * SpenceMatrix.ones = 0 := by
    ext i j
    change (∑ k, (r (k-i) : ℚ) * 1) = 0
    simp only [mul_one]
    have hcast : (∑ x, (r x : ℚ)) = 0 := by exact_mod_cast hsum
    calc
      _ = ∑ x, (r x : ℚ) := by
        exact Fintype.sum_equiv (Equiv.subRight i) _ _ (by intro x; rfl)
      _ = 0 := hcast
  have hsquare : Q * Q = (Fintype.card G : ℚ) • 1 - SpenceMatrix.ones := by
    conv_lhs => arg 2; rw [← hsymQ]
    ext i j
    dsimp only [Q]
    rw [GroupAlgebra.regular_gram_apply]
    have hi : (∑ x, (r x : ℚ) * (r (x+(i-j)) : ℚ)) =
        if i-j=0 then (Fintype.card G : ℚ)-1 else -1 := by exact_mod_cast hcorr (i-j)
    rw [hi]
    by_cases hij : i = j <;>
      simp [SpenceMatrix.ones, Matrix.smul_apply, Matrix.one_apply, hij, sub_eq_zero]
  refine ⟨hsymQ, ?_, ?_, hsquare, ?_, hright⟩
  · intro i
    change (r (i-i) : ℚ) = 0
    simp [hz]
  · intro i j hij
    obtain h | h := hsign (j-i) (sub_ne_zero.mpr hij.symm)
    · left; change (r (j-i) : ℚ) = 1; exact_mod_cast h
    · right; change (r (j-i) : ℚ) = -1; exact_mod_cast h
  · have ht := congrArg Matrix.transpose hright
    simpa only [Matrix.transpose_mul, hsymQ, Matrix.transpose_zero,
      show (SpenceMatrix.ones : Matrix G G ℚ).transpose = SpenceMatrix.ones from rfl] using ht

def checkPaleySign (p : ℕ) (a : ℕ → ℤ) : Bool :=
  (List.range p).all (fun x => decide (if x=0 then a x=0 else a x=1 ∨ a x= -1))

def checkPaleySymmetric (p : ℕ) (a : ℕ → ℤ) : Bool :=
  (List.range p).all (fun x => decide (a ((p-x)%p) = a x))

def checkPaleyCorr (p : ℕ) (a : ℕ → ℤ) : Bool :=
  (List.range p).all (fun s => decide
    (natSum p (fun x => a x * a ((x+s)%p)) = if s=0 then (p:ℤ)-1 else -1))

theorem prime_paley_checked {p : ℕ} [NeZero p] (a : ℕ → ℤ)
    (hsg : checkPaleySign p a = true) (hsy : checkPaleySymmetric p a = true)
    (hsm : natSum p a = 0) (hc : checkPaleyCorr p a = true) :
    SpenceMatrix.IsPaley (GroupAlgebra.regular (fun x : ZMod p => (a x.val : ℚ))) := by
  apply paley_of_integer_rows (fun x : ZMod p => a x.val)
  · have h := of_decide_eq_true
      ((List.all_eq_true.mp hsg) 0 (List.mem_range.mpr (NeZero.pos p)))
    simpa using h
  · intro x
    have h := of_decide_eq_true
      ((List.all_eq_true.mp hsy) x.val (List.mem_range.mpr x.val_lt))
    have he : (-x).val = (p-x.val)%p := by
      have hcast : ((p-x.val:ℕ):ZMod p) = -x := by
        rw [Nat.cast_sub x.val_lt.le]; simp
      rw [← hcast, ZMod.val_natCast]
    simpa only [he] using h
  · intro x hx
    have h := of_decide_eq_true
      ((List.all_eq_true.mp hsg) x.val (List.mem_range.mpr x.val_lt))
    have hn : x.val ≠ 0 := by
      intro hh
      apply hx
      rw [← ZMod.natCast_zmod_val x, hh]; simp
    simpa [hn] using h
  · simpa only [natSum_eq_sum_zmod] using hsm
  · intro s
    have h := of_decide_eq_true
      ((List.all_eq_true.mp hc) s.val (List.mem_range.mpr s.val_lt))
    rw [natSum_eq_sum_zmod] at h
    have hz : s.val=0 ↔ s=0 := by
      constructor
      · intro hh; rw [← ZMod.natCast_zmod_val s,hh]; simp
      · rintro rfl; simp
    simpa only [ZMod.val_add, hz, ZMod.card] using h

end ButsonTowers.Data

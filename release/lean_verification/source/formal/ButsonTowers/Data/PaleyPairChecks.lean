import ButsonTowers.Data.PaleyChecks

namespace ButsonTowers.Data
open scoped BigOperators

abbrev Paley121Group := ZMod 11 × ZMod 11

def pairEncode (x : Paley121Group) : ℕ := x.1.val + 11*x.2.val
def pairToZMod (x : Paley121Group) : ZMod 121 := pairEncode x

theorem pairEncode_lt (x : Paley121Group) : pairEncode x < 121 := by
  have h₁ := x.1.val_lt
  have h₂ := x.2.val_lt
  unfold pairEncode
  omega

theorem pairEncode_mod (x : Paley121Group) : pairEncode x % 11 = x.1.val := by
  have h := x.1.val_lt
  unfold pairEncode
  omega

theorem pairEncode_div (x : Paley121Group) : pairEncode x / 11 = x.2.val := by
  have h := x.1.val_lt
  unfold pairEncode
  omega

theorem pairEncode_zero (x : Paley121Group) : pairEncode x = 0 ↔ x = 0 := by
  constructor
  · intro h
    apply Prod.ext <;> apply ZMod.val_injective
    · change x.1.val = 0
      unfold pairEncode at h
      omega
    · change x.2.val = 0
      unfold pairEncode at h
      omega
  · rintro rfl
    rfl

theorem pairToZMod_bijective : Function.Bijective pairToZMod := by
  constructor
  · intro x y h
    have hh := congrArg ZMod.val h
    simp only [pairToZMod, ZMod.val_natCast, Nat.mod_eq_of_lt (pairEncode_lt _)] at hh
    have hm := congrArg (fun n => n%11) hh
    have hd := congrArg (fun n => n/11) hh
    rw [pairEncode_mod,pairEncode_mod] at hm
    rw [pairEncode_div,pairEncode_div] at hd
    exact Prod.ext (ZMod.val_injective 11 hm) (ZMod.val_injective 11 hd)
  · intro y
    let x : Paley121Group := ((y.val%11 : ℕ), (y.val/11 : ℕ))
    refine ⟨x, ?_⟩
    have hd : y.val/11 < 11 := by have h := y.val_lt; omega
    have he : pairEncode x = y.val := by
      simp only [pairEncode, x, ZMod.val_natCast, Nat.mod_mod, Nat.mod_eq_of_lt hd]
      omega
    change (pairEncode x : ZMod 121) = y
    rw [he, ZMod.natCast_zmod_val]

theorem natSum_eq_sum_pair (f : ℕ → ℤ) :
    natSum 121 f = ∑ x : Paley121Group, f (pairEncode x) := by
  rw [natSum_eq_sum_zmod]
  symm
  apply Fintype.sum_equiv (Equiv.ofBijective pairToZMod pairToZMod_bijective)
  intro x
  simp [pairToZMod, ZMod.val_natCast, Nat.mod_eq_of_lt (pairEncode_lt x)]

def pairAddNat (x s : ℕ) : ℕ :=
  (x%11+s%11)%11 + 11*((x/11+s/11)%11)

def pairNegNat (x : ℕ) : ℕ := (11-x%11)%11 + 11*((11-x/11)%11)

theorem pairEncode_add (x y : Paley121Group) :
    pairEncode (x+y) = pairAddNat (pairEncode x) (pairEncode y) := by
  unfold pairAddNat
  rw [pairEncode_mod,pairEncode_mod,pairEncode_div,pairEncode_div]
  simp only [pairEncode, Prod.fst_add, Prod.snd_add, ZMod.val_add]

theorem val_neg_eleven (x : ZMod 11) : (-x).val = (11-x.val)%11 := by
  have h : ((11-x.val : ℕ) : ZMod 11) = -x := by
    rw [Nat.cast_sub x.val_lt.le]
    rw [ZMod.natCast_self]
    simp
  rw [← h, ZMod.val_natCast]

theorem pairEncode_neg (x : Paley121Group) :
    pairEncode (-x) = pairNegNat (pairEncode x) := by
  unfold pairNegNat
  rw [pairEncode_mod,pairEncode_div]
  simp only [pairEncode, Prod.fst_neg, Prod.snd_neg, val_neg_eleven]

def checkPairPaleySymmetric (a : ℕ → ℤ) : Bool :=
  (List.range 121).all (fun x => decide (a (pairNegNat x) = a x))

def checkPairPaleyCorr (a : ℕ → ℤ) : Bool :=
  (List.range 121).all (fun s => decide
    (natSum 121 (fun x => a x * a (pairAddNat x s)) = if s=0 then 120 else -1))

theorem pair_paley_checked (a : ℕ → ℤ)
    (hsg : checkPaleySign 121 a = true) (hsy : checkPairPaleySymmetric a = true)
    (hsm : natSum 121 a = 0) (hc : checkPairPaleyCorr a = true) :
    SpenceMatrix.IsPaley (GroupAlgebra.regular (fun x : Paley121Group => (a (pairEncode x) : ℚ))) := by
  apply paley_of_integer_rows (fun x : Paley121Group => a (pairEncode x))
  · have h := of_decide_eq_true
      ((List.all_eq_true.mp hsg) 0 (List.mem_range.mpr (by decide : 0<121)))
    simpa [pairEncode] using h
  · intro x
    have h := of_decide_eq_true
      ((List.all_eq_true.mp hsy) (pairEncode x) (List.mem_range.mpr (pairEncode_lt x)))
    simpa only [← pairEncode_neg] using h
  · intro x hx
    have h := of_decide_eq_true
      ((List.all_eq_true.mp hsg) (pairEncode x) (List.mem_range.mpr (pairEncode_lt x)))
    simpa [pairEncode_zero, hx] using h
  · simpa only [natSum_eq_sum_pair] using hsm
  · intro s
    have h := of_decide_eq_true
      ((List.all_eq_true.mp hc) (pairEncode s) (List.mem_range.mpr (pairEncode_lt s)))
    rw [natSum_eq_sum_pair] at h
    simpa [← pairEncode_add, pairEncode_zero, Fintype.card_prod, ZMod.card] using h

end ButsonTowers.Data

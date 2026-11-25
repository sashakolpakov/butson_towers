import ButsonTowers.Data.FiniteChecks
import ButsonTowers.GroupAlgebra
import ButsonTowers.SpenceMatrix

namespace ButsonTowers.Data
open scoped BigOperators

def quadDot (a : Fin 4 → ℕ → ℤ) (x y : ℕ) : ℤ :=
  a 0 x * a 0 y + a 1 x * a 1 y + a 2 x * a 2 y + a 3 x * a 3 y

theorem quadDot_eq_sum (a : Fin 4 → ℕ → ℤ) (x y : ℕ) :
    quadDot a x y = ∑ i : Fin 4, a i x * a i y := by
  simp [quadDot, Fin.sum_univ_succ, add_assoc]

def quadCorr (p : ℕ) (a : Fin 4 → ℕ → ℤ) (s : ℕ) : ℤ :=
  natSum p (fun x => quadDot a x ((x+s)%p))

def checkQuadCorr (p : ℕ) (a : Fin 4 → ℕ → ℤ) (w : ℤ) : Bool :=
  (List.range p).all (fun s => decide (quadCorr p a s = if s=0 then w else 0))

def checkQuadSign (p : ℕ) (a : Fin 4 → ℕ → ℤ) : Bool :=
  (List.range p).all (fun x => decide (∀ i : Fin 4, a i x = 1 ∨ a i x = -1))

def checkQuadSymmetric (p : ℕ) (a : Fin 4 → ℕ → ℤ) : Bool :=
  (List.range p).all (fun x => decide (∀ i : Fin 4, a i ((p-x)%p) = a i x))

theorem val_neg_nat {p : ℕ} [NeZero p] (x : ZMod p) :
    (-x).val = (p-x.val)%p := by
  have h : ((p-x.val : ℕ) : ZMod p) = -x := by
    rw [Nat.cast_sub x.val_lt.le]
    simp
  rw [← h, ZMod.val_natCast]

theorem val_eq_zero_iff {p : ℕ} [NeZero p] (x : ZMod p) :
    x.val = 0 ↔ x = 0 := by
  constructor
  · intro h
    rw [← ZMod.natCast_zmod_val x, h]
    simp
  · rintro rfl
    simp

theorem quad_corr_sound {p : ℕ} [NeZero p] {a : Fin 4 → ℕ → ℤ} {w : ℤ}
    (hc : checkQuadCorr p a w = true) :
    ∀ s : ZMod p, (∑ i : Fin 4, ∑ x : ZMod p, a i x.val * a i (x+s).val) =
      if s=0 then w else 0 := by
  intro s
  have h := of_decide_eq_true
    ((List.all_eq_true.mp hc) s.val (List.mem_range.mpr s.val_lt))
  change quadCorr p a s.val = _ at h
  rw [quadCorr, natSum_eq_sum_zmod] at h
  rw [Finset.sum_comm]
  simpa only [ZMod.val_add, quadDot_eq_sum, val_eq_zero_iff] using h

theorem quad_sign_sound {p : ℕ} [NeZero p] {a : Fin 4 → ℕ → ℤ}
    (hc : checkQuadSign p a = true) (i : Fin 4) (x : ZMod p) :
    a i x.val = 1 ∨ a i x.val = -1 := by
  have h : ∀ i : Fin 4, a i x.val=1 ∨ a i x.val= -1 :=
    of_decide_eq_true ((List.all_eq_true.mp hc) x.val (List.mem_range.mpr x.val_lt))
  exact h i

theorem quad_symmetric_sound {p : ℕ} [NeZero p] {a : Fin 4 → ℕ → ℤ}
    (hc : checkQuadSymmetric p a = true) (i : Fin 4) (x : ZMod p) :
    a i (-x).val = a i x.val := by
  have h : ∀ i : Fin 4, a i ((p-x.val)%p) = a i x.val :=
    of_decide_eq_true ((List.all_eq_true.mp hc) x.val (List.mem_range.mpr x.val_lt))
  simpa only [val_neg_nat] using h i

def rationalRows {p : ℕ} (a : Fin 4 → ℕ → ℤ) (i : Fin 4) (x : ZMod p) : ℚ :=
  (a i x.val : ℚ)

def rowMatrices (p : ℕ) (a : Fin 4 → ℕ → ℤ) :
    Fin 4 → Matrix (ZMod p) (ZMod p) ℚ :=
  fun i => GroupAlgebra.regular (rationalRows a i)

theorem quad_matrices_sign {p : ℕ} [NeZero p] {a : Fin 4 → ℕ → ℤ}
    (h : checkQuadSign p a = true) : ∀ i, IsSignMatrix (rowMatrices p a i) := by
  intro i x y
  obtain h | h := quad_sign_sound h i (y-x)
  · left
    change ((a i (y-x).val : ℤ) : ℚ) = 1
    exact_mod_cast h
  · right
    change ((a i (y-x).val : ℤ) : ℚ) = -1
    exact_mod_cast h

theorem quad_matrices_symmetric {p : ℕ} [NeZero p] {a : Fin 4 → ℕ → ℤ}
    (h : checkQuadSymmetric p a = true) :
    ∀ i, (rowMatrices p a i).transpose = rowMatrices p a i := by
  intro i
  ext x y
  change (a i (x-y).val : ℚ) = (a i (y-x).val : ℚ)
  have hi := quad_symmetric_sound h i (y-x)
  have hh : -(y-x) = x-y := by abel
  rw [hh] at hi
  exact_mod_cast hi

theorem quad_matrices_commute {p : ℕ} [NeZero p] (a : Fin 4 → ℕ → ℤ) :
    ∀ i j, rowMatrices p a i * rowMatrices p a j = rowMatrices p a j * rowMatrices p a i := by
  intro i j
  exact GroupAlgebra.regular_mul_comm _ _

theorem quad_matrices_norm {p : ℕ} [NeZero p] {a : Fin 4 → ℕ → ℤ} {w : ℤ}
    (h : checkQuadCorr p a w = true) :
    (∑ i : Fin 4, rowMatrices p a i * (rowMatrices p a i).transpose) =
      (w : ℚ) • (1 : Matrix (ZMod p) (ZMod p) ℚ) := by
  apply GroupAlgebra.sum_regular_gram
  intro s
  change (∑ i : Fin 4, ∑ x : ZMod p,
    (a i x.val : ℚ) * (a i (x+s).val : ℚ)) = if s=0 then (w : ℚ) else 0
  have hi := quad_corr_sound h s
  exact_mod_cast hi

end ButsonTowers.Data

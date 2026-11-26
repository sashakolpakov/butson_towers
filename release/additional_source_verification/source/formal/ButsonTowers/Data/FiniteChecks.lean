import Mathlib.Data.ZMod.Basic
import Mathlib.Algebra.BigOperators.Fin
import Mathlib.Tactic.IntervalCases
import Mathlib.Tactic.NormNum

/-! Sound small checkers for literal finite rows. All reflected facts below
are closed with ordinary kernel reduction, never an external oracle. -/

namespace ButsonTowers.Data

open scoped BigOperators

def natSum (n : ℕ) (f : ℕ → ℤ) : ℤ :=
  (List.range n).foldl (fun acc i => acc + f i) 0

theorem natSum_eq_sum_range (n : ℕ) (f : ℕ → ℤ) :
    natSum n f = ∑ x ∈ Finset.range n, f x := by
  induction n with
  | zero => simp [natSum]
  | succ n ih =>
    simpa [natSum, List.range_succ, List.foldl_append, Finset.sum_range_succ]
      using congrArg (fun z => z + f n) ih

theorem natSum_eq_sum_zmod (n : ℕ) [NeZero n] (f : ℕ → ℤ) :
    natSum n f = ∑ x : ZMod n, f x.val := by
  rw [natSum_eq_sum_range]
  have h := Fintype.sum_equiv (ZMod.finEquiv n).toEquiv
    (fun x : Fin n => f x.val) (fun x : ZMod n => f x.val) (by
      intro x
      cases n with
      | zero => exact (NeZero.ne 0 rfl).elim
      | succ n => rfl)
  rw [← h, Fin.sum_univ_eq_sum_range]

def unitCode (a : ℕ) (i : Fin 4) : ℤ :=
  if a / 2 = i.val then (if a % 2 = 0 then 1 else -1) else 0

def dotCode (a b : ℕ) : ℤ :=
  if a / 2 = b / 2 then (if a % 2 = b % 2 then 1 else -1) else 0

theorem dotCode_eq_sum (a b : ℕ) (ha : a < 8) (hb : b < 8) :
    dotCode a b = ∑ i : Fin 4, unitCode a i * unitCode b i := by
  interval_cases a <;> interval_cases b <;> decide

theorem unitCode_onehot (a : ℕ) (ha : a < 8) :
    ∃ c : Fin 4, (unitCode a c = 1 ∨ unitCode a c = -1) ∧
      ∀ d : Fin 4, d ≠ c → unitCode a d = 0 := by
  interval_cases a <;> decide

def codeCorr (p : ℕ) (codes : ℕ → ℕ) (s : ℕ) : ℤ :=
  natSum p (fun x => dotCode (codes x) (codes ((x+s)%p)))

def checkCodesBound (p : ℕ) (codes : ℕ → ℕ) : Bool :=
  (List.range p).all (fun x => decide (codes x < 8))

def checkCodeCorr (p : ℕ) (codes : ℕ → ℕ) : Bool :=
  (List.range p).all (fun s =>
    decide (codeCorr p codes s = if s = 0 then (p : ℤ) else 0))

theorem codes_bound_sound {p : ℕ} {codes : ℕ → ℕ}
    (h : checkCodesBound p codes = true) (x : ℕ) (hx : x < p) :
    codes x < 8 := by
  have := (List.all_eq_true.mp h) x (List.mem_range.mpr hx)
  exact of_decide_eq_true this

theorem code_corr_sound {p : ℕ} [NeZero p] {codes : ℕ → ℕ}
    (hb : checkCodesBound p codes = true) (hc : checkCodeCorr p codes = true) :
    ∀ s : ZMod p,
      (∑ i : Fin 4, ∑ x : ZMod p,
        unitCode (codes x.val) i * unitCode (codes (x+s).val) i) =
        if s = 0 then (p : ℤ) else 0 := by
  intro s
  have h := of_decide_eq_true
    ((List.all_eq_true.mp hc) s.val (List.mem_range.mpr s.val_lt))
  change codeCorr p codes s.val = _ at h
  rw [codeCorr, natSum_eq_sum_zmod] at h
  rw [Finset.sum_comm]
  have hz : (s.val = 0) ↔ s = 0 := by
    constructor
    · intro hs
      rw [← ZMod.natCast_zmod_val s, hs]
      simp
    · rintro rfl
      simp
  simpa only [ZMod.val_add, dotCode_eq_sum _ _
    (codes_bound_sound hb _ (ZMod.val_lt _))
    (codes_bound_sound hb _ (Nat.mod_lt _ (NeZero.pos p))), hz] using h

theorem cast_int_corr {p : ℕ} [NeZero p] (b : Fin 4 → ZMod p → ℤ)
    (h : ∀ s : ZMod p, (∑ i : Fin 4, ∑ x : ZMod p, b i x * b i (x+s)) =
      if s = 0 then (p : ℤ) else 0) :
    ∀ s : ZMod p, (∑ i : Fin 4, ∑ x : ZMod p,
      (b i x : ℚ) * (b i (x+s) : ℚ)) = if s = 0 then (p : ℚ) else 0 := by
  intro s
  exact_mod_cast h s

theorem cast_int_onehot {G : Type*} (b : Fin 4 → G → ℤ)
    (h : ∀ x, ∃ c : Fin 4, (b c x = 1 ∨ b c x = -1) ∧
      ∀ d : Fin 4, d ≠ c → b d x = 0) :
    ∀ x, ∃ c : Fin 4, ((b c x : ℚ) = 1 ∨ (b c x : ℚ) = -1) ∧
      ∀ d : Fin 4, d ≠ c → (b d x : ℚ) = 0 := by
  intro x
  obtain ⟨c, hc, hd⟩ := h x
  refine ⟨c, ?_, ?_⟩
  · rcases hc with hc | hc
    · left; exact_mod_cast hc
    · right; exact_mod_cast hc
  · intro d hdc
    exact_mod_cast hd d hdc

end ButsonTowers.Data

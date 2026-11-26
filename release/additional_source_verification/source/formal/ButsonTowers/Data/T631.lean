import ButsonTowers.Data.T631Blocked

namespace ButsonTowers.Data
open scoped BigOperators
set_option maxRecDepth 1000000
set_option maxHeartbeats 0

theorem t631_codes_bound : checkCodesBound 631 t631CodesFast = true := by decide +kernel

theorem t631_codes_corr : checkCodeCorr 631 t631CodesFast = true :=
  t631_codes_corr_blocked

def t631Int (i : Fin 4) (x : ZMod 631) : ℤ := unitCode (t631CodesFast x.val) i

def t631 (i : Fin 4) (x : ZMod 631) : ℚ := (t631Int i x : ℚ)

theorem t631_int_corr : ∀ s : ZMod 631,
    (∑ i : Fin 4, ∑ x : ZMod 631, t631Int i x * t631Int i (x+s)) =
      if s=0 then (631 : ℤ) else 0 :=
  @code_corr_sound 631 inferInstance t631CodesFast t631_codes_bound t631_codes_corr

theorem t631_corr : ∀ s : ZMod 631,
    (∑ i : Fin 4, ∑ x : ZMod 631, t631 i x * t631 i (x+s)) =
      if s=0 then (631 : ℚ) else 0 :=
  @cast_int_corr 631 inferInstance t631Int t631_int_corr

theorem t631_onehot : ∀ x : ZMod 631, ∃ c : Fin 4,
    (t631 c x=1 ∨ t631 c x= -1) ∧ ∀ d : Fin 4, d≠c → t631 d x=0 := by
  apply cast_int_onehot t631Int
  intro x
  exact unitCode_onehot _
    (@codes_bound_sound 631 t631CodesFast t631_codes_bound x.val x.val_lt)

#print axioms t631_codes_corr
#print axioms t631_corr
#print axioms t631_onehot
end ButsonTowers.Data

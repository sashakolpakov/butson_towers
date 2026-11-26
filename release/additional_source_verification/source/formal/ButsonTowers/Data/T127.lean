import ButsonTowers.Data.FiniteChecks
import ButsonTowers.Data.TreeLiterals

namespace ButsonTowers.Data
open scoped BigOperators
set_option maxRecDepth 1000000
set_option maxHeartbeats 0

theorem t127_codes_bound : checkCodesBound 127 t127CodesFast = true := by decide +kernel

theorem t127_codes_corr : checkCodeCorr 127 t127CodesFast = true := by decide +kernel

def t127Int (i : Fin 4) (x : ZMod 127) : ℤ := unitCode (t127CodesFast x.val) i

def t127 (i : Fin 4) (x : ZMod 127) : ℚ := (t127Int i x : ℚ)

theorem t127_int_corr : ∀ s : ZMod 127,
    (∑ i : Fin 4, ∑ x : ZMod 127, t127Int i x * t127Int i (x+s)) =
      if s=0 then (127 : ℤ) else 0 :=
  code_corr_sound t127_codes_bound t127_codes_corr

theorem t127_corr : ∀ s : ZMod 127,
    (∑ i : Fin 4, ∑ x : ZMod 127, t127 i x * t127 i (x+s)) =
      if s=0 then (127 : ℚ) else 0 :=
  cast_int_corr t127Int t127_int_corr

theorem t127_onehot : ∀ x : ZMod 127, ∃ c : Fin 4,
    (t127 c x=1 ∨ t127 c x= -1) ∧ ∀ d : Fin 4, d≠c → t127 d x=0 := by
  apply cast_int_onehot t127Int
  intro x
  exact unitCode_onehot _ (codes_bound_sound t127_codes_bound _ x.val_lt)

#print axioms t127_codes_corr
#print axioms t127_corr
#print axioms t127_onehot
end ButsonTowers.Data

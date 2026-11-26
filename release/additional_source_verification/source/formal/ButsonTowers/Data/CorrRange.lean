import ButsonTowers.Data.FiniteChecks

namespace ButsonTowers.Data

theorem code_corr_of_pointwise {p : ℕ} {codes : ℕ → ℕ}
    (h : ∀ s, s < p → codeCorr p codes s = if s=0 then (p : ℤ) else 0) :
    checkCodeCorr p codes = true := by
  apply List.all_eq_true.mpr
  intro s hs
  exact decide_eq_true (h s (List.mem_range.mp hs))

def checkCorrRange (p : ℕ) (codes : ℕ → ℕ) (start count : ℕ) : Bool :=
  (List.range' start count).all (fun s =>
    decide (codeCorr p codes s = if s=0 then (p : ℤ) else 0))

theorem corr_range_sound {p start count : ℕ} {codes : ℕ → ℕ}
    (h : checkCorrRange p codes start count = true) (s : ℕ)
    (hlo : start ≤ s) (hhi : s < start+count) :
    codeCorr p codes s = if s=0 then (p : ℤ) else 0 := by
  apply of_decide_eq_true
  apply (List.all_eq_true.mp h) s
  apply List.mem_range'.mpr
  refine ⟨s-start, ?_, ?_⟩
  · omega
  · simp only [one_mul]
    omega

end ButsonTowers.Data

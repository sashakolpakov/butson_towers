import ButsonTowers.AdditionalSpenceFamilies

/-! Further finite specializations and the shifted q=109 family. These
existence statements have no literature-status assumptions. -/

namespace ButsonTowers

theorem hadamard_2985892 : HadamardExists 2985892 := by
  simpa using spence_family_13 2

theorem hadamard_38816596 : HadamardExists 38816596 := by
  simpa using spence_family_13 3

theorem hadamard_291145924 : HadamardExists 291145924 := by
  simpa using spence_family_61 2

theorem hadamard_15131380 : HadamardExists 15131380 := by
  simpa using spence_family_109 1

theorem new_spence_family_109 (r : ℕ) : HadamardExists (15131380*109^r) := by
  have h := spence_family_109 (r+1)
  have he : 138820*109^(r+1) = 15131380*109^r := by
    rw [pow_succ]
    ring
  rw [he] at h
  exact h

end ButsonTowers

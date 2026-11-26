import ButsonTowers.Endpoints
import ButsonTowers.Data.AdditionalSpence

/-! Six additional all-parameter existence theorems from literal strong
starters. These propositions assert existence, independently of novelty. -/

namespace ButsonTowers

theorem spence_family_5 (r : ℕ) : HadamardExists (7572*5^r) := by
  have h := spence_exists_from_t Data.t631 Data.t631_corr Data.t631_onehot
    Data.additionalSpence5 Data.additionalPaley5
    (by simpa using Data.additionalSpence5_strong)
    Data.additionalPaley5_isPaley (by norm_num) r
  simpa using h

theorem spence_family_53 (r : ℕ) : HadamardExists (68148*53^r) := by
  have h := spence_exists_from_t Data.t631 Data.t631_corr Data.t631_onehot
    Data.additionalSpence53 Data.additionalPaley53
    (by simpa using Data.additionalSpence53_strong)
    Data.additionalPaley53_isPaley (by norm_num) r
  simpa using h

theorem spence_family_97 (r : ℕ) : HadamardExists (123676*97^r) := by
  have h := spence_exists_from_t Data.t631 Data.t631_corr Data.t631_onehot
    Data.additionalSpence97 Data.additionalPaley97
    (by simpa using Data.additionalSpence97_strong)
    Data.additionalPaley97_isPaley (by norm_num) r
  simpa using h

theorem spence_family_109 (r : ℕ) : HadamardExists (138820*109^r) := by
  have h := spence_exists_from_t Data.t631 Data.t631_corr Data.t631_onehot
    Data.additionalSpence109 Data.additionalPaley109
    (by simpa using Data.additionalSpence109_strong)
    Data.additionalPaley109_isPaley (by norm_num) r
  simpa using h

theorem spence_family_113 (r : ℕ) : HadamardExists (143868*113^r) := by
  have h := spence_exists_from_t Data.t631 Data.t631_corr Data.t631_onehot
    Data.additionalSpence113 Data.additionalPaley113
    (by simpa using Data.additionalSpence113_strong)
    Data.additionalPaley113_isPaley (by norm_num) r
  simpa using h

theorem spence_family_137 (r : ℕ) : HadamardExists (174156*137^r) := by
  have h := spence_exists_from_t Data.t631 Data.t631_corr Data.t631_onehot
    Data.additionalSpence137 Data.additionalPaley137
    (by simpa using Data.additionalSpence137_strong)
    Data.additionalPaley137_isPaley (by norm_num) r
  simpa using h

#print axioms spence_family_5
#print axioms spence_family_53
#print axioms spence_family_97
#print axioms spence_family_109
#print axioms spence_family_113
#print axioms spence_family_137

end ButsonTowers

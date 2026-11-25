import ButsonTowers.GoethalsSeidelEntry
import ButsonTowers.OrthogonalSubstitution

namespace ButsonTowers.OrthogonalDesign

open scoped BigOperators
open GroupAlgebra GoethalsSeidel

variable {G : Type*} [AddCommGroup G] [Fintype G] [DecidableEq G]

def combine (x : Fin 4 → ℚ) (b : Fin 4 → G → ℚ) : Fin 4 → G → ℚ :=
  ![x 0 • b 0 - x 1 • b 1 - x 2 • b 2 - x 3 • b 3,
    x 1 • b 0 + x 0 • b 1 + x 3 • b 2 - x 2 • b 3,
    x 2 • b 0 - x 3 • b 1 + x 0 • b 2 + x 1 • b 3,
    x 3 • b 0 + x 2 • b 1 - x 1 • b 2 + x 0 • b 3]

def basis (i : Fin 4) : Fin 4 → ℚ := fun j => if j=i then 1 else 0

def coefficients (b : Fin 4 → G → ℚ) (i : Fin 4) :
    Matrix (Fin 4 × G) (Fin 4 × G) ℚ := gs (combine (basis i) b)

def combineSign (r c : Fin 4) : ℚ :=
  if (r=0 ∧ c≠0) ∨ (r=1 ∧ c=3) ∨ (r=2 ∧ c=1) ∨ (r=3 ∧ c=2)
    then -1 else 1

theorem combineSign_sign (r c : Fin 4) : combineSign r c = 1 ∨ combineSign r c = -1 := by
  unfold combineSign
  split_ifs
  · exact Or.inr rfl
  · exact Or.inl rfl

theorem combine_basis_apply (b : Fin 4 → G → ℚ) (r k : Fin 4) (u : G) :
    combine (basis k) b r u = combineSign r (channel r k) * b (channel r k) u := by
  fin_cases r <;> fin_cases k
  · change (1 : ℚ) * b 0 u - (0 : ℚ) * b 1 u - (0 : ℚ) * b 2 u - (0 : ℚ) * b 3 u = (1 : ℚ) * b 0 u
    simp only [zero_mul, one_mul, sub_zero]
  · change (0 : ℚ) * b 0 u - (1 : ℚ) * b 1 u - (0 : ℚ) * b 2 u - (0 : ℚ) * b 3 u = (-1 : ℚ) * b 1 u
    simp only [zero_mul, one_mul, sub_zero, zero_sub, neg_one_mul]
  · change (0 : ℚ) * b 0 u - (0 : ℚ) * b 1 u - (1 : ℚ) * b 2 u - (0 : ℚ) * b 3 u = (-1 : ℚ) * b 2 u
    simp only [zero_mul, one_mul, sub_zero, zero_sub, neg_one_mul]
  · change (0 : ℚ) * b 0 u - (0 : ℚ) * b 1 u - (0 : ℚ) * b 2 u - (1 : ℚ) * b 3 u = (-1 : ℚ) * b 3 u
    simp only [zero_mul, one_mul, sub_zero, zero_sub, neg_one_mul]
  · change (0 : ℚ) * b 0 u + (1 : ℚ) * b 1 u + (0 : ℚ) * b 2 u - (0 : ℚ) * b 3 u = (1 : ℚ) * b 1 u
    simp only [zero_mul, one_mul, sub_zero, zero_add, add_zero]
  · change (1 : ℚ) * b 0 u + (0 : ℚ) * b 1 u + (0 : ℚ) * b 2 u - (0 : ℚ) * b 3 u = (1 : ℚ) * b 0 u
    simp only [zero_mul, one_mul, sub_zero, add_zero]
  · change (0 : ℚ) * b 0 u + (0 : ℚ) * b 1 u + (0 : ℚ) * b 2 u - (1 : ℚ) * b 3 u = (-1 : ℚ) * b 3 u
    simp only [zero_mul, one_mul, zero_sub, add_zero, neg_one_mul]
  · change (0 : ℚ) * b 0 u + (0 : ℚ) * b 1 u + (1 : ℚ) * b 2 u - (0 : ℚ) * b 3 u = (1 : ℚ) * b 2 u
    simp only [zero_mul, one_mul, sub_zero, zero_add, add_zero]
  · change (0 : ℚ) * b 0 u - (0 : ℚ) * b 1 u + (1 : ℚ) * b 2 u + (0 : ℚ) * b 3 u = (1 : ℚ) * b 2 u
    simp only [zero_mul, one_mul, sub_zero, zero_add, add_zero]
  · change (0 : ℚ) * b 0 u - (0 : ℚ) * b 1 u + (0 : ℚ) * b 2 u + (1 : ℚ) * b 3 u = (1 : ℚ) * b 3 u
    simp only [zero_mul, one_mul, sub_zero, zero_add, add_zero]
  · change (1 : ℚ) * b 0 u - (0 : ℚ) * b 1 u + (0 : ℚ) * b 2 u + (0 : ℚ) * b 3 u = (1 : ℚ) * b 0 u
    simp only [zero_mul, one_mul, sub_zero, add_zero]
  · change (0 : ℚ) * b 0 u - (1 : ℚ) * b 1 u + (0 : ℚ) * b 2 u + (0 : ℚ) * b 3 u = (-1 : ℚ) * b 1 u
    simp only [zero_mul, one_mul, zero_sub, add_zero, neg_one_mul]
  · change (0 : ℚ) * b 0 u + (0 : ℚ) * b 1 u - (0 : ℚ) * b 2 u + (1 : ℚ) * b 3 u = (1 : ℚ) * b 3 u
    simp only [zero_mul, one_mul, sub_zero, zero_add, add_zero]
  · change (0 : ℚ) * b 0 u + (0 : ℚ) * b 1 u - (1 : ℚ) * b 2 u + (0 : ℚ) * b 3 u = (-1 : ℚ) * b 2 u
    simp only [zero_mul, one_mul, zero_sub, add_zero, neg_one_mul]
  · change (0 : ℚ) * b 0 u + (1 : ℚ) * b 1 u - (0 : ℚ) * b 2 u + (0 : ℚ) * b 3 u = (1 : ℚ) * b 1 u
    simp only [zero_mul, one_mul, sub_zero, zero_add, add_zero]
  · change (1 : ℚ) * b 0 u + (0 : ℚ) * b 1 u - (0 : ℚ) * b 2 u + (0 : ℚ) * b 3 u = (1 : ℚ) * b 0 u
    simp only [zero_mul, one_mul, sub_zero, add_zero]

theorem combine_basis_onehot (b : Fin 4 → G → ℚ)
    (hb : ∀ u, ∃ c, (b c u=1 ∨ b c u= -1) ∧ ∀ d, d≠c → b d u=0)
    (r : Fin 4) (u : G) :
    ∃ k, (combine (basis k) b r u=1 ∨ combine (basis k) b r u= -1) ∧
      ∀ l, l≠k → combine (basis l) b r u=0 := by
  obtain ⟨c,hc,hzero⟩ := hb u
  refine ⟨channel r c, ?_, ?_⟩
  · simp only [combine_basis_apply, channel_involutive r c]
    rcases combineSign_sign r c with hs | hs <;>
      rcases hc with hc | hc <;> simp [hs,hc]
  · intro l hl
    rw [combine_basis_apply]
    have hne : channel r l ≠ c := by
      intro he
      apply hl
      calc
        l = channel r (channel r l) := (channel_involutive r l).symm
        _ = channel r c := congrArg (channel r) he
    rw [hzero _ hne, mul_zero]

theorem coefficients_signed_partition (b : Fin 4 → G → ℚ)
    (hb : ∀ u, ∃ c, (b c u=1 ∨ b c u= -1) ∧ ∀ d, d≠c → b d u=0) :
    SignedPartition (coefficients b) := by
  intro ⟨r,u⟩ ⟨c,v⟩
  obtain ⟨k,hk,hzero⟩ := combine_basis_onehot b hb (channel r c) (entryArgument r c u v)
  refine ⟨k,?_,?_⟩
  · change gs (combine (basis k) b) (r,u) (c,v)=1 ∨
      gs (combine (basis k) b) (r,u) (c,v)= -1
    simp only [gs_apply]
    rcases entrySign_sign r c with hs | hs <;>
      rcases hk with hk | hk <;> simp [hs,hk]
  · intro l hl
    change gs (combine (basis l) b) (r,u) (c,v)=0
    rw [gs_apply, hzero l hl, mul_zero]

theorem combine_dot (x : Fin 4 → ℚ) (b : Fin 4 → G → ℚ) (u v : G) :
    (∑ r, combine x b r u * combine x b r v) =
      (∑ r, (x r)^2) * (∑ r, b r u * b r v) := by
  simp [combine, Fin.sum_univ_succ]
  ring

theorem combine_add (x y : Fin 4 → ℚ) (b : Fin 4 → G → ℚ) :
    combine (x+y) b = combine x b + combine y b := by
  ext i u
  fin_cases i <;> simp [combine] <;> ring

theorem combine_correlations (x : Fin 4 → ℚ) (b : Fin 4 → G → ℚ) (p : ℚ)
    (h : ∀ s, ∑ c, ∑ u, b c u * b c (u+s) = if s=0 then p else 0) :
    ∀ s, ∑ c, ∑ u, combine x b c u * combine x b c (u+s) =
      if s=0 then p*(∑ r, (x r)^2) else 0 := by
  intro s
  calc
    _ = ∑ u, (∑ r, (x r)^2) * (∑ c, b c u * b c (u+s)) := by
      rw [Finset.sum_comm]
      simp_rw [combine_dot]
    _ = (∑ r, (x r)^2) * (∑ c, ∑ u, b c u * b c (u+s)) := by
      rw [← Finset.mul_sum, Finset.sum_comm]
    _ = _ := by rw [h s]; split_ifs <;> ring

end ButsonTowers.OrthogonalDesign

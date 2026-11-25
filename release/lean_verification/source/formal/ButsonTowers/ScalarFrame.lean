import Mathlib.Algebra.Star.Rat
import Mathlib.LinearAlgebra.Matrix.Kronecker
import Mathlib.LinearAlgebra.Matrix.ConjTranspose
import Mathlib.Tactic.NoncommRing
import Mathlib.Tactic.SimpRw
import Mathlib.Tactic.Ext

/-!
Rational Euclidean frame transport.  The source may be a commutative
star algebra such as Q[G], and need not be a field or have root entries.
-/

namespace ButsonTowers

open scoped BigOperators

variable {A : Type*} [CommRing A] [StarRing A] [Algebra ℚ A]
variable {d : Type*} [Fintype d] [DecidableEq d]

theorem inverse_metric_adjoint
    (ρ : A →ₐ[ℚ] Matrix d d ℚ) (G J : Matrix d d ℚ)
    (hGJ : G * J = 1) (hJG : J * G = 1)
    (hadjoint : ∀ a, (ρ a).transpose * G = G * ρ (star a)) (a : A) :
    ρ (star a) * J = J * (ρ a).transpose := by
  calc
    ρ (star a) * J = (J * G) * ρ (star a) * J := by rw [hJG]; simp
    _ = J * (G * ρ (star a)) * J := by noncomm_ring
    _ = J * ((ρ a).transpose * G) * J := by rw [hadjoint]
    _ = J * (ρ a).transpose * (G * J) := by noncomm_ring
    _ = J * (ρ a).transpose := by rw [hGJ]; simp

omit [StarRing A] in
theorem framed_zero (ρ : A →ₐ[ℚ] Matrix d d ℚ) (F : Matrix d d ℚ) :
    ρ 0 * F = 0 := by simp

theorem frame_pair_transport
    (ρ : A →ₐ[ℚ] Matrix d d ℚ)
    (hstar : ∀ a, ρ (star a) = (ρ a).transpose)
    (F : Matrix d d ℚ) (c : ℚ)
    (hF : F * F.transpose = c • (1 : Matrix d d ℚ)) (a b : A) :
    (ρ a * F) * (ρ b * F).transpose = c • ρ (a * star b) := by
  rw [Matrix.transpose_mul]
  calc
    ρ a * F * (F.transpose * (ρ b).transpose) =
        ρ a * (F * F.transpose) * (ρ b).transpose := by noncomm_ring
    _ = ρ a * (c • (1 : Matrix d d ℚ)) * ρ (star b) := by rw [hF, hstar]
    _ = c • (ρ a * ρ (star b)) := by simp
    _ = c • ρ (a * star b) := by rw [map_mul]

theorem frame_sum_transport
    {k : Type*} [Fintype k]
    (ρ : A →ₐ[ℚ] Matrix d d ℚ)
    (hstar : ∀ a, ρ (star a) = (ρ a).transpose)
    (F : Matrix d d ℚ) (c : ℚ)
    (hF : F * F.transpose = c • (1 : Matrix d d ℚ)) (a b : k → A) :
    ∑ i, (ρ (a i) * F) * (ρ (b i) * F).transpose =
      c • ρ (∑ i, a i * star (b i)) := by
  simp_rw [frame_pair_transport ρ hstar F c hF]
  simp only [map_sum, Finset.smul_sum]

theorem weighted_pair_transport
    (ρ : A →ₐ[ℚ] Matrix d d ℚ) (J : Matrix d d ℚ)
    (hadjoint : ∀ a, ρ (star a) * J = J * (ρ a).transpose) (a b : A) :
    ρ a * J * (ρ b).transpose = ρ (a * star b) * J := by
  rw [map_mul]
  calc
    ρ a * J * (ρ b).transpose = ρ a * (J * (ρ b).transpose) := by noncomm_ring
    _ = ρ a * (ρ (star b) * J) := by rw [hadjoint]
    _ = ρ a * ρ (star b) * J := by noncomm_ring

theorem weighted_sum_transport
    {k : Type*} [Fintype k]
    (ρ : A →ₐ[ℚ] Matrix d d ℚ) (J : Matrix d d ℚ)
    (hadjoint : ∀ a, ρ (star a) * J = J * (ρ a).transpose) (a b : k → A) :
    ∑ i, ρ (a i) * J * (ρ (b i)).transpose =
      ρ (∑ i, a i * star (b i)) * J := by
  simp_rw [weighted_pair_transport ρ J hadjoint]
  simp only [map_sum, Finset.sum_mul]

theorem weighted_gram_block
    {k : Type*} [Fintype k]
    (ρ : A →ₐ[ℚ] Matrix d d ℚ) (J : Matrix d d ℚ)
    (hadjoint : ∀ a, ρ (star a) * J = J * (ρ a).transpose) (a b : k → A)
    (w : ℚ) (hsource : ∑ i, a i * star (b i) = algebraMap ℚ A w) :
    ∑ i, ρ (a i) * J * (ρ (b i)).transpose = w • J := by
  rw [weighted_sum_transport ρ J hadjoint, hsource, ρ.commutes]
  simp [Algebra.algebraMap_eq_smul_one]

theorem frame_gram_block
    {k : Type*} [Fintype k]
    (ρ : A →ₐ[ℚ] Matrix d d ℚ)
    (hstar : ∀ a, ρ (star a) = (ρ a).transpose)
    (F : Matrix d d ℚ) (c w : ℚ)
    (hF : F * F.transpose = c • (1 : Matrix d d ℚ)) (a b : k → A)
    (hsource : ∑ i, a i * star (b i) = algebraMap ℚ A w) :
    ∑ i, (ρ (a i) * F) * (ρ (b i) * F).transpose =
      (c * w) • (1 : Matrix d d ℚ) := by
  rw [frame_sum_transport ρ hstar F c hF, hsource, ρ.commutes]
  simp [Algebra.algebraMap_eq_smul_one, smul_smul]

def framedMatrix {n : Type*} (ρ : A →ₐ[ℚ] Matrix d d ℚ)
    (F : Matrix d d ℚ) (H : Matrix n n A) : Matrix (n × d) (n × d) ℚ :=
  fun i j => (ρ (H i.1 j.1) * F) i.2 j.2

def IsSignMatrix {n : Type*} (H : Matrix n n ℚ) : Prop :=
  ∀ i j, H i j = 1 ∨ H i j = -1

omit [StarRing A] in
theorem framed_matrix_flat {n : Type*}
    (ρ : A →ₐ[ℚ] Matrix d d ℚ) (F : Matrix d d ℚ) (H : Matrix n n A)
    (hentries : ∀ i j a b, (ρ (H i j) * F) a b = 1 ∨
      (ρ (H i j) * F) a b = -1) :
    IsSignMatrix (framedMatrix ρ F H) := by
  intro i j
  exact hentries i.1 j.1 i.2 j.2

theorem framed_matrix_gram
    {n : Type*} [Fintype n] [DecidableEq n]
    (ρ : A →ₐ[ℚ] Matrix d d ℚ)
    (hstar : ∀ a, ρ (star a) = (ρ a).transpose)
    (F : Matrix d d ℚ) (H : Matrix n n A) (c w : ℚ)
    (hF : F * F.transpose = c • (1 : Matrix d d ℚ))
    (hsource : ∀ i j, ∑ k, H i k * star (H j k) =
      algebraMap ℚ A (if i = j then w else 0)) :
    framedMatrix ρ F H * (framedMatrix ρ F H).transpose =
      (c * w) • (1 : Matrix (n × d) (n × d) ℚ) := by
  ext ⟨i,a⟩ ⟨j,b⟩
  have block := frame_gram_block ρ hstar F c (if i = j then w else 0)
    hF (fun k => H i k) (fun k => H j k) (hsource i j)
  have entry := congrArg (fun M : Matrix d d ℚ => M a b) block
  change (∑ k : n × d, (ρ (H i k.1) * F) a k.2 * (ρ (H j k.1) * F) b k.2) = _
  rw [Fintype.sum_prod_type]
  rw [Matrix.sum_apply] at entry
  change (∑ k, ∑ l, (ρ (H i k) * F) a l * (ρ (H j k) * F) b l) = _ at entry
  rw [entry]
  by_cases hij : i = j <;> by_cases hab : a = b <;>
    simp [hij, hab, Matrix.smul_apply]

def IsHadamard {n : Type*} [Fintype n] [DecidableEq n] (H : Matrix n n ℚ) : Prop :=
  0 < Fintype.card n ∧ IsSignMatrix H ∧
    H * H.transpose = (Fintype.card n : ℚ) • (1 : Matrix n n ℚ)

theorem framed_matrix_hadamard
    {n : Type*} [Fintype n] [DecidableEq n] [Nonempty n] [Nonempty d]
    (ρ : A →ₐ[ℚ] Matrix d d ℚ)
    (hstar : ∀ a, ρ (star a) = (ρ a).transpose)
    (F : Matrix d d ℚ) (H : Matrix n n A) (c w : ℚ)
    (hF : F * F.transpose = c • (1 : Matrix d d ℚ))
    (hsource : ∀ i j, ∑ k, H i k * star (H j k) =
      algebraMap ℚ A (if i = j then w else 0))
    (hentries : ∀ i j a b, (ρ (H i j) * F) a b = 1 ∨
      (ρ (H i j) * F) a b = -1)
    (hscale : c*w = (Fintype.card (n × d) : ℚ)) :
    IsHadamard (framedMatrix ρ F H) := by
  refine ⟨Fintype.card_pos, framed_matrix_flat ρ F H hentries, ?_⟩
  rw [framed_matrix_gram ρ hstar F H c w hF hsource, hscale]

end ButsonTowers

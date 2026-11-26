import ButsonTowers.ScalarFrame
import Mathlib.RingTheory.Trace.Defs
import Mathlib.LinearAlgebra.Matrix.ToLin
import Mathlib.Tactic.Ring

namespace ButsonTowers.TraceTowers

open Module
open scoped BigOperators

variable {A : Type*} [CommRing A] [StarRing A] [Algebra ℚ A] [StarModule ℚ A]
variable {ι κ : Type*} [Fintype ι] [DecidableEq ι] [Fintype κ] [DecidableEq κ]

noncomputable def traceGram (b : Basis ι ℚ A) (τ : A →ₗ[ℚ] ℚ) : Matrix ι ι ℚ :=
  fun i j => τ (star (b i) * b j)

omit [DecidableEq ι] in
theorem expand_left (b : Basis ι ℚ A) (τ : A →ₗ[ℚ] ℚ) (x y : A) :
    ∑ i, b.repr x i * τ (star (b i) * y) = τ (star x * y) := by
  calc
    _ = τ (star (∑ i, b.repr x i • b i) * y) := by
      simp only [star_sum, star_smul, star_trivial, Finset.sum_mul,
        smul_mul_assoc, map_sum, map_smul, smul_eq_mul]
    _ = _ := by rw [b.sum_repr]

omit [StarModule ℚ A] [DecidableEq ι] in
theorem expand_right (b : Basis ι ℚ A) (τ : A →ₗ[ℚ] ℚ) (x y : A) :
    ∑ i, τ (star x * b i) * b.repr y i = τ (star x * y) := by
  calc
    _ = τ (star x * (∑ i, b.repr y i • b i)) := by
      simp only [Finset.mul_sum, mul_smul_comm, map_sum, map_smul, smul_eq_mul]
      apply Finset.sum_congr rfl
      intro i _
      exact mul_comm _ _
    _ = _ := by rw [b.sum_repr]

theorem gram_mul_action (b : Basis ι ℚ A) (τ : A →ₗ[ℚ] ℚ) (a : A) (i j : ι) :
    (traceGram b τ * Algebra.leftMulMatrix b a) i j = τ (star (b i) * (a * b j)) := by
  simpa [Matrix.mul_apply, traceGram, Algebra.leftMulMatrix_eq_repr_mul] using
    expand_right b τ (b i) (a * b j)

theorem action_transpose_mul_gram (b : Basis ι ℚ A) (τ : A →ₗ[ℚ] ℚ) (a : A) (i j : ι) :
    ((Algebra.leftMulMatrix b a).transpose * traceGram b τ) i j =
      τ (star (a * b i) * b j) := by
  simpa [Matrix.mul_apply, traceGram, Algebra.leftMulMatrix_eq_repr_mul] using
    expand_left b τ (a * b i) (b j)

theorem trace_metric_adjoint (b : Basis ι ℚ A) (τ : A →ₗ[ℚ] ℚ) (a : A) :
    (Algebra.leftMulMatrix b a).transpose * traceGram b τ =
      traceGram b τ * Algebra.leftMulMatrix b (star a) := by
  ext i j
  rw [action_transpose_mul_gram, gram_mul_action, star_mul, mul_assoc]

omit [StarModule ℚ A] [Fintype ι] [DecidableEq ι] in
theorem trace_gram_symmetric (b : Basis ι ℚ A) (τ : A →ₗ[ℚ] ℚ)
    (hτ : ∀ x, τ (star x) = τ x) : (traceGram b τ).transpose = traceGram b τ := by
  ext i j
  change τ (star (b j) * b i) = τ (star (b i) * b j)
  simpa [star_mul] using hτ (star (b i) * b j)

theorem inverse_trace_metric_adjoint (b : Basis ι ℚ A) (τ : A →ₗ[ℚ] ℚ)
    (J : Matrix ι ι ℚ) (hGJ : traceGram b τ * J = 1) (hJG : J * traceGram b τ = 1)
    (a : A) : Algebra.leftMulMatrix b (star a) * J =
      J * (Algebra.leftMulMatrix b a).transpose :=
  ButsonTowers.inverse_metric_adjoint (Algebra.leftMulMatrix b) (traceGram b τ) J
    hGJ hJG (trace_metric_adjoint b τ) a

theorem weighted_trace_gram_transport
    {n : Type*} [Fintype n]
    (b : Basis ι ℚ A) (τ : A →ₗ[ℚ] ℚ)
    (J : Matrix ι ι ℚ) (hGJ : traceGram b τ * J = 1) (hJG : J * traceGram b τ = 1)
    (x y : n → A) (w : ℚ)
    (hsource : ∑ i, x i * star (y i) = algebraMap ℚ A w) :
    ∑ i, Algebra.leftMulMatrix b (x i) * J *
      (Algebra.leftMulMatrix b (y i)).transpose = w • J :=
  ButsonTowers.weighted_gram_block (Algebra.leftMulMatrix b) J
    (inverse_trace_metric_adjoint b τ J hGJ hJG) x y w hsource

noncomputable def coordinateMatrix (b : Basis ι ℚ A) (v : κ → A) : Matrix ι κ ℚ :=
  fun i j => b.repr (v j) i

omit [Fintype κ] [DecidableEq κ] in
theorem gram_congruence (b : Basis ι ℚ A) (τ : A →ₗ[ℚ] ℚ) (v : κ → A) :
    (fun i j => τ (star (v i) * v j)) =
      (coordinateMatrix b v).transpose * traceGram b τ * coordinateMatrix b v := by
  have right (i : ι) (j : κ) :
      (traceGram b τ * coordinateMatrix b v) i j = τ (star (b i) * v j) := by
    simpa [Matrix.mul_apply, traceGram, coordinateMatrix] using expand_right b τ (b i) (v j)
  rw [Matrix.mul_assoc]
  ext i j
  symm
  change (∑ k, b.repr (v i) k * (traceGram b τ * coordinateMatrix b v) k j) = _
  simp_rw [right]
  exact expand_left b τ (v i) (v j)

theorem basis_gram_congruence (b : Basis ι ℚ A) (c : Basis κ ℚ A) (τ : A →ₗ[ℚ] ℚ) :
    traceGram c τ = (coordinateMatrix b c).transpose * traceGram b τ * coordinateMatrix b c :=
  gram_congruence b τ c

section Towers

variable {K B : Type*} [CommRing K] [CommRing B]
  [Algebra ℚ K] [Algebra K B] [Algebra ℚ B] [IsScalarTower ℚ K B]

omit [DecidableEq ι] [DecidableEq κ] in
theorem finite_free_trace_transitivity (b : Basis ι ℚ K) (c : Basis κ K B) (x : B) :
    Algebra.trace ℚ K (Algebra.trace K B x) = Algebra.trace ℚ B x :=
  Algebra.trace_trace_of_basis b c x

theorem tower_action_coordinates (b : Basis ι ℚ K) (c : Basis κ K B)
    (x : B) (ik jk : ι × κ) :
    Algebra.leftMulMatrix (b.smulTower c) x ik jk =
      Algebra.leftMulMatrix b (Algebra.leftMulMatrix c x ik.2 jk.2) ik.1 jk.1 :=
  Algebra.smulTower_leftMulMatrix b c x ik jk

variable [StarRing K] [StarRing B] [StarModule ℚ K] [StarModule K B]

omit [DecidableEq κ] in
theorem tower_trace_metric (b : Basis ι ℚ K) (c : Basis κ K B)
    (i j : ι) (k l : κ) :
    traceGram (b.smulTower c) (Algebra.trace ℚ B) (i,k) (j,l) =
      (traceGram b (Algebra.trace ℚ K) *
        Algebra.leftMulMatrix b (Algebra.trace K B (star (c k) * c l))) i j := by
  rw [gram_mul_action]
  change Algebra.trace ℚ B
    (star ((b.smulTower c) (i,k)) * (b.smulTower c) (j,l)) = _
  rw [← Algebra.trace_trace_of_basis b c]
  congr 1
  simp only [Basis.smulTower_apply, star_smul, smul_mul_assoc, mul_smul_comm,
    map_smul, smul_eq_mul]
  ring

end Towers

end ButsonTowers.TraceTowers

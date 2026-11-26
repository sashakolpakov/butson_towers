import ButsonTowers.OrthogonalSubstitution
import ButsonTowers.SpenceMatrix

namespace ButsonTowers.Unreal
open scoped BigOperators Kronecker

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

def J : Matrix (Fin 3) (Fin 3) ℚ := fun _ _ => 1
def C : Matrix (Fin 3) (Fin 3) ℚ := !![0,1,-1; -1,0,1; 1,-1,0]

theorem transpose_J : J.transpose = J := by rfl
theorem transpose_C : C.transpose = -C := by decide
theorem JJ : J*J = (3:ℚ) • J := by
  ext i j; fin_cases i <;> fin_cases j <;> norm_num [J,Matrix.mul_apply,Fin.sum_univ_succ]
theorem CC : C*C.transpose = (3:ℚ) • (1 : Matrix (Fin 3) (Fin 3) ℚ) - J := by
  ext i j; fin_cases i <;> fin_cases j <;> norm_num [C,J,Matrix.mul_apply,Fin.sum_univ_succ]
theorem JC : J*C = 0 := by
  ext i j; fin_cases i <;> fin_cases j <;> norm_num [C,J,Matrix.mul_apply,Fin.sum_univ_succ]
theorem CJ : C*J = 0 := by
  ext i j; fin_cases i <;> fin_cases j <;> norm_num [C,J,Matrix.mul_apply,Fin.sum_univ_succ]

def nextA (_A B : Matrix ι ι ℚ) := B ⊗ₖ J
def nextB (A B : Matrix ι ι ℚ) := A ⊗ₖ (1 : Matrix (Fin 3) (Fin 3) ℚ) + B ⊗ₖ C

omit [DecidableEq ι] in
theorem next_amicable (A B : Matrix ι ι ℚ)
    (hab : A*B.transpose = B*A.transpose) :
    nextA A B * (nextB A B).transpose = nextB A B * (nextA A B).transpose := by
  simp only [nextA,nextB,Matrix.transpose_add,transpose_tensor,Matrix.transpose_one,
    transpose_J,Matrix.mul_add,Matrix.add_mul,← Matrix.mul_kronecker_mul,
    Matrix.mul_one,Matrix.one_mul]
  have hjc : J*C.transpose = 0 := by rw [transpose_C]; simp [JC]
  rw [hjc,CJ,hab]

theorem next_norm (A B : Matrix ι ι ℚ) (n : ℚ)
    (hab : A*B.transpose = B*A.transpose)
    (hnorm : A*A.transpose + (3:ℚ) • (B*B.transpose) = (4*n) • (1 : Matrix ι ι ℚ)) :
    nextA A B * (nextA A B).transpose +
      (3:ℚ) • (nextB A B * (nextB A B).transpose) =
      (4*(3*n)) • (1 : Matrix (ι × Fin 3) (ι × Fin 3) ℚ) := by
  have hcross : (A*B.transpose) ⊗ₖ C.transpose + (B*A.transpose) ⊗ₖ C = 0 := by
    rw [transpose_C,hab]
    ext ⟨i,a⟩ ⟨j,b⟩
    simp
  have hB : nextB A B * (nextB A B).transpose =
      (A*A.transpose) ⊗ₖ (1 : Matrix (Fin 3) (Fin 3) ℚ) +
      (B*B.transpose) ⊗ₖ ((3:ℚ) • (1 : Matrix (Fin 3) (Fin 3) ℚ) - J) := by
    simp only [nextB,Matrix.transpose_add,transpose_tensor,Matrix.transpose_one,
      Matrix.add_mul,Matrix.mul_add,← Matrix.mul_kronecker_mul,Matrix.mul_one,Matrix.one_mul,CC]
    calc
      _ = (A*A.transpose) ⊗ₖ (1 : Matrix (Fin 3) (Fin 3) ℚ) +
          (B*B.transpose) ⊗ₖ ((3:ℚ) • (1 : Matrix (Fin 3) (Fin 3) ℚ) - J) +
          ((A*B.transpose) ⊗ₖ C.transpose + (B*A.transpose) ⊗ₖ C) := by abel
      _ = _ := by rw [hcross]; simp
  have hA : nextA A B * (nextA A B).transpose =
      (B*B.transpose) ⊗ₖ ((3:ℚ) • J) := by
    simp only [nextA,transpose_tensor,transpose_J,← Matrix.mul_kronecker_mul,JJ]
  rw [hA,hB]
  have ht := congrArg (fun M : Matrix ι ι ℚ =>
    M ⊗ₖ (1 : Matrix (Fin 3) (Fin 3) ℚ)) hnorm
  simp only [Matrix.add_kronecker, Matrix.smul_kronecker,Matrix.one_kronecker_one] at ht
  ext ⟨i,a⟩ ⟨j,b⟩
  have he := congrArg (fun M : Matrix (ι × Fin 3) (ι × Fin 3) ℚ => M (i,a) (j,b)) ht
  simp only [Matrix.add_apply,Matrix.smul_apply,smul_eq_mul,Matrix.kroneckerMap_apply] at he ⊢
  simp only [Matrix.sub_apply, Matrix.smul_apply, smul_eq_mul]
  linear_combination 3*he

omit [Fintype ι] [DecidableEq ι] in
theorem nextA_sign (A B : Matrix ι ι ℚ) (hB : IsSignMatrix B) :
    IsSignMatrix (nextA A B) := by
  intro ⟨i,a⟩ ⟨j,b⟩
  simpa [nextA,J] using hB i j

omit [Fintype ι] [DecidableEq ι] in
theorem nextB_sign (A B : Matrix ι ι ℚ)
    (hA : IsSignMatrix A) (hB : IsSignMatrix B) :
    IsSignMatrix (nextB A B) := by
  intro ⟨i,a⟩ ⟨j,b⟩
  have ha := hA i j
  have hb := hB i j
  fin_cases a <;> fin_cases b <;> rcases ha with ha | ha <;>
    rcases hb with hb | hb <;> simp [nextB,C,ha,hb]

structure Pair (ι : Type) [Fintype ι] [DecidableEq ι] (n : ℚ) where
  A : Matrix ι ι ℚ
  B : Matrix ι ι ℚ
  signA : IsSignMatrix A
  signB : IsSignMatrix B
  amicable : A*B.transpose = B*A.transpose
  norm : A*A.transpose + (3:ℚ) • (B*B.transpose) = (4*n) • 1

def tripled {ι : Type} [Fintype ι] [DecidableEq ι] {n : ℚ}
    (P : Pair ι n) : Pair (ι × Fin 3) (3*n) where
  A := nextA P.A P.B
  B := nextB P.A P.B
  signA := nextA_sign P.A P.B P.signB
  signB := nextB_sign P.A P.B P.signA P.signB
  amicable := next_amicable P.A P.B P.amicable
  norm := next_norm P.A P.B n P.amicable P.norm

def tower {ι : Type} [Fintype ι] [DecidableEq ι] {n : ℚ}
    (P : Pair ι n) : (r : ℕ) → Pair (SpenceMatrix.FamilyIndex ι (Fin 3) r) (n*3^r)
  | 0 => by
    have hF : SpenceMatrix.familyFintype ι (Fin 3) 0 = (inferInstance : Fintype ι) :=
      Subsingleton.elim _ _
    have hD : SpenceMatrix.familyDecidableEq ι (Fin 3) 0 = (inferInstance : DecidableEq ι) :=
      Subsingleton.elim _ _
    simp only [pow_zero, mul_one]
    change @Pair ι (SpenceMatrix.familyFintype ι (Fin 3) 0)
      (SpenceMatrix.familyDecidableEq ι (Fin 3) 0) n
    rw [hF,hD]
    exact P
  | r+1 => by
    have h := tripled (tower P r)
    have he : 3*(n*3^r) = n*3^(r+1) := by rw [pow_succ]; ring
    rw [he] at h
    exact h

def quad {ι : Type} [Fintype ι] [DecidableEq ι] {n : ℚ}
    (P : Pair ι n) : Fin 4 → Matrix ι ι ℚ := ![P.A,P.B,P.B,P.B]

theorem quad_amicable {ι : Type} [Fintype ι] [DecidableEq ι] {n : ℚ}
    (P : Pair ι n) : AmicableQuad (quad P) (4*n) := by
  refine ⟨?_, ?_, ?_⟩
  · intro i
    fin_cases i <;> first | exact P.signA | exact P.signB
  · intro i j
    fin_cases i <;> fin_cases j <;> first | rfl | exact P.amicable | exact P.amicable.symm
  · have hthree : (3:ℚ) • (P.B*P.B.transpose) =
        P.B*P.B.transpose + (P.B*P.B.transpose + P.B*P.B.transpose) := by
      rw [show (3:ℚ)=1+1+1 by norm_num, add_smul, add_smul, one_smul]
      abel
    simpa [quad, Fin.sum_univ_succ, hthree] using P.norm

end ButsonTowers.Unreal

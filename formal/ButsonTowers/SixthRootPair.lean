import Mathlib.LinearAlgebra.Matrix.NonsingularInverse
import Mathlib.Tactic.NoncommRing
import Mathlib.Tactic.Module

/-!
The rational sixth-root matrix action in the manuscript.  The assumptions
are exactly the two equal row norms and the symmetric mixed product of
the matrices P and Q, with nonzero common scale d.  This module does not
identify a saved packed matrix or introduce a cyclotomic field.
-/

namespace ButsonTowers.SixthRootPair

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

def action (P Q : Matrix ι ι ℚ) (d : ℚ) : Matrix ι ι ℚ :=
  d⁻¹ • (Q * P.transpose)

theorem transpose_mul_of_mul_transpose
    (P : Matrix ι ι ℚ) (d : ℚ) (hd : d ≠ 0)
    (hP : P * P.transpose = d • (1 : Matrix ι ι ℚ)) :
    P.transpose * P = d • (1 : Matrix ι ι ℚ) := by
  have hright : P * (d⁻¹ • P.transpose) = 1 := by
    rw [Matrix.mul_smul, hP, smul_smul, inv_mul_cancel₀ hd, one_smul]
  have hleft := mul_eq_one_comm.mp hright
  have h := congrArg (fun M : Matrix ι ι ℚ => d • M) hleft
  simpa [Matrix.smul_mul, smul_smul, hd] using h

theorem action_mul_transpose
    (P Q : Matrix ι ι ℚ) (d : ℚ) (hd : d ≠ 0)
    (hP : P * P.transpose = d • (1 : Matrix ι ι ℚ))
    (hQ : Q * Q.transpose = d • (1 : Matrix ι ι ℚ)) :
    action P Q d * (action P Q d).transpose = 1 := by
  have hc := transpose_mul_of_mul_transpose P d hd hP
  have hm : (Q * P.transpose) * (P * Q.transpose) =
      (d*d) • (1 : Matrix ι ι ℚ) := by
    calc
      _ = Q * (P.transpose * P) * Q.transpose := by noncomm_ring
      _ = Q * (d • (1 : Matrix ι ι ℚ)) * Q.transpose := by rw [hc]
      _ = (d*d) • (1 : Matrix ι ι ℚ) := by
        rw [Matrix.mul_smul, Matrix.mul_one, Matrix.smul_mul, hQ, smul_smul]
  simp only [action, Matrix.transpose_smul, Matrix.transpose_mul,
    Matrix.transpose_transpose, Matrix.smul_mul, Matrix.mul_smul, smul_smul]
  rw [hm, smul_smul]
  have hs : d⁻¹ * d⁻¹ * (d*d) = (1 : ℚ) := by field_simp
  rw [hs, one_smul]

theorem action_add_transpose
    (P Q : Matrix ι ι ℚ) (d : ℚ) (hd : d ≠ 0)
    (hPQ : P * Q.transpose + Q * P.transpose = d • (1 : Matrix ι ι ℚ)) :
    action P Q d + (action P Q d).transpose = 1 := by
  have hm : Q * P.transpose + P * Q.transpose = d • (1 : Matrix ι ι ℚ) := by
    simpa only [add_comm] using hPQ
  simp only [action, Matrix.transpose_smul, Matrix.transpose_mul,
    Matrix.transpose_transpose]
  rw [← smul_add, hm, smul_smul, inv_mul_cancel₀ hd, one_smul]

theorem action_transpose
    (P Q : Matrix ι ι ℚ) (d : ℚ) (hd : d ≠ 0)
    (hPQ : P * Q.transpose + Q * P.transpose = d • (1 : Matrix ι ι ℚ)) :
    (action P Q d).transpose = 1 - action P Q d := by
  have h := action_add_transpose P Q d hd hPQ
  calc
    _ = (action P Q d + (action P Q d).transpose) - action P Q d := by abel
    _ = _ := by rw [h]

theorem action_polynomial
    (P Q : Matrix ι ι ℚ) (d : ℚ) (hd : d ≠ 0)
    (hP : P * P.transpose = d • (1 : Matrix ι ι ℚ))
    (hQ : Q * Q.transpose = d • (1 : Matrix ι ι ℚ))
    (hPQ : P * Q.transpose + Q * P.transpose = d • (1 : Matrix ι ι ℚ)) :
    action P Q d * action P Q d - action P Q d + 1 = 0 := by
  have hn := action_mul_transpose P Q d hd hP hQ
  have ht := action_transpose P Q d hd hPQ
  calc
    _ = action P Q d * action P Q d - action P Q d +
        action P Q d * (action P Q d).transpose := by rw [hn]
    _ = 0 := by rw [ht]; noncomm_ring

theorem action_mul_left
    (P Q : Matrix ι ι ℚ) (d : ℚ) (hd : d ≠ 0)
    (hP : P * P.transpose = d • (1 : Matrix ι ι ℚ)) :
    action P Q d * P = Q := by
  have hc := transpose_mul_of_mul_transpose P d hd hP
  calc
    _ = d⁻¹ • (Q * (P.transpose * P)) := by
      simp only [action, Matrix.smul_mul, Matrix.mul_assoc]
    _ = Q := by
      rw [hc, Matrix.mul_smul, Matrix.mul_one, smul_smul,
        inv_mul_cancel₀ hd, one_smul]

def imaginaryAction (P Q : Matrix ι ι ℚ) (d : ℚ) : Matrix ι ι ℚ :=
  (2 : ℚ) • action P Q d - 1

theorem imaginaryAction_transpose
    (P Q : Matrix ι ι ℚ) (d : ℚ) (hd : d ≠ 0)
    (hPQ : P * Q.transpose + Q * P.transpose = d • (1 : Matrix ι ι ℚ)) :
    (imaginaryAction P Q d).transpose = -imaginaryAction P Q d := by
  simp only [imaginaryAction, Matrix.transpose_sub, Matrix.transpose_smul,
    Matrix.transpose_one, action_transpose P Q d hd hPQ]
  module

theorem imaginaryAction_square
    (P Q : Matrix ι ι ℚ) (d : ℚ) (hd : d ≠ 0)
    (hP : P * P.transpose = d • (1 : Matrix ι ι ℚ))
    (hQ : Q * Q.transpose = d • (1 : Matrix ι ι ℚ))
    (hPQ : P * Q.transpose + Q * P.transpose = d • (1 : Matrix ι ι ℚ)) :
    imaginaryAction P Q d * imaginaryAction P Q d =
      (-3 : ℚ) • (1 : Matrix ι ι ℚ) := by
  have hp := action_polynomial P Q d hd hP hQ hPQ
  have hs : action P Q d * action P Q d = action P Q d - 1 := by
    calc
      _ = (action P Q d * action P Q d - action P Q d + 1) +
          (action P Q d - 1) := by abel
      _ = _ := by rw [hp]; simp
  simp only [imaginaryAction, sub_mul, mul_sub, Matrix.smul_mul,
    Matrix.mul_smul, Matrix.one_mul, Matrix.mul_one, smul_smul, hs]
  module

/-- Coordinates a+b*zeta, where zeta^2-zeta+1=0. -/
def coordinateAction (R : Matrix ι ι ℚ) (a b : ℚ) : Matrix ι ι ℚ :=
  a • 1 + b • R

theorem coordinateAction_mul
    (R : Matrix ι ι ℚ) (hR : R*R-R+1=0) (a b c e : ℚ) :
    coordinateAction R a b * coordinateAction R c e =
      coordinateAction R (a*c-b*e) (a*e+b*c+b*e) := by
  have hs : R*R=R-1 := by
    calc
      _ = (R*R-R+1)+(R-1) := by abel
      _ = _ := by rw [hR]; simp
  simp only [coordinateAction, add_mul, mul_add, Matrix.smul_mul,
    Matrix.mul_smul, Matrix.one_mul, Matrix.mul_one, smul_smul, hs]
  module

omit [Fintype ι] in
theorem coordinateAction_conjugate
    (R : Matrix ι ι ℚ) (hR : R.transpose=1-R) (a b : ℚ) :
    coordinateAction R (a+b) (-b) = (coordinateAction R a b).transpose := by
  simp only [coordinateAction, Matrix.transpose_add, Matrix.transpose_smul,
    Matrix.transpose_one, hR]
  module

end ButsonTowers.SixthRootPair

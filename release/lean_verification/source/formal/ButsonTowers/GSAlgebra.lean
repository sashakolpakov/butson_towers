import Mathlib.LinearAlgebra.Matrix.ConjTranspose
import Mathlib.LinearAlgebra.Matrix.Notation
import Mathlib.Tactic.FinCases
import Mathlib.Tactic.NoncommRing
import Mathlib.Tactic.Abel

namespace ButsonTowers.GSAlgebra
open scoped BigOperators

variable {R : Type*} [Ring R] [StarRing R]

def blocks (a : Fin 4 → R) (J : R) : Matrix (Fin 4) (Fin 4) R :=
  !![a 0, a 1*J, a 2*J, a 3*J;
    -(a 1*J), a 0, -(star (a 3)*J), star (a 2)*J;
    -(a 2*J), star (a 3)*J, a 0, -(star (a 1)*J);
    -(a 3*J), -(star (a 2)*J), star (a 1)*J, a 0]

theorem blocks_gram (a : Fin 4 → R) (J : R)
    (hJ : star J=J) (hJJ : J*J=1)
    (hJa : ∀ i, J*star (a i)=a i*J)
    (hcomm : ∀ i j, a i*a j=a j*a i)
    (hstarcomm : ∀ i j, star (a i)*a j=a j*star (a i)) :
    blocks a J * (blocks a J).conjTranspose =
      Matrix.diagonal (fun _ : Fin 4 => ∑ i, a i*star (a i)) := by
  have hJa' (i : Fin 4) : J*a i=star (a i)*J := by
    have h := congrArg (fun x : R => J*x*J) (hJa i)
    have hl : J*(J*star (a i))*J=star (a i)*J := by
      rw [← mul_assoc J J, hJJ, one_mul]
    have hr : J*(a i*J)*J=J*a i := by
      rw [mul_assoc, mul_assoc, hJJ, mul_one]
    rw [hl,hr] at h
    exact h.symm
  have hstarstar (i j : Fin 4) : star (a i)*star (a j)=star (a j)*star (a i) := by
    simpa only [star_mul] using congrArg star (hcomm j i)
  have hjjx (X : R) : J*(J*X)=X := by rw [← mul_assoc,hJJ,one_mul]
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [blocks, Matrix.mul_apply, Matrix.conjTranspose_apply,
      Fin.sum_univ_succ, star_mul, hJ]
  all_goals noncomm_ring
  all_goals try simp only [hjjx]
  all_goals try simp only [hJa,hJa']
  all_goals try simp only [← mul_assoc,hcomm,hstarcomm,hstarstar]
  all_goals abel

end ButsonTowers.GSAlgebra

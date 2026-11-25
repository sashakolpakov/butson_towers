import ButsonTowers.GoethalsSeidelCore
import ButsonTowers.GSAlgebra

namespace ButsonTowers.GoethalsSeidel
open scoped BigOperators
open GroupAlgebra
variable {G : Type*} [AddCommGroup G] [Fintype G] [DecidableEq G]

theorem blocks_gram (a : Fin 4 → G → ℚ) :
    blocks a * blockAdjoint (blocks a) =
      Matrix.diagonal (fun _ : Fin 4 => ∑ c : Fin 4,
        regular (a c) * (regular (a c)).transpose) := by
  have hs (M : Matrix G G ℚ) : star M = M.transpose := by
    rw [Matrix.star_eq_conjTranspose, Matrix.conjTranspose_eq_transpose_of_trivial]
  have hb : blocks a = GSAlgebra.blocks (fun i => regular (a i))
      (reflection : Matrix G G ℚ) := by
    simp [blocks, GSAlgebra.blocks, hs]
  have ht (B : Matrix (Fin 4) (Fin 4) (Matrix G G ℚ)) :
      blockAdjoint B = B.conjTranspose := by
    apply Matrix.ext
    intro i j
    simp [blockAdjoint, Matrix.conjTranspose_apply, hs]
  have h := GSAlgebra.blocks_gram (fun i => regular (a i))
    (reflection : Matrix G G ℚ)
    (by rw [hs, reflection_transpose]) reflection_sq
    (by
      intro i
      rw [hs, ← regular_involution, reflection_regular, involution_involution])
    (by intro i j; exact regular_mul_comm _ _)
    (by
      intro i j
      rw [hs, ← regular_involution]
      exact regular_mul_comm _ _)
  simpa only [← hb, ← ht, hs] using h

theorem gs_gram (a : Fin 4 → G → ℚ) (p : ℚ)
    (h : (∑ c, regular (a c) * (regular (a c)).transpose) =
      p • (1 : Matrix G G ℚ)) :
    gs a * (gs a).transpose = p • (1 : Matrix (Fin 4 × G) (Fin 4 × G) ℚ) := by
  change flatten (blocks a) * (flatten (blocks a)).transpose = _
  rw [← flatten_adjoint, ← flatten_mul, blocks_gram, h]
  ext ⟨i,u⟩ ⟨j,v⟩
  by_cases hij : i = j <;> by_cases huv : u = v <;>
    simp [flatten, Matrix.smul_apply, hij, huv]


end ButsonTowers.GoethalsSeidel

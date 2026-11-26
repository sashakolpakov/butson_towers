import ButsonTowers.Data.StrongChecks
import ButsonTowers.Data.PaleyChecks
import ButsonTowers.Data.LiteralTree

/-! Additional strong starter for q=113; width 57.
Channel order 0,1,2,3; channel signs (1, 1, 1, 1).
The literal Paley row is also checked in full: no primality oracle is a proof premise. -/

namespace ButsonTowers.Data
set_option maxRecDepth 1000000
set_option maxHeartbeats 0

def additionalSpence113Trees : Fin 4 → LiteralTree ℤ :=
  ![(.node 28 (.node 14 (.node 7 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (1)))))) (.node 7 (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (1))))))) (.node 14 (.node 7 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (1)))))) (.node 7 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1)))))))),
    (.node 28 (.node 14 (.node 7 (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (1)))))) (.node 7 (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (1))))))) (.node 14 (.node 7 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (1)))))) (.node 7 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1)))))))),
    (.node 28 (.node 14 (.node 7 (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (-1)))))) (.node 7 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1))))))) (.node 14 (.node 7 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (1)))))) (.node 7 (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (-1)))))))),
    (.node 28 (.node 14 (.node 7 (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (-1)))))) (.node 7 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1))))))) (.node 14 (.node 7 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (1)))))) (.node 7 (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (-1))))))))]

def additionalSpence113Row (i : Fin 4) (x : ℕ) : ℤ :=
  treeAt (additionalSpence113Trees i) x

theorem additionalSpence113_sign_checked :
    checkQuadSign 57 additionalSpence113Row = true := by decide +kernel

theorem additionalSpence113_symmetric_checked :
    checkQuadSymmetric 57 additionalSpence113Row = true := by decide +kernel

theorem additionalSpence113_minus_checked :
    checkQuadCorr 57 (pairRows additionalSpence113Row false) 4 = true := by decide +kernel

theorem additionalSpence113_plus_checked :
    checkQuadCorr 57 (pairRows additionalSpence113Row true) (4*(113 : ℤ)) = true := by decide +kernel

def additionalSpence113 : Fin 4 → Matrix (ZMod 57) (ZMod 57) ℚ :=
  rowMatrices 57 additionalSpence113Row

theorem additionalSpence113_strong : SpenceMatrix.Strong 113 0 additionalSpence113 :=
  quad_matrices_strong additionalSpence113_sign_checked additionalSpence113_symmetric_checked
    additionalSpence113_minus_checked additionalSpence113_plus_checked

def additionalPaley113Tree : LiteralTree ℤ :=
  (.node 56 (.node 28 (.node 14 (.node 7 (.node 3 (.node 1 (.leaf (0)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1)))))) (.node 7 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (-1))))))) (.node 14 (.node 7 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (1)))))) (.node 7 (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (-1)))))))) (.node 28 (.node 14 (.node 7 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (1)))))) (.node 7 (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (1))))))) (.node 14 (.node 7 (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1)))))) (.node 7 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (1)))))))))

def additionalPaley113Row : ℕ → ℤ := treeAt additionalPaley113Tree

theorem additionalPaley113_sign_checked :
    checkPaleySign 113 additionalPaley113Row = true := by decide +kernel
theorem additionalPaley113_symmetric_checked :
    checkPaleySymmetric 113 additionalPaley113Row = true := by decide +kernel
theorem additionalPaley113_sum_checked :
    natSum 113 additionalPaley113Row = 0 := by decide +kernel
theorem additionalPaley113_corr_checked :
    checkPaleyCorr 113 additionalPaley113Row = true := by decide +kernel

def additionalPaley113 : Matrix (ZMod 113) (ZMod 113) ℚ :=
  GroupAlgebra.regular (fun x => (additionalPaley113Row x.val : ℚ))

theorem additionalPaley113_isPaley : SpenceMatrix.IsPaley additionalPaley113 :=
  prime_paley_checked additionalPaley113Row additionalPaley113_sign_checked
    additionalPaley113_symmetric_checked additionalPaley113_sum_checked
    additionalPaley113_corr_checked

#print axioms additionalSpence113_strong
#print axioms additionalPaley113_isPaley

end ButsonTowers.Data

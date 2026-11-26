import ButsonTowers.Data.StrongChecks
import ButsonTowers.Data.PaleyChecks
import ButsonTowers.Data.LiteralTree

/-! Additional strong starter for q=109; width 55.
Channel order 0,1,2,3; channel signs (1, 1, 1, 1).
The literal Paley row is also checked in full: no primality oracle is a proof premise. -/

namespace ButsonTowers.Data
set_option maxRecDepth 1000000
set_option maxHeartbeats 0

def additionalSpence109Trees : Fin 4 → LiteralTree ℤ :=
  ![(.node 27 (.node 13 (.node 6 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (-1)))))) (.node 7 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (1))))))) (.node 14 (.node 7 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (1)))))) (.node 7 (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (-1)))))))),
    (.node 27 (.node 13 (.node 6 (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (-1)))))) (.node 7 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (1))))))) (.node 14 (.node 7 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (1)))))) (.node 7 (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (-1)))))))),
    (.node 27 (.node 13 (.node 6 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (-1)))))) (.node 7 (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (1))))))) (.node 14 (.node 7 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (1)))))) (.node 7 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (1)))))))),
    (.node 27 (.node 13 (.node 6 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (-1)))))) (.node 7 (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (1))))))) (.node 14 (.node 7 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (1)))))) (.node 7 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (1))))))))]

def additionalSpence109Row (i : Fin 4) (x : ℕ) : ℤ :=
  treeAt (additionalSpence109Trees i) x

theorem additionalSpence109_sign_checked :
    checkQuadSign 55 additionalSpence109Row = true := by decide +kernel

theorem additionalSpence109_symmetric_checked :
    checkQuadSymmetric 55 additionalSpence109Row = true := by decide +kernel

theorem additionalSpence109_minus_checked :
    checkQuadCorr 55 (pairRows additionalSpence109Row false) 4 = true := by decide +kernel

theorem additionalSpence109_plus_checked :
    checkQuadCorr 55 (pairRows additionalSpence109Row true) (4*(109 : ℤ)) = true := by decide +kernel

def additionalSpence109 : Fin 4 → Matrix (ZMod 55) (ZMod 55) ℚ :=
  rowMatrices 55 additionalSpence109Row

theorem additionalSpence109_strong : SpenceMatrix.Strong 109 0 additionalSpence109 :=
  quad_matrices_strong additionalSpence109_sign_checked additionalSpence109_symmetric_checked
    additionalSpence109_minus_checked additionalSpence109_plus_checked

def additionalPaley109Tree : LiteralTree ℤ :=
  (.node 54 (.node 27 (.node 13 (.node 6 (.node 3 (.node 1 (.leaf (0)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (1)))))) (.node 7 (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (1))))))) (.node 13 (.node 6 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (-1)))))) (.node 7 (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1)))))))) (.node 27 (.node 13 (.node 6 (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1)))))) (.node 7 (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1))))))) (.node 14 (.node 7 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (1)))))) (.node 7 (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1)))))))))

def additionalPaley109Row : ℕ → ℤ := treeAt additionalPaley109Tree

theorem additionalPaley109_sign_checked :
    checkPaleySign 109 additionalPaley109Row = true := by decide +kernel
theorem additionalPaley109_symmetric_checked :
    checkPaleySymmetric 109 additionalPaley109Row = true := by decide +kernel
theorem additionalPaley109_sum_checked :
    natSum 109 additionalPaley109Row = 0 := by decide +kernel
theorem additionalPaley109_corr_checked :
    checkPaleyCorr 109 additionalPaley109Row = true := by decide +kernel

def additionalPaley109 : Matrix (ZMod 109) (ZMod 109) ℚ :=
  GroupAlgebra.regular (fun x => (additionalPaley109Row x.val : ℚ))

theorem additionalPaley109_isPaley : SpenceMatrix.IsPaley additionalPaley109 :=
  prime_paley_checked additionalPaley109Row additionalPaley109_sign_checked
    additionalPaley109_symmetric_checked additionalPaley109_sum_checked
    additionalPaley109_corr_checked

#print axioms additionalSpence109_strong
#print axioms additionalPaley109_isPaley

end ButsonTowers.Data

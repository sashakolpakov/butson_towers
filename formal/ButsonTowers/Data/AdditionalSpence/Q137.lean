import ButsonTowers.Data.StrongChecks
import ButsonTowers.Data.PaleyChecks
import ButsonTowers.Data.LiteralTree

/-! Additional strong starter for q=137; width 69.
Channel order 0,1,2,3; channel signs (1, -1, 1, 1).
The literal Paley row is also checked in full: no primality oracle is a proof premise. -/

namespace ButsonTowers.Data
set_option maxRecDepth 1000000
set_option maxHeartbeats 0

def additionalSpence137Trees : Fin 4 → LiteralTree ℤ :=
  ![(.node 34 (.node 17 (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1))))))) (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1)))))))) (.node 17 (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1))))))) (.node 9 (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1))))))))),
    (.node 34 (.node 17 (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1))))))) (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1)))))))) (.node 17 (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1))))))) (.node 9 (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1))))))))),
    (.node 34 (.node 17 (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1))))))) (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))))))) (.node 17 (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1))))))) (.node 9 (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1)))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1))))))))),
    (.node 34 (.node 17 (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1))))))) (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))))))) (.node 17 (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1))))))) (.node 9 (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1)))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1)))))))))]

def additionalSpence137Row (i : Fin 4) (x : ℕ) : ℤ :=
  treeAt (additionalSpence137Trees i) x

theorem additionalSpence137_sign_checked :
    checkQuadSign 69 additionalSpence137Row = true := by decide +kernel

theorem additionalSpence137_symmetric_checked :
    checkQuadSymmetric 69 additionalSpence137Row = true := by decide +kernel

theorem additionalSpence137_minus_checked :
    checkQuadCorr 69 (pairRows additionalSpence137Row false) 4 = true := by decide +kernel

theorem additionalSpence137_plus_checked :
    checkQuadCorr 69 (pairRows additionalSpence137Row true) (4*(137 : ℤ)) = true := by decide +kernel

def additionalSpence137 : Fin 4 → Matrix (ZMod 69) (ZMod 69) ℚ :=
  rowMatrices 69 additionalSpence137Row

theorem additionalSpence137_strong : SpenceMatrix.Strong 137 0 additionalSpence137 :=
  quad_matrices_strong additionalSpence137_sign_checked additionalSpence137_symmetric_checked
    additionalSpence137_minus_checked additionalSpence137_plus_checked

def additionalPaley137Tree : LiteralTree ℤ :=
  (.node 68 (.node 34 (.node 17 (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (0)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1))))))) (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1)))))))) (.node 17 (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1))))))) (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1))))))))) (.node 34 (.node 17 (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1))))))) (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1)))))))) (.node 17 (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1))))))) (.node 9 (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1)))))) (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1))))))))))

def additionalPaley137Row : ℕ → ℤ := treeAt additionalPaley137Tree

theorem additionalPaley137_sign_checked :
    checkPaleySign 137 additionalPaley137Row = true := by decide +kernel
theorem additionalPaley137_symmetric_checked :
    checkPaleySymmetric 137 additionalPaley137Row = true := by decide +kernel
theorem additionalPaley137_sum_checked :
    natSum 137 additionalPaley137Row = 0 := by decide +kernel
theorem additionalPaley137_corr_checked :
    checkPaleyCorr 137 additionalPaley137Row = true := by decide +kernel

def additionalPaley137 : Matrix (ZMod 137) (ZMod 137) ℚ :=
  GroupAlgebra.regular (fun x => (additionalPaley137Row x.val : ℚ))

theorem additionalPaley137_isPaley : SpenceMatrix.IsPaley additionalPaley137 :=
  prime_paley_checked additionalPaley137Row additionalPaley137_sign_checked
    additionalPaley137_symmetric_checked additionalPaley137_sum_checked
    additionalPaley137_corr_checked

#print axioms additionalSpence137_strong
#print axioms additionalPaley137_isPaley

end ButsonTowers.Data

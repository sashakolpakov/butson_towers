import ButsonTowers.Data.StrongChecks
import ButsonTowers.Data.PaleyChecks
import ButsonTowers.Data.LiteralTree

/-! Additional strong starter for q=97; width 49.
Channel order 0,1,2,3; channel signs (1, 1, 1, -1).
The literal Paley row is also checked in full: no primality oracle is a proof premise. -/

namespace ButsonTowers.Data
set_option maxRecDepth 1000000
set_option maxHeartbeats 0

def additionalSpence97Trees : Fin 4 → LiteralTree ℤ :=
  ![(.node 24 (.node 12 (.node 6 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1)))))) (.node 6 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1))))))) (.node 12 (.node 6 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))))) (.node 6 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (1)))))))),
    (.node 24 (.node 12 (.node 6 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1)))))) (.node 6 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1))))))) (.node 12 (.node 6 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))))) (.node 6 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (1)))))))),
    (.node 24 (.node 12 (.node 6 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))))) (.node 6 (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1))))))) (.node 12 (.node 6 (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))))) (.node 6 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1)))))))),
    (.node 24 (.node 12 (.node 6 (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))))) (.node 6 (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1))))))) (.node 12 (.node 6 (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))))) (.node 6 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1))))))))]

def additionalSpence97Row (i : Fin 4) (x : ℕ) : ℤ :=
  treeAt (additionalSpence97Trees i) x

theorem additionalSpence97_sign_checked :
    checkQuadSign 49 additionalSpence97Row = true := by decide +kernel

theorem additionalSpence97_symmetric_checked :
    checkQuadSymmetric 49 additionalSpence97Row = true := by decide +kernel

theorem additionalSpence97_minus_checked :
    checkQuadCorr 49 (pairRows additionalSpence97Row false) 4 = true := by decide +kernel

theorem additionalSpence97_plus_checked :
    checkQuadCorr 49 (pairRows additionalSpence97Row true) (4*(97 : ℤ)) = true := by decide +kernel

def additionalSpence97 : Fin 4 → Matrix (ZMod 49) (ZMod 49) ℚ :=
  rowMatrices 49 additionalSpence97Row

theorem additionalSpence97_strong : SpenceMatrix.Strong 97 0 additionalSpence97 :=
  quad_matrices_strong additionalSpence97_sign_checked additionalSpence97_symmetric_checked
    additionalSpence97_minus_checked additionalSpence97_plus_checked

def additionalPaley97Tree : LiteralTree ℤ :=
  (.node 48 (.node 24 (.node 12 (.node 6 (.node 3 (.node 1 (.leaf (0)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))))) (.node 6 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1))))))) (.node 12 (.node 6 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))))) (.node 6 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1)))))))) (.node 24 (.node 12 (.node 6 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1)))))) (.node 6 (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1))))))) (.node 12 (.node 6 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1)))))) (.node 6 (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (1)))))))))

def additionalPaley97Row : ℕ → ℤ := treeAt additionalPaley97Tree

theorem additionalPaley97_sign_checked :
    checkPaleySign 97 additionalPaley97Row = true := by decide +kernel
theorem additionalPaley97_symmetric_checked :
    checkPaleySymmetric 97 additionalPaley97Row = true := by decide +kernel
theorem additionalPaley97_sum_checked :
    natSum 97 additionalPaley97Row = 0 := by decide +kernel
theorem additionalPaley97_corr_checked :
    checkPaleyCorr 97 additionalPaley97Row = true := by decide +kernel

def additionalPaley97 : Matrix (ZMod 97) (ZMod 97) ℚ :=
  GroupAlgebra.regular (fun x => (additionalPaley97Row x.val : ℚ))

theorem additionalPaley97_isPaley : SpenceMatrix.IsPaley additionalPaley97 :=
  prime_paley_checked additionalPaley97Row additionalPaley97_sign_checked
    additionalPaley97_symmetric_checked additionalPaley97_sum_checked
    additionalPaley97_corr_checked

#print axioms additionalSpence97_strong
#print axioms additionalPaley97_isPaley

end ButsonTowers.Data

import ButsonTowers.Data.StrongChecks
import ButsonTowers.Data.PaleyChecks
import ButsonTowers.Data.LiteralTree

/-! Additional strong starter for q=5; width 3.
Channel order 0,1,2,3; channel signs (1, -1, 1, 1).
The literal Paley row is also checked in full: no primality oracle is a proof premise. -/

namespace ButsonTowers.Data
set_option maxRecDepth 1000000
set_option maxHeartbeats 0

def additionalSpence5Trees : Fin 4 → LiteralTree ℤ :=
  ![(.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1)))),
    (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1)))),
    (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1)))),
    (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1))))]

def additionalSpence5Row (i : Fin 4) (x : ℕ) : ℤ :=
  treeAt (additionalSpence5Trees i) x

theorem additionalSpence5_sign_checked :
    checkQuadSign 3 additionalSpence5Row = true := by decide +kernel

theorem additionalSpence5_symmetric_checked :
    checkQuadSymmetric 3 additionalSpence5Row = true := by decide +kernel

theorem additionalSpence5_minus_checked :
    checkQuadCorr 3 (pairRows additionalSpence5Row false) 4 = true := by decide +kernel

theorem additionalSpence5_plus_checked :
    checkQuadCorr 3 (pairRows additionalSpence5Row true) (4*(5 : ℤ)) = true := by decide +kernel

def additionalSpence5 : Fin 4 → Matrix (ZMod 3) (ZMod 3) ℚ :=
  rowMatrices 3 additionalSpence5Row

theorem additionalSpence5_strong : SpenceMatrix.Strong 5 0 additionalSpence5 :=
  quad_matrices_strong additionalSpence5_sign_checked additionalSpence5_symmetric_checked
    additionalSpence5_minus_checked additionalSpence5_plus_checked

def additionalPaley5Tree : LiteralTree ℤ :=
  (.node 2 (.node 1 (.leaf (0)) (.leaf (1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1)))))

def additionalPaley5Row : ℕ → ℤ := treeAt additionalPaley5Tree

theorem additionalPaley5_sign_checked :
    checkPaleySign 5 additionalPaley5Row = true := by decide +kernel
theorem additionalPaley5_symmetric_checked :
    checkPaleySymmetric 5 additionalPaley5Row = true := by decide +kernel
theorem additionalPaley5_sum_checked :
    natSum 5 additionalPaley5Row = 0 := by decide +kernel
theorem additionalPaley5_corr_checked :
    checkPaleyCorr 5 additionalPaley5Row = true := by decide +kernel

def additionalPaley5 : Matrix (ZMod 5) (ZMod 5) ℚ :=
  GroupAlgebra.regular (fun x => (additionalPaley5Row x.val : ℚ))

theorem additionalPaley5_isPaley : SpenceMatrix.IsPaley additionalPaley5 :=
  prime_paley_checked additionalPaley5Row additionalPaley5_sign_checked
    additionalPaley5_symmetric_checked additionalPaley5_sum_checked
    additionalPaley5_corr_checked

#print axioms additionalSpence5_strong
#print axioms additionalPaley5_isPaley

end ButsonTowers.Data

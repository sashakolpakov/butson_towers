import ButsonTowers.Data.StrongChecks
import ButsonTowers.Data.PaleyChecks
import ButsonTowers.Data.LiteralTree

/-! Additional strong starter for q=53; width 27.
Channel order 0,1,2,3; channel signs (1, -1, 1, 1).
The literal Paley row is also checked in full: no primality oracle is a proof premise. -/

namespace ButsonTowers.Data
set_option maxRecDepth 1000000
set_option maxHeartbeats 0

def additionalSpence53Trees : Fin 4 → LiteralTree ℤ :=
  ![(.node 13 (.node 6 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (-1)))))) (.node 7 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (-1))))))),
    (.node 13 (.node 6 (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (-1)))))) (.node 7 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (-1))))))),
    (.node 13 (.node 6 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1)))))) (.node 7 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1))))))),
    (.node 13 (.node 6 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1)))))) (.node 7 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1)))))))]

def additionalSpence53Row (i : Fin 4) (x : ℕ) : ℤ :=
  treeAt (additionalSpence53Trees i) x

theorem additionalSpence53_sign_checked :
    checkQuadSign 27 additionalSpence53Row = true := by decide +kernel

theorem additionalSpence53_symmetric_checked :
    checkQuadSymmetric 27 additionalSpence53Row = true := by decide +kernel

theorem additionalSpence53_minus_checked :
    checkQuadCorr 27 (pairRows additionalSpence53Row false) 4 = true := by decide +kernel

theorem additionalSpence53_plus_checked :
    checkQuadCorr 27 (pairRows additionalSpence53Row true) (4*(53 : ℤ)) = true := by decide +kernel

def additionalSpence53 : Fin 4 → Matrix (ZMod 27) (ZMod 27) ℚ :=
  rowMatrices 27 additionalSpence53Row

theorem additionalSpence53_strong : SpenceMatrix.Strong 53 0 additionalSpence53 :=
  quad_matrices_strong additionalSpence53_sign_checked additionalSpence53_symmetric_checked
    additionalSpence53_minus_checked additionalSpence53_plus_checked

def additionalPaley53Tree : LiteralTree ℤ :=
  (.node 26 (.node 13 (.node 6 (.node 3 (.node 1 (.leaf (0)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (-1)))))) (.node 6 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (1))))))) (.node 13 (.node 6 (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (1)))))) (.node 7 (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (1))))))))

def additionalPaley53Row : ℕ → ℤ := treeAt additionalPaley53Tree

theorem additionalPaley53_sign_checked :
    checkPaleySign 53 additionalPaley53Row = true := by decide +kernel
theorem additionalPaley53_symmetric_checked :
    checkPaleySymmetric 53 additionalPaley53Row = true := by decide +kernel
theorem additionalPaley53_sum_checked :
    natSum 53 additionalPaley53Row = 0 := by decide +kernel
theorem additionalPaley53_corr_checked :
    checkPaleyCorr 53 additionalPaley53Row = true := by decide +kernel

def additionalPaley53 : Matrix (ZMod 53) (ZMod 53) ℚ :=
  GroupAlgebra.regular (fun x => (additionalPaley53Row x.val : ℚ))

theorem additionalPaley53_isPaley : SpenceMatrix.IsPaley additionalPaley53 :=
  prime_paley_checked additionalPaley53Row additionalPaley53_sign_checked
    additionalPaley53_symmetric_checked additionalPaley53_sum_checked
    additionalPaley53_corr_checked

#print axioms additionalSpence53_strong
#print axioms additionalPaley53_isPaley

end ButsonTowers.Data

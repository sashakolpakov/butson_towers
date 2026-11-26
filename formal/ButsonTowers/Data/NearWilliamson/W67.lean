import ButsonTowers.Data.NearWilliamson.Core

/-! Literal width-67 source. Primary data: https://www.cs.uleth.ca/~zaitsevv/nw/NW67.txt
Raw data SHA256: 429c347b25e2f4eaf0924b5761fff114ce489794ee86da2f30cf5e322e1733a2.
Literal four-row SHA256(int8,row-major): edb93a64a3430aeac83ad5567c75284319c04de695774b1df56e63fe6e33faf3.
First matrix column-reversed: true.
The following closed finite checks are reduced by the Lean kernel. -/

namespace ButsonTowers.Data
open scoped BigOperators
set_option maxRecDepth 1000000
set_option maxHeartbeats 0

def nearWilliamson67Trees : Fin 4 → LiteralTree ℤ :=
  ![(.node 33 (.node 16 (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (1)))))) (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1)))))))) (.node 17 (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1))))))) (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1))))))))),
    (.node 33 (.node 16 (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (1)))))) (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1)))))))) (.node 17 (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1))))))) (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1))))))))),
    (.node 33 (.node 16 (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (-1)))))) (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1)))))))) (.node 17 (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1))))))) (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1))))))))),
    (.node 33 (.node 16 (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (1)))))) (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))))))) (.node 17 (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1))))))) (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1)))))))))]

def nearWilliamson67Row (i : Fin 4) (x : ℕ) : ℤ :=
  treeAt (nearWilliamson67Trees i) x

theorem nearWilliamson67_sign_checked :
    checkQuadSign 67 nearWilliamson67Row = true := by decide +kernel

theorem nearWilliamson67_symmetric_checked :
    checkQuadSymmetric 67 (tailRows nearWilliamson67Row) = true := by decide +kernel

theorem nearWilliamson67_norm_checked :
    checkQuadCorr 67 nearWilliamson67Row 268 = true := by decide +kernel

def nearWilliamson67 : Fin 4 → Matrix (ZMod 67) (ZMod 67) ℚ :=
  nearMatrices 67 nearWilliamson67Row

theorem nearWilliamson67_quad : AmicableQuad nearWilliamson67 (4*(67 : ℚ)) :=
  near_matrices_quad nearWilliamson67_sign_checked
    nearWilliamson67_symmetric_checked nearWilliamson67_norm_checked

#print axioms nearWilliamson67_quad

end ButsonTowers.Data

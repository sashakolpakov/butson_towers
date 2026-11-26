import ButsonTowers.Data.NearWilliamson.Core

/-! Literal width-59 source. Primary data: https://www.cs.uleth.ca/~zaitsevv/nw/NW59.txt
Raw data SHA256: edcc7413d9e6d72bb516bf758871c1d45f4221db62abb4518638f2cd5e5d58a3.
Literal four-row SHA256(int8,row-major): c8942f16d5384560bd22e984fffb372b6ea482e54e0460f5fe7614275ed7439d.
First matrix column-reversed: true.
The following closed finite checks are reduced by the Lean kernel. -/

namespace ButsonTowers.Data
open scoped BigOperators
set_option maxRecDepth 1000000
set_option maxHeartbeats 0

def nearWilliamson59Trees : Fin 4 → LiteralTree ℤ :=
  ![(.node 29 (.node 14 (.node 7 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (1)))))) (.node 7 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (1))))))) (.node 15 (.node 7 (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (-1)))))) (.node 7 (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1)))))))),
    (.node 29 (.node 14 (.node 7 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1)))))) (.node 7 (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1))))))) (.node 15 (.node 7 (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (1)))))) (.node 7 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1)))))))),
    (.node 29 (.node 14 (.node 7 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1)))))) (.node 7 (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (1))))))) (.node 15 (.node 7 (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (-1)))))) (.node 7 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (-1)))))))),
    (.node 29 (.node 14 (.node 7 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (-1)))))) (.node 7 (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (1))))))) (.node 15 (.node 7 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (1)))))) (.node 7 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1))))))))]

def nearWilliamson59Row (i : Fin 4) (x : ℕ) : ℤ :=
  treeAt (nearWilliamson59Trees i) x

theorem nearWilliamson59_sign_checked :
    checkQuadSign 59 nearWilliamson59Row = true := by decide +kernel

theorem nearWilliamson59_symmetric_checked :
    checkQuadSymmetric 59 (tailRows nearWilliamson59Row) = true := by decide +kernel

theorem nearWilliamson59_norm_checked :
    checkQuadCorr 59 nearWilliamson59Row 236 = true := by decide +kernel

def nearWilliamson59 : Fin 4 → Matrix (ZMod 59) (ZMod 59) ℚ :=
  nearMatrices 59 nearWilliamson59Row

theorem nearWilliamson59_quad : AmicableQuad nearWilliamson59 (4*(59 : ℚ)) :=
  near_matrices_quad nearWilliamson59_sign_checked
    nearWilliamson59_symmetric_checked nearWilliamson59_norm_checked

#print axioms nearWilliamson59_quad

end ButsonTowers.Data

import ButsonTowers.Data.NearWilliamson.Core

/-! Literal width-73 source. Primary data: https://www.cs.uleth.ca/~zaitsevv/nw/NW73.txt
Raw data SHA256: 94b4bf2b6c0ab6aa51a110c994dcef93a2ce5beeba9af3f5535b07580ec44ab7.
Literal four-row SHA256(int8,row-major): b4b8cdc94e7c60c1de4728e916af1012960fd836a1de46ff2cbe60ffecda1973.
First matrix column-reversed: true.
The following closed finite checks are reduced by the Lean kernel. -/

namespace ButsonTowers.Data
open scoped BigOperators
set_option maxRecDepth 1000000
set_option maxHeartbeats 0

def nearWilliamson73Trees : Fin 4 → LiteralTree ℤ :=
  ![(.node 36 (.node 18 (.node 9 (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1)))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1))))))) (.node 9 (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1)))))) (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))))))) (.node 18 (.node 9 (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1)))))) (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1))))))) (.node 9 (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))))) (.node 5 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1))))))))),
    (.node 36 (.node 18 (.node 9 (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1))))))) (.node 9 (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1)))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1)))))))) (.node 18 (.node 9 (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1)))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1))))))) (.node 9 (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))))) (.node 5 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1))))))))),
    (.node 36 (.node 18 (.node 9 (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1)))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1))))))) (.node 9 (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1)))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))))))) (.node 18 (.node 9 (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1)))))) (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1))))))) (.node 9 (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1)))))) (.node 5 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1))))))))),
    (.node 36 (.node 18 (.node 9 (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1)))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1))))))) (.node 9 (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1)))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))))))) (.node 18 (.node 9 (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1)))))) (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1))))))) (.node 9 (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))))) (.node 5 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1)))))))))]

def nearWilliamson73Row (i : Fin 4) (x : ℕ) : ℤ :=
  treeAt (nearWilliamson73Trees i) x

theorem nearWilliamson73_sign_checked :
    checkQuadSign 73 nearWilliamson73Row = true := by decide +kernel

theorem nearWilliamson73_symmetric_checked :
    checkQuadSymmetric 73 (tailRows nearWilliamson73Row) = true := by decide +kernel

theorem nearWilliamson73_norm_checked :
    checkQuadCorr 73 nearWilliamson73Row 292 = true := by decide +kernel

def nearWilliamson73 : Fin 4 → Matrix (ZMod 73) (ZMod 73) ℚ :=
  nearMatrices 73 nearWilliamson73Row

theorem nearWilliamson73_quad : AmicableQuad nearWilliamson73 (4*(73 : ℚ)) :=
  near_matrices_quad nearWilliamson73_sign_checked
    nearWilliamson73_symmetric_checked nearWilliamson73_norm_checked

#print axioms nearWilliamson73_quad

end ButsonTowers.Data

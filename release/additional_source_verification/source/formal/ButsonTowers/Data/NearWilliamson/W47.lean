import ButsonTowers.Data.NearWilliamson.Core

/-! Literal width-47 source. Primary data: https://www.cs.uleth.ca/~zaitsevv/nw/db/NW47_first1000.txt
Raw data SHA256: 552ec602f53d469408b47a9ce3bd15d8519b3d5d61cc91a8a888926b7b8b407a.
Literal four-row SHA256(int8,row-major): 4b78981450b3cdd263e16fdccc064797433b897998762341d76cf40cc72d2eda.
First matrix column-reversed: true.
The following closed finite checks are reduced by the Lean kernel. -/

namespace ButsonTowers.Data
open scoped BigOperators
set_option maxRecDepth 1000000
set_option maxHeartbeats 0

def nearWilliamson47Trees : Fin 4 → LiteralTree ℤ :=
  ![(.node 23 (.node 11 (.node 5 (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1)))))) (.node 6 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1))))))) (.node 12 (.node 6 (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))))) (.node 6 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1)))))))),
    (.node 23 (.node 11 (.node 5 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))))) (.node 6 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1))))))) (.node 12 (.node 6 (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1)))))) (.node 6 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1)))))))),
    (.node 23 (.node 11 (.node 5 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1)))))) (.node 6 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1))))))) (.node 12 (.node 6 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))))) (.node 6 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))))))),
    (.node 23 (.node 11 (.node 5 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1)))))) (.node 6 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1))))))) (.node 12 (.node 6 (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 3 (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1)))))) (.node 6 (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 3 (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1))))))))]

def nearWilliamson47Row (i : Fin 4) (x : ℕ) : ℤ :=
  treeAt (nearWilliamson47Trees i) x

theorem nearWilliamson47_sign_checked :
    checkQuadSign 47 nearWilliamson47Row = true := by decide +kernel

theorem nearWilliamson47_symmetric_checked :
    checkQuadSymmetric 47 (tailRows nearWilliamson47Row) = true := by decide +kernel

theorem nearWilliamson47_norm_checked :
    checkQuadCorr 47 nearWilliamson47Row 188 = true := by decide +kernel

def nearWilliamson47 : Fin 4 → Matrix (ZMod 47) (ZMod 47) ℚ :=
  nearMatrices 47 nearWilliamson47Row

theorem nearWilliamson47_quad : AmicableQuad nearWilliamson47 (4*(47 : ℚ)) :=
  near_matrices_quad nearWilliamson47_sign_checked
    nearWilliamson47_symmetric_checked nearWilliamson47_norm_checked

#print axioms nearWilliamson47_quad

end ButsonTowers.Data

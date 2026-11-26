import ButsonTowers.Data.NearWilliamson.Core

/-! Literal width-69 source. Primary data: https://zenodo.org/records/1195693/files/williamson-69.txt?download=1
Raw data SHA256: bfc68b4172d52e881fcd417d97b4a9b779b2dca57767fa54276fff0173f9f9f3.
Literal four-row SHA256(int8,row-major): f2283a3ee736bcca9e52cc61c366a004c667fcbf198e0e78a8f6f3b713362f15.
First matrix column-reversed: false.
The following closed finite checks are reduced by the Lean kernel. -/

namespace ButsonTowers.Data
open scoped BigOperators
set_option maxRecDepth 1000000
set_option maxHeartbeats 0

def nearWilliamson69Trees : Fin 4 → LiteralTree ℤ :=
  ![(.node 34 (.node 17 (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1))))))) (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (-1)))))))) (.node 17 (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1))))))) (.node 9 (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1)))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1))))))))),
    (.node 34 (.node 17 (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1))))))) (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))))))) (.node 17 (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1))))))) (.node 9 (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (1)))))) (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (1))))))))),
    (.node 34 (.node 17 (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1))))))) (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))))))) (.node 17 (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1))))))) (.node 9 (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1)))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1))))))))),
    (.node 34 (.node 17 (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (1)) (.leaf (-1))))))) (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (1)) (.leaf (-1))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (1)))))))) (.node 17 (.node 8 (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1))))) (.node 4 (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (-1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (-1)) (.leaf (-1))))))) (.node 9 (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (1))) (.node 1 (.leaf (1)) (.node 1 (.leaf (-1)) (.leaf (-1)))))) (.node 4 (.node 2 (.node 1 (.leaf (-1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.leaf (1)))) (.node 2 (.node 1 (.leaf (1)) (.leaf (-1))) (.node 1 (.leaf (-1)) (.node 1 (.leaf (1)) (.leaf (1)))))))))]

def nearWilliamson69Row (i : Fin 4) (x : ℕ) : ℤ :=
  treeAt (nearWilliamson69Trees i) x

theorem nearWilliamson69_sign_checked :
    checkQuadSign 69 nearWilliamson69Row = true := by decide +kernel

theorem nearWilliamson69_symmetric_checked :
    checkQuadSymmetric 69 nearWilliamson69Row = true := by decide +kernel

theorem nearWilliamson69_norm_checked :
    checkQuadCorr 69 nearWilliamson69Row 276 = true := by decide +kernel

def nearWilliamson69 : Fin 4 → Matrix (ZMod 69) (ZMod 69) ℚ :=
  rowMatrices 69 nearWilliamson69Row

theorem nearWilliamson69_quad : AmicableQuad nearWilliamson69 (4*(69 : ℚ)) :=
  ordinary_matrices_quad nearWilliamson69_sign_checked
    nearWilliamson69_symmetric_checked nearWilliamson69_norm_checked

#print axioms nearWilliamson69_quad

end ButsonTowers.Data

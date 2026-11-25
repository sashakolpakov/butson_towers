import ButsonTowers.Data.StrongChecks
import ButsonTowers.Data.TreeLiterals
namespace ButsonTowers.Data
set_option maxRecDepth 1000000
set_option maxHeartbeats 0

def spence13Row (i : Fin 4) (x : ℕ) : ℤ := treeAt (spence13Trees i) x

theorem spence13_sign_checked : checkQuadSign 7 spence13Row = true := by decide +kernel

theorem spence13_symmetric_checked : checkQuadSymmetric 7 spence13Row = true := by decide +kernel

theorem spence13_minus_checked : checkQuadCorr 7 (pairRows spence13Row false) 4 = true := by decide +kernel

theorem spence13_plus_checked : checkQuadCorr 7 (pairRows spence13Row true) (4*(13 : ℤ)) = true := by decide +kernel

def spence13 : Fin 4 → Matrix (ZMod 7) (ZMod 7) ℚ := rowMatrices 7 spence13Row

theorem spence13_strong : SpenceMatrix.Strong 13 0 spence13 :=
  quad_matrices_strong spence13_sign_checked spence13_symmetric_checked
    spence13_minus_checked spence13_plus_checked



def spence17Row (i : Fin 4) (x : ℕ) : ℤ := treeAt (spence17Trees i) x

theorem spence17_sign_checked : checkQuadSign 9 spence17Row = true := by decide +kernel

theorem spence17_symmetric_checked : checkQuadSymmetric 9 spence17Row = true := by decide +kernel

theorem spence17_minus_checked : checkQuadCorr 9 (pairRows spence17Row false) 4 = true := by decide +kernel

theorem spence17_plus_checked : checkQuadCorr 9 (pairRows spence17Row true) (4*(17 : ℤ)) = true := by decide +kernel

def spence17 : Fin 4 → Matrix (ZMod 9) (ZMod 9) ℚ := rowMatrices 9 spence17Row

theorem spence17_strong : SpenceMatrix.Strong 17 0 spence17 :=
  quad_matrices_strong spence17_sign_checked spence17_symmetric_checked
    spence17_minus_checked spence17_plus_checked



def spence61Row (i : Fin 4) (x : ℕ) : ℤ := treeAt (spence61Trees i) x

theorem spence61_sign_checked : checkQuadSign 31 spence61Row = true := by decide +kernel

theorem spence61_symmetric_checked : checkQuadSymmetric 31 spence61Row = true := by decide +kernel

theorem spence61_minus_checked : checkQuadCorr 31 (pairRows spence61Row false) 4 = true := by decide +kernel

theorem spence61_plus_checked : checkQuadCorr 31 (pairRows spence61Row true) (4*(61 : ℤ)) = true := by decide +kernel

def spence61 : Fin 4 → Matrix (ZMod 31) (ZMod 31) ℚ := rowMatrices 31 spence61Row

theorem spence61_strong : SpenceMatrix.Strong 61 0 spence61 :=
  quad_matrices_strong spence61_sign_checked spence61_symmetric_checked
    spence61_minus_checked spence61_plus_checked



def spence101Row (i : Fin 4) (x : ℕ) : ℤ := treeAt (spence101Trees i) x

theorem spence101_sign_checked : checkQuadSign 51 spence101Row = true := by decide +kernel

theorem spence101_symmetric_checked : checkQuadSymmetric 51 spence101Row = true := by decide +kernel

theorem spence101_minus_checked : checkQuadCorr 51 (pairRows spence101Row false) 4 = true := by decide +kernel

theorem spence101_plus_checked : checkQuadCorr 51 (pairRows spence101Row true) (4*(101 : ℤ)) = true := by decide +kernel

def spence101 : Fin 4 → Matrix (ZMod 51) (ZMod 51) ℚ := rowMatrices 51 spence101Row

theorem spence101_strong : SpenceMatrix.Strong 101 0 spence101 :=
  quad_matrices_strong spence101_sign_checked spence101_symmetric_checked
    spence101_minus_checked spence101_plus_checked



def spence121Row (i : Fin 4) (x : ℕ) : ℤ := treeAt (spence121Trees i) x

theorem spence121_sign_checked : checkQuadSign 61 spence121Row = true := by decide +kernel

theorem spence121_symmetric_checked : checkQuadSymmetric 61 spence121Row = true := by decide +kernel

theorem spence121_minus_checked : checkQuadCorr 61 (pairRows spence121Row false) 4 = true := by decide +kernel

theorem spence121_plus_checked : checkQuadCorr 61 (pairRows spence121Row true) (4*(121 : ℤ)) = true := by decide +kernel

def spence121 : Fin 4 → Matrix (ZMod 61) (ZMod 61) ℚ := rowMatrices 61 spence121Row

theorem spence121_strong : SpenceMatrix.Strong 121 0 spence121 :=
  quad_matrices_strong spence121_sign_checked spence121_symmetric_checked
    spence121_minus_checked spence121_plus_checked



#print axioms spence121_strong
end ButsonTowers.Data

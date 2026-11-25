import ButsonTowers.Data.RowChecks
import ButsonTowers.Data.TreeLiterals

namespace ButsonTowers.Data
open scoped BigOperators
set_option maxRecDepth 1000000
set_option maxHeartbeats 0

def williamson1Row (i : Fin 4) (x : ℕ) : ℤ := treeAt (williamson1Trees i) x

theorem williamson1_sign_checked : checkQuadSign 1 williamson1Row = true := by decide +kernel

theorem williamson1_symmetric_checked : checkQuadSymmetric 1 williamson1Row = true := by decide +kernel

theorem williamson1_norm_checked : checkQuadCorr 1 williamson1Row 4 = true := by decide +kernel

def williamson1 : Fin 4 → Matrix (ZMod 1) (ZMod 1) ℚ := rowMatrices 1 williamson1Row

theorem williamson1_sign : ∀ i, IsSignMatrix (williamson1 i) :=
  quad_matrices_sign williamson1_sign_checked

theorem williamson1_symmetric : ∀ i, (williamson1 i).transpose = williamson1 i :=
  quad_matrices_symmetric williamson1_symmetric_checked

theorem williamson1_commute : ∀ i j, williamson1 i * williamson1 j = williamson1 j * williamson1 i :=
  quad_matrices_commute williamson1Row

theorem williamson1_amicable : ∀ i j, williamson1 i * (williamson1 j).transpose =
    williamson1 j * (williamson1 i).transpose := by
  intro i j
  rw [williamson1_symmetric i, williamson1_symmetric j]
  exact williamson1_commute i j

theorem williamson1_norm : (∑ i : Fin 4, williamson1 i * (williamson1 i).transpose) =
    (4 * (1 : ℚ)) • (1 : Matrix (ZMod 1) (ZMod 1) ℚ) := by
  have h := quad_matrices_norm williamson1_norm_checked
  norm_num at h ⊢
  exact h



def williamson3Row (i : Fin 4) (x : ℕ) : ℤ := treeAt (williamson3Trees i) x

theorem williamson3_sign_checked : checkQuadSign 3 williamson3Row = true := by decide +kernel

theorem williamson3_symmetric_checked : checkQuadSymmetric 3 williamson3Row = true := by decide +kernel

theorem williamson3_norm_checked : checkQuadCorr 3 williamson3Row 12 = true := by decide +kernel

def williamson3 : Fin 4 → Matrix (ZMod 3) (ZMod 3) ℚ := rowMatrices 3 williamson3Row

theorem williamson3_sign : ∀ i, IsSignMatrix (williamson3 i) :=
  quad_matrices_sign williamson3_sign_checked

theorem williamson3_symmetric : ∀ i, (williamson3 i).transpose = williamson3 i :=
  quad_matrices_symmetric williamson3_symmetric_checked

theorem williamson3_commute : ∀ i j, williamson3 i * williamson3 j = williamson3 j * williamson3 i :=
  quad_matrices_commute williamson3Row

theorem williamson3_amicable : ∀ i j, williamson3 i * (williamson3 j).transpose =
    williamson3 j * (williamson3 i).transpose := by
  intro i j
  rw [williamson3_symmetric i, williamson3_symmetric j]
  exact williamson3_commute i j

theorem williamson3_norm : (∑ i : Fin 4, williamson3 i * (williamson3 i).transpose) =
    (4 * (3 : ℚ)) • (1 : Matrix (ZMod 3) (ZMod 3) ℚ) := by
  have h := quad_matrices_norm williamson3_norm_checked
  norm_num at h ⊢
  exact h



def williamson5Row (i : Fin 4) (x : ℕ) : ℤ := treeAt (williamson5Trees i) x

theorem williamson5_sign_checked : checkQuadSign 5 williamson5Row = true := by decide +kernel

theorem williamson5_symmetric_checked : checkQuadSymmetric 5 williamson5Row = true := by decide +kernel

theorem williamson5_norm_checked : checkQuadCorr 5 williamson5Row 20 = true := by decide +kernel

def williamson5 : Fin 4 → Matrix (ZMod 5) (ZMod 5) ℚ := rowMatrices 5 williamson5Row

theorem williamson5_sign : ∀ i, IsSignMatrix (williamson5 i) :=
  quad_matrices_sign williamson5_sign_checked

theorem williamson5_symmetric : ∀ i, (williamson5 i).transpose = williamson5 i :=
  quad_matrices_symmetric williamson5_symmetric_checked

theorem williamson5_commute : ∀ i j, williamson5 i * williamson5 j = williamson5 j * williamson5 i :=
  quad_matrices_commute williamson5Row

theorem williamson5_amicable : ∀ i j, williamson5 i * (williamson5 j).transpose =
    williamson5 j * (williamson5 i).transpose := by
  intro i j
  rw [williamson5_symmetric i, williamson5_symmetric j]
  exact williamson5_commute i j

theorem williamson5_norm : (∑ i : Fin 4, williamson5 i * (williamson5 i).transpose) =
    (4 * (5 : ℚ)) • (1 : Matrix (ZMod 5) (ZMod 5) ℚ) := by
  have h := quad_matrices_norm williamson5_norm_checked
  norm_num at h ⊢
  exact h



def williamson7Row (i : Fin 4) (x : ℕ) : ℤ := treeAt (williamson7Trees i) x

theorem williamson7_sign_checked : checkQuadSign 7 williamson7Row = true := by decide +kernel

theorem williamson7_symmetric_checked : checkQuadSymmetric 7 williamson7Row = true := by decide +kernel

theorem williamson7_norm_checked : checkQuadCorr 7 williamson7Row 28 = true := by decide +kernel

def williamson7 : Fin 4 → Matrix (ZMod 7) (ZMod 7) ℚ := rowMatrices 7 williamson7Row

theorem williamson7_sign : ∀ i, IsSignMatrix (williamson7 i) :=
  quad_matrices_sign williamson7_sign_checked

theorem williamson7_symmetric : ∀ i, (williamson7 i).transpose = williamson7 i :=
  quad_matrices_symmetric williamson7_symmetric_checked

theorem williamson7_commute : ∀ i j, williamson7 i * williamson7 j = williamson7 j * williamson7 i :=
  quad_matrices_commute williamson7Row

theorem williamson7_amicable : ∀ i j, williamson7 i * (williamson7 j).transpose =
    williamson7 j * (williamson7 i).transpose := by
  intro i j
  rw [williamson7_symmetric i, williamson7_symmetric j]
  exact williamson7_commute i j

theorem williamson7_norm : (∑ i : Fin 4, williamson7 i * (williamson7 i).transpose) =
    (4 * (7 : ℚ)) • (1 : Matrix (ZMod 7) (ZMod 7) ℚ) := by
  have h := quad_matrices_norm williamson7_norm_checked
  norm_num at h ⊢
  exact h



def williamson9Row (i : Fin 4) (x : ℕ) : ℤ := treeAt (williamson9Trees i) x

theorem williamson9_sign_checked : checkQuadSign 9 williamson9Row = true := by decide +kernel

theorem williamson9_symmetric_checked : checkQuadSymmetric 9 williamson9Row = true := by decide +kernel

theorem williamson9_norm_checked : checkQuadCorr 9 williamson9Row 36 = true := by decide +kernel

def williamson9 : Fin 4 → Matrix (ZMod 9) (ZMod 9) ℚ := rowMatrices 9 williamson9Row

theorem williamson9_sign : ∀ i, IsSignMatrix (williamson9 i) :=
  quad_matrices_sign williamson9_sign_checked

theorem williamson9_symmetric : ∀ i, (williamson9 i).transpose = williamson9 i :=
  quad_matrices_symmetric williamson9_symmetric_checked

theorem williamson9_commute : ∀ i j, williamson9 i * williamson9 j = williamson9 j * williamson9 i :=
  quad_matrices_commute williamson9Row

theorem williamson9_amicable : ∀ i j, williamson9 i * (williamson9 j).transpose =
    williamson9 j * (williamson9 i).transpose := by
  intro i j
  rw [williamson9_symmetric i, williamson9_symmetric j]
  exact williamson9_commute i j

theorem williamson9_norm : (∑ i : Fin 4, williamson9 i * (williamson9 i).transpose) =
    (4 * (9 : ℚ)) • (1 : Matrix (ZMod 9) (ZMod 9) ℚ) := by
  have h := quad_matrices_norm williamson9_norm_checked
  norm_num at h ⊢
  exact h



def williamson11Row (i : Fin 4) (x : ℕ) : ℤ := treeAt (williamson11Trees i) x

theorem williamson11_sign_checked : checkQuadSign 11 williamson11Row = true := by decide +kernel

theorem williamson11_symmetric_checked : checkQuadSymmetric 11 williamson11Row = true := by decide +kernel

theorem williamson11_norm_checked : checkQuadCorr 11 williamson11Row 44 = true := by decide +kernel

def williamson11 : Fin 4 → Matrix (ZMod 11) (ZMod 11) ℚ := rowMatrices 11 williamson11Row

theorem williamson11_sign : ∀ i, IsSignMatrix (williamson11 i) :=
  quad_matrices_sign williamson11_sign_checked

theorem williamson11_symmetric : ∀ i, (williamson11 i).transpose = williamson11 i :=
  quad_matrices_symmetric williamson11_symmetric_checked

theorem williamson11_commute : ∀ i j, williamson11 i * williamson11 j = williamson11 j * williamson11 i :=
  quad_matrices_commute williamson11Row

theorem williamson11_amicable : ∀ i j, williamson11 i * (williamson11 j).transpose =
    williamson11 j * (williamson11 i).transpose := by
  intro i j
  rw [williamson11_symmetric i, williamson11_symmetric j]
  exact williamson11_commute i j

theorem williamson11_norm : (∑ i : Fin 4, williamson11 i * (williamson11 i).transpose) =
    (4 * (11 : ℚ)) • (1 : Matrix (ZMod 11) (ZMod 11) ℚ) := by
  have h := quad_matrices_norm williamson11_norm_checked
  norm_num at h ⊢
  exact h



def williamson13Row (i : Fin 4) (x : ℕ) : ℤ := treeAt (williamson13Trees i) x

theorem williamson13_sign_checked : checkQuadSign 13 williamson13Row = true := by decide +kernel

theorem williamson13_symmetric_checked : checkQuadSymmetric 13 williamson13Row = true := by decide +kernel

theorem williamson13_norm_checked : checkQuadCorr 13 williamson13Row 52 = true := by decide +kernel

def williamson13 : Fin 4 → Matrix (ZMod 13) (ZMod 13) ℚ := rowMatrices 13 williamson13Row

theorem williamson13_sign : ∀ i, IsSignMatrix (williamson13 i) :=
  quad_matrices_sign williamson13_sign_checked

theorem williamson13_symmetric : ∀ i, (williamson13 i).transpose = williamson13 i :=
  quad_matrices_symmetric williamson13_symmetric_checked

theorem williamson13_commute : ∀ i j, williamson13 i * williamson13 j = williamson13 j * williamson13 i :=
  quad_matrices_commute williamson13Row

theorem williamson13_amicable : ∀ i j, williamson13 i * (williamson13 j).transpose =
    williamson13 j * (williamson13 i).transpose := by
  intro i j
  rw [williamson13_symmetric i, williamson13_symmetric j]
  exact williamson13_commute i j

theorem williamson13_norm : (∑ i : Fin 4, williamson13 i * (williamson13 i).transpose) =
    (4 * (13 : ℚ)) • (1 : Matrix (ZMod 13) (ZMod 13) ℚ) := by
  have h := quad_matrices_norm williamson13_norm_checked
  norm_num at h ⊢
  exact h



def williamson15Row (i : Fin 4) (x : ℕ) : ℤ := treeAt (williamson15Trees i) x

theorem williamson15_sign_checked : checkQuadSign 15 williamson15Row = true := by decide +kernel

theorem williamson15_symmetric_checked : checkQuadSymmetric 15 williamson15Row = true := by decide +kernel

theorem williamson15_norm_checked : checkQuadCorr 15 williamson15Row 60 = true := by decide +kernel

def williamson15 : Fin 4 → Matrix (ZMod 15) (ZMod 15) ℚ := rowMatrices 15 williamson15Row

theorem williamson15_sign : ∀ i, IsSignMatrix (williamson15 i) :=
  quad_matrices_sign williamson15_sign_checked

theorem williamson15_symmetric : ∀ i, (williamson15 i).transpose = williamson15 i :=
  quad_matrices_symmetric williamson15_symmetric_checked

theorem williamson15_commute : ∀ i j, williamson15 i * williamson15 j = williamson15 j * williamson15 i :=
  quad_matrices_commute williamson15Row

theorem williamson15_amicable : ∀ i j, williamson15 i * (williamson15 j).transpose =
    williamson15 j * (williamson15 i).transpose := by
  intro i j
  rw [williamson15_symmetric i, williamson15_symmetric j]
  exact williamson15_commute i j

theorem williamson15_norm : (∑ i : Fin 4, williamson15 i * (williamson15 i).transpose) =
    (4 * (15 : ℚ)) • (1 : Matrix (ZMod 15) (ZMod 15) ℚ) := by
  have h := quad_matrices_norm williamson15_norm_checked
  norm_num at h ⊢
  exact h



def williamson17Row (i : Fin 4) (x : ℕ) : ℤ := treeAt (williamson17Trees i) x

theorem williamson17_sign_checked : checkQuadSign 17 williamson17Row = true := by decide +kernel

theorem williamson17_symmetric_checked : checkQuadSymmetric 17 williamson17Row = true := by decide +kernel

theorem williamson17_norm_checked : checkQuadCorr 17 williamson17Row 68 = true := by decide +kernel

def williamson17 : Fin 4 → Matrix (ZMod 17) (ZMod 17) ℚ := rowMatrices 17 williamson17Row

theorem williamson17_sign : ∀ i, IsSignMatrix (williamson17 i) :=
  quad_matrices_sign williamson17_sign_checked

theorem williamson17_symmetric : ∀ i, (williamson17 i).transpose = williamson17 i :=
  quad_matrices_symmetric williamson17_symmetric_checked

theorem williamson17_commute : ∀ i j, williamson17 i * williamson17 j = williamson17 j * williamson17 i :=
  quad_matrices_commute williamson17Row

theorem williamson17_amicable : ∀ i j, williamson17 i * (williamson17 j).transpose =
    williamson17 j * (williamson17 i).transpose := by
  intro i j
  rw [williamson17_symmetric i, williamson17_symmetric j]
  exact williamson17_commute i j

theorem williamson17_norm : (∑ i : Fin 4, williamson17 i * (williamson17 i).transpose) =
    (4 * (17 : ℚ)) • (1 : Matrix (ZMod 17) (ZMod 17) ℚ) := by
  have h := quad_matrices_norm williamson17_norm_checked
  norm_num at h ⊢
  exact h



#print axioms williamson1_norm
#print axioms williamson3_norm
#print axioms williamson5_norm
#print axioms williamson7_norm
#print axioms williamson9_norm
#print axioms williamson11_norm
#print axioms williamson13_norm
#print axioms williamson15_norm
#print axioms williamson17_norm
end ButsonTowers.Data

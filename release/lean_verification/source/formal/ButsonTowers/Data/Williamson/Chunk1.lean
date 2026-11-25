import ButsonTowers.Data.RowChecks
import ButsonTowers.Data.TreeLiterals

namespace ButsonTowers.Data
open scoped BigOperators
set_option maxRecDepth 1000000
set_option maxHeartbeats 0

def williamson19Row (i : Fin 4) (x : ℕ) : ℤ := treeAt (williamson19Trees i) x

theorem williamson19_sign_checked : checkQuadSign 19 williamson19Row = true := by decide +kernel

theorem williamson19_symmetric_checked : checkQuadSymmetric 19 williamson19Row = true := by decide +kernel

theorem williamson19_norm_checked : checkQuadCorr 19 williamson19Row 76 = true := by decide +kernel

def williamson19 : Fin 4 → Matrix (ZMod 19) (ZMod 19) ℚ := rowMatrices 19 williamson19Row

theorem williamson19_sign : ∀ i, IsSignMatrix (williamson19 i) :=
  quad_matrices_sign williamson19_sign_checked

theorem williamson19_symmetric : ∀ i, (williamson19 i).transpose = williamson19 i :=
  quad_matrices_symmetric williamson19_symmetric_checked

theorem williamson19_commute : ∀ i j, williamson19 i * williamson19 j = williamson19 j * williamson19 i :=
  quad_matrices_commute williamson19Row

theorem williamson19_amicable : ∀ i j, williamson19 i * (williamson19 j).transpose =
    williamson19 j * (williamson19 i).transpose := by
  intro i j
  rw [williamson19_symmetric i, williamson19_symmetric j]
  exact williamson19_commute i j

theorem williamson19_norm : (∑ i : Fin 4, williamson19 i * (williamson19 i).transpose) =
    (4 * (19 : ℚ)) • (1 : Matrix (ZMod 19) (ZMod 19) ℚ) := by
  have h := quad_matrices_norm williamson19_norm_checked
  norm_num at h ⊢
  exact h



def williamson21Row (i : Fin 4) (x : ℕ) : ℤ := treeAt (williamson21Trees i) x

theorem williamson21_sign_checked : checkQuadSign 21 williamson21Row = true := by decide +kernel

theorem williamson21_symmetric_checked : checkQuadSymmetric 21 williamson21Row = true := by decide +kernel

theorem williamson21_norm_checked : checkQuadCorr 21 williamson21Row 84 = true := by decide +kernel

def williamson21 : Fin 4 → Matrix (ZMod 21) (ZMod 21) ℚ := rowMatrices 21 williamson21Row

theorem williamson21_sign : ∀ i, IsSignMatrix (williamson21 i) :=
  quad_matrices_sign williamson21_sign_checked

theorem williamson21_symmetric : ∀ i, (williamson21 i).transpose = williamson21 i :=
  quad_matrices_symmetric williamson21_symmetric_checked

theorem williamson21_commute : ∀ i j, williamson21 i * williamson21 j = williamson21 j * williamson21 i :=
  quad_matrices_commute williamson21Row

theorem williamson21_amicable : ∀ i j, williamson21 i * (williamson21 j).transpose =
    williamson21 j * (williamson21 i).transpose := by
  intro i j
  rw [williamson21_symmetric i, williamson21_symmetric j]
  exact williamson21_commute i j

theorem williamson21_norm : (∑ i : Fin 4, williamson21 i * (williamson21 i).transpose) =
    (4 * (21 : ℚ)) • (1 : Matrix (ZMod 21) (ZMod 21) ℚ) := by
  have h := quad_matrices_norm williamson21_norm_checked
  norm_num at h ⊢
  exact h



def williamson23Row (i : Fin 4) (x : ℕ) : ℤ := treeAt (williamson23Trees i) x

theorem williamson23_sign_checked : checkQuadSign 23 williamson23Row = true := by decide +kernel

theorem williamson23_symmetric_checked : checkQuadSymmetric 23 williamson23Row = true := by decide +kernel

theorem williamson23_norm_checked : checkQuadCorr 23 williamson23Row 92 = true := by decide +kernel

def williamson23 : Fin 4 → Matrix (ZMod 23) (ZMod 23) ℚ := rowMatrices 23 williamson23Row

theorem williamson23_sign : ∀ i, IsSignMatrix (williamson23 i) :=
  quad_matrices_sign williamson23_sign_checked

theorem williamson23_symmetric : ∀ i, (williamson23 i).transpose = williamson23 i :=
  quad_matrices_symmetric williamson23_symmetric_checked

theorem williamson23_commute : ∀ i j, williamson23 i * williamson23 j = williamson23 j * williamson23 i :=
  quad_matrices_commute williamson23Row

theorem williamson23_amicable : ∀ i j, williamson23 i * (williamson23 j).transpose =
    williamson23 j * (williamson23 i).transpose := by
  intro i j
  rw [williamson23_symmetric i, williamson23_symmetric j]
  exact williamson23_commute i j

theorem williamson23_norm : (∑ i : Fin 4, williamson23 i * (williamson23 i).transpose) =
    (4 * (23 : ℚ)) • (1 : Matrix (ZMod 23) (ZMod 23) ℚ) := by
  have h := quad_matrices_norm williamson23_norm_checked
  norm_num at h ⊢
  exact h



def williamson25Row (i : Fin 4) (x : ℕ) : ℤ := treeAt (williamson25Trees i) x

theorem williamson25_sign_checked : checkQuadSign 25 williamson25Row = true := by decide +kernel

theorem williamson25_symmetric_checked : checkQuadSymmetric 25 williamson25Row = true := by decide +kernel

theorem williamson25_norm_checked : checkQuadCorr 25 williamson25Row 100 = true := by decide +kernel

def williamson25 : Fin 4 → Matrix (ZMod 25) (ZMod 25) ℚ := rowMatrices 25 williamson25Row

theorem williamson25_sign : ∀ i, IsSignMatrix (williamson25 i) :=
  quad_matrices_sign williamson25_sign_checked

theorem williamson25_symmetric : ∀ i, (williamson25 i).transpose = williamson25 i :=
  quad_matrices_symmetric williamson25_symmetric_checked

theorem williamson25_commute : ∀ i j, williamson25 i * williamson25 j = williamson25 j * williamson25 i :=
  quad_matrices_commute williamson25Row

theorem williamson25_amicable : ∀ i j, williamson25 i * (williamson25 j).transpose =
    williamson25 j * (williamson25 i).transpose := by
  intro i j
  rw [williamson25_symmetric i, williamson25_symmetric j]
  exact williamson25_commute i j

theorem williamson25_norm : (∑ i : Fin 4, williamson25 i * (williamson25 i).transpose) =
    (4 * (25 : ℚ)) • (1 : Matrix (ZMod 25) (ZMod 25) ℚ) := by
  have h := quad_matrices_norm williamson25_norm_checked
  norm_num at h ⊢
  exact h



def williamson27Row (i : Fin 4) (x : ℕ) : ℤ := treeAt (williamson27Trees i) x

theorem williamson27_sign_checked : checkQuadSign 27 williamson27Row = true := by decide +kernel

theorem williamson27_symmetric_checked : checkQuadSymmetric 27 williamson27Row = true := by decide +kernel

theorem williamson27_norm_checked : checkQuadCorr 27 williamson27Row 108 = true := by decide +kernel

def williamson27 : Fin 4 → Matrix (ZMod 27) (ZMod 27) ℚ := rowMatrices 27 williamson27Row

theorem williamson27_sign : ∀ i, IsSignMatrix (williamson27 i) :=
  quad_matrices_sign williamson27_sign_checked

theorem williamson27_symmetric : ∀ i, (williamson27 i).transpose = williamson27 i :=
  quad_matrices_symmetric williamson27_symmetric_checked

theorem williamson27_commute : ∀ i j, williamson27 i * williamson27 j = williamson27 j * williamson27 i :=
  quad_matrices_commute williamson27Row

theorem williamson27_amicable : ∀ i j, williamson27 i * (williamson27 j).transpose =
    williamson27 j * (williamson27 i).transpose := by
  intro i j
  rw [williamson27_symmetric i, williamson27_symmetric j]
  exact williamson27_commute i j

theorem williamson27_norm : (∑ i : Fin 4, williamson27 i * (williamson27 i).transpose) =
    (4 * (27 : ℚ)) • (1 : Matrix (ZMod 27) (ZMod 27) ℚ) := by
  have h := quad_matrices_norm williamson27_norm_checked
  norm_num at h ⊢
  exact h



def williamson29Row (i : Fin 4) (x : ℕ) : ℤ := treeAt (williamson29Trees i) x

theorem williamson29_sign_checked : checkQuadSign 29 williamson29Row = true := by decide +kernel

theorem williamson29_symmetric_checked : checkQuadSymmetric 29 williamson29Row = true := by decide +kernel

theorem williamson29_norm_checked : checkQuadCorr 29 williamson29Row 116 = true := by decide +kernel

def williamson29 : Fin 4 → Matrix (ZMod 29) (ZMod 29) ℚ := rowMatrices 29 williamson29Row

theorem williamson29_sign : ∀ i, IsSignMatrix (williamson29 i) :=
  quad_matrices_sign williamson29_sign_checked

theorem williamson29_symmetric : ∀ i, (williamson29 i).transpose = williamson29 i :=
  quad_matrices_symmetric williamson29_symmetric_checked

theorem williamson29_commute : ∀ i j, williamson29 i * williamson29 j = williamson29 j * williamson29 i :=
  quad_matrices_commute williamson29Row

theorem williamson29_amicable : ∀ i j, williamson29 i * (williamson29 j).transpose =
    williamson29 j * (williamson29 i).transpose := by
  intro i j
  rw [williamson29_symmetric i, williamson29_symmetric j]
  exact williamson29_commute i j

theorem williamson29_norm : (∑ i : Fin 4, williamson29 i * (williamson29 i).transpose) =
    (4 * (29 : ℚ)) • (1 : Matrix (ZMod 29) (ZMod 29) ℚ) := by
  have h := quad_matrices_norm williamson29_norm_checked
  norm_num at h ⊢
  exact h



def williamson31Row (i : Fin 4) (x : ℕ) : ℤ := treeAt (williamson31Trees i) x

theorem williamson31_sign_checked : checkQuadSign 31 williamson31Row = true := by decide +kernel

theorem williamson31_symmetric_checked : checkQuadSymmetric 31 williamson31Row = true := by decide +kernel

theorem williamson31_norm_checked : checkQuadCorr 31 williamson31Row 124 = true := by decide +kernel

def williamson31 : Fin 4 → Matrix (ZMod 31) (ZMod 31) ℚ := rowMatrices 31 williamson31Row

theorem williamson31_sign : ∀ i, IsSignMatrix (williamson31 i) :=
  quad_matrices_sign williamson31_sign_checked

theorem williamson31_symmetric : ∀ i, (williamson31 i).transpose = williamson31 i :=
  quad_matrices_symmetric williamson31_symmetric_checked

theorem williamson31_commute : ∀ i j, williamson31 i * williamson31 j = williamson31 j * williamson31 i :=
  quad_matrices_commute williamson31Row

theorem williamson31_amicable : ∀ i j, williamson31 i * (williamson31 j).transpose =
    williamson31 j * (williamson31 i).transpose := by
  intro i j
  rw [williamson31_symmetric i, williamson31_symmetric j]
  exact williamson31_commute i j

theorem williamson31_norm : (∑ i : Fin 4, williamson31 i * (williamson31 i).transpose) =
    (4 * (31 : ℚ)) • (1 : Matrix (ZMod 31) (ZMod 31) ℚ) := by
  have h := quad_matrices_norm williamson31_norm_checked
  norm_num at h ⊢
  exact h



def williamson33Row (i : Fin 4) (x : ℕ) : ℤ := treeAt (williamson33Trees i) x

theorem williamson33_sign_checked : checkQuadSign 33 williamson33Row = true := by decide +kernel

theorem williamson33_symmetric_checked : checkQuadSymmetric 33 williamson33Row = true := by decide +kernel

theorem williamson33_norm_checked : checkQuadCorr 33 williamson33Row 132 = true := by decide +kernel

def williamson33 : Fin 4 → Matrix (ZMod 33) (ZMod 33) ℚ := rowMatrices 33 williamson33Row

theorem williamson33_sign : ∀ i, IsSignMatrix (williamson33 i) :=
  quad_matrices_sign williamson33_sign_checked

theorem williamson33_symmetric : ∀ i, (williamson33 i).transpose = williamson33 i :=
  quad_matrices_symmetric williamson33_symmetric_checked

theorem williamson33_commute : ∀ i j, williamson33 i * williamson33 j = williamson33 j * williamson33 i :=
  quad_matrices_commute williamson33Row

theorem williamson33_amicable : ∀ i j, williamson33 i * (williamson33 j).transpose =
    williamson33 j * (williamson33 i).transpose := by
  intro i j
  rw [williamson33_symmetric i, williamson33_symmetric j]
  exact williamson33_commute i j

theorem williamson33_norm : (∑ i : Fin 4, williamson33 i * (williamson33 i).transpose) =
    (4 * (33 : ℚ)) • (1 : Matrix (ZMod 33) (ZMod 33) ℚ) := by
  have h := quad_matrices_norm williamson33_norm_checked
  norm_num at h ⊢
  exact h



def williamson37Row (i : Fin 4) (x : ℕ) : ℤ := treeAt (williamson37Trees i) x

theorem williamson37_sign_checked : checkQuadSign 37 williamson37Row = true := by decide +kernel

theorem williamson37_symmetric_checked : checkQuadSymmetric 37 williamson37Row = true := by decide +kernel

theorem williamson37_norm_checked : checkQuadCorr 37 williamson37Row 148 = true := by decide +kernel

def williamson37 : Fin 4 → Matrix (ZMod 37) (ZMod 37) ℚ := rowMatrices 37 williamson37Row

theorem williamson37_sign : ∀ i, IsSignMatrix (williamson37 i) :=
  quad_matrices_sign williamson37_sign_checked

theorem williamson37_symmetric : ∀ i, (williamson37 i).transpose = williamson37 i :=
  quad_matrices_symmetric williamson37_symmetric_checked

theorem williamson37_commute : ∀ i j, williamson37 i * williamson37 j = williamson37 j * williamson37 i :=
  quad_matrices_commute williamson37Row

theorem williamson37_amicable : ∀ i j, williamson37 i * (williamson37 j).transpose =
    williamson37 j * (williamson37 i).transpose := by
  intro i j
  rw [williamson37_symmetric i, williamson37_symmetric j]
  exact williamson37_commute i j

theorem williamson37_norm : (∑ i : Fin 4, williamson37 i * (williamson37 i).transpose) =
    (4 * (37 : ℚ)) • (1 : Matrix (ZMod 37) (ZMod 37) ℚ) := by
  have h := quad_matrices_norm williamson37_norm_checked
  norm_num at h ⊢
  exact h



#print axioms williamson19_norm
#print axioms williamson21_norm
#print axioms williamson23_norm
#print axioms williamson25_norm
#print axioms williamson27_norm
#print axioms williamson29_norm
#print axioms williamson31_norm
#print axioms williamson33_norm
#print axioms williamson37_norm
end ButsonTowers.Data

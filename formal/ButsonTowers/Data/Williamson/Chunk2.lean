import ButsonTowers.Data.RowChecks
import ButsonTowers.Data.TreeLiterals

namespace ButsonTowers.Data
open scoped BigOperators
set_option maxRecDepth 1000000
set_option maxHeartbeats 0

def williamson39Row (i : Fin 4) (x : ℕ) : ℤ := treeAt (williamson39Trees i) x

theorem williamson39_sign_checked : checkQuadSign 39 williamson39Row = true := by decide +kernel

theorem williamson39_symmetric_checked : checkQuadSymmetric 39 williamson39Row = true := by decide +kernel

theorem williamson39_norm_checked : checkQuadCorr 39 williamson39Row 156 = true := by decide +kernel

def williamson39 : Fin 4 → Matrix (ZMod 39) (ZMod 39) ℚ := rowMatrices 39 williamson39Row

theorem williamson39_sign : ∀ i, IsSignMatrix (williamson39 i) :=
  quad_matrices_sign williamson39_sign_checked

theorem williamson39_symmetric : ∀ i, (williamson39 i).transpose = williamson39 i :=
  quad_matrices_symmetric williamson39_symmetric_checked

theorem williamson39_commute : ∀ i j, williamson39 i * williamson39 j = williamson39 j * williamson39 i :=
  quad_matrices_commute williamson39Row

theorem williamson39_amicable : ∀ i j, williamson39 i * (williamson39 j).transpose =
    williamson39 j * (williamson39 i).transpose := by
  intro i j
  rw [williamson39_symmetric i, williamson39_symmetric j]
  exact williamson39_commute i j

theorem williamson39_norm : (∑ i : Fin 4, williamson39 i * (williamson39 i).transpose) =
    (4 * (39 : ℚ)) • (1 : Matrix (ZMod 39) (ZMod 39) ℚ) := by
  have h := quad_matrices_norm williamson39_norm_checked
  norm_num at h ⊢
  exact h



def williamson41Row (i : Fin 4) (x : ℕ) : ℤ := treeAt (williamson41Trees i) x

theorem williamson41_sign_checked : checkQuadSign 41 williamson41Row = true := by decide +kernel

theorem williamson41_symmetric_checked : checkQuadSymmetric 41 williamson41Row = true := by decide +kernel

theorem williamson41_norm_checked : checkQuadCorr 41 williamson41Row 164 = true := by decide +kernel

def williamson41 : Fin 4 → Matrix (ZMod 41) (ZMod 41) ℚ := rowMatrices 41 williamson41Row

theorem williamson41_sign : ∀ i, IsSignMatrix (williamson41 i) :=
  quad_matrices_sign williamson41_sign_checked

theorem williamson41_symmetric : ∀ i, (williamson41 i).transpose = williamson41 i :=
  quad_matrices_symmetric williamson41_symmetric_checked

theorem williamson41_commute : ∀ i j, williamson41 i * williamson41 j = williamson41 j * williamson41 i :=
  quad_matrices_commute williamson41Row

theorem williamson41_amicable : ∀ i j, williamson41 i * (williamson41 j).transpose =
    williamson41 j * (williamson41 i).transpose := by
  intro i j
  rw [williamson41_symmetric i, williamson41_symmetric j]
  exact williamson41_commute i j

theorem williamson41_norm : (∑ i : Fin 4, williamson41 i * (williamson41 i).transpose) =
    (4 * (41 : ℚ)) • (1 : Matrix (ZMod 41) (ZMod 41) ℚ) := by
  have h := quad_matrices_norm williamson41_norm_checked
  norm_num at h ⊢
  exact h



def williamson43Row (i : Fin 4) (x : ℕ) : ℤ := treeAt (williamson43Trees i) x

theorem williamson43_sign_checked : checkQuadSign 43 williamson43Row = true := by decide +kernel

theorem williamson43_symmetric_checked : checkQuadSymmetric 43 williamson43Row = true := by decide +kernel

theorem williamson43_norm_checked : checkQuadCorr 43 williamson43Row 172 = true := by decide +kernel

def williamson43 : Fin 4 → Matrix (ZMod 43) (ZMod 43) ℚ := rowMatrices 43 williamson43Row

theorem williamson43_sign : ∀ i, IsSignMatrix (williamson43 i) :=
  quad_matrices_sign williamson43_sign_checked

theorem williamson43_symmetric : ∀ i, (williamson43 i).transpose = williamson43 i :=
  quad_matrices_symmetric williamson43_symmetric_checked

theorem williamson43_commute : ∀ i j, williamson43 i * williamson43 j = williamson43 j * williamson43 i :=
  quad_matrices_commute williamson43Row

theorem williamson43_amicable : ∀ i j, williamson43 i * (williamson43 j).transpose =
    williamson43 j * (williamson43 i).transpose := by
  intro i j
  rw [williamson43_symmetric i, williamson43_symmetric j]
  exact williamson43_commute i j

theorem williamson43_norm : (∑ i : Fin 4, williamson43 i * (williamson43 i).transpose) =
    (4 * (43 : ℚ)) • (1 : Matrix (ZMod 43) (ZMod 43) ℚ) := by
  have h := quad_matrices_norm williamson43_norm_checked
  norm_num at h ⊢
  exact h



def williamson49Row (i : Fin 4) (x : ℕ) : ℤ := treeAt (williamson49Trees i) x

theorem williamson49_sign_checked : checkQuadSign 49 williamson49Row = true := by decide +kernel

theorem williamson49_symmetric_checked : checkQuadSymmetric 49 williamson49Row = true := by decide +kernel

theorem williamson49_norm_checked : checkQuadCorr 49 williamson49Row 196 = true := by decide +kernel

def williamson49 : Fin 4 → Matrix (ZMod 49) (ZMod 49) ℚ := rowMatrices 49 williamson49Row

theorem williamson49_sign : ∀ i, IsSignMatrix (williamson49 i) :=
  quad_matrices_sign williamson49_sign_checked

theorem williamson49_symmetric : ∀ i, (williamson49 i).transpose = williamson49 i :=
  quad_matrices_symmetric williamson49_symmetric_checked

theorem williamson49_commute : ∀ i j, williamson49 i * williamson49 j = williamson49 j * williamson49 i :=
  quad_matrices_commute williamson49Row

theorem williamson49_amicable : ∀ i j, williamson49 i * (williamson49 j).transpose =
    williamson49 j * (williamson49 i).transpose := by
  intro i j
  rw [williamson49_symmetric i, williamson49_symmetric j]
  exact williamson49_commute i j

theorem williamson49_norm : (∑ i : Fin 4, williamson49 i * (williamson49 i).transpose) =
    (4 * (49 : ℚ)) • (1 : Matrix (ZMod 49) (ZMod 49) ℚ) := by
  have h := quad_matrices_norm williamson49_norm_checked
  norm_num at h ⊢
  exact h



def williamson51Row (i : Fin 4) (x : ℕ) : ℤ := treeAt (williamson51Trees i) x

theorem williamson51_sign_checked : checkQuadSign 51 williamson51Row = true := by decide +kernel

theorem williamson51_symmetric_checked : checkQuadSymmetric 51 williamson51Row = true := by decide +kernel

theorem williamson51_norm_checked : checkQuadCorr 51 williamson51Row 204 = true := by decide +kernel

def williamson51 : Fin 4 → Matrix (ZMod 51) (ZMod 51) ℚ := rowMatrices 51 williamson51Row

theorem williamson51_sign : ∀ i, IsSignMatrix (williamson51 i) :=
  quad_matrices_sign williamson51_sign_checked

theorem williamson51_symmetric : ∀ i, (williamson51 i).transpose = williamson51 i :=
  quad_matrices_symmetric williamson51_symmetric_checked

theorem williamson51_commute : ∀ i j, williamson51 i * williamson51 j = williamson51 j * williamson51 i :=
  quad_matrices_commute williamson51Row

theorem williamson51_amicable : ∀ i j, williamson51 i * (williamson51 j).transpose =
    williamson51 j * (williamson51 i).transpose := by
  intro i j
  rw [williamson51_symmetric i, williamson51_symmetric j]
  exact williamson51_commute i j

theorem williamson51_norm : (∑ i : Fin 4, williamson51 i * (williamson51 i).transpose) =
    (4 * (51 : ℚ)) • (1 : Matrix (ZMod 51) (ZMod 51) ℚ) := by
  have h := quad_matrices_norm williamson51_norm_checked
  norm_num at h ⊢
  exact h



def williamson55Row (i : Fin 4) (x : ℕ) : ℤ := treeAt (williamson55Trees i) x

theorem williamson55_sign_checked : checkQuadSign 55 williamson55Row = true := by decide +kernel

theorem williamson55_symmetric_checked : checkQuadSymmetric 55 williamson55Row = true := by decide +kernel

theorem williamson55_norm_checked : checkQuadCorr 55 williamson55Row 220 = true := by decide +kernel

def williamson55 : Fin 4 → Matrix (ZMod 55) (ZMod 55) ℚ := rowMatrices 55 williamson55Row

theorem williamson55_sign : ∀ i, IsSignMatrix (williamson55 i) :=
  quad_matrices_sign williamson55_sign_checked

theorem williamson55_symmetric : ∀ i, (williamson55 i).transpose = williamson55 i :=
  quad_matrices_symmetric williamson55_symmetric_checked

theorem williamson55_commute : ∀ i j, williamson55 i * williamson55 j = williamson55 j * williamson55 i :=
  quad_matrices_commute williamson55Row

theorem williamson55_amicable : ∀ i j, williamson55 i * (williamson55 j).transpose =
    williamson55 j * (williamson55 i).transpose := by
  intro i j
  rw [williamson55_symmetric i, williamson55_symmetric j]
  exact williamson55_commute i j

theorem williamson55_norm : (∑ i : Fin 4, williamson55 i * (williamson55 i).transpose) =
    (4 * (55 : ℚ)) • (1 : Matrix (ZMod 55) (ZMod 55) ℚ) := by
  have h := quad_matrices_norm williamson55_norm_checked
  norm_num at h ⊢
  exact h



def williamson57Row (i : Fin 4) (x : ℕ) : ℤ := treeAt (williamson57Trees i) x

theorem williamson57_sign_checked : checkQuadSign 57 williamson57Row = true := by decide +kernel

theorem williamson57_symmetric_checked : checkQuadSymmetric 57 williamson57Row = true := by decide +kernel

theorem williamson57_norm_checked : checkQuadCorr 57 williamson57Row 228 = true := by decide +kernel

def williamson57 : Fin 4 → Matrix (ZMod 57) (ZMod 57) ℚ := rowMatrices 57 williamson57Row

theorem williamson57_sign : ∀ i, IsSignMatrix (williamson57 i) :=
  quad_matrices_sign williamson57_sign_checked

theorem williamson57_symmetric : ∀ i, (williamson57 i).transpose = williamson57 i :=
  quad_matrices_symmetric williamson57_symmetric_checked

theorem williamson57_commute : ∀ i j, williamson57 i * williamson57 j = williamson57 j * williamson57 i :=
  quad_matrices_commute williamson57Row

theorem williamson57_amicable : ∀ i j, williamson57 i * (williamson57 j).transpose =
    williamson57 j * (williamson57 i).transpose := by
  intro i j
  rw [williamson57_symmetric i, williamson57_symmetric j]
  exact williamson57_commute i j

theorem williamson57_norm : (∑ i : Fin 4, williamson57 i * (williamson57 i).transpose) =
    (4 * (57 : ℚ)) • (1 : Matrix (ZMod 57) (ZMod 57) ℚ) := by
  have h := quad_matrices_norm williamson57_norm_checked
  norm_num at h ⊢
  exact h



def williamson61Row (i : Fin 4) (x : ℕ) : ℤ := treeAt (williamson61Trees i) x

theorem williamson61_sign_checked : checkQuadSign 61 williamson61Row = true := by decide +kernel

theorem williamson61_symmetric_checked : checkQuadSymmetric 61 williamson61Row = true := by decide +kernel

theorem williamson61_norm_checked : checkQuadCorr 61 williamson61Row 244 = true := by decide +kernel

def williamson61 : Fin 4 → Matrix (ZMod 61) (ZMod 61) ℚ := rowMatrices 61 williamson61Row

theorem williamson61_sign : ∀ i, IsSignMatrix (williamson61 i) :=
  quad_matrices_sign williamson61_sign_checked

theorem williamson61_symmetric : ∀ i, (williamson61 i).transpose = williamson61 i :=
  quad_matrices_symmetric williamson61_symmetric_checked

theorem williamson61_commute : ∀ i j, williamson61 i * williamson61 j = williamson61 j * williamson61 i :=
  quad_matrices_commute williamson61Row

theorem williamson61_amicable : ∀ i j, williamson61 i * (williamson61 j).transpose =
    williamson61 j * (williamson61 i).transpose := by
  intro i j
  rw [williamson61_symmetric i, williamson61_symmetric j]
  exact williamson61_commute i j

theorem williamson61_norm : (∑ i : Fin 4, williamson61 i * (williamson61 i).transpose) =
    (4 * (61 : ℚ)) • (1 : Matrix (ZMod 61) (ZMod 61) ℚ) := by
  have h := quad_matrices_norm williamson61_norm_checked
  norm_num at h ⊢
  exact h



def williamson63Row (i : Fin 4) (x : ℕ) : ℤ := treeAt (williamson63Trees i) x

theorem williamson63_sign_checked : checkQuadSign 63 williamson63Row = true := by decide +kernel

theorem williamson63_symmetric_checked : checkQuadSymmetric 63 williamson63Row = true := by decide +kernel

theorem williamson63_norm_checked : checkQuadCorr 63 williamson63Row 252 = true := by decide +kernel

def williamson63 : Fin 4 → Matrix (ZMod 63) (ZMod 63) ℚ := rowMatrices 63 williamson63Row

theorem williamson63_sign : ∀ i, IsSignMatrix (williamson63 i) :=
  quad_matrices_sign williamson63_sign_checked

theorem williamson63_symmetric : ∀ i, (williamson63 i).transpose = williamson63 i :=
  quad_matrices_symmetric williamson63_symmetric_checked

theorem williamson63_commute : ∀ i j, williamson63 i * williamson63 j = williamson63 j * williamson63 i :=
  quad_matrices_commute williamson63Row

theorem williamson63_amicable : ∀ i j, williamson63 i * (williamson63 j).transpose =
    williamson63 j * (williamson63 i).transpose := by
  intro i j
  rw [williamson63_symmetric i, williamson63_symmetric j]
  exact williamson63_commute i j

theorem williamson63_norm : (∑ i : Fin 4, williamson63 i * (williamson63 i).transpose) =
    (4 * (63 : ℚ)) • (1 : Matrix (ZMod 63) (ZMod 63) ℚ) := by
  have h := quad_matrices_norm williamson63_norm_checked
  norm_num at h ⊢
  exact h



#print axioms williamson39_norm
#print axioms williamson41_norm
#print axioms williamson43_norm
#print axioms williamson49_norm
#print axioms williamson51_norm
#print axioms williamson55_norm
#print axioms williamson57_norm
#print axioms williamson61_norm
#print axioms williamson63_norm
end ButsonTowers.Data

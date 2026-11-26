import ButsonTowers.Data.PaleyPairChecks
import ButsonTowers.Data.TreeLiterals

namespace ButsonTowers.Data

set_option maxRecDepth 1000000
set_option maxHeartbeats 0

def paley13Row : ℕ → ℤ := treeAt paley13Tree

theorem paley13_sign_checked : checkPaleySign 13 paley13Row = true := by decide +kernel
theorem paley13_symmetric_checked : checkPaleySymmetric 13 paley13Row = true := by decide +kernel
theorem paley13_sum_checked : natSum 13 paley13Row = 0 := by decide +kernel
theorem paley13_corr_checked : checkPaleyCorr 13 paley13Row = true := by decide +kernel

def paley13 : Matrix (ZMod 13) (ZMod 13) ℚ :=
  GroupAlgebra.regular (fun x => (paley13Row x.val : ℚ))

theorem paley13_isPaley : SpenceMatrix.IsPaley paley13 :=
  prime_paley_checked paley13Row paley13_sign_checked paley13_symmetric_checked
    paley13_sum_checked paley13_corr_checked

def paley17Row : ℕ → ℤ := treeAt paley17Tree

theorem paley17_sign_checked : checkPaleySign 17 paley17Row = true := by decide +kernel
theorem paley17_symmetric_checked : checkPaleySymmetric 17 paley17Row = true := by decide +kernel
theorem paley17_sum_checked : natSum 17 paley17Row = 0 := by decide +kernel
theorem paley17_corr_checked : checkPaleyCorr 17 paley17Row = true := by decide +kernel

def paley17 : Matrix (ZMod 17) (ZMod 17) ℚ :=
  GroupAlgebra.regular (fun x => (paley17Row x.val : ℚ))

theorem paley17_isPaley : SpenceMatrix.IsPaley paley17 :=
  prime_paley_checked paley17Row paley17_sign_checked paley17_symmetric_checked
    paley17_sum_checked paley17_corr_checked

def paley61Row : ℕ → ℤ := treeAt paley61Tree

theorem paley61_sign_checked : checkPaleySign 61 paley61Row = true := by decide +kernel
theorem paley61_symmetric_checked : checkPaleySymmetric 61 paley61Row = true := by decide +kernel
theorem paley61_sum_checked : natSum 61 paley61Row = 0 := by decide +kernel
theorem paley61_corr_checked : checkPaleyCorr 61 paley61Row = true := by decide +kernel

def paley61 : Matrix (ZMod 61) (ZMod 61) ℚ :=
  GroupAlgebra.regular (fun x => (paley61Row x.val : ℚ))

theorem paley61_isPaley : SpenceMatrix.IsPaley paley61 :=
  prime_paley_checked paley61Row paley61_sign_checked paley61_symmetric_checked
    paley61_sum_checked paley61_corr_checked

def paley101Row : ℕ → ℤ := treeAt paley101Tree

theorem paley101_sign_checked : checkPaleySign 101 paley101Row = true := by decide +kernel
theorem paley101_symmetric_checked : checkPaleySymmetric 101 paley101Row = true := by decide +kernel
theorem paley101_sum_checked : natSum 101 paley101Row = 0 := by decide +kernel
theorem paley101_corr_checked : checkPaleyCorr 101 paley101Row = true := by decide +kernel

def paley101 : Matrix (ZMod 101) (ZMod 101) ℚ :=
  GroupAlgebra.regular (fun x => (paley101Row x.val : ℚ))

theorem paley101_isPaley : SpenceMatrix.IsPaley paley101 :=
  prime_paley_checked paley101Row paley101_sign_checked paley101_symmetric_checked
    paley101_sum_checked paley101_corr_checked

def paley121Row : ℕ → ℤ := treeAt paley121Tree

theorem paley121_sign_checked : checkPaleySign 121 paley121Row = true := by decide +kernel
theorem paley121_symmetric_checked : checkPairPaleySymmetric paley121Row = true := by decide +kernel
theorem paley121_sum_checked : natSum 121 paley121Row = 0 := by decide +kernel
theorem paley121_corr_checked : checkPairPaleyCorr paley121Row = true := by decide +kernel

def paley121 : Matrix Paley121Group Paley121Group ℚ :=
  GroupAlgebra.regular (fun x => (paley121Row (pairEncode x) : ℚ))

theorem paley121_isPaley : SpenceMatrix.IsPaley paley121 :=
  pair_paley_checked paley121Row paley121_sign_checked paley121_symmetric_checked
    paley121_sum_checked paley121_corr_checked

end ButsonTowers.Data


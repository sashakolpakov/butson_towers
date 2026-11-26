import ButsonTowers.Data.CorrRange
import ButsonTowers.Data.TreeLiterals
namespace ButsonTowers.Data
set_option maxRecDepth 1000000
set_option maxHeartbeats 0

theorem t631_corr_block_00 : checkCorrRange 631 t631CodesFast 0 20 = true := by decide +kernel

theorem t631_corr_block_01 : checkCorrRange 631 t631CodesFast 20 20 = true := by decide +kernel

theorem t631_corr_block_02 : checkCorrRange 631 t631CodesFast 40 20 = true := by decide +kernel

theorem t631_corr_block_03 : checkCorrRange 631 t631CodesFast 60 20 = true := by decide +kernel

theorem t631_corr_block_04 : checkCorrRange 631 t631CodesFast 80 20 = true := by decide +kernel

theorem t631_corr_block_05 : checkCorrRange 631 t631CodesFast 100 20 = true := by decide +kernel

theorem t631_corr_block_06 : checkCorrRange 631 t631CodesFast 120 20 = true := by decide +kernel

theorem t631_corr_block_07 : checkCorrRange 631 t631CodesFast 140 20 = true := by decide +kernel

#print axioms t631_corr_block_07
end ButsonTowers.Data

import ButsonTowers.Data.CorrRange
import ButsonTowers.Data.TreeLiterals
namespace ButsonTowers.Data
set_option maxRecDepth 1000000
set_option maxHeartbeats 0

theorem t631_corr_block_16 : checkCorrRange 631 t631CodesFast 320 20 = true := by decide +kernel

theorem t631_corr_block_17 : checkCorrRange 631 t631CodesFast 340 20 = true := by decide +kernel

theorem t631_corr_block_18 : checkCorrRange 631 t631CodesFast 360 20 = true := by decide +kernel

theorem t631_corr_block_19 : checkCorrRange 631 t631CodesFast 380 20 = true := by decide +kernel

theorem t631_corr_block_20 : checkCorrRange 631 t631CodesFast 400 20 = true := by decide +kernel

theorem t631_corr_block_21 : checkCorrRange 631 t631CodesFast 420 20 = true := by decide +kernel

theorem t631_corr_block_22 : checkCorrRange 631 t631CodesFast 440 20 = true := by decide +kernel

theorem t631_corr_block_23 : checkCorrRange 631 t631CodesFast 460 20 = true := by decide +kernel

#print axioms t631_corr_block_23
end ButsonTowers.Data

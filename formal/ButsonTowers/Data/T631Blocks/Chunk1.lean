import ButsonTowers.Data.CorrRange
import ButsonTowers.Data.TreeLiterals
namespace ButsonTowers.Data
set_option maxRecDepth 1000000
set_option maxHeartbeats 0

theorem t631_corr_block_08 : checkCorrRange 631 t631CodesFast 160 20 = true := by decide +kernel

theorem t631_corr_block_09 : checkCorrRange 631 t631CodesFast 180 20 = true := by decide +kernel

theorem t631_corr_block_10 : checkCorrRange 631 t631CodesFast 200 20 = true := by decide +kernel

theorem t631_corr_block_11 : checkCorrRange 631 t631CodesFast 220 20 = true := by decide +kernel

theorem t631_corr_block_12 : checkCorrRange 631 t631CodesFast 240 20 = true := by decide +kernel

theorem t631_corr_block_13 : checkCorrRange 631 t631CodesFast 260 20 = true := by decide +kernel

theorem t631_corr_block_14 : checkCorrRange 631 t631CodesFast 280 20 = true := by decide +kernel

theorem t631_corr_block_15 : checkCorrRange 631 t631CodesFast 300 20 = true := by decide +kernel

#print axioms t631_corr_block_15
end ButsonTowers.Data

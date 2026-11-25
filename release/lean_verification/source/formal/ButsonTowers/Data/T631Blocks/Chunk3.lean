import ButsonTowers.Data.CorrRange
import ButsonTowers.Data.TreeLiterals
namespace ButsonTowers.Data
set_option maxRecDepth 1000000
set_option maxHeartbeats 0

theorem t631_corr_block_24 : checkCorrRange 631 t631CodesFast 480 20 = true := by decide +kernel

theorem t631_corr_block_25 : checkCorrRange 631 t631CodesFast 500 20 = true := by decide +kernel

theorem t631_corr_block_26 : checkCorrRange 631 t631CodesFast 520 20 = true := by decide +kernel

theorem t631_corr_block_27 : checkCorrRange 631 t631CodesFast 540 20 = true := by decide +kernel

theorem t631_corr_block_28 : checkCorrRange 631 t631CodesFast 560 20 = true := by decide +kernel

theorem t631_corr_block_29 : checkCorrRange 631 t631CodesFast 580 20 = true := by decide +kernel

theorem t631_corr_block_30 : checkCorrRange 631 t631CodesFast 600 20 = true := by decide +kernel

theorem t631_corr_block_31 : checkCorrRange 631 t631CodesFast 620 11 = true := by decide +kernel

#print axioms t631_corr_block_31
end ButsonTowers.Data

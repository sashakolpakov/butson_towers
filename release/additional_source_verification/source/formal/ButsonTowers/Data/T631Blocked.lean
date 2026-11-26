import ButsonTowers.Data.T631Blocks.Chunk0
import ButsonTowers.Data.T631Blocks.Chunk1
import ButsonTowers.Data.T631Blocks.Chunk2
import ButsonTowers.Data.T631Blocks.Chunk3

namespace ButsonTowers.Data
set_option maxHeartbeats 0
theorem t631_codes_corr_blocked : checkCodeCorr 631 t631CodesFast = true := by
  apply @code_corr_of_pointwise 631 t631CodesFast
  intro s hs631

  by_cases h0 : s < 20
  · exact @corr_range_sound 631 0 20 t631CodesFast t631_corr_block_00 s (by omega) (by omega)

  by_cases h1 : s < 40
  · exact @corr_range_sound 631 20 20 t631CodesFast t631_corr_block_01 s (by omega) (by omega)

  by_cases h2 : s < 60
  · exact @corr_range_sound 631 40 20 t631CodesFast t631_corr_block_02 s (by omega) (by omega)

  by_cases h3 : s < 80
  · exact @corr_range_sound 631 60 20 t631CodesFast t631_corr_block_03 s (by omega) (by omega)

  by_cases h4 : s < 100
  · exact @corr_range_sound 631 80 20 t631CodesFast t631_corr_block_04 s (by omega) (by omega)

  by_cases h5 : s < 120
  · exact @corr_range_sound 631 100 20 t631CodesFast t631_corr_block_05 s (by omega) (by omega)

  by_cases h6 : s < 140
  · exact @corr_range_sound 631 120 20 t631CodesFast t631_corr_block_06 s (by omega) (by omega)

  by_cases h7 : s < 160
  · exact @corr_range_sound 631 140 20 t631CodesFast t631_corr_block_07 s (by omega) (by omega)

  by_cases h8 : s < 180
  · exact @corr_range_sound 631 160 20 t631CodesFast t631_corr_block_08 s (by omega) (by omega)

  by_cases h9 : s < 200
  · exact @corr_range_sound 631 180 20 t631CodesFast t631_corr_block_09 s (by omega) (by omega)

  by_cases h10 : s < 220
  · exact @corr_range_sound 631 200 20 t631CodesFast t631_corr_block_10 s (by omega) (by omega)

  by_cases h11 : s < 240
  · exact @corr_range_sound 631 220 20 t631CodesFast t631_corr_block_11 s (by omega) (by omega)

  by_cases h12 : s < 260
  · exact @corr_range_sound 631 240 20 t631CodesFast t631_corr_block_12 s (by omega) (by omega)

  by_cases h13 : s < 280
  · exact @corr_range_sound 631 260 20 t631CodesFast t631_corr_block_13 s (by omega) (by omega)

  by_cases h14 : s < 300
  · exact @corr_range_sound 631 280 20 t631CodesFast t631_corr_block_14 s (by omega) (by omega)

  by_cases h15 : s < 320
  · exact @corr_range_sound 631 300 20 t631CodesFast t631_corr_block_15 s (by omega) (by omega)

  by_cases h16 : s < 340
  · exact @corr_range_sound 631 320 20 t631CodesFast t631_corr_block_16 s (by omega) (by omega)

  by_cases h17 : s < 360
  · exact @corr_range_sound 631 340 20 t631CodesFast t631_corr_block_17 s (by omega) (by omega)

  by_cases h18 : s < 380
  · exact @corr_range_sound 631 360 20 t631CodesFast t631_corr_block_18 s (by omega) (by omega)

  by_cases h19 : s < 400
  · exact @corr_range_sound 631 380 20 t631CodesFast t631_corr_block_19 s (by omega) (by omega)

  by_cases h20 : s < 420
  · exact @corr_range_sound 631 400 20 t631CodesFast t631_corr_block_20 s (by omega) (by omega)

  by_cases h21 : s < 440
  · exact @corr_range_sound 631 420 20 t631CodesFast t631_corr_block_21 s (by omega) (by omega)

  by_cases h22 : s < 460
  · exact @corr_range_sound 631 440 20 t631CodesFast t631_corr_block_22 s (by omega) (by omega)

  by_cases h23 : s < 480
  · exact @corr_range_sound 631 460 20 t631CodesFast t631_corr_block_23 s (by omega) (by omega)

  by_cases h24 : s < 500
  · exact @corr_range_sound 631 480 20 t631CodesFast t631_corr_block_24 s (by omega) (by omega)

  by_cases h25 : s < 520
  · exact @corr_range_sound 631 500 20 t631CodesFast t631_corr_block_25 s (by omega) (by omega)

  by_cases h26 : s < 540
  · exact @corr_range_sound 631 520 20 t631CodesFast t631_corr_block_26 s (by omega) (by omega)

  by_cases h27 : s < 560
  · exact @corr_range_sound 631 540 20 t631CodesFast t631_corr_block_27 s (by omega) (by omega)

  by_cases h28 : s < 580
  · exact @corr_range_sound 631 560 20 t631CodesFast t631_corr_block_28 s (by omega) (by omega)

  by_cases h29 : s < 600
  · exact @corr_range_sound 631 580 20 t631CodesFast t631_corr_block_29 s (by omega) (by omega)

  by_cases h30 : s < 620
  · exact @corr_range_sound 631 600 20 t631CodesFast t631_corr_block_30 s (by omega) (by omega)

  exact @corr_range_sound 631 620 11 t631CodesFast t631_corr_block_31 s (by omega) (by omega)

#print axioms t631_codes_corr_blocked
end ButsonTowers.Data

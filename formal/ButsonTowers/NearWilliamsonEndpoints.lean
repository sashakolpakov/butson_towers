import ButsonTowers.Endpoints
import ButsonTowers.Data.NearWilliamson

/-! Exact existence from the literal amicable quadruples. Novelty labels
are literature claims and are not assumptions of these theorems. -/

namespace ButsonTowers

theorem hadamard_118628 : HadamardExists 118628 := by
  simpa using t631_exists Data.nearWilliamson47 Data.nearWilliamson47_quad

theorem hadamard_133772 : HadamardExists 133772 := by
  simpa using t631_exists Data.nearWilliamson53 Data.nearWilliamson53_quad

theorem hadamard_148916 : HadamardExists 148916 := by
  simpa using t631_exists Data.nearWilliamson59 Data.nearWilliamson59_quad

theorem hadamard_169108 : HadamardExists 169108 := by
  simpa using t631_exists Data.nearWilliamson67 Data.nearWilliamson67_quad

theorem hadamard_174156 : HadamardExists 174156 := by
  simpa using t631_exists Data.nearWilliamson69 Data.nearWilliamson69_quad

theorem hadamard_184252 : HadamardExists 184252 := by
  simpa using t631_exists Data.nearWilliamson73 Data.nearWilliamson73_quad

theorem hadamard_199396 : HadamardExists 199396 := by
  simpa using t631_exists Data.nearWilliamson79 Data.nearWilliamson79_quad

theorem hadamard_234732 : HadamardExists 234732 := by
  simpa using t631_exists Data.nearWilliamson93 Data.nearWilliamson93_quad

theorem hadamard_244828 : HadamardExists 244828 := by
  simpa using t631_exists Data.nearWilliamson97 Data.nearWilliamson97_quad

theorem hadamard_259972 : HadamardExists 259972 := by
  simpa using t631_exists Data.nearWilliamson103 Data.nearWilliamson103_quad

theorem hadamard_275116 : HadamardExists 275116 := by
  simpa using t631_exists Data.nearWilliamson109 Data.nearWilliamson109_quad

theorem hadamard_320548 : HadamardExists 320548 := by
  simpa using t631_exists Data.nearWilliamson127 Data.nearWilliamson127_quad

def verifiedNearWilliamsonWidths : List ℕ := [47,53,59,67,69,73,79,93,97,103,109,127]

theorem near_williamson_table_endpoints : ∀ w ∈ verifiedNearWilliamsonWidths,
    HadamardExists (2524*w) := by
  intro w hw
  simp only [verifiedNearWilliamsonWidths, List.mem_cons, List.not_mem_nil, or_false] at hw
  rcases hw with h0 | h1 | h2 | h3 | h4 | h5 | h6 | h7 | h8 | h9 | h10 | h11
  · subst w; exact hadamard_118628
  · subst w; exact hadamard_133772
  · subst w; exact hadamard_148916
  · subst w; exact hadamard_169108
  · subst w; exact hadamard_174156
  · subst w; exact hadamard_184252
  · subst w; exact hadamard_199396
  · subst w; exact hadamard_234732
  · subst w; exact hadamard_244828
  · subst w; exact hadamard_259972
  · subst w; exact hadamard_275116
  · subst w; exact hadamard_320548

#print axioms hadamard_118628
#print axioms hadamard_133772
#print axioms hadamard_148916
#print axioms hadamard_169108
#print axioms hadamard_174156
#print axioms hadamard_184252
#print axioms hadamard_199396
#print axioms hadamard_234732
#print axioms hadamard_244828
#print axioms hadamard_259972
#print axioms hadamard_275116
#print axioms hadamard_320548
#print axioms near_williamson_table_endpoints

end ButsonTowers

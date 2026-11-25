import ButsonTowers.Construction
import ButsonTowers.UnrealTripling
import ButsonTowers.Data.T127
import ButsonTowers.Data.T631
import ButsonTowers.Data.Williamson
import ButsonTowers.Data.StrongStarters
import ButsonTowers.Data.Paley
import ButsonTowers.Data.Unreal13

namespace ButsonTowers
open scoped BigOperators

/-- The announcement's complete list of eighteen author-knowledge new orders.
The formal proposition asserts existence; literature status is external. -/
def announcedOrders : List ℕ :=
  [7572,17668,22716,32812,42908,58052,73196,78244,93388,98436,
   108532,128724,153964,229684,295308,386172,4772884,13001124]

theorem t631_exists {ι : Type} [Fintype ι] [DecidableEq ι] [Nonempty ι]
    (A : Fin 4 → Matrix ι ι ℚ) (hA : AmicableQuad A (4*(Fintype.card ι : ℚ))) :
    HadamardExists (2524*Fintype.card ι) := by
  simpa using exists_from_t Data.t631 Data.t631_corr Data.t631_onehot A hA

theorem t127_exists {ι : Type} [Fintype ι] [DecidableEq ι] [Nonempty ι]
    (A : Fin 4 → Matrix ι ι ℚ) (hA : AmicableQuad A (4*(Fintype.card ι : ℚ))) :
    HadamardExists (508*Fintype.card ι) := by
  simpa using exists_from_t Data.t127 Data.t127_corr Data.t127_onehot A hA

theorem williamson1_quad : AmicableQuad Data.williamson1 (4*(Fintype.card (ZMod 1) : ℚ)) := by
  simpa using (show AmicableQuad Data.williamson1 (4*(1:ℚ)) from
    ⟨Data.williamson1_sign,Data.williamson1_amicable,Data.williamson1_norm⟩)

theorem williamson3_quad : AmicableQuad Data.williamson3 (4*(Fintype.card (ZMod 3) : ℚ)) := by
  simpa using (show AmicableQuad Data.williamson3 (4*(3:ℚ)) from
    ⟨Data.williamson3_sign,Data.williamson3_amicable,Data.williamson3_norm⟩)

theorem williamson5_quad : AmicableQuad Data.williamson5 (4*(Fintype.card (ZMod 5) : ℚ)) := by
  simpa using (show AmicableQuad Data.williamson5 (4*(5:ℚ)) from
    ⟨Data.williamson5_sign,Data.williamson5_amicable,Data.williamson5_norm⟩)

theorem williamson7_quad : AmicableQuad Data.williamson7 (4*(Fintype.card (ZMod 7) : ℚ)) := by
  simpa using (show AmicableQuad Data.williamson7 (4*(7:ℚ)) from
    ⟨Data.williamson7_sign,Data.williamson7_amicable,Data.williamson7_norm⟩)

theorem williamson9_quad : AmicableQuad Data.williamson9 (4*(Fintype.card (ZMod 9) : ℚ)) := by
  simpa using (show AmicableQuad Data.williamson9 (4*(9:ℚ)) from
    ⟨Data.williamson9_sign,Data.williamson9_amicable,Data.williamson9_norm⟩)

theorem williamson11_quad : AmicableQuad Data.williamson11 (4*(Fintype.card (ZMod 11) : ℚ)) := by
  simpa using (show AmicableQuad Data.williamson11 (4*(11:ℚ)) from
    ⟨Data.williamson11_sign,Data.williamson11_amicable,Data.williamson11_norm⟩)

theorem williamson13_quad : AmicableQuad Data.williamson13 (4*(Fintype.card (ZMod 13) : ℚ)) := by
  simpa using (show AmicableQuad Data.williamson13 (4*(13:ℚ)) from
    ⟨Data.williamson13_sign,Data.williamson13_amicable,Data.williamson13_norm⟩)

theorem williamson15_quad : AmicableQuad Data.williamson15 (4*(Fintype.card (ZMod 15) : ℚ)) := by
  simpa using (show AmicableQuad Data.williamson15 (4*(15:ℚ)) from
    ⟨Data.williamson15_sign,Data.williamson15_amicable,Data.williamson15_norm⟩)

theorem williamson17_quad : AmicableQuad Data.williamson17 (4*(Fintype.card (ZMod 17) : ℚ)) := by
  simpa using (show AmicableQuad Data.williamson17 (4*(17:ℚ)) from
    ⟨Data.williamson17_sign,Data.williamson17_amicable,Data.williamson17_norm⟩)

theorem williamson19_quad : AmicableQuad Data.williamson19 (4*(Fintype.card (ZMod 19) : ℚ)) := by
  simpa using (show AmicableQuad Data.williamson19 (4*(19:ℚ)) from
    ⟨Data.williamson19_sign,Data.williamson19_amicable,Data.williamson19_norm⟩)

theorem williamson21_quad : AmicableQuad Data.williamson21 (4*(Fintype.card (ZMod 21) : ℚ)) := by
  simpa using (show AmicableQuad Data.williamson21 (4*(21:ℚ)) from
    ⟨Data.williamson21_sign,Data.williamson21_amicable,Data.williamson21_norm⟩)

theorem williamson23_quad : AmicableQuad Data.williamson23 (4*(Fintype.card (ZMod 23) : ℚ)) := by
  simpa using (show AmicableQuad Data.williamson23 (4*(23:ℚ)) from
    ⟨Data.williamson23_sign,Data.williamson23_amicable,Data.williamson23_norm⟩)

theorem williamson25_quad : AmicableQuad Data.williamson25 (4*(Fintype.card (ZMod 25) : ℚ)) := by
  simpa using (show AmicableQuad Data.williamson25 (4*(25:ℚ)) from
    ⟨Data.williamson25_sign,Data.williamson25_amicable,Data.williamson25_norm⟩)

theorem williamson27_quad : AmicableQuad Data.williamson27 (4*(Fintype.card (ZMod 27) : ℚ)) := by
  simpa using (show AmicableQuad Data.williamson27 (4*(27:ℚ)) from
    ⟨Data.williamson27_sign,Data.williamson27_amicable,Data.williamson27_norm⟩)

theorem williamson29_quad : AmicableQuad Data.williamson29 (4*(Fintype.card (ZMod 29) : ℚ)) := by
  simpa using (show AmicableQuad Data.williamson29 (4*(29:ℚ)) from
    ⟨Data.williamson29_sign,Data.williamson29_amicable,Data.williamson29_norm⟩)

theorem williamson31_quad : AmicableQuad Data.williamson31 (4*(Fintype.card (ZMod 31) : ℚ)) := by
  simpa using (show AmicableQuad Data.williamson31 (4*(31:ℚ)) from
    ⟨Data.williamson31_sign,Data.williamson31_amicable,Data.williamson31_norm⟩)

theorem williamson33_quad : AmicableQuad Data.williamson33 (4*(Fintype.card (ZMod 33) : ℚ)) := by
  simpa using (show AmicableQuad Data.williamson33 (4*(33:ℚ)) from
    ⟨Data.williamson33_sign,Data.williamson33_amicable,Data.williamson33_norm⟩)

theorem williamson37_quad : AmicableQuad Data.williamson37 (4*(Fintype.card (ZMod 37) : ℚ)) := by
  simpa using (show AmicableQuad Data.williamson37 (4*(37:ℚ)) from
    ⟨Data.williamson37_sign,Data.williamson37_amicable,Data.williamson37_norm⟩)

theorem williamson39_quad : AmicableQuad Data.williamson39 (4*(Fintype.card (ZMod 39) : ℚ)) := by
  simpa using (show AmicableQuad Data.williamson39 (4*(39:ℚ)) from
    ⟨Data.williamson39_sign,Data.williamson39_amicable,Data.williamson39_norm⟩)

theorem williamson41_quad : AmicableQuad Data.williamson41 (4*(Fintype.card (ZMod 41) : ℚ)) := by
  simpa using (show AmicableQuad Data.williamson41 (4*(41:ℚ)) from
    ⟨Data.williamson41_sign,Data.williamson41_amicable,Data.williamson41_norm⟩)

theorem williamson43_quad : AmicableQuad Data.williamson43 (4*(Fintype.card (ZMod 43) : ℚ)) := by
  simpa using (show AmicableQuad Data.williamson43 (4*(43:ℚ)) from
    ⟨Data.williamson43_sign,Data.williamson43_amicable,Data.williamson43_norm⟩)

theorem williamson49_quad : AmicableQuad Data.williamson49 (4*(Fintype.card (ZMod 49) : ℚ)) := by
  simpa using (show AmicableQuad Data.williamson49 (4*(49:ℚ)) from
    ⟨Data.williamson49_sign,Data.williamson49_amicable,Data.williamson49_norm⟩)

theorem williamson51_quad : AmicableQuad Data.williamson51 (4*(Fintype.card (ZMod 51) : ℚ)) := by
  simpa using (show AmicableQuad Data.williamson51 (4*(51:ℚ)) from
    ⟨Data.williamson51_sign,Data.williamson51_amicable,Data.williamson51_norm⟩)

theorem williamson55_quad : AmicableQuad Data.williamson55 (4*(Fintype.card (ZMod 55) : ℚ)) := by
  simpa using (show AmicableQuad Data.williamson55 (4*(55:ℚ)) from
    ⟨Data.williamson55_sign,Data.williamson55_amicable,Data.williamson55_norm⟩)

theorem williamson57_quad : AmicableQuad Data.williamson57 (4*(Fintype.card (ZMod 57) : ℚ)) := by
  simpa using (show AmicableQuad Data.williamson57 (4*(57:ℚ)) from
    ⟨Data.williamson57_sign,Data.williamson57_amicable,Data.williamson57_norm⟩)

theorem williamson61_quad : AmicableQuad Data.williamson61 (4*(Fintype.card (ZMod 61) : ℚ)) := by
  simpa using (show AmicableQuad Data.williamson61 (4*(61:ℚ)) from
    ⟨Data.williamson61_sign,Data.williamson61_amicable,Data.williamson61_norm⟩)

theorem williamson63_quad : AmicableQuad Data.williamson63 (4*(Fintype.card (ZMod 63) : ℚ)) := by
  simpa using (show AmicableQuad Data.williamson63 (4*(63:ℚ)) from
    ⟨Data.williamson63_sign,Data.williamson63_amicable,Data.williamson63_norm⟩)

theorem hadamard_7572 : HadamardExists 7572 := by
  simpa using t631_exists Data.williamson3 williamson3_quad

theorem hadamard_17668 : HadamardExists 17668 := by
  simpa using t631_exists Data.williamson7 williamson7_quad

theorem hadamard_22716 : HadamardExists 22716 := by
  simpa using t631_exists Data.williamson9 williamson9_quad

theorem hadamard_32812 : HadamardExists 32812 := by
  simpa using t631_exists Data.williamson13 williamson13_quad

theorem hadamard_42908 : HadamardExists 42908 := by
  simpa using t631_exists Data.williamson17 williamson17_quad

theorem hadamard_58052 : HadamardExists 58052 := by
  simpa using t631_exists Data.williamson23 williamson23_quad

theorem hadamard_73196 : HadamardExists 73196 := by
  simpa using t631_exists Data.williamson29 williamson29_quad

theorem hadamard_78244 : HadamardExists 78244 := by
  simpa using t631_exists Data.williamson31 williamson31_quad

theorem hadamard_93388 : HadamardExists 93388 := by
  simpa using t631_exists Data.williamson37 williamson37_quad

theorem hadamard_98436 : HadamardExists 98436 := by
  simpa using t631_exists Data.williamson39 williamson39_quad

theorem hadamard_108532 : HadamardExists 108532 := by
  simpa using t631_exists Data.williamson43 williamson43_quad

theorem hadamard_128724 : HadamardExists 128724 := by
  simpa using t631_exists Data.williamson51 williamson51_quad

theorem hadamard_153964 : HadamardExists 153964 := by
  simpa using t631_exists Data.williamson61 williamson61_quad

theorem spence_family_13 (r : ℕ) : HadamardExists (17668*13^r) := by
  have h := spence_exists_from_t Data.t631 Data.t631_corr Data.t631_onehot
    Data.spence13 Data.paley13 (by simpa using Data.spence13_strong)
    Data.paley13_isPaley (by norm_num) r
  simpa using h

theorem spence_family_17 (r : ℕ) : HadamardExists (22716*17^r) := by
  have h := spence_exists_from_t Data.t631 Data.t631_corr Data.t631_onehot
    Data.spence17 Data.paley17 (by simpa using Data.spence17_strong)
    Data.paley17_isPaley (by norm_num) r
  simpa using h

theorem spence_family_61 (r : ℕ) : HadamardExists (78244*61^r) := by
  have h := spence_exists_from_t Data.t631 Data.t631_corr Data.t631_onehot
    Data.spence61 Data.paley61 (by simpa using Data.spence61_strong)
    Data.paley61_isPaley (by norm_num) r
  simpa using h

theorem spence_family_101 (r : ℕ) : HadamardExists (128724*101^r) := by
  have h := spence_exists_from_t Data.t631 Data.t631_corr Data.t631_onehot
    Data.spence101 Data.paley101 (by simpa using Data.spence101_strong)
    Data.paley101_isPaley (by norm_num) r
  simpa using h

theorem spence_family_121 (r : ℕ) : HadamardExists (153964*121^r) := by
  have h := spence_exists_from_t Data.t631 Data.t631_corr Data.t631_onehot
    Data.spence121 Data.paley121 (by simpa using Data.spence121_strong)
    Data.paley121_isPaley (by norm_num) r
  simpa using h

theorem hadamard_229684 : HadamardExists 229684 := by
  simpa using spence_family_13 1

theorem hadamard_386172 : HadamardExists 386172 := by
  simpa using spence_family_17 1

theorem hadamard_4772884 : HadamardExists 4772884 := by
  simpa using spence_family_61 1

theorem hadamard_13001124 : HadamardExists 13001124 := by
  simpa using spence_family_101 1

def unreal13Pair : Unreal.Pair (Fin 13) 13 where
  A := Data.unreal13A
  B := Data.unreal13B
  signA := Data.unreal13A_sign
  signB := Data.unreal13B_sign
  amicable := Data.unreal13_amicable
  norm := Data.unreal13_norm

theorem unreal_cubic_family (r : ℕ) : HadamardExists (32812*3^r) := by
  let P := Unreal.tower unreal13Pair r
  have hq := Unreal.quad_amicable P
  have hc : Fintype.card (SpenceMatrix.FamilyIndex (Fin 13) (Fin 3) r) = 13*3^r := by
    simp [SpenceMatrix.family_card]
  have hA : AmicableQuad (Unreal.quad P)
      (4*(Fintype.card (SpenceMatrix.FamilyIndex (Fin 13) (Fin 3) r) : ℚ)) := by
    simpa [P, hc, Nat.cast_mul, Nat.cast_pow] using hq
  have h := t631_exists (Unreal.quad P) hA
  simpa [hc, ← Nat.mul_assoc] using h

theorem hadamard_295308 : HadamardExists 295308 := by
  simpa using unreal_cubic_family 2

theorem announced_orders_exist : ∀ n ∈ announcedOrders, HadamardExists n := by
  intro n hn
  simp only [announcedOrders,List.mem_cons,List.not_mem_nil,or_false] at hn
  rcases hn with h0 | h1 | h2 | h3 | h4 | h5 | h6 | h7 | h8 | h9 | h10 | h11 | h12 | h13 | h14 | h15 | h16 | h17
  · subst n; exact hadamard_7572
  · subst n; exact hadamard_17668
  · subst n; exact hadamard_22716
  · subst n; exact hadamard_32812
  · subst n; exact hadamard_42908
  · subst n; exact hadamard_58052
  · subst n; exact hadamard_73196
  · subst n; exact hadamard_78244
  · subst n; exact hadamard_93388
  · subst n; exact hadamard_98436
  · subst n; exact hadamard_108532
  · subst n; exact hadamard_128724
  · subst n; exact hadamard_153964
  · subst n; exact hadamard_229684
  · subst n; exact hadamard_295308
  · subst n; exact hadamard_386172
  · subst n; exact hadamard_4772884
  · subst n; exact hadamard_13001124

def verifiedWilliamsonWidths : List ℕ := [1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,37,39,41,43,49,51,55,57,61,63]

theorem williamson_table_endpoints : ∀ w ∈ verifiedWilliamsonWidths,
    HadamardExists (2524*w) ∧ HadamardExists (508*w) := by
  intro w hw
  simp only [verifiedWilliamsonWidths,List.mem_cons,List.not_mem_nil,or_false] at hw
  rcases hw with h0 | h1 | h2 | h3 | h4 | h5 | h6 | h7 | h8 | h9 | h10 | h11 | h12 | h13 | h14 | h15 | h16 | h17 | h18 | h19 | h20 | h21 | h22 | h23 | h24 | h25 | h26
  · subst w
    exact ⟨by simpa using t631_exists Data.williamson1 williamson1_quad,
      by simpa using t127_exists Data.williamson1 williamson1_quad⟩
  · subst w
    exact ⟨by simpa using t631_exists Data.williamson3 williamson3_quad,
      by simpa using t127_exists Data.williamson3 williamson3_quad⟩
  · subst w
    exact ⟨by simpa using t631_exists Data.williamson5 williamson5_quad,
      by simpa using t127_exists Data.williamson5 williamson5_quad⟩
  · subst w
    exact ⟨by simpa using t631_exists Data.williamson7 williamson7_quad,
      by simpa using t127_exists Data.williamson7 williamson7_quad⟩
  · subst w
    exact ⟨by simpa using t631_exists Data.williamson9 williamson9_quad,
      by simpa using t127_exists Data.williamson9 williamson9_quad⟩
  · subst w
    exact ⟨by simpa using t631_exists Data.williamson11 williamson11_quad,
      by simpa using t127_exists Data.williamson11 williamson11_quad⟩
  · subst w
    exact ⟨by simpa using t631_exists Data.williamson13 williamson13_quad,
      by simpa using t127_exists Data.williamson13 williamson13_quad⟩
  · subst w
    exact ⟨by simpa using t631_exists Data.williamson15 williamson15_quad,
      by simpa using t127_exists Data.williamson15 williamson15_quad⟩
  · subst w
    exact ⟨by simpa using t631_exists Data.williamson17 williamson17_quad,
      by simpa using t127_exists Data.williamson17 williamson17_quad⟩
  · subst w
    exact ⟨by simpa using t631_exists Data.williamson19 williamson19_quad,
      by simpa using t127_exists Data.williamson19 williamson19_quad⟩
  · subst w
    exact ⟨by simpa using t631_exists Data.williamson21 williamson21_quad,
      by simpa using t127_exists Data.williamson21 williamson21_quad⟩
  · subst w
    exact ⟨by simpa using t631_exists Data.williamson23 williamson23_quad,
      by simpa using t127_exists Data.williamson23 williamson23_quad⟩
  · subst w
    exact ⟨by simpa using t631_exists Data.williamson25 williamson25_quad,
      by simpa using t127_exists Data.williamson25 williamson25_quad⟩
  · subst w
    exact ⟨by simpa using t631_exists Data.williamson27 williamson27_quad,
      by simpa using t127_exists Data.williamson27 williamson27_quad⟩
  · subst w
    exact ⟨by simpa using t631_exists Data.williamson29 williamson29_quad,
      by simpa using t127_exists Data.williamson29 williamson29_quad⟩
  · subst w
    exact ⟨by simpa using t631_exists Data.williamson31 williamson31_quad,
      by simpa using t127_exists Data.williamson31 williamson31_quad⟩
  · subst w
    exact ⟨by simpa using t631_exists Data.williamson33 williamson33_quad,
      by simpa using t127_exists Data.williamson33 williamson33_quad⟩
  · subst w
    exact ⟨by simpa using t631_exists Data.williamson37 williamson37_quad,
      by simpa using t127_exists Data.williamson37 williamson37_quad⟩
  · subst w
    exact ⟨by simpa using t631_exists Data.williamson39 williamson39_quad,
      by simpa using t127_exists Data.williamson39 williamson39_quad⟩
  · subst w
    exact ⟨by simpa using t631_exists Data.williamson41 williamson41_quad,
      by simpa using t127_exists Data.williamson41 williamson41_quad⟩
  · subst w
    exact ⟨by simpa using t631_exists Data.williamson43 williamson43_quad,
      by simpa using t127_exists Data.williamson43 williamson43_quad⟩
  · subst w
    exact ⟨by simpa using t631_exists Data.williamson49 williamson49_quad,
      by simpa using t127_exists Data.williamson49 williamson49_quad⟩
  · subst w
    exact ⟨by simpa using t631_exists Data.williamson51 williamson51_quad,
      by simpa using t127_exists Data.williamson51 williamson51_quad⟩
  · subst w
    exact ⟨by simpa using t631_exists Data.williamson55 williamson55_quad,
      by simpa using t127_exists Data.williamson55 williamson55_quad⟩
  · subst w
    exact ⟨by simpa using t631_exists Data.williamson57 williamson57_quad,
      by simpa using t127_exists Data.williamson57 williamson57_quad⟩
  · subst w
    exact ⟨by simpa using t631_exists Data.williamson61 williamson61_quad,
      by simpa using t127_exists Data.williamson61 williamson61_quad⟩
  · subst w
    exact ⟨by simpa using t631_exists Data.williamson63 williamson63_quad,
      by simpa using t127_exists Data.williamson63 williamson63_quad⟩

end ButsonTowers

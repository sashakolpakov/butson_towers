import Mathlib

/-!
The exact algebraic core of Spence's recurrence.  The carrier is a
commutative ring, not a field.  Doubled outputs avoid requiring division
by two and apply in particular to rational finite group algebras.
-/

namespace ButsonTowers.Spence

variable {R : Type*} [CommRing R]

def minusNorm (d e f g : R) : R := (d - e) ^ 2 + (f - g) ^ 2
def plusNorm (d e f g : R) : R := (d + e) ^ 2 + (f + g) ^ 2

def doubledPlus (a b j t : R) : R := (a - b) * j + (a + b) * (t + 1)
def doubledMinus (a b j t : R) : R := (a - b) * j + (a + b) * (t - 1)

theorem doubled_minus_norm (d e f g j t : R) :
    minusNorm (doubledPlus d e j t) (doubledMinus d e j t)
      (doubledPlus f g j t) (doubledMinus f g j t) =
      4 * plusNorm d e f g := by
  unfold minusNorm plusNorm doubledPlus doubledMinus
  ring

theorem doubled_plus_norm (d e f g j t q : R)
    (hj : j ^ 2 = q * j) (ht : t ^ 2 = q - j) (hjt : j * t = 0) :
    plusNorm (doubledPlus d e j t) (doubledMinus d e j t)
      (doubledPlus f g j t) (doubledMinus f g j t) =
      4 * (q * plusNorm d e f g +
        (q * minusNorm d e f g - plusNorm d e f g) * j) := by
  unfold minusNorm plusNorm doubledPlus doubledMinus
  linear_combination
    4 * ((d-e)^2 + (f-g)^2) * hj +
    4 * ((d+e)^2 + (f+g)^2) * ht +
    8 * ((d-e)*(d+e) + (f-g)*(f+g)) * hjt

theorem sum_squares_from_strong_norms (d e f g m p : R)
    (hm : minusNorm d e f g = m) (hp : plusNorm d e f g = p) :
    2 * (d ^ 2 + e ^ 2 + f ^ 2 + g ^ 2) = m + p := by
  rw [← hm, ← hp]
  unfold minusNorm plusNorm
  ring

theorem strong_norm_induction (q : R) (M P J : ℕ → R)
    (hm0 : M 0 = 4) (hp0 : P 0 = 4*q)
    (hm : ∀ r, M (r+1) = P r)
    (hp : ∀ r, P (r+1) = q * P r + (q * M r - P r) * J r) :
    ∀ r, M r = 4*q^r ∧ P r = 4*q^(r+1) := by
  intro r
  induction r with
  | zero => simpa using And.intro hm0 hp0
  | succ r ih =>
    constructor
    · simpa [hm r] using ih.2
    · rw [hp r, ih.1, ih.2]
      simp only [pow_succ]
      ring

def entryStep (a b t : ℚ) (diagonal plus : Bool) : ℚ :=
  (a-b)/2 + (a+b)/2 * (t + if diagonal then (if plus then 1 else -1) else 0)

theorem entry_step_sign (a b t : ℚ) (diagonal plus : Bool)
    (ha : a = 1 ∨ a = -1) (hb : b = 1 ∨ b = -1)
    (ht : if diagonal then t = 0 else t = 1 ∨ t = -1) :
    entryStep a b t diagonal plus = 1 ∨ entryStep a b t diagonal plus = -1 := by
  cases diagonal with
  | false =>
    have ht' : t = 1 ∨ t = -1 := ht
    rcases ht' with rfl | rfl <;> rcases ha with rfl | rfl <;>
      rcases hb with rfl | rfl <;> cases plus <;> norm_num [entryStep]
  | true =>
    have ht' : t = 0 := ht
    subst t
    rcases ha with rfl | rfl <;> rcases hb with rfl | rfl <;>
      cases plus <;> norm_num [entryStep]

end ButsonTowers.Spence

import Mathlib.Data.Int.Basic
import Mathlib.Data.Fin.VecNotation

/-! Balanced literal lookup has logarithmic kernel reduction depth; Array
lookup alone reduces through Lists in the kernel. No native evaluator is used. -/

namespace ButsonTowers.Data

inductive LiteralTree (α : Type) where
  | leaf : α → LiteralTree α
  | node : ℕ → LiteralTree α → LiteralTree α → LiteralTree α

def treeAt {α : Type} : LiteralTree α → ℕ → α
  | .leaf value, _ => value
  | .node split left right, x =>
      if x < split then treeAt left x else treeAt right (x-split)

end ButsonTowers.Data

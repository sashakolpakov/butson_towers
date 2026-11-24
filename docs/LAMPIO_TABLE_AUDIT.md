# Audit of the Lampio–Östergård–Szöllősi tables

This note distinguishes three different claims that are easy to conflate:

1. a matrix is known to exist;
2. an exact number of monomial-equivalence classes is published; and
3. that exact number has a compact, independently replayed completeness
   certificate in this repository.

The source paper is P. H. J. Lampio, P. R. J. Östergård, and F. Szöllősi,
*Orderly generation of Butson Hadamard matrices*, Mathematics of Computation
89 (2020), 313–331:

<https://arxiv.org/abs/1707.02287>

The current electronic supplement is:

<https://wiki.aalto.fi/spaces/Butson/pages/120482304/Matrices+up+to+monomial+equivalence>

The supplement reports that it was last updated on 2025-11-17. Our local
snapshot and audit were made on 2026-07-26.

## Meaning of the table entries

The table covers `2 <= n <= 21` and `2 <= q <= 17`, hence 320 parameter
pairs.

- A number is the claimed exact count up to **monomial equivalence**.
- A white/blank cell is exact nonexistence from the Lam–Leung obstruction.
- `E` means existence is known but full classification is not.
- `U` means existence is unknown.

Hadamard equivalence in the paper is coarser: it additionally permits one
global automorphism of `Z/qZ`. The headline table is not counted under that
coarser relation.

## What the 2017 paper had settled

In the arXiv source of the original table:

| Status | Cells |
|---|---:|
| Numerical exact count, including displayed zeroes | 102 |
| Lam–Leung-obstructed blank cells | 162 |
| Existence only (`E`) | 38 |
| Unknown (`U`) | 18 |

Thus 264 of the 320 displayed parameter pairs had a completely determined
answer, if theoretical obstruction cells are counted as classifications.
Only 102 had a numerical table entry produced by classification or a separate
nonexistence result.

The paper makes a stronger rectangular statement for its core range: every
case with `n <= 11` and `q <= 17` was settled. In that 160-cell rectangle,
74 cells have numerical exact counts and 86 are Lam–Leung obstructions.

Beyond that rectangle, the paper reports selected classifications. Its main
new large enumerations include:

- `BH(16,4)`: 1,786,763 monomial-equivalence classes;
- `BH(14,6)`: 167,776 classes;
- `BH(21,3)`: 72 classes.

The narrative also records complete classification for the real case through
order 32, the cubic case through order 21, and the quaternary case through
order 16, drawing on both this paper and earlier literature.

## What the current supplement says

The 2025 supplement has:

| Status | Cells |
|---|---:|
| Numerical exact count, including displayed zeroes | 113 |
| Lam–Leung-obstructed cells | 162 |
| Existence only (`E`) | 40 |
| Unknown (`U`) | 5 |

It therefore gives a determined answer for 275 of 320 cells, or 85.94%.
Forty additional cells have examples but not a full classification.

The five cells still marked unknown are:

```text
BH(15,12), BH(16,15), BH(17,12), BH(17,15), BH(19,15).
```

Compared with the paper's original table, the supplement adds the exact count
`BH(12,12) = 4,497,733`, changes ten former unknown cells to exact
nonexistence, and changes three former unknown cells to existence-only:

```text
New exact count:
  BH(12,12) = 4,497,733

New exact nonexistence:
  BH(13,14), BH(13,15), BH(14,15), BH(15,10), BH(15,14),
  BH(17,10), BH(17,14), BH(19,10), BH(19,14), BH(21,10)

Unknown -> existence:
  BH(18,10), BH(18,14), BH(18,15)
```

The machine-readable snapshot used by our regression tests is
`tests/data/lampio_table_2025.csv`. Every one of its 162 white cells is checked
against an independent implementation of the Lam–Leung semigroup criterion.

## What “certified” means here

The paper's numerical entries are peer-reviewed exact classification claims
based on orderly generation. The supplement provides representative files for
many positive cells and search traces for many zero cells. These are
substantially stronger evidence than an existence table.

They are not, however, one uniform proof-certificate format from which a small
independent checker can reconstruct exhaustiveness, inequivalence, and all
search pruning decisions. In particular, this repository has not independently
repeated the searches producing millions of classes. We treat those cells as
published exact classifications, not as locally recertified results.

Order four is different. The supplement has 24 representatives across
`q = 2,4,...,16`. We snapshot those eight official `BH1` files under
`tests/data/lampio_order4/` and check all of the following:

1. every record decodes correctly, including punctuation-valued exponents;
2. every matrix satisfies all row and column inner products by exact
   cyclotomic arithmetic;
3. every representative maps by explicit row/column permutation and
   re-dephasing to the certified standard family;
4. the representatives hit every orbit under `t -> -t` and `t -> t^(-1)`
   exactly once; and
5. the independent Laurent-saturated minimal-prime certificate proves that
   there are no other order-four complex or Butson families.

Consequently the complete order-four row is independently confirmed here, and
the formula

```text
# BH(4,q) / monomial equivalence = floor(q/4) + 1  for even q
```

extends that confirmation to every positive `q`, not only the eight columns
printed in the table.

# Broad input and scalar-restriction batch

Selected published inputs and completed scalar-restriction trees are retained
as reproducible construction controls. Their source data, exact matrix recipes
and verification reports are listed below.

## Inputs and target arithmetic

All524 matrices in the selected published catalogue cells were acquired and
their full Grams checked: `BH(7,6)`, `BH(7,12)`, `BH(9,6)`, `BH(9,12)` and
`BH(13,6)`. Complete root row/column gauge tests retain33 representatives in
their original alphabets. This covers the listed matrices, not all possible
Butson matrices at these orders.

All455 sixth-root catalogue inputs were then embedded into twelfth roots
and checked under every twelfth-root row/column gauge. This **rescues80**
inputs excluded in the smaller alphabet: four at order9 and76 at order13.
With14 already admissible cases, the provider supplies94 fixed larger-field
representatives. They remain monomially equivalent to their published
inputs; no new projective Butson classes are claimed.

The expanded phase-independent preflight checks **649 input bindings**,
representing629 distinct literal matrices modulo global phase, against all
156 orders in the attributed frontier snapshot. It covers1,216 distinct
anchor alphabets and2,908 divisible source/anchor/target cases. Exactly1,374
pass the checked necessary conditions, covering22 target orders.

The declared input phases are **2,4,6,8,12,16,20,24,40,80,120**. The code
computes the actual entry field rather than treating a declared phase as
minimal: the fifth-order input labelled phase20 actually lies in Q(zeta5),
and one phase24-labelled input reduces further after a global anchor.
The sole retained target above20,000 is23960; its candidate routes still
need unknown frames of orders2396 or4792. Necessary-condition survivors
are not completed matrices.

## Completed actual matrix trees

Every node below is a completed Butson/sign matrix with actual entries,
saved source data and exact scalar-edge proofs. Normalization is represented
by explicitly different input/output coordinate bases. Only scalar
restrictions and exact literal duplicate pruning are child operations.

| Run | Fixed input bindings | Completed nodes | Duplicate branches pruned |
|---|---:|---:|---:|
| Full published/registry pool |555|765|30|
| Larger-root gauge representatives |94|282|0|
| Partial third/fifth and repeated-prime restrictions |90|464|89|

These are separate saved runs, with overlaps between them; the row counts
must not be added as global matrix-equivalence classes. They use about27.6MB
of literal bundles and checkpoints in total. Every declared queue is empty.
The first two runs'765 and282 bundles were all replayed independently.

The partial-prime tree includes32 fifth-root and33 third-root edges,
all admissible input-frame positions, both repeated-prime conventions,
and both odd-prime removal orders where available. Its53 real endpoint
matrices have independently checked full sign Grams. Their orders are
28,36,40,104,312,480,640,936,1120 and1280, all known. Every one of its
464 saved nodes replays, and `SavedNode.as_verified_source()` exposes its
actual matrix for further work.

## Replays

```sh
PYTHONPATH=src:. python -m experiments.scalar_lift.large_butson_factory --database search-state/broad-root-flat-all-inputs-v1.sqlite3 --all-verified --seconds 1 --report certificates/scalar-lift/broad-root-flat-all-inputs-tree.json
PYTHONPATH=src:. python -m experiments.scalar_lift.large_butson_factory --database search-state/expanded-gauge-root-flat-v1.sqlite3 --expanded-only --seconds 1 --report certificates/scalar-lift/expanded-gauge-root-flat-tree.json
```

These commands resume completed checkpoints without redoing matrix nodes.
The source-acquisition and partial-descent guides give the component replay
commands and their exact scopes.

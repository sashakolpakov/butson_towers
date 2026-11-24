# Scalar restrictions with explicit coordinate matrices and metrics

The tree applies only scalar restriction to already verified matrices.
Its first completed inventory starts with the fixed character matrices
\(BH(5329,24)\) and \(BH(58081,120)\). It descends through all proper
**full cyclotomic divisor fields**, retaining matrices even when their
entries are sums of roots or rational numbers rather than roots themselves.
Other abelian subfields, such as the real quadratic subfield of
\(\mathbb Q(\zeta_5)\), are outside this implementation.

For a basis \(\beta_0,\ldots,\beta_{d-1}\) of \(L/K\), the saved action is
the literal multiplication matrix \(\rho_{L/K}(z)\). The Hermitian metric is

\[
G_{ij}=\operatorname{Tr}_{L/K}(\overline{\beta_i}\beta_j).
\]

If the starting matrix has \(HH^*=nI\), its restricted matrix \(A\), of
order \(nd\), satisfies the exact identity

\[
A\,(I_n\otimes G^{-1})\,A^*
=n\,(I_n\otimes G^{-1}).
\]

The scale remains the original source order \(n\), not the expanded order
\(nd\). No Fourier frame, source modification, or Euclidean normalization
is inserted by this tree.

Each scalar edge is checked coefficientwise against the full original-root
action table. The new metric is also independently checked by trace
transitivity. The basis change to the direct canonical restriction is saved
and checked for both scalar conjugacy and Hermitian metric congruence.
This equivalence is **not** used to discard a coordinate candidate: a general
basis change need not preserve literal signs. Only identical used action
tables and trace metrics, with the same source and lower field, are merged.

The completed two-source run saved **104 distinct coordinate matrices** in
about 2.1 MB, with no exact duplicates and no pending paths. Its depth
counts are 2, 16, 38, 36 and 12. All 52 rational endpoints are integral:
eight have order 42,632 and 44 have order 1,858,592. Every one has zero
entries, so none is a literal sign matrix. These are completed weighted
orthogonal matrices, not new real Hadamard matrices.

The existing exact odd-relation certificate is replayed after constructing
the inventory. For both character sources, the actual scalar alphabet
contains \(1,\zeta_6,\zeta_6^{-1}\), with
\(1-\zeta_6-\zeta_6^{-1}=0\). It excludes a uniform rational-linear map
from that fixed alphabet to equally sized sign blocks, including arbitrary
uniform scalar towers and basis choices. The weighted matrices themselves
are retained. This flag does not claim to exclude source modifications or
nonuniform maps.

`--retained-registry` uses a separate checkpoint and the 18 **existing
verified** registry sources not excluded by that saved alphabet test.
It performs the same scalar restrictions, without constructing new sources.
Passing the alphabet test is necessary information, not a proof of a sign
endpoint. Dense-source entry counts use the actual exponent histogram;
character-source counts also account for group development.

That retained-source inventory is also complete: **368 distinct coordinate
matrices**, no exact duplicates, no pending tasks, and about 8.2 MB of saved
data. Its depth counts are 18, 82, 138, 102 and 28. All 184 rational endpoints
are integral and contain zeros; none is a literal sign matrix. Their distinct
orders are 6, 26, 40, 104, 312, 480, 640, 936, 1120, 1280, 2808, 4000, 8000,
8424 and 16000. Together the two completed inventories contain 472 weighted
coordinate matrices and 236 integral rational endpoints. These counts cover
the declared divisor-field power-basis paths only, not arbitrary basis/frame
choices or every abelian subfield.

An additional completed run uses four fixed inputs: the acquired catalogue
`P7` sixth-root and twelfth-root representatives, `D14`, and the previously
saved literal unreal `BH(9,6)`. No source construction is a tree edge.
`--additional-inputs` selects them, and repeated `--source SOURCE_ID` options
select a subset. Their [completed report](../certificates/scalar-lift/weighted-additional-inputs-tree.json)
contains **12 literal matrix nodes**, no duplicates and no pending tasks,
using 105,467 bytes. All saved nodes passed exact replay. The rational
weighted endpoints have orders 14, 18 and 28 and contain zeros.
The acquired inputs' larger ambient fields are handled separately in
[the period-field runs](AMBIENT_FIELD_CANDIDATES.md); the small orders in this
divisor-field inventory do not cap their possible ambient degrees.

```sh
PYTHONPATH=src:. ./.native-env/bin/python -m experiments.scalar_lift.weighted_scalar_tree --database search-state/weighted-additional-inputs-v1.sqlite3 --additional-inputs --seconds 30 --report certificates/scalar-lift/weighted-additional-inputs-tree.json
```

There is a separate, stronger conclusion about **new orders from these fixed
inputs**, rather than about all matrix coordinates. Restriction to the
rationals always has order \(n[L:\mathbb Q]\), independently of the
intermediate fields, because relative degrees multiply. In the current
20-input registry, 16 of these orders already have verified real
constructions, two are the impossible real orders 6 and 26, and the two
character inputs have the exact uniform-linearity obstruction above.
The [fixed-input order audit](../certificates/scalar-lift/restriction-fixed-input-order-audit.json)
therefore closes additional real orders from this particular input/ambient
field/rank combination. It does not exclude other starting matrices, larger
ambient fields, representation multiplicities, nonuniform maps, or
inequivalent matrices at an existing order. It also does not equate the
weighted Gram scale with the expanded matrix order.

The reports are
[`weighted-scalar-only-tree.json`](../certificates/scalar-lift/weighted-scalar-only-tree.json)
and
[`weighted-retained-registry-tree.json`](../certificates/scalar-lift/weighted-retained-registry-tree.json).

Each NPZ contains the original source coefficient array, every scalar action,
the metric, its inverse, and reconstruction/proof metadata.
`SavedWeightedMatrix.entry(i,j)` returns an exact lower-field coefficient
vector and denominator. `weight_entry(i,j)` returns the corresponding entry
of \(I_n\otimes G^{-1}\). `replay()` reconstructs the source and checks all
scalar edges, action tables and metrics. Rational endpoint histograms count
all entries of the entire implicit matrix exactly.

The SQLite queue resumes breadth first and writes a completed-node bundle
before marking that node done. Its report records field paths, basis data,
explicit equivalence certificates, duplicate edges and unfinished tasks.
The original root-preserving 21-matrix inventory is a separate artifact:
those steps insert Fourier frames and must not be identified with these
pure restrictions merely because they eventually reach the same field.

Tests include complete literal weighted Gram checks of orders 100 and 6,
distinct coordinate bases with verified metric congruence, saved-bundle
replay, the exact field-path count, locking and checkpoint resume.

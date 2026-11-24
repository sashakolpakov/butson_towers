# Finite phase-independent restriction preflight

The [executed certificate](../certificates/scalar-lift/phase-independent-frontier.json)
covers **555 verified input bindings and all156 target orders** in the
attributed frontier snapshot. This is a complete finite preflight of the
declared inputs, not a construction of the surviving target matrices.

The inputs include all20 original registry entries; three acquired P7/D14
representatives; two saved raw M13 candidates; the saved order9 source;
the fixed `BH(10,4)` and `BH(12,4)` controls; three published real668/716
inputs; and all524 newly acquired catalogue entries. The new catalogue
breakdown is2 `BH(7,6)`,4 `BH(7,12)`,17 `BH(9,6)`,65 `BH(9,12)`, and436
`BH(13,6)` representatives. Their fixed acquisition gauges and exact source
proofs are replayed by their provider. Inputs excluded by an alphabet
relation are retained in this inventory so that exclusions remain visible.

Declared phases checked are **2,4,6,8,12,16,20,24,40,80,120**. No sextic
pair identity is imposed on the other phases.

## What was executed

There are86,580 source/target divisibility checks,841 distinct raw or
used-entry-anchored scalar alphabets, and3,542 explicit used-entry anchor
recipes. These produce1,823 divisible anchor/target records. Of these,
289 survive the necessary tests, from41 source bindings, covering22
distinct target orders. All declared cases finished; no solver or source
construction search ran.

Exact comparison after removing a single global phase identifies543
literal input classes among the555 bindings, and40 such classes among the
41 surviving bindings. In particular, the phase24-labelled order13 seed
at tripling level0 is a global-phase version of the existing sextic seed,
not an independent new construction. This equivalence is reported rather
than confused with a literal equality of matrices or of weighted bases.
Arbitrary row/column monomial equivalence is not used for pruning: it can
change the uniform scalar alphabet and the corresponding feasibility
conditions.

Every record is bound to the verified source hash and proof hash. For each
alphabet U, every anchor `e0 in U` is recorded as the explicit global gauge
`U-e0`. Identical alphabets share necessary tests, but no literal matrix
nodes are merged. Galois equivalence is not used to discard cases.

The twelfth-root gauge exclusions were independently checked against
**all nine** missing-even/missing-odd antipodal-class patterns, on all69
published twelfth-root inputs, using a transpose-first constraint search.
Every result agrees with the acquisition provider:19 feasible classes.
Its three patterns `(0,1),(0,3),(0,5)` cover all nine because a global
even shift sends any missing pair `(e,o)` to `(0,o-e)`. Setting the first
row shift to zero is an independent opposite row/column gauge that leaves
entries unchanged. The [independent audit](../certificates/scalar-lift/phase-independent-q12-gauge-audit.json)
records all69 comparisons. Thus the provider does not restrict the
twelfth-root check to only three translated two-coset palettes.

## Exact field and sign gates

For source order n and target N, the scalar width is `d=N/n`. The model
uses uniform scalar input/output bases with each used scalar block having
Gram `dI`. Its necessary tests are:

1. The actual entry-field degree divides d.
2. A sign block has possible order: `d=1`, `d=2`, or `4|d`.
3. The actual used roots have no proved odd integer zero relation.
4. The existing determinant and dyadic necessary tests for a standard
   orthogonal scalar representation pass.

The first test uses **actual roots**, not the phase label. If the declared
phase is q and the used exponents have gcd g with q, their field is
`Q(zeta_(q/g))`, of degree `phi(q/g)`. Each anchor may reduce this field.
The certificate records both the original exponents and the effective
exponents divided by g, so subsequent executors can construct the correct
scalar action.

For example, `fifth-5-q20` actually uses exponents `{4,12,16}` modulo20.
Its field is `Q(zeta5)` of degree4. Its three anchored palettes remain
fifth-root palettes; they are not converted into a sextic problem.
The phase24-labelled cubic seed has raw field `Q(zeta12)` of degree4;
its anchored field is `Q(zeta6)` of degree2.

Failure of these tests excludes only the stated uniform-block model and
the displayed input/gauge. Arbitrary nonuniform source-coordinate maps
and unsearched row/column gauges are outside the exclusions. Passing is
not a flat-frame theorem, an ambient-field construction, or a new Hadamard
matrix.

## Surviving target orders

The necessary-test survivors cover

`2676,4676,5028,5748,5892,6852,6876,7572,7716,7924,8684,9308,9716,10052,10284,10356,11364,11460,11596,11652,11892,23960`.

The only surviving target above20,000 in this finite inventory is
**23960**, with two different source routes:

| Fixed source | Actual field | Required width | Status |
|---|---|---:|---|
| `fifth-5-q20` | `Q(zeta5)`, degree4 | 4792 | All four raw/anchored alphabet cases pass necessary tests |
| `saved-modified-gs10-q4` | `Q(i)`, degree2 | 2396 | Full quaternary alphabet passes necessary tests |

Neither route has a constructed sign frame. In particular, a quaternary
source does not supply the required order2396 Hadamard frame by itself.
The larger fifth tensor sources do not supply additional routes here:
their widths fail the exact field-degree or sign-block-order tests.

## Replay and API

```sh
PYTHONPATH=src:. ./.native-env/bin/python -m experiments.scalar_lift.phase_independent_frontier --check certificates/scalar-lift/phase-independent-frontier.json
```

The [implementation](../experiments/scalar_lift/phase_independent_frontier.py)
exposes `all_sources()`, `anchor_cases(source)`, and
`route_records(include_excluded=False, source_ids=None)`. The last returns
the necessary-test survivors by default; it does not return constructed
matrices. The saved providers continue to expose exact entry accessors.

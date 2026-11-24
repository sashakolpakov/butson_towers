# Additional fixed inputs for scalar restriction

Three complete Butson matrices have been added as starting data. After choosing
each fixed input, this module performs no source construction, deformation, or
order-changing composition. All subsequent tree operations remain scalar
restrictions. These inputs supply no new real Hadamard order.

The Bruzda–Tadej–Życzkowski catalogue was retrieved at commit
`34328f8ff1083607e058bed04e3f55e09a1f9bc1`. Literal exponent data, immutable
source links, and source-file hashes are saved under
`experiments/scalar_lift/catalogue_inputs/`.

| Fixed input | Exact provenance | Actual alphabet | Matching recorded targets |
|---|---|---|---|
| `catalogue-p7-q6` | [Published P7](https://github.com/matrix-toolbox/chm/blob/34328f8ff1083607e058bed04e3f55e09a1f9bc1/CHM/P7.m), fixed row and column shifts `(0,0,0,1,1,1,0)` | `{0,1,3,4}` modulo 6 | 4676, 7924, 9716, 10052 |
| `catalogue-p7-q12` | Existing exact order-seven input, independently matched by the full monomial canonical key to the published `P7(1/12)` | `mu_12` minus `mu_4` | Same four targets |
| `catalogue-d14-q4` | [Published D14](https://github.com/matrix-toolbox/chm/blob/34328f8ff1083607e058bed04e3f55e09a1f9bc1/CHM/D14.m), all five parameters zero | All fourth roots | Same four targets |

The seventh-order twelfth-root input is the literal class-3 matrix already saved
in `certificates/scalar-lift/inverse-dense-lifts.json`. The acquisition replay
checks its exact equivalence to the catalogue specialization, both full Grams,
and a cross ratio of order 12. Thus it cannot be reduced to sixth-root entries
by row and column phase gauges. No new parameter-family search was run.

The full cyclotomic Gram is checked independently for every retained matrix.
Each passes the complete odd-integer-relation test for its actual scalar
alphabet. Passing that necessary test does **not** prove the existence of a
sign frame. The order-seven inputs require total widths 668, 1132, 1388, and
1436; D14 requires 334, 566, 694, and 718. For a uniform orthogonal scalar
action with a normalized frame `FF^T=dI`, each used block must itself satisfy
`B_e B_e^T=dI`. An all-sign block therefore requires `d=1,2` or `4|d`.
All four D14 widths fail this necessary condition, so D14 supplies no surviving
route to these targets within that model. D14 remains a valid published input
for weighted restriction; the width test does not exclude general weighted or
nonuniform coordinate maps. Suitable ambient fields and final sign frames for
the two seventh-order inputs remain separate work. None of these source orders
divides the twenty recorded targets above 20,000.

Two other literal catalogue inputs were acquired and checked:
[M13B](https://github.com/matrix-toolbox/chm/blob/34328f8ff1083607e058bed04e3f55e09a1f9bc1/CHM/M13B.m)
and [M13C](https://github.com/matrix-toolbox/chm/blob/34328f8ff1083607e058bed04e3f55e09a1f9bc1/CHM/M13C.m).
The catalogue credits Nicoara, LaClair, Geist, and Wintenberg (May 2016).
Their full Grams pass, but their raw alphabets have an odd scalar relation.
Moreover, all sixth-root row/column gauges fail to omit a cubic coset. The
complete ternary row-prefix audit terminates after six and eight visited
prefixes respectively: a column whose three cubic colours have all appeared
can no longer be repaired by any column shift. Fixing the first row shift to
zero loses no solutions. Signs and row/column permutations do not affect this
criterion. Every such gauge therefore uses all three classes modulo 3.
One representative from each class, with suitable integer coefficient signs,
gives a rotated `1+zeta3+zeta3^2=0`, an odd three-term relation. This excludes
these sixth-root gauge orbits from uniform scalar-to-sign maps;
it is not a nonexistence result for other order-thirteen inputs or nonuniform
coordinate maps. Gauges by roots in larger fields were not tested.

The frontier comparison uses the dated, attributed 156-order snapshot in
`certificates/scalar-lift/forward-novelty-frontier-2026-11-22.json`. A match is
an arithmetic target, not a construction or an independent novelty proof.

Replay:

```sh
PYTHONPATH=src:. ./.native-env/bin/python -m experiments.scalar_lift.catalogue_source_acquisition --check certificates/scalar-lift/catalogue-source-acquisition.json
```

The API `load_source(source_id)` returns the existing `VerifiedSource` protocol,
including the literal matrix, exact proof, source hash, and entry accessor.
The old twenty-source registry and completed restriction checkpoints are
unchanged.

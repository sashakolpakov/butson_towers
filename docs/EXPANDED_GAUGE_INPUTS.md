# Larger initial root alphabets rescue 80 published inputs

All **455** published sixth-root classes of orders 7, 9 and 13 were
embedded into twelfth roots by `e -> 2e`, then audited over **all
twelfth-root row/column gauges**. This is a separate acquisition of
fixed starting representatives; the existing registry, restriction
checkpoints and phase-independent frontier report are unchanged.

| Order | Originally admissible over sixth roots | Rescued over twelfth roots | Still excluded |
|---:|---:|---:|---:|
| 7 | 2 | 0 | 0 |
| 9 | 2 | 4 | 11 |
| 13 | 10 | 76 | 350 |
| Total | **14** | **80** | **361** |

The exhaustive audit finishes all 455 cases, with **76,182 row-prefix
visits**, at most 544 per input, and no undecided cases. Every original
and selected gauged matrix passes its full exact cyclotomic Gram.
The [complete certificate](../certificates/scalar-lift/expanded-gauge-inputs.json)
stores all witnesses and exclusions, source hashes and dataset provenance.

For example, the published `BH(9,6)` class 001 has no admissible sixth-root
row/column gauge. Its explicit twelfth-root gauge gives the actual palette
`{0,1,2,3,6,7,8,9}` modulo twelve, which passes the exact odd-relation test.
This supplies another concrete source toward `6876=9*764`. The analogous
order-thirteen rescues supply sources toward 8684, 9308 and 11596.

All 80 rescues have literal entry root order twelve. They remain
monomially equivalent to their original sixth-root classes when larger
row/column phases are allowed; their projective root invariant is
unchanged. Thus they are useful new fixed representatives and entry-field
choices, not 80 previously unknown projective Butson classes.

The gauge CSP uses all nine maximal antipodal palettes modulo an explicit
global-phase symmetry. Its exact completeness argument and independent
palette/brute-force checks are described in
[BROAD_PUBLISHED_INPUTS.md](BROAD_PUBLISHED_INPUTS.md). This extension covers
all twelfth-root gauges; no assertion is made about still larger phases.

The separate [denominator gate](NORMALIZATION_DENOMINATOR_GATE.md) applies
to these inputs as well: dyadic orthogonal mixing cannot give the required
normalized sign images at widths four times an odd integer greater than
one. Suitable rational bases must introduce denominator primes of the
odd width factor. Passing the initial alphabet test does not supply such
a basis or establish a new Hadamard order.

API in [expanded_gauge_inputs.py](../experiments/scalar_lift/expanded_gauge_inputs.py):

* `sources()` lists all 455 audited descriptors.
* `admissible_sources()` lists the 94 retained fixed representatives.
* `load_source('expanded-9-001-q12')` returns the existing immutable
  `VerifiedSource` protocol with exact proofs and the explicit gauge.

```sh
PYTHONPATH=src:. ./.native-env/bin/python -m experiments.scalar_lift.expanded_gauge_inputs --check certificates/scalar-lift/expanded-gauge-inputs.json
```

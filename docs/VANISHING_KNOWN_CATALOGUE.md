# Complete scalar vanishing atoms at phases 30 and 120

**There are exactly 844 literal minimal positive vanishing root sets in \(\mu_{120}\), in nine rotation classes. Their largest weight is nine.** A search for new irreducible phase-120 sums of larger weight would only find decomposable sums. The useful matrix search is the compatible placement and joint Gram contribution of these known atoms.

This statement allows repeated terms in the definition of a vanishing sum; the theorem implies that every minimal phase-120 sum in fact has multiplicity one. Equal root values at different matrix columns still represent different indexed choices.

## Why the catalogue is complete at every weight

[Lam and Leung, Corollary 3.2](https://arxiv.org/pdf/math/9511209), reduces a minimal sum by rotation to square-free relative order. At phase 120 that order divides 30. [Christie, Dykema and Klep, Theorem 3.2](https://arxiv.org/html/2008.11268v2#S3), classifies all minimal sums of relative order dividing \(2pq\), without a weight cutoff. Taking \(p=3,q=5\) gives the prime polygons \(R_2,R_3,R_5\) and the following explicit sums.

Let \(\omega=\zeta_3\), \(\theta=\zeta_5\). For every proper, nonempty subset \(J\subset\{0,1,2,3,4\}\), take

\[
 \sum_{j\notin J}\theta^j
 +\sum_{j\in J}(-\omega\theta^j-\omega^2\theta^j)=0.
\]

Their weights are \(5+|J|\), namely six through nine. These formulas, their rotations, and the three prime polygons exhaust the catalogue. The theorem is stronger here than the repository's general Poonen–Rubinstein catalogue through weight twelve; no claim is made about replacing that catalogue at other conductors.

The nine rotation classes have weights \(2,3,5,6,7,7,8,8,9\). The two classes at weight seven correspond to the two cyclic arrangements of a two-element subset of five vertices; their complements give the two weight-eight classes. If Galois automorphisms are also identified, these become seven types. A Galois quotient is not silently used when matching fixed matrix entries.

| Weight | Literal phase-30 profiles | Literal phase-120 profiles |
| --- | ---: | ---: |
| 2 | 15 | 60 |
| 3 | 10 | 40 |
| 5 | 6 | 24 |
| 6 | 30 | 120 |
| 7 | 60 | 240 |
| 8 | 60 | 240 |
| 9 | 30 | 120 |
| Total | 211 | 844 |

The phase-120 instances occupy the four cosets of \(\mu_{30}\): embed a phase-30 profile by multiplying its exponents by four, then add offsets zero through three. All internal rotations are already included in the 211 lower profiles.

## Exact implementation and indexed matching

[`vanishing_known_catalogue.py`](../experiments/scalar_lift/vanishing_known_catalogue.py) exports:

- `minimal_profiles(30 or 120)`: all literal profiles, with no weight truncation.
- `rotation_representatives(30 or 120)`: the nine canonical root multisets.
- `iter_indexed_blocks(residues, modulus=120)`: lazy matches retaining every original column label.
- `root_vectors` and `exact_sum`: pure integer power-basis arithmetic.

The arithmetic uses \(\Phi_{30}(x)=x^8+x^7-x^5-x^4-x^3+x+1\) and \(\Phi_{120}(x)=\Phi_{30}(x^4)\). Every literal sum is checked. All 1,370 proper nonempty subsums of the nine representatives are independently checked nonzero; rotations and field embeddings preserve this minimality. Tests compare the root tables against the repository's independent `CyclotomicField` implementation.

The indexed matcher is deliberately an iterator. Matching a small root profile against a long row with many repeated values can produce many different column subsets. A cap on consumption must be reported as a partial indexed search; it does not change the complete scalar profile catalogue.

## What this supplies to the missing order-1499 source

The admitted phase-120 entry set has two cubic colors, four nonzero fifth colors, and all eight dyadic phases. Ratios of two admitted entries cover the full \(\mu_{120}\), so the catalogue applies to pairwise row-inner-product terms. A zero inner product of length 1499 decomposes into these small atoms, but that decomposition need not be unique. It does not enforce simultaneous orthogonality to other rows, nor compatibility of the ratios at shared positions.

For the 585-to-1499 packet construction, the productive collision key is a **joint matrix or sector-Gram signature**, including both row and column metrics, admitted coefficient domains, and conjugation. A scalar sum collision alone does not complete those coupled equations. Balanced packet replacements can themselves be decomposable as scalar sums; rejecting them merely because they are not irreducible would discard useful matrix constructions.

```sh
PYTHONPATH=src:. ./.native-env/bin/python -m experiments.scalar_lift.vanishing_known_catalogue --check certificates/scalar-lift/vanishing-known-catalogue.json
PYTHONPATH=src:. ./.native-env/bin/python -m pytest -q tests/test_vanishing_known_catalogue.py
```

This module is a complete scalar substrate, not a Butson source or a new Hadamard order.

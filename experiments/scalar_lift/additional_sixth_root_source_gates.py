"""Complete source Grams and finite palette obstructions for fixed BH sources.

These are gates for the fixed T631 P/Q images of exponents0,1,3,4 only.
They exclude neither other representatives at the same orders nor other maps.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np

from experiments.scalar_lift.catalogue_source_acquisition import cubic_gauge
from experiments.scalar_lift.tripled_unreal_scalar_endpoints import source_gram

W19_ROWS = (
    "3011005435321135430", "0013314245151433150", "0014242432413314510",
    "1242124424503511340", "2543204201424153130", "0354504531345341310",
    "5435323011001135430", "4245150013311433150", "2432410014243314510",
    "4424501242123511340", "4201422543204153130", "4531340354505341310",
    "5533215533210000113", "5231535231530013010", "3355123355121300010",
    "1325351325350101041", "2513532513531041100", "3152353152351010401",
    "0000000000004110100",
)


def obstruction():
    """Independently exhaust all3^7 row gauges of the eight-by-eight pattern."""
    d = np.array(((1, 2), (2, 1)), dtype=np.int64)
    e = np.kron(np.eye(4, dtype=np.int64), d)
    count = 0
    first_bad = [0] * 8
    for tail in itertools.product(range(3), repeat=7):
        rows = np.array((0, *tail), dtype=np.int64)
        colors = (e + rows[:, None]) % 3
        bad = [j for j in range(8) if set(map(int, colors[:, j])) == {0, 1, 2}]
        if not bad:
            raise ValueError("claimed obstruction has a compatible column gauge")
        first_bad[bad[0]] += 1
        count += 1
    assert count == 2187
    return {"cubic_exponents": e.tolist(), "all_fixed_first_row_gauges_checked": count,
            "compatible_gauges": 0, "first_failed_column_histogram": first_bad,
            "proof": "Fix row0 shift0 by global row/column compensation. For every other ternary row gauge, some column contains all3 residues; no column shift can omit a residue.",
            "scope": "all sixth-root row/column phase gauges, including permutations of this fixed pattern"}


def square_source(p):
    if p not in (3, 5, 7, 11):
        raise ValueError("finite source checks cover odd primes3,5,7,11")
    residues = {x * x % p for x in range(1, p)}
    paley = np.array([[0 if i == j else (1 if (j - i) % p in residues else -1)
                       for j in range(p)] for i in range(p)], dtype=np.int64)
    i = np.eye(p, dtype=np.int64)
    j = np.ones((p, p), dtype=np.int64)
    assert np.array_equal(paley @ paley.T, p * i - j)
    assert not np.any(paley @ j) and not np.any(j @ paley)
    a, b, c = np.kron(paley, paley), np.kron(j, i), np.kron(i, j)
    # H=a+b+omega*c, omega=zeta6^2. The diagonal1+omega is zeta6.
    e = np.where(a < 0, 3, 0).astype(np.int64)
    e[(c == 1) & (b == 0)] = 2
    e[(c == 1) & (b == 1)] = 1
    assert np.all((a != 0) | (b != 0) | (c != 0))
    return e


def certificate():
    small = obstruction()
    pattern = np.asarray(small["cubic_exponents"], dtype=np.int64)
    records = []
    for p in (5, 7, 11):
        e = square_source(p)
        proof = source_gram(e)
        indices = [p * outer + inner for outer in range(4) for inner in range(2)]
        assert np.array_equal(e[np.ix_(indices, indices)] % 3, pattern)
        records.append({"prime": p, "source_order": p * p, "source_full_Gram": proof,
                        "actual_exponent_alphabet": sorted(map(int, np.unique(e))),
                        "obstruction_row_and_column_indices": indices,
                        "compatible_sixth_root_gauge_exists": False,
                        "rejected_only_for_fixed_T631_palette": True})
    control = square_source(3)
    control_proof = source_gram(control)
    control_gauge = cubic_gauge(control)
    assert control_gauge["exists"]
    gauged = (control + np.asarray(control_gauge["row_shifts"])[:, None]
              + np.asarray(control_gauge["column_shifts"])[None, :]) % 6
    assert set(map(int, np.unique(gauged))) <= {0, 1, 3, 4}
    source_gram(gauged)
    w19 = np.array([[int(x) for x in row] for row in W19_ROWS], dtype=np.int64)
    assert w19.shape == (19, 19)
    w19_proof, w19_gauge = source_gram(w19), cubic_gauge(w19)
    assert not w19_gauge["exists"]
    return {"schema": "complete-fixed-sixth-root-source-palette-gates-v1",
            "target_scalar_frame": {"field": "Q(zeta_6)/Q", "field_degree": 2,
                                    "module_rank": 1262, "rational_dimension": 2524,
                                    "entry_images": {"0": "P", "1": "Q", "3": "-P", "4": "-Q"},
                                    "excluded_images": {"2": "Q-P", "5": "P-Q"},
                                    "frame_gram_scale": 2524,
                                    "trace_metric": [[2, 1], [1, 2]]},
            "prime_square_primary_source": {
                "author": "Ferenc Szollosi",
                "title": "Construction, classification and parametrization of complex Hadamard matrices",
                "url": "https://arxiv.org/pdf/1110.5590", "location": "Theorem1.4.41, printedp30",
                "attribution_in_primary_source": "Craigen and Szollosi",
                "formula": "H=P tensor P+J tensor I+omega*I tensor J, P the Paley character matrix, omega=zeta6^2",
                "actual_palette": [0, 1, 2, 3],
                "cubic_color_pattern": "p diagonal blocks of sizep with diagonal1,offdiagonal2; off-block colors0",
                "general_gate": "For every odd prime p>=5, select4 outer blocks and2 inner indices. This contains the checked8x8 obstruction, so no sixth-root monomial gauge places this fixed source in the admitted palette.",
                "different_BH_representatives_at_same_orders_excluded": False,
                "general_nonuniform_scalar_maps_excluded": False},
            "independent_finite_obstruction": small,
            "complete_prime_square_source_checks": records,
            "positive_control_order9": {"full_Gram": control_proof, "admitted_gauge": control_gauge},
            "fixed_order19_source": {
                "author": "Ferenc Szollosi", "title": "A note on the existence of BH(19,6) matrices",
                "primary_url": "https://ajc.maths.uq.edu.au/pdf/55/ajc_v55_p031.pdf",
                "location": "printedp33, literal W19 exponent matrix",
                "row_strings": list(W19_ROWS), "full_Gram": w19_proof,
                "complete_sixth_root_gauge_test": w19_gauge,
                "scope": "This fixed published W19 only; other order19 Butson matrices are not excluded"},
            "new_endpoint_constructed": False, "new_Hadamard_order_claimed": False,
            "implementation_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = certificate()
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"complete_source_orders_checked": [19, 25, 49, 121],
                      "obstruction_row_gauges_checked": 2187,
                      "order9_positive_control": True,
                      "all_exclusions_are_fixed_source_palette_gates": True}))


if __name__ == "__main__":
    main()

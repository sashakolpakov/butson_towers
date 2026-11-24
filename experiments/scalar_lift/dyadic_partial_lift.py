"""Exact partial sign lifts for 2/3 of phases 24, 48, and 96.

The source is still unknown at the large target orders. This file supplies
both the faithful rational action and its flat frame, with exact checks of
every admitted and excluded phase. See docs/DYADIC_PARTIAL_LIFT.md.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from itertools import product
from pathlib import Path

import numpy as np

from experiments.scalar_lift.dyadic_cyclotomic_action import action_numerator, verify_action
from experiments.scalar_lift.inverse_blocks import known_order_thirteen_source
from hadamard.cgw import is_cgw_exponent_matrix


def local_frame():
    """Enumerate all local sign columns up to negation; they form a Hadamard."""
    matrix = action_numerator(4)
    columns = np.array([
        (1, *tail) for tail in product((-1, 1), repeat=7)
        if np.all(np.abs(matrix @ np.array((1, *tail))) == 2)
    ], dtype=np.int64).T
    if columns.shape != (8, 8) or not np.array_equal(columns @ columns.T, 8 * np.eye(8)):
        raise ValueError("the exhaustive local columns did not form the claimed frame")
    return columns


def sylvester(order):
    if order < 1 or order & (order - 1):
        raise ValueError("a positive power of two is required")
    result = np.ones((1, 1), dtype=np.int64)
    while len(result) < order:
        result = np.block([[result, result], [result, -result]])
    return result


def sign_frame(negacyclic_size):
    """Tensor the local frame in the residue-class decomposition of U."""
    n = negacyclic_size
    if n not in (4, 8, 16):
        raise ValueError("the certified negacyclic sizes are 4, 8, and 16")
    copies = n // 4
    indices = [offset + r + copies * t
               for r in range(copies) for offset in (0, n) for t in range(4)]
    frame = np.empty((2 * n, 2 * n), dtype=np.int64)
    frame[indices, :] = np.kron(sylvester(copies), local_frame())
    return frame


def phase_blocks(negacyclic_size):
    """Generate every S^j P using integer division with a divisibility check."""
    n = negacyclic_size
    matrix = action_numerator(n)
    frame = sign_frame(n)
    blocks = []
    current = frame.copy()
    for _ in range(6 * n):
        blocks.append(current.copy())
        numerator = matrix @ current
        if np.any(numerator % 2):
            raise ValueError("a phase block unexpectedly lost integrality")
        current = numerator // 2
    if not np.array_equal(current, frame):
        raise ValueError("the phase orbit did not close")
    return blocks


def verify_dyadic_source_factorization(exponents, phase):
    """Check the source equations over Z[zeta_d], using x^(d/2)=-1."""
    if phase not in (24, 48, 96):
        raise ValueError("the certified phases are 24, 48, and 96")
    d = phase // 3
    degree = d // 2
    m = len(exponents)
    if m == 0:
        raise ValueError("a nonempty source is required")
    omega = next(e for e in (d, 2 * d) if e % 3 == 1)
    a = [[None] * m for _ in range(m)]
    c = [[None] * m for _ in range(m)]
    for i, row in enumerate(exponents):
        if len(row) != m:
            raise ValueError("the source must be square")
        for j, exponent in enumerate(row):
            if isinstance(exponent, (bool, np.bool_)) or not isinstance(exponent, (int, np.integer)):
                raise ValueError("source exponents must be exact integers")
            e = int(exponent) % phase
            if e % 3 == 0:
                a[i][j] = e // 3
            elif e % 3 == 1:
                c[i][j] = ((e - omega) // 3) % d
            else:
                raise ValueError("the source uses an excluded phase")

    def product_coefficients(first, second, i, j):
        coefficients = np.zeros(degree, dtype=np.int64)
        for left, right in zip(first[i], second[j], strict=True):
            if left is not None and right is not None:
                exponent = (left - right) % d
                coefficients[exponent % degree] += 1 if exponent < degree else -1
        return coefficients

    for i in range(m):
        for j in range(m):
            ac = product_coefficients(a, c, i, j)
            ca = product_coefficients(c, a, i, j)
            if np.any(ca - ac):
                raise ValueError("the dyadic amicability equation failed")
            gram = (product_coefficients(a, a, i, j)
                    + product_coefficients(c, c, i, j) - ac)
            gram[0] -= m if i == j else 0
            if np.any(gram):
                raise ValueError("the dyadic norm equation failed")
    return {
        "dyadic_phase": d,
        "cube_root_exponent_in_full_phase": omega,
        "A_nonzero_entries": sum(e is not None for row in a for e in row),
        "C_nonzero_entries": sum(e is not None for row in c for e in row),
        "CA_star_equals_AC_star": True,
        "AA_star_plus_CC_star_minus_AC_star_equals_mI": True,
        "verification_ring": f"Z[x]/(x^{degree}+1)",
    }


def verify_partial_lift(negacyclic_size):
    n = negacyclic_size
    phase, dimension = 6 * n, 2 * n
    action = verify_action(action_numerator(n), phase)
    blocks = phase_blocks(n)
    identity = np.eye(dimension, dtype=np.int64)
    if not all(np.array_equal(block @ block.T, dimension * identity) for block in blocks):
        raise ValueError("a phase block lost its exact Gram identity")
    admitted = [j for j, block in enumerate(blocks) if np.all(np.abs(block) == 1)]
    expected = [j for j in range(phase) if j % 3 in (0, 1)]
    if admitted != expected:
        raise ValueError("the flat alphabet was different from the claimed two cosets")
    if any(not np.all(blocks[j] % 2 == 0) for j in range(2, phase, 3)):
        raise ValueError("an excluded phase was not entrywise even")
    # A 3-by-3 Butson positive control: diagonal omega, off-diagonal 1,
    # where omega is a primitive cube root in the admitted coset.
    root_exponent = next(j for j in (phase // 3, 2 * phase // 3) if j % 3 == 1)
    control = np.block([[blocks[root_exponent] if i == j else blocks[0]
                         for j in range(3)] for i in range(3)])
    if not np.all(np.abs(control) == 1) or not np.array_equal(
        control @ control.T, 3 * dimension * np.eye(3 * dimension, dtype=np.int64)
    ):
        raise ValueError("the nontrivial source-to-sign control failed")
    # Embed the independently verified catalogue BH(13,6), already gauged
    # into its unreal alphabet. A global phase +2 places its two mod-3
    # cosets in the flat alphabet for every q considered here.
    source13 = [[(n * e + 2) % phase for e in row]
                for row in known_order_thirteen_source()]
    if not is_cgw_exponent_matrix(source13, phase, 13):
        raise ValueError("the embedded order-thirteen source failed exact verification")
    if any(e not in admitted for row in source13 for e in row):
        raise ValueError("the embedded order-thirteen source left the flat alphabet")
    control13 = np.block([[blocks[e] for e in row] for row in source13])
    if not np.all(np.abs(control13) == 1) or not np.array_equal(
        control13 @ control13.T,
        13 * dimension * np.eye(13 * dimension, dtype=np.int64),
    ):
        raise ValueError("the order-thirteen source-to-sign control failed")
    return {
        "phase": phase,
        "block_order": dimension,
        "field_degree": action["faithful_field_degree"],
        "frame": blocks[0].tolist(),
        "frame_gram": dimension,
        "flat_exponents": admitted,
        "flat_exponent_rule": "j mod 3 is 0 or 1",
        "flat_alphabet_size": len(admitted),
        "all_phase_blocks_have_exact_gram": True,
        "excluded_phase_blocks_are_even": True,
        "source_order_three_control": {
            "diagonal_phase_exponent": root_exponent,
            "off_diagonal_phase_exponent": 0,
            "real_output_order": 3 * dimension,
            "signs_and_gram_verified": True,
        },
        "source_order_thirteen_control": {
            "embedded_source_exponents": source13,
            "real_output_order": 13 * dimension,
            "signs_and_gram_verified": True,
            "dyadic_source_factorization": verify_dyadic_source_factorization(source13, phase),
            "binary_sha256": hashlib.sha256(
                ((control13 + 1) // 2).astype(np.uint8).tobytes()
            ).hexdigest(),
        },
    }


def build_report():
    return {
        "scope": "Explicit partial scalar morphisms with flat frames; no large Butson source",
        "lifts": [verify_partial_lift(n) for n in (4, 8, 16)],
        "conditional_targets": [
            {"source_order": 1499, "phase": 96, "real_order": 47968},
            {"source_order": 2999, "phase": 48, "real_order": 47984},
            {"source_order": 2999, "phase": 96, "real_order": 95968},
        ],
        "large_source_constructed": False,
        "new_real_Hadamard_order_constructed": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", type=Path)
    args = parser.parse_args()
    report = build_report()
    if args.check and json.loads(args.check.read_text()) != report:
        raise SystemExit("the partial-lift certificate differs from its exact replay")
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"verified_partial_lifts": [
        {key: item[key] for key in ("phase", "block_order", "flat_alphabet_size")}
        for item in report["lifts"]
    ], "new_real_Hadamard_order_constructed": False}))


if __name__ == "__main__":
    main()

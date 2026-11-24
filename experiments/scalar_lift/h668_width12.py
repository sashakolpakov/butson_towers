"""Rebuild the printed width-12 order-668 example with exact checks.

This is an additional existing-order input, not a new order or an
independent certificate of the paper's inequivalence claim.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from itertools import product
from pathlib import Path

import numpy as np

from hadamard.inverse_lift import is_sign_hadamard

SEED = Path(__file__).with_name("h668_width12_seed.json")


def decode_signs(word, length):
    if len(word) * 4 < length:
        raise ValueError("the hexadecimal word is too short")
    bits = f"{int(word, 16):0{4 * len(word)}b}"
    if any(bit != "0" for bit in bits[length:]):
        raise ValueError("the unused trailing bits must be zero padding")
    return np.array([1 - 2 * int(bit) for bit in bits[:length]], dtype=np.int64)


def group_indices(moduli):
    elements = list(product(*(range(modulus) for modulus in moduli)))
    position = {element: i for i, element in enumerate(elements)}
    difference = np.array([[position[tuple((b - a) % modulus for a, b, modulus
                                          in zip(x, y, moduli, strict=True))]
                            for y in elements] for x in elements], dtype=np.int64)
    negative = [position[tuple(-value % modulus for value, modulus
                               in zip(x, moduli, strict=True))] for x in elements]
    return elements, difference, negative


def goethals_seidel(sequences, moduli):
    _elements, difference, negative = group_indices(moduli)
    a, b, c, d = (sequence[difference] for sequence in sequences)
    return np.block([[a, b[:, negative], c[:, negative], d[:, negative]],
                     [-b[:, negative], a, d.T[:, negative], -c.T[:, negative]],
                     [-c[:, negative], -d.T[:, negative], a, b.T[:, negative]],
                     [-d[:, negative], c.T[:, negative], -b.T[:, negative], a]])


def known_h668_width12():
    """Return the independently checked dense sign matrix and a compact report."""
    seed = json.loads(SEED.read_text())
    core_order, quotient_order, border = 164, 4, 12
    sequences = np.array([decode_signs(word, core_order) for word in seed["core_hex"]])
    row_sums = sequences.sum(axis=1).tolist()
    if row_sums != seed["expected_core_row_sums"]:
        raise ValueError("the transcribed sequences have incorrect row sums")
    beta = sequences.reshape(4, 4, 41).sum(axis=2)
    dq = np.array([decode_signs(word, 16) for word in seed["quotient_frame_columns_D_hex"]]).T
    eq = np.array([decode_signs(word, 16) for word in seed["quotient_frame_rows_E_hex"]])
    quotient_core = goethals_seidel(beta, seed["quotient_moduli"])
    corner_numerator = -(eq @ quotient_core.T @ dq)
    if np.any(corner_numerator % (4 * quotient_order)):
        raise ValueError("the corner formula is not integral")
    corner = corner_numerator // (4 * quotient_order)
    if not np.all(abs(corner) == 1) or not np.array_equal(corner @ corner.T, border * np.eye(border, dtype=np.int64)):
        raise ValueError("the computed corner is not a sign Hadamard matrix")
    pulled_indices = np.array([f * quotient_order + x // 41 for f in range(4) for x in range(core_order)])
    core = goethals_seidel(sequences, seed["group_moduli"])
    matrix = np.block([[corner, eq[:, pulled_indices]], [dq[pulled_indices], core]])
    if matrix.shape != (668, 668) or not np.all(abs(matrix) == 1):
        raise ValueError("the assembled matrix has invalid dimensions or entries")
    if not np.array_equal(matrix @ matrix.T, 668 * np.eye(668, dtype=np.int64)):
        raise ValueError("the full integer Gram identity failed")
    # A separate bit-distance validator does not share the array multiplication.
    if not is_sign_hadamard(matrix.tolist()):
        raise ValueError("the independent bit-distance Gram check failed")
    report = {
        "scope": "Exact reconstruction of the published width-12 order-668 seed; no new order",
        "source_url": seed["source_url"], "source_identifier": seed["source_identifier"],
        "source_author_as_printed": seed["source_author_as_printed"],
        "seed_sha256": hashlib.sha256(SEED.read_bytes()).hexdigest(),
        "order": 668, "group_moduli": seed["group_moduli"], "border_width": border,
        "core_row_sums": row_sums, "quotient_coset_sums": beta.tolist(),
        "derived_corner": corner.tolist(),
        "matrix_binary_sha256": hashlib.sha256(((matrix + 1) // 2).astype(np.uint8).tobytes()).hexdigest(),
        "full_integer_and_independent_bit_distance_Gram_checks_passed": True,
        "inequivalence_independently_verified": False,
        "scalar_lift_pair_or_new_order_constructed": False,
    }
    return matrix, report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", type=Path)
    args = parser.parse_args()
    _matrix, report = known_h668_width12()
    if args.check is not None and json.loads(args.check.read_text()) != report:
        raise SystemExit("the width-12 order-668 report differs")
    if args.output is not None:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"order": 668, "border_width": 12,
                      "matrix_binary_sha256": report["matrix_binary_sha256"], "exact_Gram_verified": True}))


if __name__ == "__main__":
    main()

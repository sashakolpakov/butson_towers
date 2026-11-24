"""Exact finite-field audit of the fifth-root repeated-border seed choices.

    Both the initial and one-rotation seeds are now excluded by joint
    constraints. Their rowwise feasibility facts below remain valid.
"""

from __future__ import annotations

import argparse
import json
from itertools import product
from pathlib import Path

import numpy as np


def field_multiply(a, b):
    """F16 in the power basis modulo x^4+x^3+x^2+x+1."""
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & 16:
            a ^= 31
    return result


def field_power(a, exponent):
    result = 1
    for _ in range(exponent):
        result = field_multiply(result, a)
    return result


MULTIPLY = np.array([[field_multiply(a, b) for b in range(16)] for a in range(16)],
                    dtype=np.uint8)
CONJUGATE = np.array([field_power(a, 4) for a in range(16)], dtype=np.uint8)
TRACE = np.array([a ^ field_power(a, 2) ^ field_power(a, 4) ^ field_power(a, 8)
                  for a in range(16)], dtype=np.uint8)
FIFTH_ROOTS = [field_power(2, k) for k in range(5)]
SUM_ALPHABET = sorted({a ^ b for a in FIFTH_ROOTS[1:] for b in FIFTH_ROOTS[1:]})
TRACE_KERNEL = [a for a in range(16) if TRACE[a] == 0]


def matrix_multiply(a, b):
    a, b = np.asarray(a, dtype=np.uint8), np.asarray(b, dtype=np.uint8)
    result = np.zeros((a.shape[0], b.shape[1]), dtype=np.uint8)
    for k in range(a.shape[1]):
        result ^= MULTIPLY[a[:, k, None], b[None, k, :]]
    return result


def dagger(matrix):
    return CONJUGATE[np.asarray(matrix, dtype=np.uint8).T]


def chirped_five():
    return np.array([[FIFTH_ROOTS[(3 * (i + j) ** 2 + 1) % 5]
                      for j in range(5)] for i in range(5)], dtype=np.uint8)


def binary_rank(vectors):
    pivots = {}
    for value in vectors:
        value = int(value)
        while value:
            pivot = value.bit_length() - 1
            if pivot in pivots:
                value ^= pivots[pivot]
            else:
                pivots[pivot] = value
                break
    return len(pivots)


def packed_row(row):
    return sum(int(value) << (4 * j) for j, value in enumerate(row))


def boundary_image_audit(rows):
    rows = np.asarray(rows, dtype=np.uint8)
    # These three allowed sums are an F2 basis of ker(absolute trace).
    basis = [6, 7, 10]
    generators = [packed_row(MULTIPLY[b, row]) for row in rows for b in basis]
    rank = binary_rank(generators)
    targets = matrix_multiply(rows, dagger(chirped_five()))
    compatible = [binary_rank(generators + [packed_row(target)]) == rank for target in targets]
    return {"F2_image_rank": rank, "F2_codomain_dimension": 20,
            "representative_row_count": len(rows),
            "compatible_representative_targets": sum(compatible),
            "all_representative_targets_compatible": all(compatible)}


def scalar_lookup():
    """Exact sums, with all possible Hermitian norms, using three columns."""
    lookup = {d: {} for d in range(16)}
    for a, b, c in product(SUM_ALPHABET, repeat=3):
        value = a ^ b ^ int(MULTIPLY[FIFTH_ROOTS[3], c])
        norm = field_power(a, 5) ^ field_power(b, 5) ^ field_power(c, 5)
        lookup[value].setdefault(norm, (a, b, c))
    if any(len(options) < 3 for options in lookup.values()):
        raise ValueError("some scalar target has too few possible norms")
    return lookup


def row_witness(target, desired_norm=1):
    """A finite lookup/4-state DP, not a simultaneous matrix solver."""
    d = chirped_five()
    coordinates = matrix_multiply(np.asarray(target)[None, :], dagger(d))[0]
    lookup = scalar_lookup()
    paths = {0: []}
    for coordinate in coordinates:
        next_paths = {}
        for current_norm, path in paths.items():
            for norm, triple in lookup[int(coordinate)].items():
                next_paths.setdefault(current_norm ^ norm, path + [triple])
        paths = next_paths
    if desired_norm not in paths:
        raise ValueError("row norm could not be attained")
    triples = paths[desired_norm]
    row = np.array([triple[k] for k in range(3) for triple in triples], dtype=np.uint8)
    active_t = np.vstack([d, d, MULTIPLY[FIFTH_ROOTS[3], d]])
    if (not np.array_equal(matrix_multiply(row[None, :], active_t)[0], target)
            or matrix_multiply(row[None, :], dagger(row[None, :]))[0, 0] != desired_norm
            or any(int(a) not in SUM_ALPHABET for a in row)):
        raise ValueError("row-local witness failed")
    return row.tolist()


def odd_unitary_control():
    lam = next(a for a in range(2, 16) if field_power(a, 4) == a)
    c = np.array([[int((j - i) % 7 in (1, 2, 4)) for j in range(7)]
                  for i in range(7)], dtype=np.uint8)
    e7 = MULTIPLY[lam, c] ^ MULTIPLY[field_power(lam, 2), c.T]
    e2 = np.array([[lam, field_power(lam, 2)], [field_power(lam, 2), lam]], dtype=np.uint8)
    for matrix in (e2, e7):
        if (not np.array_equal(matrix_multiply(matrix, dagger(matrix)),
                               np.eye(len(matrix), dtype=np.uint8))
                or any(int(value) not in SUM_ALPHABET for value in matrix.flat)):
            raise ValueError("finite-field unitary control failed")
    return {"E2": e2.tolist(), "E7": e7.tolist(), "full_unitary_Grams_verified": True,
            "all_entries_in_exact_sum_alphabet": True,
            "odd_orders_747_and_1497_available_by_E7_plus_E2_direct_sums": True,
            "compatibility_with_prescribed_T_claimed": False}


def build_report():
    if ([a for a in TRACE_KERNEL if a != 1] != SUM_ALPHABET
            or any(TRACE[a] != 1 for a in FIFTH_ROOTS[1:])):
        raise ValueError("fifth-root trace alphabet failed")
    d = chirped_five()
    if (not np.array_equal(matrix_multiply(d, dagger(d)), np.eye(5, dtype=np.uint8))
            or not np.all(np.bitwise_xor.reduce(d, axis=1) == FIFTH_ROOTS[1])):
        raise ValueError("chirped D reduction failed")
    ones = np.ones((1, 5), dtype=np.uint8)
    initial_rows = np.vstack([d, MULTIPLY[FIFTH_ROOTS[1], ones]])
    w = matrix_multiply(initial_rows, dagger(d))
    if np.any(w > 1) or not np.all(matrix_multiply(w, np.ones((5, 5), dtype=np.uint8)) == 1):
        raise ValueError("the unrotated binary-descent contradiction failed")
    rotated_rows = np.vstack([d, MULTIPLY[FIFTH_ROOTS[3], d],
                              MULTIPLY[FIFTH_ROOTS[1], ones]])
    b_rotation_rows = np.vstack([initial_rows, MULTIPLY[FIFTH_ROOTS[2], ones]])
    initial = boundary_image_audit(initial_rows)
    rotated = boundary_image_audit(rotated_rows)
    b_rotated = boundary_image_audit(b_rotation_rows)
    if (initial["F2_image_rank"] != 15 or initial["compatible_representative_targets"]
            or rotated["F2_image_rank"] != 20
            or not rotated["all_representative_targets_compatible"]
            or b_rotated["F2_image_rank"] != 16
            or not b_rotated["all_representative_targets_compatible"]):
        raise ValueError("variant linear preflight failed")
    lookup = scalar_lookup()
    targets = matrix_multiply(rotated_rows, dagger(d))
    witnesses = [row_witness(target) for target in targets]
    return {
        "scope": "Exact residue-field necessity and rowwise feasibility only; no source completion",
        "field": "F2[x]/(x^4+x^3+x^2+x+1), hexadecimal/power-basis integer encoding",
        "fifth_roots": FIFTH_ROOTS, "absolute_trace_kernel": TRACE_KERNEL,
        "exact_X_plus_Y_alphabet": SUM_ALPHABET,
        "D_mod_two": d.tolist(), "D_unitary_verified": True,
        "D_row_and_column_sum": FIFTH_ROOTS[1],
        "initial_seed": {**initial, "W_equals_T_D_star_is_binary": True,
                         "W_J_is_all_ones": True, "excluded_by_trace": True},
        "one_D_block_rotated_by_zeta5_cubed": {
            **rotated, "exact_sum_alphabet_reaches_all_scalar_targets": True,
            "minimum_number_of_attainable_scalar_norms": min(map(len, lookup.values())),
            "scalar_sum_and_norm_lookup": {str(a): {str(n): list(t) for n, t in choices.items()}
                                           for a, choices in lookup.items()},
            "rowwise_boundary_and_unit_norm_witnesses": witnesses,
            "representative_targets": targets.tolist(),
            "all_11_representative_rows_verified": True,
            "simultaneous_row_orthogonality_verified": False,
            "trace_obstruction_removed": True,
            "later_joint_trace_unitarity_obstruction_excludes_this_seed": True,
            "later_proof": "docs/FIFTH_GLOBAL_RESIDUE.md",
        },
        "one_B_block_rotated_from_coset1_to_coset2": {
            **b_rotated, "later_joint_trace_unitarity_obstruction_excludes_this_seed": True},
        "odd_unitary_counterexample": odd_unitary_control(),
        "applies_to_Dita_dyadic_B8_phases": [8, 16],
        "dyadic_B8_reduction": "all entries equal1, independently of its primitive cross ratio",
        "high_source_orders": [1499, 2999],
        "all_fifth_root_bordered_sources_excluded": False,
        "new_Hadamard_matrix_constructed": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", type=Path)
    args = parser.parse_args()
    report = build_report()
    if args.check and json.loads(args.check.read_text()) != report:
        raise SystemExit("fifth boundary audit differs from certificate")
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"initial_seed_excluded": True, "rotated_row_image_rank": 20,
                      "simultaneous_completion_constructed": False}))


if __name__ == "__main__":
    main()

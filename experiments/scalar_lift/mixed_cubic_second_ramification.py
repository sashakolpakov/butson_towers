"""Exact second-bit lifts of a repeated cubic-support template.

This module works only in F4[pi]/(pi^4). A positive result is a
higher-precision root residue, not a characteristic-zero Butson.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from experiments.scalar_lift.cubic_core_support_lift import (
    boundary_data,
    large_problem,
    phase_coordinates,
    template,
)
from experiments.scalar_lift.mixed_cubic_ramification import CONJUGATE, MULTIPLY


def product(left, right):
    """F4 multiplication via four exact integer matrix products."""
    left, right = np.asarray(left, dtype=np.int64), np.asarray(right, dtype=np.int64)
    a, b, c, d = left & 1, left >> 1, right & 1, right >> 1
    bd = b @ d
    return ((a @ c + bd) % 2) ^ (2 * ((a @ d + b @ c + bd) % 2))


def derivative(fixed, variable):
    return product(variable, CONJUGATE[fixed].T) ^ product(fixed, CONJUGATE[variable].T)


def sparse_derivative(fixed, entries):
    output = np.zeros((len(fixed), len(fixed)), dtype=np.int64)
    conjugate = CONJUGATE[fixed]
    for i, j, value in entries:
        if value:
            output[i] ^= MULTIPLY[value, conjugate[:, j]]
            output[:, i] ^= MULTIPLY[fixed[:, j], CONJUGATE[value]]
    return output


def descriptor(i, j, copies, width=10):
    stop = copies * width
    if i < stop and j < stop:
        return ("same" if i // width == j // width else "different", i % width, j % width)
    if i < stop:
        return ("to_filler", i % width)
    if j < stop:
        return ("from_filler", j % width)
    return ("filler_diagonal" if i == j else "filler_off_diagonal",)


def phase_variables(copies=2, fillers=2, width=10):
    size = copies * width + fillers
    groups = {}
    for i in range(size):
        for j in range(size):
            groups.setdefault(descriptor(i, j, copies, width), []).append((i, j))
    return list(groups), list(groups.values())


def encode_blocks(*blocks):
    values = np.concatenate([np.asarray(block, dtype=np.uint8).ravel() for block in blocks])
    bits = ((values[:, None] >> np.array([0, 1], dtype=np.uint8)) & 1).ravel()
    return int.from_bytes(np.packbits(bits, bitorder="little").tobytes(), "little")


def solve_binary_columns(columns, target):
    pivots = {}
    for index, column in enumerate(columns):
        value, combination = column, 1 << index
        while value:
            pivot = value.bit_length() - 1
            if pivot in pivots:
                basis, used = pivots[pivot]
                value ^= basis
                combination ^= used
            else:
                pivots[pivot] = (value, combination)
                break
    value, combination = target, 0
    while value:
        pivot = value.bit_length() - 1
        if pivot not in pivots:
            return None, {"rank": len(pivots), "augmented_rank": len(pivots) + 1,
                          "nonzero_reduced_target_bit": pivot}
        basis, used = pivots[pivot]
        value ^= basis
        combination ^= used
    return combination, {"rank": len(pivots), "augmented_rank": len(pivots),
                         "nullity": len(columns) - len(pivots)}


def polynomial_gram_from_roots(colors, dyadic_exponents_mod4):
    """Independent root/inverse-root expansion, including conjugating pi."""
    colors = np.asarray(colors, dtype=np.int64)
    b = np.asarray(dyadic_exponents_mod4, dtype=np.int64) % 4
    roots = np.array([1, 2, 3], dtype=np.int64)[colors]
    inverse_roots = CONJUGATE[roots]
    inverse_b = (-b) % 4
    coefficients = [roots, MULTIPLY[roots, b & 1], MULTIPLY[roots, b >> 1],
                    MULTIPLY[roots, (b & 1) * (b >> 1)]]
    inverse = [inverse_roots, MULTIPLY[inverse_roots, inverse_b & 1],
               MULTIPLY[inverse_roots, inverse_b >> 1],
               MULTIPLY[inverse_roots, (inverse_b & 1) * (inverse_b >> 1)]]
    gram = []
    for degree in range(4):
        result = np.zeros((len(colors), len(colors)), dtype=np.int64)
        for first in range(degree + 1):
            result ^= product(coefficients[first], inverse[degree - first].T)
        gram.append(result)
    return gram


def representative_problem(c0, b0, a0, border_support, border_second_bits):
    c0, b0, a0, border_support, border_second_bits = (
        np.asarray(value, dtype=np.int64) for value in (c0, b0, a0, border_support, border_second_bits))
    width, copies, fillers = len(c0), 2, 2
    s, t = copies * width + fillers, len(border_support)
    correction = np.zeros((s, s), dtype=np.int64)
    b = np.zeros_like(correction)
    for index in range(copies):
        part = slice(index * width, (index + 1) * width)
        correction[part, part] = np.ones_like(c0) ^ c0
        b[part, part] = b0
    correction[2 * width:, 2 * width:] = np.eye(fillers, dtype=np.int64)
    c = np.ones_like(correction) ^ correction
    w = b ^ (c @ b.T @ c % 2)
    active = np.concatenate([a0] * copies + [np.zeros(fillers, dtype=np.int64)])
    # The border fixed vector is passed separately through the support:
    # choose the same one used in the positive first-layer certificate.
    return {"s": s, "t": t, "width": width, "C": c, "W": w, "active": active,
            "border_support": border_support, "border_second_bits": border_second_bits}


def solve_second_bits(problem, border_vector):
    s, t, width = problem["s"], problem["t"], problem["width"]
    a = problem["active"][:, None] * np.asarray(border_vector, dtype=np.int64)[None, :]
    ones = np.ones((s, s), dtype=np.int64)
    y0 = ones ^ (3 * problem["C"])
    u0 = np.ones((s, t), dtype=np.int64) ^ (3 * a)
    d0 = np.ones((t, t), dtype=np.int64) ^ (3 * problem["border_support"])
    y1 = MULTIPLY[y0, problem["W"]]
    u1 = 2 * a
    zero = np.zeros((s, s), dtype=np.int64)
    h0 = np.block([[ones, y0, u0], [y0, ones, u0], [u0.T, u0.T, d0]])
    h1 = np.block([[zero, y1, u1], [y1, zero, u1], [u1.T, u1.T, np.zeros((t, t), dtype=np.int64)]])
    if not np.array_equal(product(h0, CONJUGATE[h0].T), np.eye(len(h0), dtype=np.int64)):
        raise ValueError("the representative zeroth Gram failed")
    if np.any(derivative(h0, h1)):
        raise ValueError("the representative first Gram failed")
    fixed_h2 = np.zeros_like(h0)
    fixed_h2[-t:, -t:] = MULTIPLY[d0, problem["border_second_bits"]]
    fixed = product(h1, CONJUGATE[h1].T) ^ product(h0, CONJUGATE[h1].T)
    target2 = fixed ^ derivative(h0, fixed_h2)
    target3 = fixed ^ derivative(h1, fixed_h2)
    names, groups = phase_variables(width=width)
    columns = []
    variable_names = []
    for core in ("X", "Y"):
        for name, positions in zip(names, groups, strict=True):
            full_positions = ([(i, j) for i, j in positions] + [(s + i, s + j) for i, j in positions]
                              if core == "X" else
                              [(i, s + j) for i, j in positions] + [(s + i, j) for i, j in positions])
            second = [(i, j, int(h0[i, j])) for i, j in full_positions]
            third = [(i, j, int(h1[i, j])) for i, j in full_positions]
            degree2 = sparse_derivative(h0, second)
            degree3 = sparse_derivative(h0, third) ^ sparse_derivative(h1, second)
            columns.append(encode_blocks(degree2, degree3))
            variable_names.append((core, *name))
    target = encode_blocks(target2, target3)
    solution, rank = solve_binary_columns(columns, target)
    digest = hashlib.sha256()
    byte_width = (4 * len(h0) ** 2 + 7) // 8
    for vector in [*columns, target]:
        digest.update(vector.to_bytes(byte_width, "little"))
    report = {"representative_core_order": s, "representative_source_order": len(h0),
              "number_second_bit_variables": len(columns), "exact_binary_system": rank,
              "binary_system_sha256": digest.hexdigest(),
              "second_bit_slice": "independent same/different template-copy and filler blocks for X and Y",
              "status": "SOLVED" if solution is not None else "INCONSISTENT_IN_THIS_EQUIVARIANT_SLICE"}
    if solution is None:
        return report, None
    assignments = {name: (solution >> index) & 1 for index, name in enumerate(variable_names)}
    x2, y2 = np.zeros((s, s), dtype=np.int64), np.zeros((s, s), dtype=np.int64)
    for name, positions in zip(names, groups, strict=True):
        for i, j in positions:
            x2[i, j] = assignments[("X", *name)]
            y2[i, j] = assignments[("Y", *name)]
    first_bits = np.block([[zero, problem["W"], a], [problem["W"], zero, a],
                           [a.T, a.T, np.zeros((t, t), dtype=np.int64)]])
    second_bits = np.block([[x2, y2, np.zeros((s, t), dtype=np.int64)],
                            [y2, x2, np.zeros((s, t), dtype=np.int64)],
                            [np.zeros((t, s), dtype=np.int64), np.zeros((t, s), dtype=np.int64),
                             problem["border_second_bits"]]])
    colors = (h0 == 2).astype(np.int64)
    grams = polynomial_gram_from_roots(colors, first_bits + 2 * second_bits)
    if not np.array_equal(grams[0], np.eye(len(h0), dtype=np.int64)) or any(np.any(g) for g in grams[1:]):
        raise ValueError("the solved second bits fail independent inverse-root polynomial Gram verification")
    report.update({"nonzero_variable_indices": [index for index in range(len(columns)) if (solution >> index) & 1],
                   "all_four_representative_Gram_coefficients_verified_independently": True,
                   "representative_X_second_bits": x2.tolist(), "representative_Y_second_bits": y2.tolist()})
    return report, assignments


def verify_actual_large_lift(problem, assignments):
    """Check complete columns for representatives of every row-pair orbit.

    Simultaneous permutations of template copies, or of inactive
    fillers, preserve every residue entry. A row pair uses at most two
    copies and two fillers. Thus the representatives cover all actual
    row-pair types, including equality and the two different core halves.
    """
    s = len(problem["C"])
    groups, fillers = problem["template_groups"], problem["inactive_fillers"]
    if len(groups) < 2 or len(fillers) < 2:
        raise ValueError("the orbit certificate needs two copies and two fillers")
    width, t = len(groups[0]), len(problem["border_C"])
    local_rows = groups[0] + groups[1] + fillers[:2]
    canonical = np.empty(s, dtype=np.int64)
    for index, group in enumerate(groups):
        canonical[group] = np.arange(index * width, (index + 1) * width)
    canonical[fillers] = np.arange(len(groups) * width, s)
    xbits = np.array([[assignments[("X", *descriptor(int(canonical[i]), int(canonical[j]), len(groups), width))]
                       for j in range(s)] for i in local_rows], dtype=np.int64)
    ybits = np.array([[assignments[("Y", *descriptor(int(canonical[i]), int(canonical[j]), len(groups), width))]
                       for j in range(s)] for i in local_rows], dtype=np.int64)
    u_colors, u_dyadic = phase_coordinates(problem["U_exponents"], problem["phase"])
    d_colors, d_dyadic = phase_coordinates(problem["D_exponents"], problem["phase"])
    if not np.array_equal(u_colors, problem["A"]) or np.any(u_dyadic % 4 != problem["A"]):
        raise ValueError("the actual boundary phases differ at precision pi^4")
    c, w = problem["C"][local_rows], problem["W"][local_rows]
    a = u_colors[local_rows]
    zeros = np.zeros_like(c)
    colors = np.block([[zeros, c, a], [c, zeros, a],
                       [u_colors.T, u_colors.T, d_colors]])
    core_second_x, core_second_y = 2 * xbits, 2 * ybits
    dyadic = np.block([[core_second_x, w + core_second_y, u_dyadic[local_rows] % 4],
                       [w + core_second_y, core_second_x, u_dyadic[local_rows] % 4],
                       [u_dyadic.T % 4, u_dyadic.T % 4, d_dyadic % 4]])
    grams = polynomial_gram_from_roots(colors, dyadic)
    size = 2 * s + t
    if not np.array_equal(grams[0], np.eye(len(colors), dtype=np.int64)) or any(np.any(g) for g in grams[1:]):
        raise ValueError("the actual large root residue fails one of its four Gram coefficients")
    # Verify the group data, rather than assuming that every row uses
    # the same local support. The construction's global C=J+T has
    # nontrivial T only on these exact template and singleton blocks.
    small = template()
    for group in groups:
        where = np.ix_(group, group)
        if (not np.array_equal(problem["C"][where], small["C"])
                or not np.array_equal(problem["W"][where], small["W"])
                or not np.array_equal(problem["active"][group], small["a"])):
            raise ValueError("a full-size template group is not a copy of the certified local template")
    if int(np.count_nonzero(problem["W"])) != len(groups) * int(np.count_nonzero(small["W"])):
        raise ValueError("the first phase bits have an undeclared off-template term")
    expected_c_zeros = len(groups) * int(np.count_nonzero(1 ^ small["C"])) + len(fillers)
    if int(np.count_nonzero(1 ^ problem["C"])) != expected_c_zeros:
        raise ValueError("the cubic support has an undeclared off-template correction")
    return {"source_order": size, "phase": problem["phase"],
            "template_copies": len(groups), "inactive_fillers": len(fillers),
            "all_actual_columns_used_in_Gram": colors.shape[1],
            "representative_source_rows": local_rows + [s + i for i in local_rows] + list(range(2 * s, size)),
            "number_ordered_representative_row_pairs_checked": len(colors) ** 2,
            "all_source_row_pairs_covered_by_exact_permutation_symmetry": size * (size + 1) // 2,
            "all_four_actual_Gram_coefficients_verified": True,
            "every_entry_is_the_reduction_of_an_admitted_root": True,
            "representative_colors_and_dyadic_bits_sha256":
            hashlib.sha256(colors.astype(np.uint8).tobytes() + dyadic.astype(np.uint8).tobytes()).hexdigest(),
            "characteristic_zero_Gram_not_claimed": True}


def certificate():
    small = template()
    reports = []
    for order, phase in ((1499, 96), (2999, 48), (2999, 96)):
        boundary = boundary_data(order, phase)
        representative = representative_problem(small["C"], small["B"], small["a"],
                                                 boundary["border_C"], boundary["border_next_bits"])
        model_report, assignments = solve_second_bits(representative, boundary["boundary_vector"])
        if assignments is None:
            raise ValueError("the previously certified positive second-bit slice no longer solves")
        expected_bit = int(order == 2999)
        if set(assignments.values()) != {expected_bit}:
            raise ValueError("the positive second-bit solution changed")
        large = large_problem(order, phase)
        actual_report = verify_actual_large_lift(large, assignments)
        reports.append({"source_order": order, "phase": phase,
                        "conditional_real_order": order * (phase // 3),
                        "core_X_and_Y_second_bits_are_constant": expected_bit,
                        "boundary_second_bits_are_zero": True,
                        "border_second_bits": boundary["border_next_bits"].tolist(),
                        "exact_linear_model": model_report, "actual_large_verification": actual_report,
                        "unitary_root_residue_modulo_pi_fourth": True,
                        "characteristic_zero_source_constructed": False})
    return {"schema": "mixed-cubic-second-ramification-v1",
            "coefficient_ring": "F4[pi]/(pi^4), pi=1-zeta_d, d in{16,32}",
            "template": {key: value.tolist() for key, value in small.items()},
            "first_layer_source": "certificates/scalar-lift/cubic-core-support-lift.json",
            "fixed_first_bits": "X:0; Y:W; upper/lower boundary:A; border:0",
            "entry_second_and_third_coefficients": "H2=H0*b1, H3=H0*b0*b1",
            "conjugation_of_pi": "pi+pi^2+pi^3 modulo pi^4",
            "degree2_equation": "H2 H0* + H1 H1* + H0(H2*+H1*)=0",
            "degree3_equation": "H3 H0* + H2 H1* + H1(H2*+H1*) + H0(H3*+H1*)=0",
            "routes": reports, "native_probe_run": False,
            "new_characteristic_zero_Butson_or_real_Hadamard_claimed": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--output", type=Path)
    action.add_argument("--check", type=Path)
    args = parser.parse_args()
    report = certificate()
    if args.check:
        if json.loads(args.check.read_text()) != report:
            raise ValueError("the second-ramification certificate failed replay")
    else:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"root_residues_verified_modulo_pi_fourth": [[item["source_order"], item["phase"]]
                                                               for item in report["routes"]],
                      "characteristic_zero_source_constructed": False, "native_probe_run": False}))


if __name__ == "__main__":
    main()

"""Joint third core/boundary phase-bit audit for the variable-frame branch."""

from __future__ import annotations

import argparse
import hashlib
import json
from itertools import product as cartesian_product
from pathlib import Path

import numpy as np

from experiments.scalar_lift.cubic_core_support_lift import boundary_data, phase_coordinates
from experiments.scalar_lift.fifth_source_assembly import exact_column_gram, root_vectors
from experiments.scalar_lift.labelled_copy_ramification import (
    complete_binary_solution_space,
    copy_descriptor,
    representative_data,
    variable_groups,
)
from experiments.scalar_lift.mixed_cubic_ramification import MULTIPLY
from experiments.scalar_lift.mixed_cubic_second_ramification import encode_blocks
from experiments.scalar_lift.mixed_cubic_third_ramification import (
    gram_coefficients,
    root_coefficients,
)
from experiments.scalar_lift.variable_cubic_frame_bits import variable_problem


def short_interval_gauges(matrix, modulus=8, width=3):
    matrix = np.asarray(matrix, dtype=np.int64) % modulus
    checked, accepted = 0, 0
    example = None
    for other_rows in cartesian_product(range(modulus), repeat=len(matrix) - 1):
        rows = np.array((0, *other_rows), dtype=np.int64)
        for outputs in cartesian_product(range(width + 1), repeat=matrix.shape[1]):
            columns = np.array(outputs) - matrix[0]
            checked += 1
            if np.all((matrix + rows[:, None] + columns[None, :]) % modulus <= width):
                accepted += 1
                if example is None:
                    example = {"row_gauges": rows.tolist(), "column_gauges": (columns % modulus).tolist()}
    return {"modulus": modulus, "interval_width": width, "normalized_gauges_checked": checked,
            "admissible_gauges": accepted, "first_admissible_gauge": example}


def variable_representative(order, phase):
    data, boundary = representative_data(order, phase), boundary_data(order, phase)
    s = data["s"]
    labels = np.array([0] * 30 + [1] * 30 + [0] * 2, dtype=np.int64)
    shift = 2 * int(order == 2999)
    x = np.full((s, s), shift, dtype=np.int64)
    y = (shift + data["W"] * (1 + 2 * labels[:, None])) % 4
    u = data["A"] * (1 + 2 * labels[:, None])
    _, db = phase_coordinates(boundary["D_exponents"], phase)
    colors = (data["H0"] == 2).astype(np.int64)
    exponents = np.block([[x, y, u], [y, x, u], [u.T, u.T, db % 8]])
    return data, boundary, colors, exponents


def joint_model(order, phase):
    data, boundary, colors, base = variable_representative(order, phase)
    s = data["s"]
    forward, inverse = root_coefficients(colors, base)
    base_grams = gram_coefficients(colors, base)
    if not np.array_equal(base_grams[0], np.eye(len(colors), dtype=np.int64)) or any(np.any(g) for g in base_grams[1:4]):
        raise ValueError("the variable-boundary branch failed its known pi^4 Gram")
    names, groups = variable_groups()
    variables = []
    for core in ("X", "Y"):
        for name, positions in zip(names, groups, strict=True):
            positions = ([(i, j) for i, j in positions] + [(s + i, s + j) for i, j in positions]
                         if core == "X" else
                         [(i, s + j) for i, j in positions] + [(s + i, j) for i, j in positions])
            variables.append(((core, *name), positions))
    for label in (0, 1):
        for pair in range(3):
            positions = []
            for copy in range(3 * label, 3 * label + 3):
                for local in (2 * pair, 2 * pair + 1):
                    row = 10 * copy + local
                    for column in np.flatnonzero(boundary["boundary_vector"]):
                        j = 2 * s + int(column)
                        positions.extend(((row, j), (s + row, j), (j, row), (j, s + row)))
            variables.append((("boundary_pair", label, pair), positions))
    columns = []
    for _, positions in variables:
        blocks = []
        for remainder in range(4):
            block = np.zeros((len(colors), len(colors)), dtype=np.int64)
            for first in range(remainder + 1):
                other = remainder - first
                for i, j in positions:
                    value, inverse_value = forward[first][i, j], inverse[first][i, j]
                    if value:
                        block[i] ^= MULTIPLY[value, inverse[other][:, j]]
                    if inverse_value:
                        block[:, i] ^= MULTIPLY[forward[other][:, j], inverse_value]
            blocks.append(block)
        columns.append(encode_blocks(*blocks))
    target = encode_blocks(*base_grams[4:])
    particular, kernel, rank = complete_binary_solution_space(columns, target)
    digest = hashlib.sha256()
    byte_width = (8 * len(colors) ** 2 + 7) // 8
    for vector in [*columns, target]:
        digest.update(vector.to_bytes(byte_width, "little"))
    report = {"source_order": order, "phase": phase, "core_third_bit_variables": 1284,
              "physical_boundary_pair_third_bit_variables": 6, "exact_binary_system": rank,
              "binary_system_sha256": digest.hexdigest(),
              "complete_kernel_basis_nonzero_variable_indices":
              [[i for i in range(len(variables)) if (word >> i) & 1] for word in kernel],
              "status": "SOLVED" if particular is not None else "INCONSISTENT_IN_THIS_JOINT_LABELLED_SLICE"}
    if particular is None:
        report["general_variable_boundary_source_excluded"] = False
        return report, None, None
    if len(kernel) > 8:
        raise ValueError("the affine space exceeds this bounded exhaustive anti-collapse audit")
    candidates, selected = [], None
    border_column = int(np.flatnonzero(boundary["boundary_vector"])[0])
    minor_indices = [0, 30, 2 * s + border_column]
    for mask in range(1 << len(kernel)):
        word = particular
        for k, vector in enumerate(kernel):
            if (mask >> k) & 1:
                word ^= vector
        exponents = base.copy()
        for index, (_, positions) in enumerate(variables):
            if (word >> index) & 1:
                for i, j in positions:
                    exponents[i, j] += 4
        exponents %= 8
        grams = gram_coefficients(colors, exponents)
        if not np.array_equal(grams[0], np.eye(len(colors), dtype=np.int64)) or any(np.any(g) for g in grams[1:]):
            raise ValueError("a joint affine candidate failed inverse-root polynomial replay")
        minor = exponents[np.ix_(minor_indices, minor_indices)]
        audit = short_interval_gauges(minor)
        candidate = {"kernel_mask": mask, "all_eight_representative_Gram_coefficients_verified": True,
                     "minor_indices": minor_indices, "dyadic_minor_mod8": minor.tolist(),
                     "minor_short_interval_gauge_audit": audit}
        candidates.append(candidate)
        if selected is None or (audit["admissible_gauges"] == 0
                                and selected[2]["minor_short_interval_gauge_audit"]["admissible_gauges"] != 0):
            selected = (word, exponents, candidate)
    report["complete_affine_candidate_audit"] = candidates
    report["selected_kernel_mask"] = selected[2]["kernel_mask"]
    report["selected_minor_excludes_every_width3_interval_gauge_mod8"] = (
        selected[2]["minor_short_interval_gauge_audit"]["admissible_gauges"] == 0)
    assignments = {name: (selected[0] >> i) & 1 for i, (name, _) in enumerate(variables)}
    return report, assignments, selected[1]


def packed_bits(matrix):
    matrix = np.asarray(matrix, dtype=np.uint8)
    return [[int.from_bytes(row.tobytes(), "little") for row in np.packbits((matrix >> bit) & 1, axis=1, bitorder="little")]
            for bit in (0, 1)]


def full_column_bitset_gram(colors, exponents):
    forward, inverse = root_coefficients(colors, exponents)
    f, g = [packed_bits(value) for value in forward], [packed_bits(value) for value in inverse]
    count = 0
    for degree in range(8):
        for i in range(len(colors)):
            for j in range(len(colors)):
                constant = omega = 0
                for k in range(degree + 1):
                    a, b = f[k][0][i], f[k][1][i]
                    c, d = g[degree - k][0][j], g[degree - k][1][j]
                    bd = (b & d).bit_count()
                    constant ^= (a & c).bit_count() ^ bd
                    omega ^= (a & d).bit_count() ^ (b & c).bit_count() ^ bd
                if (constant & 1) != int(degree == 0 and i == j) or omega & 1:
                    raise ValueError("the actual full-column bitset Gram has a nonzero coefficient")
                count += 1
    return count


def verify_actual(order, phase, assignments):
    problem = variable_problem(order, phase)
    groups, fillers, labels = problem["template_groups"], problem["inactive_fillers"], problem["copy_second_bits"]
    s, width, t = len(problem["C"]), 10, len(problem["border_C"])
    type_counts = np.bincount(labels, minlength=2).tolist()
    row_third = np.zeros(s, dtype=np.int64)
    for group, label in zip(groups, labels, strict=True):
        for pair in range(3):
            row_third[group[2 * pair:2 * pair + 2]] = assignments[("boundary_pair", int(label), pair)]
    for first, second in problem["actual_row_pairs"]:
        if row_third[first] != row_third[second]:
            raise ValueError("the new third bit disagrees within a physical row pair")
    u = (problem["U_exponents"] + 12 * row_third[:, None] * problem["A"]) % phase
    expected = np.zeros((t, t, root_vectors(phase).shape[1]), dtype=np.int64)
    expected[np.arange(t), np.arange(t), 0] = s
    if not np.array_equal(exact_column_gram(u, phase), expected):
        raise ValueError("the jointly changed actual boundary lost its exact partial-frame Gram")
    chosen = [int(i) for value in (0, 1) for i in np.flatnonzero(labels == value)[:2]]
    rows = [i for copy in chosen for i in groups[copy]] + fillers[:2]
    canonical = np.empty(s, dtype=np.int64)
    for index, group in enumerate(groups):
        canonical[group] = np.arange(10 * index, 10 * (index + 1))
    canonical[fillers] = np.arange(len(groups) * width, s)
    xb = np.array([[assignments[("X", *copy_descriptor(int(canonical[i]), int(canonical[j]), type_counts))]
                    for j in range(s)] for i in rows], dtype=np.int64)
    yb = np.array([[assignments[("Y", *copy_descriptor(int(canonical[i]), int(canonical[j]), type_counts))]
                    for j in range(s)] for i in rows], dtype=np.int64)
    c, w = problem["C"][rows], problem["W"][rows]
    uc, ub = phase_coordinates(u, phase)
    dc, db = phase_coordinates(problem["D_exponents"], phase)
    colors = np.block([[np.zeros_like(c), c, uc[rows]], [c, np.zeros_like(c), uc[rows]], [uc.T, uc.T, dc]])
    shift = problem["core_base_shift"]
    x = shift + 4 * xb
    y = ((shift + w * (1 + 2 * problem["row_copy_bits"][rows, None])) % 4) + 4 * yb
    exponents = np.block([[x, y, ub[rows] % 8], [y, x, ub[rows] % 8], [ub.T % 8, ub.T % 8, db % 8]])
    count = full_column_bitset_gram(colors, exponents)
    return {"source_order": order, "phase": phase, "actual_variable_boundary_Gram_reverified": True,
            "boundary_third_bits_by_type_and_local_physical_pair":
            [[assignments[("boundary_pair", label, pair)] for pair in range(3)] for label in (0, 1)],
            "U_exponents": u.tolist(), "all_actual_source_columns_used": colors.shape[1],
            "representative_source_rows": rows + [s + i for i in rows] + list(range(2 * s, order)),
            "eight_coefficient_ordered_row_pair_checks": count,
            "all_source_row_pair_orbits_covered_by_copy_type_and_filler_symmetry": True,
            "source_row_pairs_covered": order * (order + 1) // 2,
            "all_eight_actual_Gram_coefficients_verified": True,
            "phase_rows_sha256": hashlib.sha256(colors.astype(np.uint8).tobytes() + exponents.astype(np.uint8).tobytes()).hexdigest()}


def certificate():
    routes = []
    for order, phase in ((1499, 96), (2999, 48), (2999, 96)):
        model, assignments, _ = joint_model(order, phase)
        route = {"source_order": order, "phase": phase, "joint_linear_model": model,
                 "characteristic_zero_source_constructed": False}
        if assignments is not None:
            route["actual_large_verification"] = verify_actual(order, phase, assignments)
        routes.append(route)
    return {"schema": "variable-cubic-third-ramification-v1", "coefficient_ring": "F4[pi]/(pi^8)",
            "previous_source": "experiments/scalar_lift/variable_cubic_frame_bits.py",
            "routes": routes, "native_probe_run": False,
            "no_general_characteristic_zero_feasibility_inferred": True,
            "new_real_Hadamard_order_constructed": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--output", type=Path)
    action.add_argument("--check", type=Path)
    args = parser.parse_args()
    report = certificate()
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    elif json.loads(args.check.read_text()) != report:
        raise ValueError("the variable-boundary third-bit certificate failed replay")
    print(json.dumps({"routes": [[r["source_order"], r["phase"], r["joint_linear_model"]["status"],
                                  r["joint_linear_model"].get("selected_minor_excludes_every_width3_interval_gauge_mod8")]
                                 for r in report["routes"]], "native_probe_run": False}))


if __name__ == "__main__":
    main()

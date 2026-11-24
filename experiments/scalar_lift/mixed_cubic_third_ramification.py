"""Exact third phase-bit lifting modulo pi eighth at the actual sources."""

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
from experiments.scalar_lift.mixed_cubic_second_ramification import (
    descriptor,
    encode_blocks,
    phase_variables,
    product,
    representative_problem,
    solve_binary_columns,
)


def root_coefficients(colors, exponents, precision=8):
    colors, exponents = np.asarray(colors, dtype=np.int64), np.asarray(exponents, dtype=np.int64) % precision
    roots = np.array([1, 2, 3], dtype=np.int64)[colors]
    inverse_roots, inverse_exponents = CONJUGATE[roots], (-exponents) % precision
    forward = [MULTIPLY[roots, ((exponents & k) == k).astype(np.int64)] for k in range(precision)]
    inverse = [MULTIPLY[inverse_roots, ((inverse_exponents & k) == k).astype(np.int64)] for k in range(precision)]
    return forward, inverse


def gram_coefficients(colors, exponents, precision=8):
    forward, inverse = root_coefficients(colors, exponents, precision)
    gram = []
    for k in range(precision):
        result = np.zeros((len(colors), len(colors)), dtype=np.int64)
        for j in range(k + 1):
            result ^= product(forward[j], inverse[k - j].T)
        gram.append(result)
    return gram


def representative_roots(order, phase):
    small = template()
    boundary = boundary_data(order, phase)
    local = representative_problem(small["C"], small["B"], small["a"],
                                   boundary["border_C"], boundary["border_next_bits"])
    s, t = local["s"], local["t"]
    a = local["active"][:, None] * boundary["boundary_vector"][None, :]
    zero = np.zeros((s, s), dtype=np.int64)
    colors = np.block([[zero, local["C"], a], [local["C"], zero, a],
                       [a.T, a.T, boundary["border_C"]]])
    _, border_dyadic = phase_coordinates(boundary["D_exponents"], phase)
    core_shift = 2 * int(order == 2999)
    core_x = np.full((s, s), core_shift, dtype=np.int64)
    core_y = local["W"] + core_shift
    exponents = np.block([[core_x, core_y, a], [core_y, core_x, a],
                          [a.T, a.T, border_dyadic % 8]])
    if t not in (3, 7):
        raise ValueError("the represented border is not one of the actual targets")
    return local, boundary, colors, exponents


def solve_third_bits(order, phase):
    local, boundary, colors, exponents = representative_roots(order, phase)
    s = local["s"]
    forward, inverse = root_coefficients(colors, exponents)
    grams = gram_coefficients(colors, exponents)
    if not np.array_equal(grams[0], np.eye(len(colors), dtype=np.int64)) or any(np.any(g) for g in grams[1:4]):
        raise ValueError("the saved lower phase bits no longer give the pi^4 root residue")
    names, groups = phase_variables(width=local["width"])
    columns, variable_names = [], []
    for core in ("X", "Y"):
        for name, positions in zip(names, groups, strict=True):
            positions = ([(i, j) for i, j in positions] + [(s + i, s + j) for i, j in positions]
                         if core == "X" else
                         [(i, s + j) for i, j in positions] + [(s + i, j) for i, j in positions])
            blocks = []
            for remainder in range(4):
                block = np.zeros((len(colors), len(colors)), dtype=np.int64)
                for first in range(remainder + 1):
                    other = remainder - first
                    for i, j in positions:
                        # Adding4 to an exponent changes forward and
                        # inverse root polynomials by pi^4 times their
                        # lower-degree polynomials, modulo pi^8.
                        value, inverse_value = forward[first][i, j], inverse[first][i, j]
                        if value:
                            block[i] ^= MULTIPLY[value, inverse[other][:, j]]
                        if inverse_value:
                            block[:, i] ^= MULTIPLY[forward[other][:, j], inverse_value]
                blocks.append(block)
            columns.append(encode_blocks(*blocks))
            variable_names.append((core, *name))
    target = encode_blocks(*grams[4:])
    solution, rank = solve_binary_columns(columns, target)
    digest = hashlib.sha256()
    byte_width = (8 * len(colors) ** 2 + 7) // 8
    for vector in [*columns, target]:
        digest.update(vector.to_bytes(byte_width, "little"))
    report = {"source_order": order, "phase": phase, "number_third_bit_variables": len(columns),
              "exact_binary_system": rank, "binary_system_sha256": digest.hexdigest(),
              "border_dyadic_exponents_mod8": exponents[-local["t"]:, -local["t"]:].tolist(),
              "status": "SOLVED" if solution is not None else "INCONSISTENT_IN_THIS_EQUIVARIANT_SLICE"}
    if solution is None:
        report["general_source_or_nonequivariant_third_bits_excluded"] = False
        return report, None
    assignments = {name: (solution >> index) & 1 for index, name in enumerate(variable_names)}
    for core in ("X", "Y"):
        for name, positions in zip(names, groups, strict=True):
            for i, j in positions:
                if core == "X":
                    exponents[i, j] += 4 * assignments[(core, *name)]
                    exponents[s + i, s + j] += 4 * assignments[(core, *name)]
                else:
                    exponents[i, s + j] += 4 * assignments[(core, *name)]
                    exponents[s + i, j] += 4 * assignments[(core, *name)]
    new_grams = gram_coefficients(colors, exponents)
    if not np.array_equal(new_grams[0], np.eye(len(colors), dtype=np.int64)) or any(np.any(g) for g in new_grams[1:]):
        raise ValueError("the third-bit solution failed the independent root-polynomial Gram")
    report.update({"nonzero_variable_indices": [index for index in range(len(columns)) if (solution >> index) & 1],
                   "all_eight_representative_Gram_coefficients_verified": True})
    return report, assignments


def verify_actual_large_lift(problem, assignments):
    s, groups, fillers = len(problem["C"]), problem["template_groups"], problem["inactive_fillers"]
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
    c, w, a = problem["C"][local_rows], problem["W"][local_rows], u_colors[local_rows]
    colors = np.block([[np.zeros_like(c), c, a], [c, np.zeros_like(c), a],
                       [u_colors.T, u_colors.T, d_colors]])
    shift = 2 * int(problem["source_order"] == 2999)
    core_x, core_y = shift + 4 * xbits, w + shift + 4 * ybits
    exponents = np.block([[core_x, core_y, u_dyadic[local_rows] % 8],
                          [core_y, core_x, u_dyadic[local_rows] % 8],
                          [u_dyadic.T % 8, u_dyadic.T % 8, d_dyadic % 8]])
    grams = gram_coefficients(colors, exponents)
    if not np.array_equal(grams[0], np.eye(len(colors), dtype=np.int64)) or any(np.any(g) for g in grams[1:]):
        raise ValueError("the actual high-order root candidate failed a pi^8 Gram coefficient")
    if not np.array_equal(u_colors, problem["A"]) or np.any(u_dyadic % 8 != problem["A"]):
        raise ValueError("the actual boundary phases are not equivariant at this precision")
    small = template()
    for group in groups:
        where = np.ix_(group, group)
        if (not np.array_equal(problem["C"][where], small["C"])
                or not np.array_equal(problem["W"][where], small["W"])
                or not np.array_equal(problem["active"][group], small["a"])):
            raise ValueError("the target template symmetry failed")
    if (int(np.count_nonzero(problem["W"])) != len(groups) * int(np.count_nonzero(small["W"]))
            or int(np.count_nonzero(1 ^ problem["C"])) != len(groups) * int(np.count_nonzero(1 ^ small["C"])) + len(fillers)):
        raise ValueError("an undeclared off-template term breaks the orbit coverage")
    size = 2 * s + t
    return {"source_order": size, "phase": problem["phase"],
            "all_actual_columns_used_in_Gram": colors.shape[1],
            "template_copies": len(groups), "inactive_fillers": len(fillers),
            "representative_source_rows": local_rows + [s + i for i in local_rows] + list(range(2 * s, size)),
            "number_ordered_representative_row_pairs_checked": len(colors) ** 2,
            "all_source_row_pairs_covered_by_exact_permutation_symmetry": size * (size + 1) // 2,
            "all_eight_actual_Gram_coefficients_verified": True,
            "actual_boundary_and_border_phase_exponents_retained_mod8": True,
            "representative_colors_and_dyadic_exponents_sha256":
            hashlib.sha256(colors.astype(np.uint8).tobytes() + exponents.astype(np.uint8).tobytes()).hexdigest()}


def certificate():
    routes = []
    for order, phase in ((1499, 96), (2999, 48), (2999, 96)):
        model, assignments = solve_third_bits(order, phase)
        route = {"source_order": order, "phase": phase, "exact_linear_model": model,
                 "conditional_real_order": order * (phase // 3),
                 "characteristic_zero_source_constructed": False}
        if assignments is not None:
            route["actual_large_verification"] = verify_actual_large_lift(large_problem(order, phase), assignments)
            route["core_third_bit_values"] = sorted(set(assignments.values()))
            route["unitary_root_residue_modulo_pi_eighth"] = True
        else:
            route["unitary_root_residue_modulo_pi_eighth"] = False
        routes.append(route)
    return {"schema": "mixed-cubic-third-ramification-v1",
            "coefficient_ring": "F4[pi]/(pi^8), pi=1-zeta_d, d in{16,32}",
            "previous_certificate": "certificates/scalar-lift/mixed-cubic-second-ramification.json",
            "new_bit_enters_at_degree4_and_its_square_vanishes_modulo_pi8": True,
            "involution_implemented_by_actual_inverse_root_exponents_mod8": True,
            "routes": routes, "native_probe_run": False,
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
            raise ValueError("the third-ramification certificate failed replay")
    else:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"routes": [[item["source_order"], item["phase"], item["exact_linear_model"]["status"]]
                                 for item in report["routes"]],
                      "precision": "pi^8", "characteristic_zero_source_constructed": False}))


if __name__ == "__main__":
    main()

"""Direct pi-fourth root lifts of the nondegenerate twelve-row templates.

No projector identity is imposed. The finite linear chart has two
labelled copy types and fixed genuine boundary phase classes.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from experiments.scalar_lift.cubic_core_support_lift import boundary_data, phase_coordinates
from experiments.scalar_lift.labelled_copy_ramification import (
    complete_binary_solution_space,
    four_difference_witness,
    variable_groups,
)
from experiments.scalar_lift.mixed_cubic_ramification import CONJUGATE, MULTIPLY
from experiments.scalar_lift.mixed_cubic_second_ramification import (
    derivative,
    encode_blocks,
    polynomial_gram_from_roots,
    product,
    sparse_derivative,
)
from experiments.scalar_lift.nondegenerate_cubic_projectors import retained_templates
from experiments.scalar_lift.variable_cubic_third_ramification import short_interval_gauges


def template_problem(candidate, order, phase):
    boundary = boundary_data(order, phase)
    counts = (3, 3) if order == 1499 else (2, 2)
    width, fillers = len(candidate["C"]), 2
    copies = sum(counts)
    s, t = copies * width + fillers, len(boundary["border_C"])
    correction = np.zeros((s, s), dtype=np.int64)
    w = np.zeros_like(correction)
    for copy in range(copies):
        part = slice(copy * width, (copy + 1) * width)
        correction[part, part], w[part, part] = candidate["T"], candidate["W"]
    correction[-fillers:, -fillers:] = np.eye(fillers, dtype=np.int64)
    c = 1 ^ correction
    active = np.concatenate([candidate["a"]] * copies + [np.zeros(fillers, dtype=np.int64)])
    a = active[:, None] * boundary["boundary_vector"][None, :]
    labels = np.array([0] * (width * counts[0]) + [1] * (width * counts[1]) + [0] * fillers)
    ones, zeros = np.ones((s, s), dtype=np.int64), np.zeros((s, s), dtype=np.int64)
    y0, u0, d0 = 1 ^ (3 * c), 1 ^ (3 * a), 1 ^ (3 * boundary["border_C"])
    y1, u1 = MULTIPLY[y0, w], 2 * a
    h0 = np.block([[ones, y0, u0], [y0, ones, u0], [u0.T, u0.T, d0]])
    h1 = np.block([[zeros, y1, u1], [y1, zeros, u1], [u1.T, u1.T, np.zeros((t, t), dtype=np.int64)]])
    u2 = 2 * a * labels[:, None]
    fixed_h2 = np.zeros_like(h0)
    fixed_h3 = np.zeros_like(h0)
    for target in (fixed_h2, fixed_h3):
        target[:s, 2 * s:] = target[s:2 * s, 2 * s:] = u2
        target[2 * s:, :s] = target[2 * s:, s:2 * s] = u2.T
    fixed_h2[-t:, -t:] = MULTIPLY[d0, boundary["border_next_bits"]]
    if not np.array_equal(product(h0, CONJUGATE[h0].T), np.eye(len(h0), dtype=np.int64)):
        raise ValueError("the twelve-row template failed the representative zeroth Gram")
    if np.any(derivative(h0, h1)):
        raise ValueError("the twelve-row template failed the full first Gram derivative")
    return {"s": s, "t": t, "width": width, "counts": counts,
            "H0": h0, "H1": h1, "fixed_H2": fixed_h2, "fixed_H3": fixed_h3,
            "C": c, "W": w, "A": a, "labels": labels, "boundary": boundary}


def solve_template(candidate, order, phase):
    data = template_problem(candidate, order, phase)
    h0, h1, s, t = data["H0"], data["H1"], data["s"], data["t"]
    fixed = product(h1, CONJUGATE[h1].T) ^ product(h0, CONJUGATE[h1].T)
    target2 = fixed ^ derivative(h0, data["fixed_H2"])
    target3 = fixed ^ derivative(h1, data["fixed_H2"]) ^ derivative(h0, data["fixed_H3"])
    names, groups = variable_groups(data["counts"], width=data["width"])
    columns = []
    for core in ("X", "Y"):
        for positions in groups:
            positions = ([(i, j) for i, j in positions] + [(s + i, s + j) for i, j in positions]
                         if core == "X" else
                         [(i, s + j) for i, j in positions] + [(s + i, j) for i, j in positions])
            second = [(i, j, int(h0[i, j])) for i, j in positions]
            third = [(i, j, int(h1[i, j])) for i, j in positions]
            columns.append(encode_blocks(sparse_derivative(h0, second),
                                         sparse_derivative(h0, third) ^ sparse_derivative(h1, second)))
    target = encode_blocks(target2, target3)
    solution, kernel, rank = complete_binary_solution_space(columns, target)
    digest = hashlib.sha256()
    byte_width = (4 * len(h0) ** 2 + 7) // 8
    for vector in [*columns, target]:
        digest.update(vector.to_bytes(byte_width, "little"))
    report = {"template_sha256": candidate["template_sha256"], "source_order": order, "phase": phase,
              "template_width": data["width"], "template_active_weight": int(candidate["a"].sum()),
              "actual_copy_type_counts": [31, 31] if order == 1499 else [62, 62],
              "actual_inactive_fillers": (order - t) // 2 - (62 if order == 1499 else 124) * data["width"],
              "representative_copy_type_counts": list(data["counts"]),
              "representative_source_order": len(h0), "core_second_bit_variables": len(columns),
              "exact_binary_system": rank, "binary_system_sha256": digest.hexdigest(),
              "zeroth_and_first_Gram_verified_directly": True,
              "projector_identity_not_required": True,
              "status": "SOLVED" if solution is not None else "INCONSISTENT_IN_THIS_LABELLED_SECOND_BIT_SLICE",
              "unrestricted_template_or_source_nonexistence_claimed": False}
    if solution is None:
        return report, None
    values = [(solution >> i) & 1 for i in range(len(columns))]
    x2, y2 = np.zeros((s, s), dtype=np.int64), np.zeros((s, s), dtype=np.int64)
    for index, positions in enumerate(groups):
        for i, j in positions:
            x2[i, j], y2[i, j] = values[index], values[len(names) + index]
    u = data["A"] * (1 + 2 * data["labels"][:, None])
    _, db = phase_coordinates(data["boundary"]["D_exponents"], phase)
    exponents = np.block([[2 * x2, data["W"] + 2 * y2, u],
                          [data["W"] + 2 * y2, 2 * x2, u], [u.T, u.T, db % 4]])
    colors = (h0 == 2).astype(np.int64)
    grams = polynomial_gram_from_roots(colors, exponents)
    if not np.array_equal(grams[0], np.eye(len(h0), dtype=np.int64)) or any(np.any(g) for g in grams[1:]):
        raise ValueError("the twelve-row second lift failed independent inverse-root polynomial Gram")
    border = int(np.flatnonzero(data["boundary"]["boundary_vector"])[0])
    minor_indices = [0, data["width"] * data["counts"][0], 2 * s + border]
    minor = exponents[np.ix_(minor_indices, minor_indices)]
    minor_audit = short_interval_gauges(minor, modulus=4, width=1)
    difference = four_difference_witness(exponents)
    anti = minor_audit["admissible_gauges"] == 0 or difference is not None
    report.update({"all_four_representative_Gram_coefficients_verified_by_inverse_roots": True,
                   "second_bit_variable_descriptors_per_core": [list(name) for name in names],
                   "particular_second_bit_assignment": values,
                   "complete_affine_kernel_basis_nonzero_indices":
                   [[i for i in range(len(columns)) if (word >> i) & 1] for word in kernel],
                   "minor_indices": minor_indices, "dyadic_minor_mod4": minor.tolist(),
                   "minor_adjacent_channel_audit": minor_audit,
                   "four_difference_witness": difference,
                   "adjacent_channel_collapse_ruled_out_for_this_candidate": anti,
                   "actual_large_matrix_claimed_before_frame_reconstruction": False})
    assignments = {(core, *name): values[index + (len(names) if core == "Y" else 0)]
                   for core in ("X", "Y") for index, name in enumerate(names)}
    return report, assignments


def scan_report(*, progress=False):
    candidates = retained_templates()
    results, selected = [], {}
    for index, candidate in enumerate(candidates):
        for order, phase in ((1499, 24), (2999, 48)):
            if str(order) in selected:
                continue
            report, assignment = solve_template(candidate, order, phase)
            report["retained_candidate_index"] = index
            results.append(report)
            if assignment is not None and report["adjacent_channel_collapse_ruled_out_for_this_candidate"]:
                selected[str(order)] = index
            if progress:
                print(json.dumps({"candidate": index, "order": order, "status": report["status"],
                                  "rank": report["exact_binary_system"]["rank"],
                                  "anti_collapse": report.get("adjacent_channel_collapse_ruled_out_for_this_candidate")}), flush=True)
        if len(selected) == 2:
            break
    return {"schema": "nondegenerate-cubic-second-lift-v1",
            "input_templates": "certificates/scalar-lift/nondegenerate-cubic-first-lifts.json",
            "retained_input_count": len(candidates), "maximum_templates_examined": 15,
            "coefficient_ring": "F4[pi]/(pi^4)", "results": results,
            "selected_candidate_indices_by_source_order": selected,
            "native_probe_run": False, "new_characteristic_zero_source_claimed": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--output", type=Path)
    action.add_argument("--check", type=Path)
    parser.add_argument("--progress", action="store_true")
    args = parser.parse_args()
    report = scan_report(progress=args.progress)
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    elif json.loads(args.check.read_text()) != report:
        raise ValueError("the nondegenerate second-lift certificate failed replay")
    print(json.dumps({"exact_cases_checked": len(report["results"]),
                      "selected_candidate_indices_by_source_order": report["selected_candidate_indices_by_source_order"],
                      "native_probe_run": False}))


if __name__ == "__main__":
    main()

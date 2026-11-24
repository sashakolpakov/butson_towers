"""Joint second-layer root lifts with exact physical-pair boundary actions.

The core is the finite two-labelled-copy chart. Boundary parameters are
chosen simultaneously: odd exponents on active row pairs and even
exponents on inactive row pairs. No projector identity is required.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from experiments.scalar_lift.cubic_core_support_lift import phase_coordinates, template
from experiments.scalar_lift.labelled_copy_ramification import (
    complete_binary_solution_space,
    four_difference_witness,
    variable_groups,
)
from experiments.scalar_lift.mixed_cubic_ramification import CONJUGATE
from experiments.scalar_lift.mixed_cubic_second_ramification import (
    derivative,
    encode_blocks,
    polynomial_gram_from_roots,
    product,
    sparse_derivative,
)
from experiments.scalar_lift.nondegenerate_cubic_projectors import retained_templates
from experiments.scalar_lift.nondegenerate_cubic_second_lift import template_problem


def build_system(candidate, order, phase):
    data = template_problem(candidate, order, phase)
    h0, h1, s, t, width = (data[key] for key in ("H0", "H1", "s", "t", "width"))
    if width % 2 or any(candidate["a"][j] != candidate["a"][j + 1] for j in range(0, width, 2)):
        raise ValueError("the chosen adjacent physical pairs must have constant first phase bit")
    # The fixed part contains only the actual border; every pair boundary
    # second bit is now a joint variable, including the old fixed h slice.
    h2 = np.zeros_like(h0)
    h2[-t:, -t:] = data["fixed_H2"][-t:, -t:]
    fixed = product(h1, CONJUGATE[h1].T) ^ product(h0, CONJUGATE[h1].T)
    target = encode_blocks(fixed ^ derivative(h0, h2), fixed ^ derivative(h1, h2))
    names, groups = variable_groups(data["counts"], width=width)
    columns, descriptors = [], []
    for core in ("X", "Y"):
        for name, positions in zip(names, groups, strict=True):
            positions = ([(i, j) for i, j in positions] + [(s + i, s + j) for i, j in positions]
                         if core == "X" else
                         [(i, s + j) for i, j in positions] + [(s + i, j) for i, j in positions])
            second = [(i, j, int(h0[i, j])) for i, j in positions]
            third = [(i, j, int(h1[i, j])) for i, j in positions]
            columns.append(encode_blocks(sparse_derivative(h0, second),
                                         sparse_derivative(h0, third) ^ sparse_derivative(h1, second)))
            descriptors.append((core, *name))
    core_count = len(columns)
    for label in (0, 1):
        for pair in range(width // 2):
            active = bool(candidate["a"][2 * pair])
            second, third = [], []
            for copy in range(sum(data["counts"])):
                if int(copy >= data["counts"][0]) != label:
                    continue
                for local in (2 * pair, 2 * pair + 1):
                    row = copy * width + local
                    for border in np.flatnonzero(data["boundary"]["boundary_vector"]):
                        for i, j in ((row, 2 * s + border), (s + row, 2 * s + border),
                                     (2 * s + border, row), (2 * s + border, s + row)):
                            second.append((i, j, 2 if active else 1))
                            if active:
                                third.append((i, j, 2))
            columns.append(encode_blocks(sparse_derivative(h0, second),
                                         sparse_derivative(h0, third) ^ sparse_derivative(h1, second)))
            descriptors.append(("boundary_pair", label, pair, int(active)))
    data.update({"candidate": candidate, "phase": phase, "core_names": names, "core_groups": groups,
                 "core_variable_count": core_count, "columns": columns,
                 "variable_descriptors": descriptors, "target": target})
    return data


def root_assignment(data, assignment):
    s, width = data["s"], data["width"]
    x2, y2 = np.zeros((s, s), dtype=np.int64), np.zeros((s, s), dtype=np.int64)
    count = len(data["core_names"])
    for index, positions in enumerate(data["core_groups"]):
        for i, j in positions:
            x2[i, j] = (assignment >> index) & 1
            y2[i, j] = (assignment >> (index + count)) & 1
    u = data["A"].copy()
    for index, descriptor in enumerate(data["variable_descriptors"]):
        if descriptor[0] != "boundary_pair" or not ((assignment >> index) & 1):
            continue
        _, label, pair, _ = descriptor
        for copy in range(sum(data["counts"])):
            if int(copy >= data["counts"][0]) == label:
                rows = [copy * width + 2 * pair, copy * width + 2 * pair + 1]
                u[rows] += 2 * data["boundary"]["boundary_vector"]
    _, db = phase_coordinates(data["boundary"]["D_exponents"], data["phase"])
    db = db.copy()
    for index, descriptor in enumerate(data["variable_descriptors"]):
        if descriptor[0] == "border_second_bit" and ((assignment >> index) & 1):
            db += 2 * data["border_second_bit_masks"][descriptor[1]]
    exponents = np.block([[2 * x2, data["W"] + 2 * y2, u],
                          [data["W"] + 2 * y2, 2 * x2, u], [u.T, u.T, db % 4]])
    return (data["H0"] == 2).astype(np.int64), exponents


def check_root_assignment(data, assignment):
    colors, exponents = root_assignment(data, assignment)
    grams = polynomial_gram_from_roots(colors, exponents)
    if not np.array_equal(grams[0], np.eye(len(colors), dtype=np.int64)) or any(np.any(g) for g in grams[1:]):
        raise ValueError("the joint solution failed the independent full inverse-root Gram")
    return {"all_four_representative_Gram_coefficients_verified_by_inverse_roots": True,
            "four_difference_witness": four_difference_witness(exponents)}


def solve_system(data, *, inactive_pairs):
    descriptors = data["variable_descriptors"]
    selected = [i for i, descriptor in enumerate(descriptors)
                if descriptor[0] != "boundary_pair" or inactive_pairs or descriptor[3]]
    columns = [data["columns"][i] for i in selected]
    solution, kernel, ranks = complete_binary_solution_space(columns, data["target"])
    byte_width = (4 * len(data["H0"]) ** 2 + 7) // 8
    digest = hashlib.sha256()
    for vector in [*columns, data["target"]]:
        digest.update(vector.to_bytes(byte_width, "little"))
    report = {"boundary_variation": "all_physical_pairs" if inactive_pairs else "active_physical_pairs_only",
              "variables": len(columns), "core_variables": data["core_variable_count"],
              "boundary_variables": len(columns) - data["core_variable_count"],
              "exact_binary_system": ranks, "binary_system_sha256": digest.hexdigest(),
              "status": "SOLVED" if solution is not None else "INCONSISTENT_IN_THIS_JOINT_LABELLED_SLICE"}
    if solution is None:
        return report
    def expand(word):
        return sum(1 << old for new, old in enumerate(selected) if (word >> new) & 1)
    particular = expand(solution)
    report.update(check_root_assignment(data, particular))
    report["particular_nonzero_variable_descriptors"] = [list(descriptors[i]) for i in range(len(descriptors))
                                                         if (particular >> i) & 1]
    report["complete_affine_kernel_basis_nonzero_variable_descriptors"] = [
        [list(descriptors[i]) for i in range(len(descriptors)) if (expand(word) >> i) & 1] for word in kernel]
    # Independently replay every basis neighbor, beyond verifying the binary
    # matrix kernel. Arbitrary combinations follow from the exact affine map.
    for word in kernel:
        check_root_assignment(data, particular ^ expand(word))
    report["all_kernel_basis_neighbors_verified_by_inverse_root_Gram"] = True
    return report


def solve_template(candidate, order, phase):
    data = build_system(candidate, order, phase)
    data["phase"] = phase
    return {"template_sha256": candidate["template_sha256"], "source_order": order, "phase": phase,
            "template_width": data["width"], "active_pair_count": int(candidate["a"].sum()) // 2,
            "representative_copy_type_counts": list(data["counts"]),
            "representative_source_order": len(data["H0"]),
            "actual_copy_type_counts": [31, 31] if order == 1499 else [62, 62],
            "actual_inactive_fillers": (order - data["t"]) // 2 - (62 if order == 1499 else 124) * data["width"],
            "filler_pair_phases_fixed": True,
            "zeroth_and_first_Gram_verified_directly": True,
            "projector_identity_not_required": True,
            "slices": [solve_system(data, inactive_pairs=False), solve_system(data, inactive_pairs=True)],
            "unrestricted_template_or_source_nonexistence_claimed": False,
            "actual_large_matrix_claimed_before_frame_reconstruction": False}


def certificate(*, progress=False):
    controls, results = [], []
    positive = {**template(), "template_sha256": "known-ten-row-positive-control"}
    for order, phase in ((1499, 24), (2999, 48)):
        controls.append(solve_template(positive, order, phase))
    for index, candidate in enumerate(retained_templates()):
        for order, phase in ((1499, 24), (2999, 48)):
            report = solve_template(candidate, order, phase)
            report["retained_candidate_index"] = index
            results.append(report)
            if progress:
                print(json.dumps({"candidate": index, "order": order,
                                  "slices": [{"status": part["status"], "ranks": part["exact_binary_system"]}
                                             for part in report["slices"]]}), flush=True)
    return {"schema": "nondegenerate-cubic-joint-boundary-lift-v1",
            "input_templates": "certificates/scalar-lift/nondegenerate-cubic-first-lifts.json",
            "coefficient_ring": "F4[pi]/(pi^4)", "controls": controls, "results": results,
            "native_probe_run": False, "new_characteristic_zero_source_claimed": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--output", type=Path)
    action.add_argument("--check", type=Path)
    parser.add_argument("--progress", action="store_true")
    args = parser.parse_args()
    report = certificate(progress=args.progress)
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    elif json.loads(args.check.read_text()) != report:
        raise ValueError("the joint-boundary certificate failed replay")
    print(json.dumps({"cases": len(report["results"]), "controls": len(report["controls"]),
                      "solved_slices": sum(part["status"] == "SOLVED" for row in report["results"] for part in row["slices"])}))


if __name__ == "__main__":
    main()

"""Joint twelve-row pi-fourth lifts allowing the order-seven border to move."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from experiments.scalar_lift.cubic_core_support_lift import template
from experiments.scalar_lift.mixed_cubic_second_ramification import encode_blocks, sparse_derivative
from experiments.scalar_lift.nondegenerate_cubic_joint_boundary_lift import (
    build_system,
    solve_system,
)
from experiments.scalar_lift.nondegenerate_cubic_projectors import retained_templates


def actual_affine_masks():
    """Three individually exact Petrescu branches, in the stored D7 gauge.

    The sign matters for the exact phase family; only its parity enters
    the pi-fourth second-bit problem. Simultaneous parameters are not
    asserted to give an exact characteristic-zero family.
    """
    supports = (((0, 1), (5, 6), (3, 4), (5, 6)),
                ((0, 2), (4, 6), (2, 3), (0, 5)),
                ((1, 2), (4, 5), (2, 4), (0, 6)))
    masks = []
    for first_rows, second_rows, first_columns, second_columns in supports:
        mask = np.zeros((7, 7), dtype=np.int64)
        mask[np.ix_(first_rows, first_columns)] = 1
        mask[np.ix_(second_rows, second_columns)] = -1
        masks.append(mask)
    return masks


def with_border_directions(data, masks):
    extended = dict(data)
    extended["columns"] = list(data["columns"])
    extended["variable_descriptors"] = list(data["variable_descriptors"])
    extended["border_second_bit_masks"] = [mask % 2 for mask in masks]
    s, t = data["s"], data["t"]
    for index, mask in enumerate(masks):
        if mask.shape != (t, t):
            raise ValueError("a border phase mask has the wrong dimensions")
        entries = [(2 * s + i, 2 * s + j, int(data["H0"][2 * s + i, 2 * s + j]))
                   for i, j in zip(*np.nonzero(mask % 2), strict=True)]
        extended["columns"].append(encode_blocks(sparse_derivative(data["H0"], entries),
                                                  sparse_derivative(data["H1"], entries)))
        extended["variable_descriptors"].append(("border_second_bit", index))
    return extended


def audit_affine_masks():
    # Independent Laurent cancellation: for every row pair, collect
    # phase differences by their integer parameter exponent. Each
    # coefficient must vanish as an exact cyclotomic integer.
    from experiments.scalar_lift.cubic_core_support_lift import boundary_data
    from experiments.scalar_lift.fifth_source_assembly import root_vectors

    border = boundary_data(2999, 48)["D_exponents"]
    roots = np.asarray(root_vectors(48), dtype=np.int64)
    reports = []
    for mask in actual_affine_masks():
        coefficient_count = 0
        for i in range(7):
            for j in range(i):
                powers = mask[i] - mask[j]
                for power in np.unique(powers):
                    terms = roots[(border[i, powers == power] - border[j, powers == power]) % 48]
                    if np.any(terms.sum(axis=0)):
                        raise ValueError("a proposed border affine branch fails exact Laurent cancellation")
                    coefficient_count += 1
        reports.append({"signed_mask": mask.tolist(), "exact_Laurent_coefficients_checked": coefficient_count,
                        "all_unit_parameters_preserve_the_exact_border_Gram": True})
    return reports


def solve_candidate(candidate):
    data = build_system(candidate, 2999, 48)
    basis = [np.eye(49, dtype=np.int64)[index].reshape(7, 7) for index in range(49)]
    broad = solve_system(with_border_directions(data, basis), inactive_pairs=True)
    broad["border_directions"] = 49
    broad["border_model"] = "arbitrary_root_second_bits_in_the_fixed_first_layer"
    branches = []
    for index, mask in enumerate(actual_affine_masks()):
        branch = solve_system(with_border_directions(data, [mask]), inactive_pairs=True)
        branch.update({"border_directions": 1, "individually_exact_Petrescu_branch": index})
        branches.append(branch)
    return {"template_sha256": candidate["template_sha256"], "source_order": 2999, "phase": 48,
            "template_width": data["width"], "all_twelve_labelled_physical_pair_bits_free": True,
            "filler_pair_phases_fixed": True, "broad_border_second_bit_chart": broad,
            "individually_exact_border_branches": branches,
            "no_unrestricted_source_or_template_nonexistence_claimed": True}


def certificate(*, progress=False):
    affine = audit_affine_masks()
    control = solve_candidate({**template(), "template_sha256": "known-ten-row-positive-control"})
    results = []
    for index, candidate in enumerate(retained_templates()):
        result = solve_candidate(candidate)
        result["retained_candidate_index"] = index
        results.append(result)
        if progress:
            print(json.dumps({"candidate": index, "broad_status": result["broad_border_second_bit_chart"]["status"],
                              "ranks": result["broad_border_second_bit_chart"]["exact_binary_system"],
                              "exact_branches": [branch["status"] for branch in result["individually_exact_border_branches"]]}), flush=True)
    return {"schema": "nondegenerate-cubic-moving-border-lift-v1", "coefficient_ring": "F4[pi]/(pi^4)",
            "exact_border_affine_branch_audit": affine, "known_positive_control": control, "results": results,
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
        raise ValueError("the moving-border lift certificate failed replay")
    print(json.dumps({"candidates": len(report["results"]),
                      "broad_solutions": sum(row["broad_border_second_bit_chart"]["status"] == "SOLVED" for row in report["results"])}))


if __name__ == "__main__":
    main()

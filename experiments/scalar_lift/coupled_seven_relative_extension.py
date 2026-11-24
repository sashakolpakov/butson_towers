"""Exact coupled quadratic extensions of the known order-seven source.

These six controls deform known matrices; they do not construct either
of the missing large prime-order sources.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from experiments.scalar_lift.fifth_source_assembly import exact_column_gram, root_vectors
from experiments.scalar_lift.mixed_cubic_boundaries import seven_border
from experiments.scalar_lift.nondegenerate_cubic_border_lift import actual_affine_masks
from experiments.scalar_lift.quadratic_layer_obstructions import binary_rank
from experiments.scalar_lift.relative_source_tower import ZERO, sparse_row_product


def split_dyadic(source, upper_phase):
    """H=A+eta B, eta=zeta_upper^3, eta^2=zeta_lower^3."""
    if upper_phase not in (24, 48):
        raise ValueError("this certificate covers upper phases 24 and 48")
    source = np.asarray(source, dtype=np.int64) % upper_phase
    lower = upper_phase // 2
    a = np.where(source % 2 == 0, source // 2, ZERO)
    b = np.where(source % 2 == 1, (source - 3) // 2 % lower, ZERO)
    reconstructed = np.where(a != ZERO, 2 * a, 2 * b + 3) % upper_phase
    if not np.array_equal(reconstructed, source):
        raise ValueError("the dyadic split failed to reconstruct its source")
    return a, b


def expected_gram(order, phase):
    result = np.zeros((order, order, root_vectors(phase).shape[1]), dtype=np.int64)
    result[np.arange(order), np.arange(order), 0] = order
    return result


def extension_control(upper_phase, branch):
    lower = upper_phase // 2
    original = seven_border(upper_phase)[0]
    mask = actual_affine_masks()[branch]
    if np.any(original % 2):
        raise ValueError("the starting source does not descend to the lower field")
    if not np.array_equal(exact_column_gram(original // 2, lower), expected_gram(7, lower)):
        raise ValueError("the starting lower-field Butson check failed")

    # Check the Laurent switching law independently in this exact gauge.
    roots = root_vectors(upper_phase)
    laurent_count = 0
    for i in range(7):
        for j in range(i):
            powers = mask[i] - mask[j]
            for power in np.unique(powers):
                occupied = powers == power
                if np.any(roots[(original[i, occupied] - original[j, occupied]) % upper_phase].sum(axis=0)):
                    raise ValueError("the exact Laurent switching identity failed")
                laurent_count += 1

    source = (original + 3 * mask) % upper_phase
    if np.any(source % 3 == 2) or not np.array_equal(source % 3, original % 3):
        raise ValueError("switching left the admitted two-cubic-color alphabet")
    a, b = split_dyadic(source, upper_phase)
    if not np.all((a != ZERO).astype(int) + (b != ZERO).astype(int) == 1):
        raise ValueError("the two root layers do not partition the entry support")
    aa, bb = sparse_row_product(a, a, lower), sparse_row_product(b, b, lower)
    ab, ba = sparse_row_product(a, b, lower), sparse_row_product(b, a, lower)
    negative_carry = root_vectors(lower)[(np.arange(ab.shape[-1]) - 3) % lower]
    if not np.array_equal(aa + bb, expected_gram(7, lower)):
        raise ValueError("AA*+BB*=7I failed")
    if np.any(ba + ab @ negative_carry):
        raise ValueError("BA*+xi^-1 AB*=0 failed")
    nonzero = np.argwhere(np.any(ab != 0, axis=-1))
    if len(nonzero) != 8:
        raise ValueError("the expected genuine cross-coupling was lost")
    expected = expected_gram(7, upper_phase)
    if (not np.array_equal(sparse_row_product(source, source, upper_phase), expected)
            or not np.array_equal(exact_column_gram(source, upper_phase), expected)):
        raise ValueError("the independent complete upper-field Gram failed")

    collapsed = np.maximum(a, b)
    collapse_defect = sparse_row_product(collapsed, collapsed, lower) - expected_gram(7, lower)
    if not np.any(collapse_defect):
        raise ValueError("the fixed-coefficient specialization unexpectedly stayed Hadamard")
    parity = (source - source[:, :1] - source[:1, :] + source[0, 0]) % 2
    rank = binary_rank(parity)
    if rank != 2:
        raise ValueError("the dephased parity rank is not two")
    i, j = map(int, np.argwhere(parity)[0])
    cross_ratio = int((source[i, j] - source[i, 0] - source[0, j] + source[0, 0]) % upper_phase)
    cross_i, cross_j = map(int, nonzero[0])
    defect_i, defect_j = map(int, np.argwhere(np.any(collapse_defect != 0, axis=-1))[0])
    return {
        "order": 7, "lower_phase": lower, "upper_phase": upper_phase, "branch": branch,
        "original_lower_exponents": (original // 2).tolist(), "signed_switch_mask": mask.tolist(),
        "source_exponents": source.tolist(), "A_exponents": a.tolist(), "B_exponents": b.tolist(),
        "zero_sentinel": ZERO, "layer_support_sizes": [int(np.count_nonzero(a != ZERO)), int(np.count_nonzero(b != ZERO))],
        "supports_disjoint_exhaustive_and_source_reconstructed": True,
        "admitted_upper_exponent_residues_mod3": [0, 1], "all_cubic_colors_preserved": True,
        "exact_Laurent_coefficients_checked": laurent_count,
        "lower_input_Gram_verified": True, "upper_row_and_column_Grams_verified": True,
        "AA_dagger_plus_BB_dagger_equals_7I": True, "BA_dagger_plus_xi_inverse_AB_dagger_zero": True,
        "AB_dagger_nonzero_positions": int(len(nonzero)),
        "AB_dagger_witness": {"row": cross_i, "column": cross_j, "lower_power_basis_coefficients": ab[cross_i, cross_j].tolist()},
        "collapsed_A_plus_B_is_Hadamard": False,
        "collapse_Gram_defect_witness": {"row": defect_i, "column": defect_j,
                                        "lower_power_basis_coefficients": collapse_defect[defect_i, defect_j].tolist()},
        "dephased_parity_matrix": parity.tolist(), "dephased_parity_rank": rank,
        "odd_cross_ratio": {"rows": [0, i], "columns": [0, j], "upper_exponent": cross_ratio},
        "row_column_phase_gauge_descent_to_lower_field_excluded": True,
    }


def certificate():
    return {
        "schema": "coupled-seven-relative-extension-v1",
        "relative_generator": "eta=zeta_Q^3, xi=eta^2=zeta_(Q/2)^3",
        "relative_equations": ["AA*+BB*=7I", "BA*+xi^-1 AB*=0"],
        "controls": [extension_control(phase, branch) for phase in (24, 48) for branch in range(3)],
        "scope": "Known order-seven Laurent switch, specialized in two quadratic towers; no new large source or real Hadamard order.",
        "fixed_coefficient_warning": "The Laurent switch uses eta^-1=xi^-1 eta, so B depends on xi. It is not a fixed A+tB family valid for every unit t.",
        "large_prime_sources_constructed": False,
    }


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
        raise ValueError("the coupled order-seven certificate failed replay")
    print(json.dumps({"exact_coupled_controls": len(report["controls"]), "new_large_source": False}))


if __name__ == "__main__":
    main()

from __future__ import annotations

from itertools import product

import numpy as np

from experiments.scalar_lift.paley_sextic_operator import (
    build_paley_sextic_operator_report,
    dense_sextic_operator_numerator,
    flat_sign_images,
    necessary_regular_block_sum_audit,
)


def test_dense_sextic_operator_meets_local_requirements_but_has_no_sign_frame():
    report = build_paley_sextic_operator_report()
    assert report["rank_C_modulo_167"] == 170
    assert report["rank_stacked_B_BD_modulo_167"] == 4
    assert report["dimension_of_flat_vector_span"] == 4
    assert report["flat_vector_count"] == 8
    assert report["R_reduced_denominators_per_row"] == {"2": 1, "334": 667}
    assert report["full_sign_frame_excluded_for_this_operator"]
    assert not report["new_real_hadamard_order_claimed"]
    ranks = [row["rank_stacked_B_BD_modulo_167"] for row in report["quaternion_orientation_audit"]]
    assert ranks.count(2) == ranks.count(4) == 8
    candidate = report["compatible_opposite_chirality_candidate"]
    assert candidate["rank_C_modulo_167"] == 168
    assert candidate["C_row_sums_by_block"] == [501, -167, -167, -167]
    assert candidate["full_sign_frame_excluded"]
    assert candidate["integer_left_null_contradiction"] == 16032
    assert not candidate["full_sign_frame_constructed"]
    audit = report["known_hadamard_frame_audit"]
    assert len(audit) == 8
    assert max(max(row["native_flat_columns_by_Hall_mask"]) for row in audit) == 1
    assert max(max(row["border_grouped_flat_columns_by_Hall_mask"]) for row in audit) == 3
    candidate_audit = next(row for row in audit if row["D_upper_triangle"] == [1, 1, 1, -1, 1, -1])
    assert candidate_audit["native_flat_columns_by_Hall_mask"] == [0] * 16
    assert candidate_audit["border_grouped_flat_columns_by_Hall_mask"] == [0] * 16


def test_flat_congruence_matches_direct_images_and_binary_kernel():
    p = 167
    c = dense_sextic_operator_numerator()
    # All constant-block candidates, and perturbations with nontrivial block sums.
    for block_signs in product((-1, 1), repeat=4):
        base = np.repeat(block_signs, p)
        for perturbations in ((), (0,), (0, 1), (0, p, 2 * p, 3 * p)):
            signs = base.copy()
            signs[list(perturbations)] *= -1
            numerator = p * signs + c @ signs
            direct = bool(np.all(np.abs(numerator) == 2 * p))
            binary = (1 - signs) // 2
            criterion = bool(np.all((c @ binary) % p == 0) and np.sum(binary) % 2 == 1)
            assert flat_sign_images(c, signs) == direct == criterion


def test_extremal_frame_sums_reduce_to_four_orthogonal_directions():
    audit = necessary_regular_block_sum_audit()
    assert len(audit["integral_sigma_candidates"]) == 16
    assert len(audit["sigma_candidates_with_odd_R_image_block_sums"]) == 8
    assert audit["required_column_multiplicity_per_antipodal_direction"] == 167
    assert audit["all_16_blocks_must_have_constant_row_and_column_sums_M"]
    assert not audit["circulant_blocks_assumed"]
    assert not audit["frame_constructed"]

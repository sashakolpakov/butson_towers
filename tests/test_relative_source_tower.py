import numpy as np
import pytest

from experiments.scalar_lift.fifth_source_assembly import (
    complete_tensor_source,
    extended_eight_seed,
)
from experiments.scalar_lift.relative_source_tower import (
    cross_ratio,
    reconstruct_gram,
    reconstruct_source,
    relative_control,
    relative_gram_coefficients,
    sparse_row_product,
    split_layers,
    target_metadata,
    tower_compatibility,
)


@pytest.mark.parametrize("lower,upper", ((20, 40), (40, 80), (20, 80)))
def test_sparse_relative_layers_reconstruct_and_satisfy_every_exact_gram_coefficient(lower, upper):
    source = complete_tensor_source(upper)
    report = relative_control(source, lower, upper)
    assert report["source_order"] == 40
    assert report["all_nonconstant_Ck_zero"]
    assert report["independent_upper_row_and_column_Grams_verified"]
    for rule in report["fifth_color_channels"]:
        j = rule["layer_index"]
        assert rule["allowed_lower_exponent_residues_mod5"] == [a for a in range(5) if ((upper // lower) * a + j) % 5]


def test_nonhadamard_matrix_checks_the_negative_carry_and_every_relative_coefficient():
    source = np.array([[1, 7, 19], [38, 2, 71]], dtype=np.int64)
    for lower in (20, 40):
        layers = split_layers(source, lower, 80)
        coefficients = relative_gram_coefficients(layers, lower, 80)
        assert np.any(coefficients[1:])
        assert np.array_equal(reconstruct_source(layers, lower, 80), source)
        assert np.array_equal(reconstruct_gram(coefficients, lower, 80), sparse_row_product(source, source, 80))
    assert tower_compatibility(source)["all_four_Gram_coefficients_identical"]


def test_two_binary_splits_agree_with_four_layers_for_the_primitive_seed():
    source = extended_eight_seed(80)
    report = tower_compatibility(source)
    assert report["all_four_sparse_layers_identical"]
    assert report["same_source_reconstructed"]
    assert cross_ratio(source, 80)["root_order"] == 16


def test_prime_targets_require_the_stated_extension_without_claiming_existence():
    reports = target_metadata()
    assert [r["excluded_lower_phase"] for r in reports[:2]] == [20, 40]
    assert all(r["quadratic_A_B_dagger_must_be_nonzero"] for r in reports[:2])
    assert reports[2]["descent_to_phase40_not_excluded_by_this_test"]
    assert not reports[2]["existence_at_phase40_or80_asserted"]

import numpy as np
import pytest

from experiments.scalar_lift.quadratic_layer_obstructions import (
    binary_rank,
    even_weighing_control,
    fourier_control,
    parity_dephase,
    rank_two_block_identity,
)


@pytest.mark.parametrize("phase", (40, 80))
def test_cross_coupling_changes_under_gauge_but_rank_one_collapse_does_not(phase):
    report = fourier_control(phase)
    first, second = report["controls"]
    assert not first["relative_equations"]["ABstar_is_nonzero"]
    assert second["relative_equations"]["ABstar_is_nonzero"]
    assert len(second["nonzero_ABstar_entries"]) == 2
    assert all(item["rank_one_collapse"]["dephased_parity_rank"] == 1 for item in (first, second))


@pytest.mark.parametrize("phase", (40, 80))
def test_even_control_really_has_two_nontrivial_scalar_gram_layers(phase):
    assert even_weighing_control(phase)["both_layer_Grams_equal_identity"]


def test_parity_rank_two_is_distinguished_and_phase_gauges_preserve_the_invariant():
    # This only tests a support invariant, not a Butson source.
    source = np.array([[1, 1, 1], [1, 2, 1], [1, 1, 2]])
    _, parity = parity_dephase(source, 40)
    assert binary_rank(parity) == 2
    gauged = (source + 5 * np.array([1, 0, 1])[:, None] + 5 * np.array([0, 1, 1])[None, :]) % 40
    _, other = parity_dephase(gauged, 40)
    assert np.array_equal(other, parity)


@pytest.mark.parametrize("phase", (40, 80))
def test_three_way_coefficients_match_an_independent_nonzero_full_gram(phase):
    report = rank_two_block_identity(phase)
    assert report["all_block_Gram_formulas_match_full_exact_Gram"]
    assert report["nonzero_cross_residual_present"]

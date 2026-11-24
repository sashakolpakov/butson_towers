import numpy as np
import pytest

from experiments.scalar_lift.coupled_seven_relative_extension import (
    extension_control,
    split_dyadic,
)
from experiments.scalar_lift.fifth_source_assembly import root_vectors
from experiments.scalar_lift.relative_source_tower import ZERO, sparse_row_product


@pytest.mark.parametrize("phase", (24, 48))
@pytest.mark.parametrize("branch", range(3))
def test_complete_order_seven_extension_has_essential_relative_coupling(phase, branch):
    report = extension_control(phase, branch)
    assert report["AB_dagger_nonzero_positions"] == 8
    assert report["dephased_parity_rank"] == 2
    assert report["odd_cross_ratio"]["upper_exponent"] % 2
    assert not report["collapsed_A_plus_B_is_Hadamard"]
    assert report["all_cubic_colors_preserved"]


def test_dyadic_carry_reconstructs_a_nonhadamard_rectangular_gram():
    source = np.array([[0, 3, 23], [7, 12, 5]], dtype=np.int64)
    a, b = split_dyadic(source, 24)
    assert np.all((a != ZERO) ^ (b != ZERO))
    c0 = sparse_row_product(a, a, 12) + sparse_row_product(b, b, 12)
    ab = sparse_row_product(a, b, 12)
    carry = root_vectors(12)[(np.arange(ab.shape[-1]) - 3) % 12]
    c1 = sparse_row_product(b, a, 12) + ab @ carry
    upper_roots = root_vectors(24)
    reconstructed = c0 @ upper_roots[2 * np.arange(c0.shape[-1])]
    reconstructed += c1 @ upper_roots[(2 * np.arange(c1.shape[-1]) + 3) % 24]
    assert np.any(c1)
    assert np.array_equal(reconstructed, sparse_row_product(source, source, 24))

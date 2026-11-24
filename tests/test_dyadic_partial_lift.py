import numpy as np
import pytest

from experiments.scalar_lift.dyadic_partial_lift import (
    phase_blocks,
    verify_dyadic_source_factorization,
    verify_partial_lift,
)


def test_exact_partial_alphabets_and_source_control():
    for n in (4, 8, 16):
        report = verify_partial_lift(n)
        assert report["flat_alphabet_size"] == 4 * n
        assert report["source_order_three_control"]["signs_and_gram_verified"]
        assert report["source_order_thirteen_control"]["signs_and_gram_verified"]


def test_cross_gram_represents_scalar_quotients():
    # This verifies the multiplicative identity needed for arbitrary source
    # row inner products, not just the individual flat-block Gram checks.
    for n in (4, 8, 16):
        blocks = phase_blocks(n)
        reference = blocks[0]
        for a in range(6 * n):
            for b in range(6 * n):
                assert np.array_equal(
                    blocks[a] @ blocks[b].T,
                    blocks[(a - b) % (6 * n)] @ reference.T,
                )


def test_source_checker_rejects_lossy_exponents_and_false_gram():
    for value in (0.5, 0.0, True):
        with pytest.raises(ValueError, match="exact integers"):
            verify_dyadic_source_factorization([[value]], 48)
    with pytest.raises(ValueError, match="norm equation"):
        verify_dyadic_source_factorization([[0, 0], [0, 0]], 48)
    with pytest.raises(ValueError, match="excluded phase"):
        verify_dyadic_source_factorization([[2]], 48)

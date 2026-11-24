import numpy as np
import pytest
from sympy import Matrix, cyclotomic_poly

from experiments.scalar_lift.dyadic_cyclotomic_action import (
    action_numerator,
    build_report,
    verify_action,
)


@pytest.mark.parametrize(("phase", "size"), [(24, 4), (48, 8), (96, 16)])
def test_exact_orthogonal_actions_and_independent_characteristic_polynomial(phase, size):
    numerator = action_numerator(size)
    report = verify_action(numerator, phase)
    action = Matrix(numerator.tolist()) / 2
    characteristic = action.charpoly()
    assert characteristic.as_expr() == cyclotomic_poly(phase, characteristic.gen)
    assert action * action.T == Matrix.eye(2 * size)
    assert report["faithful_field_degree"] == 2 * size


def test_corrupt_or_nonfaithful_actions_are_rejected():
    corrupted = action_numerator(8)
    corrupted[0, 0] += 1
    with pytest.raises(ValueError, match="metric"):
        verify_action(corrupted, 48)
    with pytest.raises(ValueError, match="cyclotomic identity"):
        verify_action(2 * np.eye(16, dtype=np.int64), 48)


def test_large_orders_are_conditional_and_no_sign_witness_is_claimed():
    report = build_report()
    assert not report["source_matrices_constructed"]
    assert not report["compatible_sign_bases_constructed"]
    assert not report["new_real_Hadamard_order_constructed"]
    for route in report["conditional_target_examples"]:
        assert route["order"] == route["source_order"] * route["block_order"]

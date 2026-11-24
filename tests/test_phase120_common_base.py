import numpy as np
import pytest

from experiments.scalar_lift.phase120_common_base import (
    action_power,
    admitted_exponents,
    build_report,
    common_base,
    lift_source,
    phase_blocks,
    source_fifteen,
)


def test_complete_action_alphabet_and_genuine_phase_source_control():
    report = build_report()
    assert report["faithful_order_120_orthogonal_action"]
    assert report["dimension"] == 32
    assert all(r["all_14400_cross_Grams_verified"] for r in report["frames"])
    assert all(len(r["allowed_exponents"]) == 64 for r in report["frames"])
    assert report["control"]["primitive_cross_ratio_exponent"] == 119
    assert report["control"]["full_integer_Gram_verified"]
    assert report["coupled_control"]["real_order"] == 1120
    assert report["coupled_control"]["full_integer_Gram_verified"]
    assert not report["conditional_high_route"]["source_constructed"]


def test_relative_tensor_and_raw_frame_shift():
    data = common_base()
    for exponent, key, denominator in ((15, "dyadic_action", 1), (40, "cubic_numerator", 2)):
        numerator, actual_denominator = action_power(exponent)
        assert np.array_equal(denominator * numerator, actual_denominator * data[key])
    raw, shifted = phase_blocks(False), phase_blocks(True)
    assert all(np.array_equal(shifted[e], raw[(e + 24) % 120]) for e in range(120))
    assert admitted_exponents(False) == [e for e in range(120) if e % 3 in (0, 1) and e % 5]


def test_source_validation_does_not_turn_an_arbitrary_matrix_into_a_control():
    with pytest.raises(ValueError, match="excluded phase"):
        lift_source([[1]])
    with pytest.raises(ValueError, match="exact integers"):
        lift_source([[0.0]])
    with pytest.raises(ValueError, match="cyclotomic Gram"):
        lift_source([[0, 0], [0, 0]])
    source = source_fifteen()
    source[0][0] = (source[0][0] + 15) % 120
    with pytest.raises(ValueError, match="cyclotomic Gram"):
        lift_source(source, False)

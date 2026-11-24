import json
from pathlib import Path

import numpy as np
import pytest

from experiments.scalar_lift.fifth_source_descent import (
    certificate,
    crt_coordinates,
    descend_source,
    split_coefficients,
)


def test_full_iterated_control_and_finite_scalar_identities_replay():
    path = Path(__file__).resolve().parents[1] / "certificates/scalar-lift/fifth-source-descent.json"
    report = certificate()
    assert report == json.loads(path.read_text())
    stages = report["complete_source_control"]["stages"]
    assert [(s["order"], s["phase"]) for s in stages] == [(40, 80), (80, 40), (160, 20)]
    assert all(s["complete_row_and_column_Grams_verified"] for s in stages)
    assert not report["large_Butson_source_or_new_real_order_constructed"]


def test_every_admitted_scalar_has_one_lower_coefficient_of_the_same_fifth_color():
    for upper in (40, 80):
        lower, d = upper // 2, upper // 10
        for exponent in range(upper):
            if exponent % 5 == 0:
                continue
            color, dyadic = crt_coordinates(exponent, upper)
            a, b = split_coefficients([[exponent]], upper)
            assert (a[0, 0] >= 0) ^ (b[0, 0] >= 0)
            value = int(b[0, 0] if dyadic % 2 else a[0, 0])
            assert value % 5 != 0
            assert pow(d, -1, 5) * value % 5 == color
            assert (2 * value + 5 * (dyadic % 2)) % upper == exponent
            assert 0 <= value < lower


def test_descent_checks_the_source_gram_alphabet_and_exact_input():
    with pytest.raises(ValueError, match="excluded dyadic"):
        descend_source([[0]], 80)
    with pytest.raises(ValueError, match="exact integers"):
        descend_source([[1.5]], 80)
    with pytest.raises(ValueError, match="cyclotomic Gram"):
        descend_source([[1, 1], [1, 1]], 80)
    result = descend_source([[np.int64(1)]], 80)
    assert result.shape == (2, 2)
    assert np.all(result % 5)

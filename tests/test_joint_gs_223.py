from __future__ import annotations

from random import Random

import numpy as np
import pytest

from experiments.scalar_lift.joint_gs_223 import (
    SYMBOLS,
    goethals_seidel_array,
    orbit_data,
    row_sum_types,
)


@pytest.mark.parametrize(("prime", "multiplier", "expected_orbits", "expected_shifts", "expected_types"),
                         [(223, 3, 75, 37, 30), (491, 7, 71, 35, 8)])
def test_compressed_coefficients_equal_every_expanded_correlation(prime, multiplier, expected_orbits,
                                                               expected_shifts, expected_types):
    orbits, index, shifts, pairs, coefficients = orbit_data(prime, multiplier)
    assert (len(orbits), len(shifts)) == (expected_orbits, expected_shifts)
    assert sorted(x for orbit in orbits for x in orbit) == list(range(prime))
    random = Random(prime)
    symbols = [SYMBOLS[random.randrange(4)] for _ in orbits]
    # A numeric specialization also tests both diagonal and mixed terms.
    values = [2 * a + 5 * b for a, b in symbols]
    expanded = [values[index[x]] for x in range(prime)]
    compressed_products = np.array([values[i] * values[j] for i, j in pairs])
    actual = coefficients @ compressed_products
    for k, shift in enumerate(shifts):
        assert int(actual[k]) == sum(expanded[x] * expanded[(x + shift) % prime] for x in range(prime))
    types = row_sum_types(prime, multiplier)
    assert len(types) == expected_types
    for row in types:
        a, b = row[:4], row[4:]
        assert sum(x * x for x in a) == 3 * prime
        assert sum(x * x for x in b) == prime
        assert sum(x * y for x, y in zip(a, b, strict=True)) == 0
        assert all((x + y) % 2 == 1 for x, y in zip(a, b, strict=True))


def test_goethals_seidel_gram_and_design_to_pair_identities():
    sequences = np.array([[2, -1, 3, 0, 1], [-2, 0, 1, 3, -1],
                          [1, 1, -2, 0, 2], [3, 1, 0, -1, -2]])
    matrix = goethals_seidel_array(sequences)
    correlations = [sum(int(row[x] * row[(x + shift) % 5]) for row in sequences for x in range(5))
                    for shift in range(5)]
    expected_block = np.array([[correlations[(j - i) % 5] for j in range(5)] for i in range(5)])
    assert np.array_equal(matrix @ matrix.T, np.kron(np.eye(4, dtype=int), expected_block))
    a = np.array([[0, 1, 1], [-1, 0, 1], [-1, -1, 0], [-1, 1, -1]])
    b = np.eye(4, 3, dtype=int)
    u, v = goethals_seidel_array(a), goethals_seidel_array(b)
    unit = np.eye(12, dtype=int)
    assert np.array_equal(u @ u.T, 9 * unit)
    assert np.array_equal(v @ v.T, 3 * unit)
    assert not np.any(u @ v.T + v @ u.T)
    first, second = u - v, u + v
    assert np.all(abs(first) == 1) and np.all(abs(second) == 1)
    assert np.array_equal(first @ second.T + second @ first.T, 12 * unit)

import numpy as np
from sympy import Matrix, eye

from experiments.scalar_lift.alpoge_mixed_basis_probe import (
    binary_rank,
    canonical_fraction,
    entry_statistics,
    reflected_partner,
    restricted_entry,
    tensor_minor,
    twice_root_action,
)
from experiments.scalar_lift.forward_source_registry import load_source


def sylvester(order):
    result = np.ones((1, 1), dtype=np.int64)
    while len(result) < order:
        result = np.block([[result, result], [result, -result]])
    return result


def test_dense_reflection_formula_against_independent_rational_conjugacy():
    frame = sylvester(8)
    vector = np.array([1, -2, 0, 1, 2, -1, 1, 1], dtype=np.int64)
    numerator, denominator = reflected_partner(frame, vector)
    v = Matrix(vector.tolist())
    q = eye(8)-2*v*v.T/(v.T*v)[0]
    r = Matrix(twice_root_action(np.eye(8, dtype=np.int64)).tolist())/2
    exact = Matrix(numerator.tolist())/denominator
    f = Matrix(frame.tolist())
    assert exact == q*r*q.T*f
    assert exact*exact.T == 8*eye(8)
    assert exact*f.T+f*exact.T == 8*eye(8)
    witness = tensor_minor(vector)
    assert witness and witness['integer_minor'] != 0
    n, d = canonical_fraction(numerator, denominator)
    assert Matrix(n.tolist())/d == exact


def test_positive_control_constructs_full_known_order52_from_verified_source():
    frame = sylvester(4)
    numerator, denominator = reflected_partner(frame)
    assert entry_statistics(numerator, denominator)['all_entries_signs']
    partner = numerator//denominator
    source = load_source('unreal13-q6')
    result = np.array([[restricted_entry(source, frame, partner, 1, i, j)
                        for j in range(52)] for i in range(52)], dtype=object)
    assert np.all(np.abs(result) == 1)
    assert np.array_equal(result @ result.T, 52*np.eye(52, dtype=np.int64))


def test_rational_nonflat_entries_are_counted_and_witnessed_exactly():
    numerator = np.array([[2, -2, 0], [1, 4, -3], [0, 0, 2]])
    stats = entry_statistics(numerator, 2)
    assert stats['sign_entries'] == 3
    assert stats['zero_entries'] == 3
    assert stats['other_integer_entries'] == 1
    assert stats['noninteger_entries'] == 2
    assert stats['first_failed_entry'] == {'row': 0, 'column': 2, 'numerator': 0, 'denominator': 2}


def test_binary_rank_matches_exhaustive_span():
    rows = np.array([[1, 1, 0], [0, 1, 1], [1, 0, 1]], dtype=np.uint8)
    span = {tuple(((np.array([(mask >> i) & 1 for i in range(3)]) @ rows) % 2).tolist())
            for mask in range(8)}
    assert len(span) == 2**binary_rank(rows) == 4

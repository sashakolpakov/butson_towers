import numpy as np
import pytest
from sympy import Matrix, eye, kronecker_product

from experiments.scalar_lift.exact_norm_descent import from_butson, positive_control
from experiments.scalar_lift.relative_scalar_restriction import edge, regular_matrix
from experiments.scalar_lift.vanishing_collision_search import integral_roots


def test_actual_source_tower_keeps_norm_and_yields_known_sign_endpoint():
    report = positive_control()
    assert report["direct_and_tower_conjugacy_verified"]
    assert report["endpoint"]["sign_Hadamard_verified"]
    assert report["endpoint"]["dimension"] == 16
    assert report["terminal"]["parent"] is not None
    assert not report["new_unknown_order_constructed"]


def test_nonconstant_incoming_metric_transfers_by_trace_transitivity():
    source = from_butson([[1]], 30)
    middle = source.restrict(6)
    assert np.any(middle.metric.numerators[:, :, 1])
    terminal = middle.restrict(2)
    direct = source.restrict(2)
    exponents = [a+5*b for a in range(4) for b in range(2)]
    roots = np.asarray(integral_roots(30), dtype=object)
    path_basis = Matrix(roots[exponents].T.tolist())
    change = edge(30, 2).inverse_basis_matrix*path_basis
    assert change.T*regular_matrix(direct.metric, 2)*change == regular_matrix(terminal.metric, 2)
    assert change.inv()*regular_matrix(direct.operator, 2)*change == regular_matrix(terminal.operator, 2)


def test_normalization_and_flatness_are_independent_certificates():
    exponents = np.outer(np.arange(4), np.arange(4)) % 4
    terminal = from_butson(exponents, 4).restrict(2)
    # No frame: the weighted norm is exact but is the wrong ordinary scale.
    unframed = terminal.endpoint(eye(8), eye(8))
    assert unframed["weighted_norm_verified"]
    assert not unframed["Euclidean_norm_verified"]
    # A different exact norm frame yields the right scale but has zeros.
    frame = Matrix([[1, 1], [1, -1]])
    s = kronecker_product(eye(4), frame)
    q = Matrix([[3, -4], [4, 3]])/5
    moved = terminal.endpoint(s, kronecker_product(eye(4), q))
    assert moved["Euclidean_norm_verified"]
    assert moved["flat_entry_defects"] > 0
    assert not moved["sign_Hadamard_verified"]


def test_inexact_or_nonorthogonal_sources_are_rejected():
    with pytest.raises(ValueError, match="norm identity"):
        from_butson([[0, 0], [0, 0]], 6)
    with pytest.raises(ValueError, match="exact exponent"):
        from_butson([[0.5]], 6)

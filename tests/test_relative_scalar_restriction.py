import numpy as np
import pytest
from sympy import Matrix, Rational, diag

from experiments.scalar_lift.relative_scalar_restriction import (
    edge,
    regular_matrix,
    restriction_in_basis,
)


def test_requested_odd_removal_edges_have_exact_actions_and_trace_metrics():
    for upper, lower, width, inverse_denominator in (
        (30, 10, 2, 3), (30, 6, 4, 5), (30, 2, 8, 15),
        (120, 40, 2, 3), (120, 24, 4, 5), (120, 8, 8, 15),
    ):
        restriction = edge(upper, lower)
        assert restriction.width == width
        assert restriction.actions.numerators.shape == (upper, width, width, restriction.lower_degree)
        assert restriction.metric.denominator == 1
        assert restriction.inverse_metric.denominator == inverse_denominator
        assert restriction.basis_matrix * restriction.inverse_basis_matrix == Matrix.eye(restriction.upper_degree)
        assert restriction.proof()["all_upper_scalar_powers_verified_by_exact_basis_conjugacy"] == upper
        assert not restriction.proof()["ordinary_identity_metric_or_flat_entries_claimed"]
        assert restriction.trace_exponent(0).numerators.tolist() == [width] + [0] * (restriction.lower_degree - 1)


def test_sextic_to_rational_action_and_weighted_complete_source_Gram():
    restriction = edge(6, 2)
    assert regular_matrix(restriction.action(1), 2) == Matrix([[0, -1], [1, 1]])
    assert regular_matrix(restriction.metric, 2) == Matrix([[2, 1], [1, 2]])
    inverse = regular_matrix(restriction.inverse_metric, 2)
    assert inverse == Matrix([[Rational(2, 3), Rational(-1, 3)],
                              [Rational(-1, 3), Rational(2, 3)]])
    # The existing complete BH3 uses exponent4 on its diagonal and2 off.
    exponents = 2 + 2 * np.eye(3, dtype=np.int64)
    blocks = [[regular_matrix(restriction.action(int(e)), 2) for e in row] for row in exponents]
    source = Matrix.vstack(*(Matrix.hstack(*row) for row in blocks))
    weight = diag(inverse, inverse, inverse)
    assert source * weight * source.T == 3 * weight
    assert any(value == 0 for value in source)
    assert source * source.T != 6 * Matrix.eye(6)


def test_coordinate_distinct_path_basis_retains_action_and_metric_congruence():
    canonical = edge(30, 2)
    path = restriction_in_basis(30, 2, tuple(i + 3*j for i in range(2) for j in range(4)))
    change = canonical.inverse_basis_matrix * path.basis_matrix
    assert change != Matrix.eye(8)
    assert change.T * regular_matrix(canonical.metric, 2) * change == regular_matrix(path.metric, 2)
    for power in range(30):
        assert change.inv() * regular_matrix(canonical.action(power), 2) * change == regular_matrix(path.action(power), 2)
    assert canonical.metadata()["exact_definition_sha256"] != path.metadata()["exact_definition_sha256"]


def test_identity_and_equal_rational_fields_keep_distinct_declared_root_labels():
    identity = edge(30, 30)
    assert identity.width == 1
    assert identity.metric.numerators.tolist() == [[[1, 0, 0, 0, 0, 0, 0, 0]]]
    first, second = edge(6, 1), edge(6, 2)
    assert np.array_equal(first.actions.numerators, second.actions.numerators)
    assert np.array_equal(first.metric.numerators, second.metric.numerators)
    assert first.metadata()["declared_lower_root_label_retained"] == 1
    assert second.metadata()["declared_lower_root_label_retained"] == 2
    assert first.metadata()["lower_field_conductor"] == second.metadata()["lower_field_conductor"] == 1
    with pytest.raises(ValueError, match="dividing"):
        edge(30, 8)
    with pytest.raises(ValueError, match="singular"):
        restriction_in_basis(30, 2, (0,) * 8)

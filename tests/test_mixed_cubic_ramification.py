import numpy as np

from experiments.scalar_lift.mixed_cubic_ramification import (
    certificate,
    core_derivative_formula,
)


def test_exact_first_ramification_obstruction_at_both_actual_high_source_orders():
    report = certificate()
    assert [(target["source_order"], target["phase"]) for target in report["targets"]] == [
        (1499, 96), (2999, 48), (2999, 96)]
    for target in report["targets"]:
        assert target["no_root_phase_lift_modulo_pi_squared_for_this_fixed_boundary_and_core_support"]
        assert target["actual_border_first_bits_are_zero"]
        assert target["nonzero_A_entries"] > 0
        assert not target["general_Butson_source_order_excluded"]
    assert report["oracle"]["full_boundary_derivatives_checked"] == 20
    assert not report["native_probe_run"]


def test_sum_difference_equations_do_not_hide_a_full_gram_derivative():
    triangle = np.array([[0, 1, 1, 0], [1, 0, 1, 0],
                         [1, 1, 0, 0], [0, 0, 0, 0]])
    diagonal_0, diagonal_w, cross_0, cross_w = core_derivative_formula(triangle, triangle)
    assert np.array_equal(diagonal_0, triangle)
    assert np.array_equal(cross_0, triangle)
    assert not np.any(diagonal_w)
    assert not np.any(cross_w)

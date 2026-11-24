import numpy as np

from experiments.scalar_lift.mixed_cubic_ramification import MULTIPLY, product
from experiments.scalar_lift.nondegenerate_cubic_projectors import (
    binary_affine_solutions,
    certify_retained,
    retained_templates,
    template_checks,
)


def test_all_literal_first_lifts_pass_complete_field_checks_and_remove_isolated_coordinates():
    replay = certify_retained()
    assert replay["verified_first_lifts"] == 15
    assert replay["local_full_Gram_and_derivative_checks"] == 60
    assert not replay["projector_family_claimed"]
    for item in retained_templates():
        audit = template_checks(item["C"], item["B"], item["a"])
        assert item["C"].shape == (12, 12) and item["a"].sum() == 4
        assert audit["first_lift_acceptable"]
        assert audit["T_bipartite_connected"]
        assert min(audit["W_row_weights"]) > 0 and min(audit["W_column_weights"]) > 0


def test_affine_solver_and_extra_projector_condition_are_distinct():
    item = retained_templates()[0]
    c, b, active = item["C"], item["B"], item["a"]
    rows, columns, particular, _ = binary_affine_solutions(c, active)
    solved = np.zeros_like(c)
    solved[rows, columns] = [(particular >> j) & 1 for j in range(len(rows))]
    assert template_checks(c, solved, active)["linear_criterion_verified"]
    f = MULTIPLY[2, c @ b.T % 2] ^ MULTIPLY[3, b @ c.T % 2]
    defect = product(f, f) ^ f
    audit = template_checks(c, b, active)
    row, column, value = audit["projector_failure_entry_over_F4"]
    assert defect[row, column] == value != 0
    assert not audit["projector_F_squared_equals_F"]

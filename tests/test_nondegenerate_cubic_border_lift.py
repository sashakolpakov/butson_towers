from experiments.scalar_lift.cubic_core_support_lift import template
from experiments.scalar_lift.nondegenerate_cubic_border_lift import (
    audit_affine_masks,
    solve_candidate,
)
from experiments.scalar_lift.nondegenerate_cubic_projectors import retained_templates


def test_three_actual_petrescu_branches_have_exact_laurent_gram_cancellation():
    reports = audit_affine_masks()
    assert len(reports) == 3
    assert all(report["all_unit_parameters_preserve_the_exact_border_Gram"] for report in reports)
    assert all(report["exact_Laurent_coefficients_checked"] == 37 for report in reports)


def test_known_ten_row_control_survives_moving_border_and_independent_root_gram():
    report = solve_candidate({**template(), "template_sha256": "known-ten-row-positive-control"})
    for part in [report["broad_border_second_bit_chart"], *report["individually_exact_border_branches"]]:
        assert part["status"] == "SOLVED"
        assert part["all_four_representative_Gram_coefficients_verified_by_inverse_roots"]
        assert part["all_kernel_basis_neighbors_verified_by_inverse_root_Gram"]


def test_first_twelve_row_template_fails_even_with_all_border_second_bits():
    report = solve_candidate(retained_templates()[0])
    broad = report["broad_border_second_bit_chart"]
    assert broad["variables"] == 1889
    assert broad["exact_binary_system"]["rank"] == 1885
    assert broad["exact_binary_system"]["augmented_rank"] == 1886
    assert broad["status"] == "INCONSISTENT_IN_THIS_JOINT_LABELLED_SLICE"
    assert all(part["status"] == "INCONSISTENT_IN_THIS_JOINT_LABELLED_SLICE"
               for part in report["individually_exact_border_branches"])

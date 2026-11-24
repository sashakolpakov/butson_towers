import pytest

from experiments.scalar_lift.cubic_core_support_lift import template
from experiments.scalar_lift.nondegenerate_cubic_projectors import retained_templates
from experiments.scalar_lift.nondegenerate_cubic_second_lift import solve_template


@pytest.mark.parametrize("order,phase", ((1499, 24), (2999, 48)))
def test_known_ten_row_family_passes_the_same_complete_second_bit_operator(order, phase):
    control = {**template(), "template_sha256": "known-ten-row-positive-control"}
    report, assignment = solve_template(control, order, phase)
    assert assignment is not None
    assert report["status"] == "SOLVED"
    assert report["all_four_representative_Gram_coefficients_verified_by_inverse_roots"]
    assert report["adjacent_channel_collapse_ruled_out_for_this_candidate"]


@pytest.mark.parametrize("order,phase", ((1499, 24), (2999, 48)))
def test_twelve_row_first_lift_can_fail_this_specific_second_bit_slice(order, phase):
    report, assignment = solve_template(retained_templates()[0], order, phase)
    assert assignment is None
    assert report["zeroth_and_first_Gram_verified_directly"]
    assert report["projector_identity_not_required"]
    assert report["exact_binary_system"]["rank"] == 1828
    assert report["exact_binary_system"]["augmented_rank"] == 1829
    assert report["status"] == "INCONSISTENT_IN_THIS_LABELLED_SECOND_BIT_SLICE"
    assert not report["unrestricted_template_or_source_nonexistence_claimed"]

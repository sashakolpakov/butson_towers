import pytest

from experiments.scalar_lift.cubic_core_support_lift import template
from experiments.scalar_lift.mixed_cubic_ramification import CONJUGATE
from experiments.scalar_lift.mixed_cubic_second_ramification import (
    derivative,
    encode_blocks,
    product,
)
from experiments.scalar_lift.nondegenerate_cubic_joint_boundary_lift import (
    build_system,
    solve_template,
)
from experiments.scalar_lift.nondegenerate_cubic_projectors import retained_templates


@pytest.mark.parametrize("order,phase", ((1499, 24), (2999, 48)))
def test_known_ten_row_positive_family_passes_joint_active_and_inactive_pair_model(order, phase):
    control = {**template(), "template_sha256": "known-ten-row-positive-control"}
    report = solve_template(control, order, phase)
    for part in report["slices"]:
        assert part["status"] == "SOLVED"
        assert part["exact_binary_system"]["nullity"] == 2
        assert part["all_four_representative_Gram_coefficients_verified_by_inverse_roots"]
        assert part["all_kernel_basis_neighbors_verified_by_inverse_root_Gram"]


def test_joint_model_contains_the_previous_fixed_boundary_slice_exactly():
    data = build_system(retained_templates()[0], 1499, 24)
    h0, h1 = data["H0"], data["H1"]
    fixed = product(h1, CONJUGATE[h1].T) ^ product(h0, CONJUGATE[h1].T)
    old_target = encode_blocks(fixed ^ derivative(h0, data["fixed_H2"]),
                               fixed ^ derivative(h1, data["fixed_H2"]) ^ derivative(h0, data["fixed_H3"]))
    new_target = data["target"]
    for descriptor, column in zip(data["variable_descriptors"], data["columns"], strict=True):
        if descriptor[0] == "boundary_pair" and descriptor[1] == 1 and descriptor[3]:
            new_target ^= column
    assert old_target == new_target


@pytest.mark.parametrize("order,phase", ((1499, 24), (2999, 48)))
def test_first_width_twelve_candidate_fails_only_the_reported_joint_slice(order, phase):
    report = solve_template(retained_templates()[0], order, phase)
    assert report["projector_identity_not_required"]
    assert report["filler_pair_phases_fixed"]
    assert not report["unrestricted_template_or_source_nonexistence_claimed"]
    for part in report["slices"]:
        assert part["status"] == "INCONSISTENT_IN_THIS_JOINT_LABELLED_SLICE"
        assert part["exact_binary_system"]["augmented_rank"] == part["exact_binary_system"]["rank"] + 1
    assert report["slices"][1]["variables"] == 1840
    assert report["slices"][1]["boundary_variables"] == 12

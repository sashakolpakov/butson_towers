import pytest

from experiments.scalar_lift.mixed_cubic_boundaries import (
    certificate,
    dense_certificate,
    metric_control,
    pair_action_audit,
)


@pytest.mark.parametrize("phase", (24, 48, 96))
def test_pair_mixing_is_exactly_unitary_and_root_valued_on_all_sign_pairs(phase):
    report = pair_action_audit(phase)
    assert report["twice_pair_action_Gram_equals_4I"]
    assert report["all_four_sign_pair_images_verified"]


def test_actual_large_mixed_boundary_and_full1499_residue():
    report = certificate()
    assert report["route"]["active_rows"] == 372
    assert report["route"]["adjoint_boundary_U_equals_V_is_false"]
    assert report["route"]["cross_ratio"]["primitive_full_phase_order"] == 96
    assert report["full_source_residue"]["complete_F4_source_row_Gram_equals_I"]
    assert report["full_source_residue"]["all_row_pairs_including_diagonal_checked"] == 1499 * 1500 // 2
    assert not report["characteristic_zero_large_source_constructed"]


def test_metric_factorization_uses_the_correct_transposes_and_conjugates():
    report = metric_control()
    assert report["E_conjugate_U_equals_minus_U_Dstar_verified"]
    assert report["E_Estar_equals_11I_minus_2U_Ustar_verified"]
    assert not report["F_or_BH11_constructed"]


def test_actual2999_frames_dense_residues_and_characteristic_zero_metric_seed():
    report = dense_certificate()
    routes = report["routes"]
    assert [route["source_order_sought"] for route in routes] == [1499, 2999, 2999]
    assert [route["conditional_real_lift_order"] for route in routes[1:]] == [47984, 95968]
    for route in routes[1:]:
        assert route["active_rows"] == 748
        assert route["even_fixed_vector"] == [1, 1, 0, 1, 1, 0, 0]
        assert route["border_column_permutation"] == [1, 2, 0, 4, 5, 3, 6]
        assert route["Ustar_U_and_Vstar_V_equal_1496I"]
        assert route["Dstar_D_and_D_Dstar_equal_7I"]
        assert route["cross_ratio"]["primitive_full_phase_order"] == route["phase"]
    for residue in report["dense_full_source_residues"]:
        order = residue["source_residue_order"]
        assert residue["all_row_pairs_including_diagonal_checked"] == order * (order + 1) // 2
        assert residue["complete_F4_source_row_Gram_equals_I"]
        assert residue["core_cubic_support_is_connected"]
    assert report["historical_v1"]["old_PX_zero_PY_identity_support_cannot_lift_at_any_dyadic_level"]
    assert all(control["z_conjugate_z_equals_2999"] for control in report["exact_norm_2999_controls"])
    assert not report["characteristic_zero_large_source_constructed"]
    assert not report["new_real_Hadamard_constructed"]

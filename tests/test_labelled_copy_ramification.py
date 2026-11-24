from experiments.scalar_lift.labelled_copy_ramification import certificate


def test_two_labelled_copy_types_have_no_new_second_bit_direction():
    report = certificate()
    assert len(report["routes"]) == 3
    for route in report["routes"]:
        assert route["number_second_bit_variables"] == 1284
        assert route["exact_binary_system"] == {"rank": 1284, "augmented_rank": 1284, "nullity": 0}
        assert route["full_affine_kernel_basis_nonzero_variable_indices"] == []
        assert route["all_four_representative_root_Gram_coefficients_verified"]
        assert route["entire_affine_space_has_explicit_adjacent_channel_gauge"]
        assert route["no_new_anti_collapse_candidate_in_this_labelled_slice"]
        assert route["remaining_dyadic_channels_mod4"] == [0, 1]
    assert not report["native_probe_run"]
    assert not report["higher_bit_solving_performed"]
    assert not report["new_large_residue_or_characteristic_zero_source_claimed"]

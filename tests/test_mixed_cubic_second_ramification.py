from experiments.scalar_lift.mixed_cubic_second_ramification import certificate


def test_actual_high_order_roots_pass_all_four_ramification_gram_coefficients():
    report = certificate()
    assert [(route["source_order"], route["phase"]) for route in report["routes"]] == [
        (1499, 96), (2999, 48), (2999, 96)]
    for route in report["routes"]:
        expected_bit = int(route["source_order"] == 2999)
        assert route["core_X_and_Y_second_bits_are_constant"] == expected_bit
        model = route["exact_linear_model"]
        assert model["status"] == "SOLVED"
        assert model["number_second_bit_variables"] == 444
        assert model["exact_binary_system"]["rank"] == 444
        assert model["all_four_representative_Gram_coefficients_verified_independently"]
        actual = route["actual_large_verification"]
        assert actual["all_actual_columns_used_in_Gram"] == route["source_order"]
        assert actual["all_four_actual_Gram_coefficients_verified"]
        assert actual["all_source_row_pairs_covered_by_exact_permutation_symmetry"] == (
            route["source_order"] * (route["source_order"] + 1) // 2)
        assert not route["characteristic_zero_source_constructed"]
    assert not report["native_probe_run"]
    assert not report["new_characteristic_zero_Butson_or_real_Hadamard_claimed"]

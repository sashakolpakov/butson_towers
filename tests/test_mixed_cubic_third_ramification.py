from experiments.scalar_lift.mixed_cubic_third_ramification import certificate


def test_actual_high_order_roots_pass_all_eight_ramification_gram_coefficients():
    report = certificate()
    assert [(route["source_order"], route["phase"]) for route in report["routes"]] == [
        (1499, 96), (2999, 48), (2999, 96)]
    assert [route["core_third_bit_values"] for route in report["routes"]] == [[0], [0], [1]]
    for route in report["routes"]:
        model = route["exact_linear_model"]
        assert model["status"] == "SOLVED"
        assert model["number_third_bit_variables"] == 444
        assert model["exact_binary_system"]["rank"] == 444
        assert model["all_eight_representative_Gram_coefficients_verified"]
        actual = route["actual_large_verification"]
        assert actual["all_actual_columns_used_in_Gram"] == route["source_order"]
        assert actual["all_eight_actual_Gram_coefficients_verified"]
        assert actual["actual_boundary_and_border_phase_exponents_retained_mod8"]
        assert actual["all_source_row_pairs_covered_by_exact_permutation_symmetry"] == (
            route["source_order"] * (route["source_order"] + 1) // 2)
        assert route["unitary_root_residue_modulo_pi_eighth"]
        assert not route["characteristic_zero_source_constructed"]
    assert not report["native_probe_run"]
    assert not report["new_characteristic_zero_Butson_or_real_Hadamard_claimed"]

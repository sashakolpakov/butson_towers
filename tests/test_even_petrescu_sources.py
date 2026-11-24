from experiments.scalar_lift.even_petrescu_sources import (
    build_report,
    corner_formula_controls,
    symmetric_three_phase_corner_audit,
    two_phase_corner_audit,
)


def test_known_inputs_and_exact_norm_witnesses_without_claiming_completion():
    report = build_report()
    assert report["Paley_base_complete_check"]["order"] == 500
    for record in report["known_input_certificates"]:
        assert record["T_row_Gram_inherited_from_full_verified_input"]
        assert record["T_column_Gram_all_off_diagonal_pairs_checked"] > 0
    for route in report["routes"]:
        assert route["corner_row_sum_witness"]["full_row_sum_and_norm_verified"]
        assert route["corner_determinant_norm_witness"]["factor_norms_verified_exactly"]
        assert not route["previous_odd_border_obstruction_applies"]
        assert not route["large_source_constructed"]


def test_two_phase_corner_exclusion_includes_arbitrary_unit_phases():
    audits = [two_phase_corner_audit(order) for order in (1000, 2000)]
    assert [audit["maximum_k_minus_lambda"] for audit in audits] == [209, 1]
    for audit in audits:
        assert audit["all_two_global_unit_phase_corners_excluded"]
        assert audit["phases_need_not_be_roots_of_unity_for_this_exclusion"]
        assert not audit["corners_with_three_or_more_phases_excluded"]


def test_symmetric_corner_census_and_independent_gram_expansion():
    audits = [symmetric_three_phase_corner_audit(t, q) for t, q
              in ((1000, 24), (1000, 48), (2000, 48))]
    assert [audit["integral_alpha_beta_pair_count"] for audit in audits] == [14, 14, 16]
    for audit in audits:
        assert not audit["valid_strongly_regular_parameter_sets"]
        assert not audit["all_three_phase_corners_excluded"]
    controls = corner_formula_controls()
    assert controls["regular_graph_expansion"]["complete_power_basis_Gram_expansion_verified"]
    assert controls["regular_graph_expansion"]["nonzero_Petrescu_Gram_defect_verified"]
    assert controls["published_Petrescu_corner"]["both_Grams_equal_5I_plus_2J"]

from experiments.scalar_lift.fifth_gauged_seed_escape import build_report


def test_actual_small_butson_gauges_surject_onto_the_compression_codomain():
    report = build_report()
    for family in report["explicit_passing_gauge_family"]:
        assert family["M_exact_BH5_Gram_verified"]
        assert family["M_avoids_fifth_root_one"]
        assert family["compression_F2_rank"] == 100
        assert family["full_compression_equation_verified"]
        assert not family["witness_realized_as_relative_trace_of_unitary_E"]
        symmetry = family["symmetric_completion_preflight"]
        assert symmetry["symmetric_E_completion_excluded"]
        assert symmetry["dimension_after_adjoining_T_D_plus_D_star"] == 16
        assert symmetry["dual_vanishes_on_every_T_row"]
        assert not symmetry["nonsymmetric_E_completion_excluded"]
    assert report["pure_column_permutation_control"]["permutations_checked"] == 120
    assert not report["simultaneous_unitary_residue_completion_constructed"]

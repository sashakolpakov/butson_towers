from experiments.scalar_lift.even_petrescu_conference_family import build_report


def test_exact_complete_conference_boundary_sources():
    report = build_report()
    for formal in report["formal_full_Gram_certificates"]:
        assert formal["complete_Gram_equals_twice_weight_plus_one_I"]
    for source in report["complete_known_sources"]:
        assert source["normalized_first_row_and_column_all_one"]
        assert source["cross_ratio_order"] == source["phase"]
        assert sum(source["root_exponent_counts"]) == source["order"] ** 2
        assert not source["higher_Petrescu_source_completed"]
    assert all(control["complete_entrywise_cyclotomic_Gram_checked"]
               for control in report["independent_small_direct_Gram_controls"])


def test_partial_corner_has_only_the_claimed_sector_and_lattice():
    corner = build_report()["partial_corner"]
    assert corner["difference_sector_Gram"] == "998 I_500"
    assert corner["sum_sector_Gram"] == "1002 I_500"
    assert corner["all_blocks_preserve_old_Paley_lattice_modulo_125"]
    assert not corner["desired_Petrescu_corner_constructed"]


def test_variable_boundary_has_a_new_coupled_graph_lattice():
    for kernel in build_report()["new_boundary_kernel_certificates"]:
        assert kernel["kernel_rank"] == 500
        assert kernel["each_Paley_summand_contributes_rank"] == 250
        assert kernel["different_from_two_copies_of_old_Paley_kernel"]
        assert kernel["graph_is_totally_isotropic_modulo_125"]
        assert all(channel["rank_over_cyclotomic_ring_mod_125"] == 1
                   and channel["determinant_zero_exactly"] for channel in kernel["channels"])

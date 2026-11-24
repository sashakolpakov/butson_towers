from experiments.scalar_lift.even_petrescu_lattice import build_report


def test_full_paley_lattice_and_dense_positive_negative_controls():
    report = build_report()
    assert report["P_squared_equals_P"]
    assert report["Pt_equals_I_minus_P"]
    assert report["rank_P_modulo_5_by_full_elimination"] == 250
    controls = report["dense_controls"]
    assert controls["both_lattice_projections_zero"]
    assert controls["nonzero_negative_projection_entries"] == 500 * 500
    assert not controls["Petrescu_corner_Gram_claimed"]


def test_walsh_factorization_both_directions_and_scoped_support_bound():
    report = build_report()
    for control in report["direct_small_Walsh_factor_controls"]:
        assert control["forward_factorization_verified"]
        assert control["inverse_recovers_all_block_projections"]
        assert control["both_zero_congruences_hold"] != control["one_entry_flipped"]
    for audit in report["root_support_audits"]:
        assert audit["minimum_support_of_nonzero_single_row_change"] == 63
        assert audit["minimum_support_of_nonzero_single_row_sign_flips"] == 125
        assert not audit["multirow_trade_lower_bound_claimed"]

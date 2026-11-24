from experiments.scalar_lift.fifth_source_boundary_audit import build_report


def test_initial_obstruction_and_rotated_exact_row_preflight():
    report = build_report()
    assert report["initial_seed"]["excluded_by_trace"]
    rotated = report["one_D_block_rotated_by_zeta5_cubed"]
    assert rotated["F2_image_rank"] == 20
    assert rotated["minimum_number_of_attainable_scalar_norms"] == 3
    assert rotated["all_11_representative_rows_verified"]
    assert not rotated["simultaneous_row_orthogonality_verified"]
    assert rotated["later_joint_trace_unitarity_obstruction_excludes_this_seed"]


def test_odd_unitary_counterexample_limits_obstruction_scope():
    report = build_report()
    control = report["odd_unitary_counterexample"]
    assert control["full_unitary_Grams_verified"]
    assert control["all_entries_in_exact_sum_alphabet"]
    assert not control["compatibility_with_prescribed_T_claimed"]
    assert not report["all_fifth_root_bordered_sources_excluded"]

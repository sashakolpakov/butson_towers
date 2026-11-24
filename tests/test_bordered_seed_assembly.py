import numpy as np

from experiments.scalar_lift.bordered_seed_assembly import assemble_exponents, build_report
from hadamard.cgw import is_cgw_exponent_matrix


def test_verified_seed_inputs_have_required_orders_and_no_completion_claim():
    report = build_report()
    assert [item["all_row_pairs_checked"] for item in report["full_real_seed_checks"]] == [279378, 1118260]
    assert [route["real_target_order"] for route in report["routes"]] == [47968, 47984, 95968]
    assert not report["new_Butson_sources_constructed"]


def test_bordered_assembly_on_an_exact_full_alphabet_control():
    # Known cubic H9 from s=4,t=1. It tests the block/adjoint placement;
    # it deliberately uses all three cosets, so is not a partial-lift source.
    skew = np.array([[0, 1, 1, 1], [-1, 0, 1, -1],
                     [-1, -1, 0, 1], [-1, 1, -1, 0]])
    x = (skew == 1).astype(np.int64)
    y = (skew == -1).astype(np.int64)
    np.fill_diagonal(x, 2)
    np.fill_diagonal(y, 2)
    source = assemble_exponents(x, y, np.zeros((4, 1), dtype=np.int64), [[1]], 3)
    assert is_cgw_exponent_matrix(source.tolist(), 3, 9)
    assert set(source.ravel()) == {0, 1, 2}

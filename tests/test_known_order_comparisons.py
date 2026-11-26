import json

import numpy as np
import pytest

from experiments.scalar_lift import known_order_comparisons as comparison


def test_known_order_comparison_replays_printed_pair_certificate():
    path = comparison.ROOT / "certificates/scalar-lift/known-order-paley-comparison.json"
    report = comparison.audit()
    assert report == json.loads(path.read_text())
    assert report["order"] == 885924 and report["order_status"] == "known"
    assert report["printed_pair_source"]["rank"] == 195
    assert report["printed_pair_source"]["matrix_mod3_uint8_sha256"] == (
        "f19f073b02e065754dcfb08a9f44b5af7726ab9f7dac6d041617cf40ad6fe39a"
    )
    assert report["endpoint_rank"] == 246090
    assert report["Paley_I"]["rank_mod3"] == 442962
    assert [(c["control_endpoint_order"], c["direct_endpoint_rank"])
            for c in report["small_literal_frame_controls"]] == [
        (52, 52), (156, 78), (468, 178), (84, 42), (364, 364)
    ]
    assert report["largest_matrix_dimension_used"] == 702
    assert report["new_Hadamard_order_claimed"] is False
    assert report["global_new_equivalence_class_claimed"] is False
    assert report["rank_comparison_kernel_checked_in_Lean"] is False


def test_rank_certificate_proves_both_bounds_and_detects_corruption():
    matrix = np.array([[1, 2, 0, 1], [2, 1, 0, 2], [0, 0, 1, 1]])
    u, v, left, pivots = comparison.rank_factorization(matrix)
    assert comparison.verify_rank_factorization(matrix, u, v, left, pivots) == 2
    assert pivots == [0, 2]
    assert np.array_equal(u, matrix[:, pivots])
    for name, value in (("U", u), ("V", v), ("L", left)):
        bad = value.copy()
        bad.flat[0] = (bad.flat[0]+1) % 3
        args = [matrix, u, v, left, pivots]
        args[{"U": 1, "V": 2, "L": 3}[name]] = bad
        with pytest.raises(ValueError):
            comparison.verify_rank_factorization(*args)
    with pytest.raises(ValueError, match="integer residues"):
        comparison.verify_rank_factorization(matrix, u.astype(float), v, left, pivots)
    # The factorization and left inverse alone are insufficient unless the
    # selected-column binding is checked too.
    with pytest.raises(ValueError, match="selected-column identity"):
        comparison.verify_rank_factorization(matrix, 2*u % 3, 2*v % 3, 2*left % 3, pivots)


def test_sixth_root_action_and_printed_pair_gauge_agree_entrywise():
    for e in range(6):
        for f in range(6):
            left = comparison.dual_number_matrix(np.array([[e]]))
            right = comparison.dual_number_matrix(np.array([[f]]))
            product = comparison.dual_number_matrix(np.array([[(e+f) % 6]]))
            assert np.array_equal(left @ right % 3, product)
    e = np.array([[0, 1], [3, 4]])
    printed, a, b = comparison.printed_pair_matrix(e)
    identity, zero = np.eye(2, dtype=np.int16), np.zeros((2, 2), dtype=np.int16)
    f = np.block([[identity, identity], [zero, identity]])
    c = np.block([[2*identity, zero], [zero, identity]])
    assert set(np.unique(a)) == {-1, 1} and set(np.unique(b)) == {-1, 1}
    assert np.array_equal(f @ c @ printed @ c % 3, comparison.dual_number_matrix(e))
    with pytest.raises(ValueError, match="admitted frame palette"):
        comparison.printed_pair_matrix(np.array([[2]]))


def test_numpy_only_rank_fallback(monkeypatch):
    monkeypatch.setattr(comparison, "nmod_mat", None)
    assert comparison.checked_rank(np.array([[1, 2], [2, 1]])) == 1
    assert comparison.checked_rank(np.zeros((2, 3), dtype=np.int16)) == 0

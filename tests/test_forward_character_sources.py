from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from experiments.scalar_lift import forward_character_sources as forward


def test_completed_large_source_Grams_and_primitive_phases():
    for prime, order, phase in ((73, 5329, 24), (241, 58081, 120)):
        source = forward.character_source(prime)
        assert source.order == order and source.phase == phase
        assert source.coefficients.shape == (prime, prime)
        assert source.coefficients.nbytes == prime * prime * 8
        proof = forward.formal_gram_check(source)
        assert proof["all_matrix_Gram_entries_proved_by_group_ring_identity"]
        assert proof["number_of_Gram_entries_covered"] == order**2
        assert forward.primitive_cross_ratio(source)["cross_ratio_order"] == phase
        assert len(forward.direct_correlation_check(source)["literal_checks"]) == 6
        row, column = (prime - 1, 2), (1, prime - 1)
        assert source.entry(row, column) == source.coefficients[2, prime - 3]


def test_gamma_freedom_and_raw_signed_triangle_obstruction():
    # Changing gamma is an actual completed-source parameter, not a guessed
    # modification of its Gram.  Test the small full correlation independently.
    for gamma in range(12):
        source = forward.character_source(5, 4, gamma=gamma)
        result = forward.direct_correlation_check(source, all_shifts=True)
        assert len(result["literal_checks"]) == 25
        assert forward.uniform_sign_obstruction(source)[
            "uniform_Q_linear_final_sign_map_excluded_for_raw_source"]
    with pytest.raises(ValueError, match="nontrivial"):
        forward.character_source(73, 5)
    with pytest.raises(ValueError, match="prime"):
        forward.character_source(75, 2)


def test_tiny_gauge_obstructions_and_nonvacuous_positive_controls():
    records = forward.cubic_gauge_census()
    assert [row["feasible_row_gauges"] for row in records] == [0, 0, 2, 0, 2, 0]
    assert all(row["row_gauges_exhausted"] == 2187 for row in records)
    assert records[2]["witness"] and records[4]["witness"]
    for prime in (73, 241):
        for gamma in range(3):
            source = forward.character_source(prime, gamma=gamma)
            gate = forward.literal_gauge_template_check(source)
            assert gate["arbitrary_root_row_column_gauge_to_current_two_cubic_cosets"] == "EXCLUDED"
            assert gate["other_partial_frames_or_nonuniform_descent_maps_excluded"] is False
            # All exponents are genuine roots, including the fixed axis values.
            assert np.all((source.coefficients >= 0) & (source.coefficients < source.phase))


def test_saved_forward_source_certificate():
    path = Path(__file__).resolve().parents[1] / "certificates/scalar-lift/forward-character-sources.json"
    report = json.loads(json.dumps(forward.certificate()))
    assert json.loads(path.read_text()) == report
    assert report["small_control"]["full_dense_order25_Gram_verified"]
    assert report["new_real_Hadamard_order_claimed"] is False

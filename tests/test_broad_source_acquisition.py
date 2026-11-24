from itertools import product

import numpy as np

from experiments.scalar_lift import broad_source_acquisition as acquisition
from experiments.scalar_lift.scalar_restriction_obstructions import alphabet_obstruction


def test_antipodal_palette_criterion_matches_all_exact_root_relation_tests():
    for phase in (6, 12):
        for mask in range(1, 1 << phase):
            used = [e for e in range(phase) if (mask >> e) & 1]
            antipodal = 0
            for exponent in used:
                antipodal |= 1 << (exponent % (phase//2))
            excluded = ((antipodal & 7) == 7 if phase == 6 else
                        (antipodal & 21) == 21 or (antipodal & 42) == 42)
            gate = alphabet_obstruction(phase, used)
            assert gate["complete_for_odd_integer_relations"]
            assert gate["excluded"] == excluded


def _brute_sixth_gauge_exists(matrix):
    n = len(matrix)
    columns = np.array(list(product(range(3), repeat=n)))
    for tail in product(range(3), repeat=n-1):
        rows = np.array((0, *tail))
        values = (matrix[None, :, :]+rows[None, :, None]+columns[:, None, :]) % 3
        masks = np.bitwise_or.reduce((1 << values).reshape(len(columns), -1), axis=1)
        if np.any(masks != 7):
            return True
    return False


def test_gauge_csp_matches_unpruned_exhaustive_gauges_with_acceptances_and_rejections():
    rng = np.random.default_rng(20261123)
    outcomes = set()
    for _ in range(10):
        matrix = rng.integers(6, size=(5, 5))
        actual = acquisition.root_gauge(matrix, 6)["exists"]
        assert actual == _brute_sixth_gauge_exists(matrix)
        outcomes.add(actual)
    assert outcomes == {False, True}


def test_published_dataset_counts_and_complete_exact_gauge_inventory():
    report = acquisition.certificate()
    assert report["source_count"] == 524
    assert report["retained_fixed_representative_count"] == 33
    assert report["gauge_orbits_excluded_in_their_own_root_alphabet"] == 491
    assert report["raw_data_bytes"] == 84288
    for row in report["sources"]:
        assert row["proof"]["full_Gram_verified"]
        assert row["proof"]["acquisition"]["raw_source_full_Gram"]["full_Gram_verified"]
        assert row["proof"]["monomial_gauge_audit"]["complete_over_all_row_column_root_gauges"]
        assert not row["proof"]["acquisition"]["tree_edge_added"]
        assert 0 in row["used_alphabet"]


def test_genuine_twelfth_root_literal_survives_with_exact_bound_proof():
    source = acquisition.load_source("aalto-9-12-064")
    assert not source.proof["alphabet_gate"]["excluded"]
    assert any(exponent % 2 for exponent in source.used_alphabet)
    assert source.proof["exponent_sha256"] == source.source_hash
    assert source.entry(0, 0) == 0

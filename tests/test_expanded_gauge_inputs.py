import numpy as np

from experiments.scalar_lift import expanded_gauge_inputs as expanded
from experiments.scalar_lift.broad_source_acquisition import dataset, load_source
from experiments.scalar_lift.forward_butson_descent import verify_source


def test_larger_alphabet_rescues_a_known_excluded_source_with_an_explicit_gauge():
    previous = load_source("aalto-9-6-001")
    assert previous.proof["alphabet_gate"]["excluded"]
    assert not previous.proof["monomial_gauge_audit"]["exists"]
    source = expanded.load_source("expanded-9-001-q12")
    proof = source.proof
    assert not proof["alphabet_gate"]["excluded"]
    assert proof["rescued_by_larger_initial_root_gauge"]
    assert proof["effective_literal_root_order"] == 12
    assert proof["projective_root_order"] == 6
    gauge = proof["twelfth_root_gauge_audit"]
    raw = dataset(9, 6)[1]
    expected = (2*raw+np.array(gauge["row_shifts"])[:, None]
                +np.array(gauge["column_shifts"])[None, :]
                -proof["acquisition"]["global_anchor_subtracted_after_gauge"]) % 12
    assert np.array_equal(source.exponents, expected)
    assert verify_source(expected, 12)["exponent_sha256"] == source.source_hash


def test_all_declared_extended_gauges_finish_with_exact_counts_and_no_missing_inputs():
    report = expanded.certificate()
    assert report["original_published_input_count"] == 455
    assert report["originally_admissible_count"] == 14
    assert report["rescued_input_count"] == 80
    assert report["still_excluded_count"] == 361
    assert report["total_admissible_twelfth_root_representatives"] == 94
    assert report["unknown_cases"] == 0
    assert report["all_rescues_have_literal_entry_root_order12"]
    assert len(expanded.admissible_sources()) == 94
    for record in report["sources"]:
        assert record["proof"]["full_Gram_verified"]
        assert not record["proof"]["acquisition"]["scalar_restriction_tree_edge_added"]


def test_original_raw_input_is_unchanged_by_expanded_acquisition():
    raw = dataset(13, 6)[4].copy()
    expanded.load_source("expanded-13-004-q12")
    assert np.array_equal(raw, dataset(13, 6)[4])
    assert expanded.load_source("expanded-13-004-q12").proof["rescued_by_larger_initial_root_gauge"]

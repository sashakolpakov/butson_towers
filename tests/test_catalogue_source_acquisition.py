import numpy as np
import pytest

from experiments.scalar_lift.catalogue_source_acquisition import (
    SOURCE_IDS,
    block_width_gate,
    certificate,
    cubic_gauge,
    load_source,
    twelfth_root_provenance,
)
from experiments.scalar_lift.forward_butson_descent import verify_source


def test_all_fixed_inputs_have_exact_gram_and_pass_alphabet_gate():
    for source_id in SOURCE_IDS:
        source = load_source(source_id)
        assert verify_source(source.exponents, source.phase)["full_Gram_verified"]
        assert not source.proof["alphabet_gate"]["excluded"]
        assert all(source.entry(i, j) == source.exponents[i, j]
                   for i, j in np.ndindex(source.exponents.shape))
    with pytest.raises(KeyError):
        load_source("unverified")


def test_twelfth_root_input_has_published_provenance_and_primitive_cross_ratio():
    result = twelfth_root_provenance()
    assert result["exact_monomial_equivalence_verified"]
    assert result["primitive_cross_ratio"]["root_order"] == 12
    assert result["not_equivalent_to_a_sixth_root_matrix"]


def test_gauge_exclusions_and_frontier_routes_are_replayed():
    report = certificate()
    assert [r["two_cubic_coset_gauge"]["visited_row_prefixes"]
            for r in report["rejected_uniform_alphabet_candidates"]] == [6, 8]
    assert report["largest_matching_frontier_order"] == 10052
    routes = report["retained_inputs"][2]["frontier_factor_routes"]
    assert [r["required_total_width"] for r in routes] == [334, 566, 694, 718]
    assert not any(r["normalized_uniform_sign_block_gate"]["passes"] for r in routes)
    assert all(block_width_gate(d)["passes"] for d in (1, 2, 4, 668))
    assert not block_width_gate(6)["passes"]
    assert not report["new_Hadamard_order_constructed"]
    # The Fourier3 matrix really does admit a two-coset gauge.
    assert cubic_gauge([[0, 0, 0], [0, 1, 2], [0, 2, 1]])["exists"]

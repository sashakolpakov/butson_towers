import json
from pathlib import Path

import numpy as np
import pytest

from experiments.scalar_lift.ambient_field_candidates import load_ambient_source, period_field
from experiments.scalar_lift.mixed_basis_restriction import (
    MixedRestriction,
    independent_small_verification,
    load_saved_node,
    recipes,
)


def test_genuinely_mixed_basis_matches_independent_full_restriction_and_Gram():
    field = period_field(3, 31, (2, 3, 1))
    source = load_ambient_source("catalogue-p7-q12")
    base = MixedRestriction(field, source)
    mixed = MixedRestriction(field, source, recipes()[-1])
    assert mixed.literal_hash() != base.literal_hash()
    assert mixed.census() != base.census()
    assert mixed.replay_inverse_shears()
    proof = independent_small_verification(mixed)
    assert proof["complete_weighted_Gram_checked"]
    assert proof["complete_matrix_determinant_checked"]
    assert sum(mixed.census().values()) == 84**2


def test_inverse_shears_are_literal_duplicates_including_the_metric():
    field = period_field(3, 31, (2, 3, 1))
    source = load_ambient_source("catalogue-p7-q12")
    base = MixedRestriction(field, source)
    returned = MixedRestriction(field, source, [(0, 5, 2), (0, 5, -2)])
    assert returned.literal_hash() == base.literal_hash()
    assert np.array_equal(returned.metric, base.metric)
    with pytest.raises(ValueError, match="both tensor factors"):
        MixedRestriction(field, source, [(0, 1, 1)])
    with pytest.raises(ValueError, match="outside"):
        base.entry(base.order, 0)


def test_executed_large_batch_counts_actual_entries_and_preserves_odd_determinant():
    path = Path("certificates/scalar-lift/mixed-basis-restriction.json")
    report = json.loads(path.read_text())
    assert report["node_visits"] == 90
    assert report["distinct_nodes"] == 70
    assert report["literal_duplicates_pruned"] == 20
    assert report["pending_recipes"] == 0
    assert [case["matrix_order"] for case in report["cases"]] == [84, 9308]
    for case in report["cases"]:
        assert case["determinant_is_odd"]
        assert not case["unequal_domain_codomain_bases_excluded"]
        for record in case["nodes"]:
            if "literal_duplicate_of" in record:
                assert (case["nodes"][record["literal_duplicate_of"]]["action_metric_sha256"]
                        == record["action_metric_sha256"])
                continue
            assert sum(record["complete_entry_counts"].values()) == case["matrix_order"]**2
            assert not record["all_entries_signs"]
    large = load_saved_node(path, case=1)
    assert large.order == 9308
    assert large.census() == report["cases"][1]["nodes"][-1]["complete_entry_counts"]

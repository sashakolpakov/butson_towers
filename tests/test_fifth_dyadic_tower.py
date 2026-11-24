from math import gcd

import numpy as np
import pytest

from experiments.scalar_lift.fifth_dyadic_tower import (
    action_power,
    build_report,
    lift_source,
    phase_blocks,
    tower_level,
)


def test_complete_finite_audit_and_source_to_sign_controls():
    report = build_report()
    audit = report["complex_structure_audit"]
    assert audit["total_complex_structures_tested"] == 1680
    assert audit["commuting_count"] == 6
    viable = {(r["source_order"], r["source_phase"]) for r in report["high_source_route_audit"]
              if r["passes_checked_necessary_source_filters"]}
    assert viable == {(1499, 40), (1499, 80), (2999, 80)}
    for level in report["levels"]:
        phase = level["phase"]
        assert level["dimension"] == level["cyclotomic_degree"]
        assert level["allowed_exponents_sparse_frame"] == [e for e in range(phase) if e % 5]
        assert level["control"]["both_exact_integer_Grams_verified"]
        assert any(gcd(e, phase) == 1 for row in level["control"]["source_exponents"] for e in row)


def test_tower_really_extends_the_dyadic_scalar():
    for phase in (40, 80):
        upper = tower_level(phase)
        lower = tower_level(phase // 2)
        assert np.array_equal(upper["dyadic_action"] @ upper["dyadic_action"],
                              np.kron(np.eye(2, dtype=np.int64), lower["dyadic_action"]))
        numerator, denominator = action_power(upper, 5)
        assert np.array_equal(numerator, denominator * upper["dyadic_action"])


def test_sparse_frame_is_never_mistaken_for_an_admitted_sign_block():
    for phase in (20, 40, 80):
        assert not np.all(np.abs(phase_blocks(phase)[0]) == 1)
        with pytest.raises(ValueError, match="excluded fifth-root"):
            lift_source([[0]], phase)
        with pytest.raises(ValueError, match="exact integers"):
            lift_source([[1.0]], phase)
        with pytest.raises(ValueError, match="Gram"):
            lift_source([[1, 1], [1, 1]], phase)

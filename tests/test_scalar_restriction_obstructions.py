import json
from pathlib import Path

import pytest

from experiments.scalar_lift.scalar_restriction_obstructions import (
    alphabet_obstruction,
    certificate,
    registry_certificate,
    verify_relation,
)


def test_fixed_source_certificate_and_literal_char_occurrences():
    report = certificate()
    saved = Path(__file__).resolve().parents[1] / "certificates/scalar-lift/scalar-restriction-obstructions.json"
    assert report == json.loads(saved.read_text())
    for row in report["sources"][:2]:
        assert row["alphabet_gate"]["excluded"]
        assert all(x["row"] == 0 for x in row["actual_root_occurrences"])
        assert len(row["actual_root_occurrences"]) == 3


def test_signed_relations_are_characteristic_zero_and_not_mod_two():
    gate = alphabet_obstruction(24, [0, 4, 20])
    assert gate["excluded"]
    assert gate["relation"]["coefficient_sum"] % 2
    assert set(gate["relation"]["zero_integral_vector"]) == {0}
    # A nonzero even integral value is not a characteristic-zero relation.
    with pytest.raises(ValueError, match="not an exact odd integer relation"):
        verify_relation(8, [0], [{"exponent": 0, "coefficient": 2}])
    with pytest.raises(ValueError, match="does not occur"):
        verify_relation(24, [0, 4], gate["relation"]["terms"])


def test_complete_fifth_and_high_dyadic_phase_checks():
    for phase in (5, 10, 20, 40, 80, 120):
        used = [(phase // 5) * j for j in range(5)]
        gate = alphabet_obstruction(phase, used)
        assert gate["excluded"]
        assert gate["relation"]["signed_weight"] == 5
        assert gate["complete_for_odd_integer_relations"]
    for phase in (8, 16, 32, 64):
        gate = alphabet_obstruction(phase, range(phase))
        assert gate["status"] == "SKIP"
        assert gate["complete_for_odd_integer_relations"]
        assert not gate["sign_flattenability_proved"]


def test_unchecked_odd_prime_and_empty_alphabet_do_not_get_excluded():
    # The genuine seven-term relation is outside the reused q30 catalogue.
    gate = alphabet_obstruction(7, range(7))
    assert gate["status"] == "SKIP"
    assert not gate["complete_for_odd_integer_relations"]
    assert alphabet_obstruction(120, [])["status"] == "SKIP"


def test_entire_fixed_registry_audit_replays_without_constructive_claims():
    report = registry_certificate()
    saved = Path(__file__).resolve().parents[1] / "certificates/scalar-lift/scalar-restriction-all-sources.json"
    assert report == json.loads(saved.read_text())
    assert report["source_count"] == 20
    assert report["counts"] == {"excluded": 2, "retained": 18}
    assert report["excluded_uniform_scalar_branches"] == ["character-73", "character-241"]
    assert report["retained_does_not_assert_sign_flattenability"]

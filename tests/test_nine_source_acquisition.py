import json
from pathlib import Path

import numpy as np
import pytest

from experiments.scalar_lift.nine_source_acquisition import (
    SOURCE_HASH,
    certificate,
    load_source,
    period_preflight,
)
from hadamard.cgw import is_cgw_exponent_matrix


def test_saved_literal_input_has_independent_exact_gram_and_is_read_only():
    source = load_source("saved-unreal9-q6")
    assert (source.order, source.phase) == (9, 6)
    assert source.source_hash == SOURCE_HASH
    assert source.used_alphabet == (1, 2, 4, 5)
    assert not source.exponents.flags.writeable
    assert is_cgw_exponent_matrix(source.exponents.tolist(), 6, 9)
    assert source.proof["odd_integer_relation_gate"]["complete_for_odd_integer_relations"]
    assert not source.proof["odd_integer_relation_gate"]["excluded"]
    assert all(source.entry(i, j) == source.exponents[i, j] for i in range(9) for j in range(9))
    with pytest.raises(ValueError, match="unknown saved"):
        load_source("unavailable")


def test_source_acquisition_never_executes_tripling(monkeypatch):
    import experiments.scalar_lift.inverse_blocks as inverse

    def prohibited(*args, **kwargs):
        raise AssertionError("a source construction was called during literal acquisition")

    monkeypatch.setattr(inverse, "triple_unreal_source", prohibited)
    load_source.cache_clear()
    source = load_source()
    assert source.proof["source_construction_performed"] is False


def test_certificate_replays_and_period_candidates_are_exact_arithmetic_only():
    report = certificate()
    saved = Path(__file__).resolve().parents[1] / "certificates/scalar-lift/nine-source-acquisition.json"
    assert report == json.loads(saved.read_text())
    assert report["projective_field_degree"] == 2
    assert report["phase_field_degree_divides764"]
    assert not report["real_Hadamard_of_order6876_constructed"]
    for row in period_preflight():
        ell, h = row["prime_conductor"], row["subgroup_order"]
        a, b, t = row["norm_witness_a_b_t"]
        assert ell == 191*h + 1
        assert a*a + 191*b*b == ell*t*t
        assert (a % 2, b % 2, t % 2) == (1, 0, 1)
        assert pow(2, h, ell) == 1
        assert row["rational_source_matrix_order"] == 6876
        assert not row["sign_flattenability_proved"]
    assert np.unique(load_source().exponents).tolist() == [1, 2, 4, 5]

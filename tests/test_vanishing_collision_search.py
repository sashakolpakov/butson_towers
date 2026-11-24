"""Exactness, completeness on a small control, and checkpoint integrity."""

from itertools import combinations_with_replacement

import pytest

from experiments.scalar_lift.vanishing_collision_search import (
    CollisionSearch,
    Fingerprint,
    canonical_relation,
    exact_sum,
    integral_roots,
    minimality_certificate,
    read_status,
)
from hadamard.coefficients import CyclotomicField


def test_integer_basis_and_modular_false_positive():
    for q in (6, 30, 120, 210):
        field = CyclotomicField(q)
        for exponent in (0, 1, q // 2, q - 1):
            assert integral_roots(q)[exponent] == tuple(
                int(c) for c in (field.zeta() ** exponent).coefficients)
    tiny_hash = Fingerprint(3, primes=(7, 7))
    assert tiny_hash.value((0,) * 7) == (0, 0)
    assert not minimality_certificate((0,) * 7, 3)["exact_zero"]


def test_minimality_mixed_relation_and_hidden_decomposition():
    atom = (0, 6, 12, 18, 19, 29)
    assert minimality_certificate(atom, 30)["minimal"]
    duplicated = minimality_certificate(atom + atom, 30)
    assert duplicated["exact_zero"] and not duplicated["minimal"]
    assert not any(exact_sum(duplicated["proper_subsum"], 30))


def test_small_complete_search_matches_independent_exhaustion(tmp_path):
    expected = set()
    for weight in range(2, 5):
        for terms in combinations_with_replacement(range(6), weight):
            if any(exact_sum(terms, 6)):
                continue
            # Deliberately use full subset enumeration instead of the worker's MITM.
            if all(any(exact_sum(tuple(e for i, e in enumerate(terms) if mask >> i & 1), 6))
                   for mask in range(1, (1 << weight) - 1)):
                expected.add(canonical_relation(terms, 6))
    search = CollisionSearch(tmp_path / "small.sqlite3", 6, 2)
    try:
        report = search.run(seconds=10)
        assert report["state"]["finished"]
        assert {tuple(row["exponents"]) for row in report["examples"]} == expected
        assert report["state"]["halves_visited"] == 27
    finally:
        search.close()


def test_resume_inside_bucket_and_single_writer(tmp_path):
    path = tmp_path / "resume.sqlite3"
    search = CollisionSearch(path, 6, 3)
    search.run(seconds=10, steps=12)
    assert search.state["pending"] is not None
    with pytest.raises(RuntimeError, match="another worker"):
        CollisionSearch(path, 6, 3)
    assert read_status(path)["state"] == search.state
    search.close()
    resumed = CollisionSearch(path, 6, 3)
    direct = CollisionSearch(tmp_path / "direct.sqlite3", 6, 3)
    try:
        resumed.run(seconds=10)
        direct.run(seconds=10)
        assert resumed.state == direct.state
        assert resumed.report()["examples"] == direct.report()["examples"]
    finally:
        resumed.close()
        direct.close()
    with pytest.raises(ValueError, match="different search"):
        CollisionSearch(path, 30, 3)


def test_interrupted_transaction_does_not_skip_candidates(tmp_path):
    search = CollisionSearch(tmp_path / "rollback.sqlite3", 6, 2)
    before = read_status(search.path)["state"]
    step = search.step

    def interrupted():
        step()
        raise KeyboardInterrupt

    search.step = interrupted
    try:
        with pytest.raises(KeyboardInterrupt):
            search.run(seconds=10)
        assert search.state == before == read_status(search.path)["state"]
        assert search.db.execute("SELECT COUNT(*) FROM halves").fetchone()[0] == 0
    finally:
        search.close()

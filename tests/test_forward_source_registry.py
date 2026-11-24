import json

import numpy as np
import pytest

from experiments.scalar_lift.forward_butson_descent import (
    ForwardTower,
    exponent_hash,
    verify_source,
)
from experiments.scalar_lift.forward_source_registry import load_source, sources


def test_all_twenty_sources_are_complete_and_replay_bounded_or_implicit_proofs():
    descriptors = sources()
    assert len(descriptors) == len({row["id"] for row in descriptors}) == 20
    json.dumps(descriptors)
    for descriptor in descriptors:
        source = load_source(descriptor)
        assert (source.order, source.phase) == (descriptor["order"], descriptor["phase"])
        assert source.proof["full_Gram_verified"] is True
        assert source.source_hash == source.hash
        assert source is load_source(source.id)
        assert source.entry(0, source.order - 1) in source.used_alphabet
        assert source.entry(source.order - 1, 0) in source.used_alphabet
        json.dumps(source.metadata())
        if source.exponents is not None:
            assert source.order <= 2000
            assert source.source_hash == exponent_hash(source.exponents)
            assert tuple(map(int, np.unique(source.exponents))) == source.used_alphabet
            assert not source.exponents.flags.writeable
        else:
            assert descriptor["family"] == "character"
            assert "implicit" in source.proof["hash_kind"]


def test_actual_dense_seed_and_tripling_gram_match_inductive_source_certificate():
    for identifier in ("cubic13-q24-t0", "cubic13-q24-t1", "cubic13-q24-t2"):
        source = load_source(identifier)
        assert verify_source(source.exponents, source.phase)["exponent_sha256"] == source.source_hash
        assert all(e % 3 in (0, 1) for e in source.used_alphabet)
    largest = load_source("cubic13-q24-t4")
    assert largest.order == 1053
    assert [row["order"] for row in largest.proof["stages"]] == [13, 39, 117, 351, 1053]
    assert largest.proof["all_literal_tensor_expansions_and_root_coefficients_verified"]


def test_shifted_complete_sextic_sources_recover_original_unreal_inputs():
    from experiments.scalar_lift.inverse_blocks import known_order_thirteen_source

    for n in (3, 13):
        source = load_source(f"unreal{n}-q6")
        assert source.order == n and source.phase == 6
        assert set(source.used_alphabet).issubset({0, 1, 3, 4})
        assert verify_source(source.exponents, 6)["exponent_sha256"] == source.source_hash
        original = (2 + 2 * np.eye(3, dtype=np.int64) if n == 3
                    else np.asarray(known_order_thirteen_source(), dtype=np.int64))
        assert np.array_equal((source.exponents + 1) % 6, original)
        assert source.proof["inverse_shift_for_original_unreal_substitution"] == 1


def test_implicit_linear_entry_accessor_matches_translation_invariant_source():
    from experiments.scalar_lift.forward_character_sources import character_source

    for prime in (73, 241):
        source = load_source(f"character-{prime}")
        original = character_source(prime)
        assert source.exponents is None
        for i, j in ((0, 1), (prime - 1, prime + 2), (source.order - 1, source.order - 2)):
            assert source.entry(i, j) == original.entry(divmod(i, prime), divmod(j, prime))
        assert source.proof["formal_Gram"]["number_of_Gram_entries_covered"] == source.order**2


def test_sources_feed_forward_executor_and_unfinished_targets_are_rejected():
    source = load_source("conference-1000-q8")
    tower = ForwardTower(source.exponents, source.phase)
    assert tower.proof(source.proof)["terminal_order"] == 4000
    assert tower.sign_entry(3999, 3998) in (-1, 1)
    with pytest.raises(ValueError, match="unfinished targets"):
        load_source("missing-1499-q120")
    changed = dict(sources()[0], order=1499)
    with pytest.raises(ValueError, match="differs"):
        load_source(changed)
    with pytest.raises(ValueError, match="outside"):
        source.entry(0.5, 0)

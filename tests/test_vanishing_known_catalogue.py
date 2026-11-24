import pytest

from experiments.scalar_lift.vanishing_known_catalogue import (
    certificate,
    exact_sum,
    iter_indexed_blocks,
    minimal_profiles,
    root_vectors,
    rotation_representatives,
)
from hadamard.coefficients import CyclotomicField


def test_complete_counts_and_minimality_proof_replay():
    report = certificate()
    assert report["proper_subsums_independently_checked"] == 1370
    assert [row["literal_profile_count"] for row in report["catalogues"]] == [211, 844]
    assert len(rotation_representatives()) == 9
    assert report["catalogues"][1]["weight_histogram"] == {
        "2": 60, "3": 40, "5": 24, "6": 120, "7": 240, "8": 240, "9": 120}


@pytest.mark.parametrize("modulus", (30, 120))
def test_integer_basis_against_independent_repository_field(modulus):
    field = CyclotomicField(modulus)
    zeta = field.zeta()
    table = root_vectors(modulus)
    assert all(tuple((zeta ** exponent).coefficients) == row for exponent, row in enumerate(table))
    assert any(exact_sum((0, 1), modulus))
    assert not any(exact_sum((0, modulus // 2), modulus))


def test_indexed_matcher_retains_duplicate_column_choices():
    assert sorted(iter_indexed_blocks([0, 60, 0])) == [(0, 1), (1, 2)]
    assert list(iter_indexed_blocks([0, 40, 80])) == [(0, 1, 2)]
    assert list(iter_indexed_blocks([])) == []
    with pytest.raises(ValueError):
        minimal_profiles(210)
    with pytest.raises(ValueError):
        list(iter_indexed_blocks([0, 1.0]))

import json

import numpy as np
import pytest
from flint import fmpz_mat

from experiments.scalar_lift.ambient_field_candidates import (
    family,
    first_split_conductors,
    period_field,
    source_coordinate_view,
    target_preflight,
)


def test_period_multiplication_against_independent_prime_cyclotomic_convolution():
    field = period_field(3, 31, (2, 3, 1))
    periods = np.zeros((3, 31), dtype=np.int64)
    for k, rep in enumerate(field.coset_representatives):
        periods[k, rep*field.subgroup % 31] = 1
    for j in range(3):
        direct = np.zeros(31, dtype=np.int64)
        for x in np.flatnonzero(periods[0]):
            for y in np.flatnonzero(periods[j]):
                direct[(x+y) % 31] += 1
        represented = field.multiplication[:, j] @ periods
        # The only rational relation on the 31 roots is Phi_31=sum roots.
        assert len(set(map(int, direct-represented))) == 1
    assert first_split_conductors(3, count=1, maximum_h=10)[0]['conductor'] == 31


def test_new_fields_have_exact_odd_trace_frames_and_faithful_period_action():
    for r, conductor in ((179, 268501), (223, 189997)):
        field = period_field(r)
        assert field.conductor == conductor
        assert field.multiplication.shape == (r, r)
        assert field.basis_denominator % 2 == 1
        assert np.array_equal(field.basis_numerator % 2, np.eye(r, dtype=np.int64))
        G, S, M = (fmpz_mat(x.tolist()) for x in
                   (field.trace_gram, field.basis_numerator, field.multiplication))
        unit = fmpz_mat(np.eye(r, dtype=np.int64).tolist())
        assert S.transpose()*G*S == field.basis_denominator**2*unit
        assert M.transpose()*G == G*M
        assert field.proof()['compositum_degree'] == 4*r
    with pytest.raises(ValueError, match='split-at-two'):
        period_field(3, 7, (2, 1, 1))


def test_actual_source_coordinate_views_have_complete_Gram_but_are_nonflat(tmp_path):
    for r, N in ((179, 9308), (223, 11596)):
        field = period_field(r)
        compact, metadata, _ = source_coordinate_view(field)
        assert compact.shape == (52, 52)
        assert np.array_equal(compact@compact.T, N*np.eye(52, dtype=np.int64))
        assert 0 in metadata['integer_entry_alphabet']
        assert max(abs(x) for x in metadata['integer_entry_alphabet']) > 1
        path = tmp_path/f'{r}.npz'
        report = family(r, bundle=path)
        with np.load(path, allow_pickle=False) as bundle:
            assert np.array_equal(bundle['compact_integer_matrix'], compact)
            assert np.array_equal(bundle['period_multiplication'], field.multiplication)
            assert json.loads(str(bundle['metadata_json'])) == report
        assert report['target_listed_unresolved_in_used_snapshot']
        assert not report['new_Hadamard_order_claimed']


def test_target_preflight_does_not_mistake_high_order_or_phase120_input_for_novelty():
    report = target_preflight()
    assert report['fixed_source_routes_above_20000'] == []
    source13 = next(x for x in report['source_routes'] if x['source_order'] == 13)
    assert source13['field_degree_compatible_targets'] == [8684, 9308, 11596]
    source15 = next(x for x in report['source_routes'] if x['source_order'] == 15)
    assert source15['order_divisible_targets'] == [11460]
    assert source15['field_degree_compatible_targets'] == []

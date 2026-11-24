import numpy as np
from sympy import Matrix

from experiments.scalar_lift.fifth_source_assembly import root_vectors
from experiments.scalar_lift.forward_source_registry import load_source
from experiments.scalar_lift.phase_independent_frontier import (
    _route,
    anchor_cases,
    field_from_alphabet,
    source_global_phase_hash,
)


def test_literal_fifth_root_field_is_not_the_declared_phase20_field():
    source = load_source('fifth-5-q20')
    actual = field_from_alphabet(source.phase, source.used_alphabet)
    assert actual['actual_root_generator_order'] == 5
    assert actual['entry_field_degree'] == 4
    roots = root_vectors(20)
    assert Matrix(roots[[0, 4, 8, 12, 16]].tolist()).rank() == 4
    for case in anchor_cases(source):
        row = _route(source, case, 23960)
        assert row['necessary_tests_pass']
        assert row['required_total_width'] == 4792


def test_global_anchor_reduces_cubic_phase24_without_claiming_a_new_source():
    source = load_source('cubic13-q24-t0')
    cases = anchor_cases(source)
    assert {row['entry_field_degree'] for row in cases} == {2, 4}
    assert {row['actual_root_generator_order'] for row in cases} == {6, 12}
    assert source_global_phase_hash(source) == source_global_phase_hash(load_source('unreal13-q6'))
    assert all(_route(source, case, 8684)['necessary_tests_pass'] for case in cases)


def test_odd_relation_and_width_are_separate_exact_rejections():
    from experiments.scalar_lift.catalogue_source_acquisition import data
    from experiments.scalar_lift.catalogue_source_acquisition import load_source as catalogue
    from experiments.scalar_lift.phase_independent_frontier import _literal_source

    d14 = catalogue('catalogue-d14-q4')
    assert all(any(e['kind'] == 'sign_block_width' for e in _route(d14, c, 4676)['exclusions'])
               for c in anchor_cases(d14))
    record = data('M13B')
    m13 = _literal_source('test-m13b', np.asarray(record['exponents']), 6, {'fixed_literal': True})
    assert all(any(e['kind'] == 'odd_integer_scalar_relation' for e in _route(m13, c, 8684)['exclusions'])
               for c in anchor_cases(m13))


def test_optional_higher_root_inputs_do_not_change_the_fixed_default_catalogue():
    from experiments.scalar_lift.forward_butson_descent import exponent_hash
    from experiments.scalar_lift.phase_independent_frontier import all_sources

    original = all_sources()
    expanded = all_sources(include_expanded=True)
    assert len(original) == 555
    assert len(expanded) == 649
    assert tuple(s.id for s in expanded[:555]) == tuple(s.id for s in original)
    assert all_sources() is original
    assert len({s.id for s in expanded}) == len(expanded)
    for source in expanded[555:]:
        assert source.phase == 12
        assert source.proof['full_Gram_verified']
        assert source.source_hash == exponent_hash(source.exponents)

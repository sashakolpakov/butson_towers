from contextlib import suppress

import numpy as np
import pytest

from experiments.scalar_lift.forward_butson_descent import ForwardTower, verify_source
from experiments.scalar_lift.forward_source_registry import load_source
from experiments.scalar_lift.partial_fifth_descent import (
    SavedNode,
    block_table,
    build_report,
    load_input,
    restrict,
    restrict_third,
    scalar_certificate,
    third_certificate,
)
from hadamard.inverse_lift import is_sign_hadamard


def test_relative_fifth_basis_and_actual_small_Butson_child():
    for q in (20, 40, 60, 80, 120):
        certificate = scalar_certificate(q)
        assert certificate['relative_field_degree'] == 4
        assert certificate['faithful_field_basis_exactly_identified']
        table = block_table(q)
        assert [e for e in range(q) if np.all(table[e] >= 0)] == [e for e in range(q) if e % 5]
    source = load_source('fifth-5-q20')
    child = restrict(source.exponents, 20)
    assert child.shape == (20, 20)
    proof = verify_source(child, 4)
    tower = ForwardTower(child, 4)
    assert tower.proof(proof)['terminal_order'] == 40
    assert is_sign_hadamard((1-2*tower.materialize()).tolist())


def test_bad_raw_common35_domain_is_rescued_by_an_explicit_input_basis():
    source = load_source('common120-35')
    with pytest.raises(ValueError, match='forbidden fifth coset'):
        restrict(source.exponents, 120)
    child = restrict(source.exponents, 120, 2)
    assert child.shape == (140, 140)
    assert verify_source(child, 24)['full_Gram_verified']


def test_partial_third_and_fifth_compose_to_a_full_sign_matrix():
    for q in (12, 24, 60, 120):
        assert third_certificate(q)['faithful_field_basis_exactly_identified']
    source = load_source('common120-15')
    child = restrict_third(source.exponents, 120)
    assert child.shape == (30, 30)
    assert verify_source(child, 40)['full_Gram_verified']
    candidates = []
    for shift in range(5):
        with suppress(ValueError):
            candidates.append(restrict(child, 40, shift))
    assert candidates
    terminal = ForwardTower(candidates[0], 8).materialize()
    assert terminal.shape == (480, 480)
    assert is_sign_hadamard((1-2*terminal).tolist())


def test_a_targeted_fixed_twelfth_root_gauge_and_recursive_saved_nodes(tmp_path):
    source = load_input('gauge3-aalto-9-12-000')
    assert source.proof['full_Gram_verified']
    assert len({e % 3 for e in source.used_alphabet}) < 3
    report = build_report(tmp_path, source_ids=('fifth-5-q20',))
    assert report['real_orders'] == [40]
    assert report['pending_nodes'] == 0
    assert report['exact_literal_duplicates_pruned'] > 0
    for record in report['nodes']:
        saved = SavedNode(record['bundle'])
        assert all(saved.replay().values())
        verified = saved.as_verified_source()
        assert verified.entry(0, 0) == saved.entry(0, 0)
        assert not verified.exponents.flags.writeable

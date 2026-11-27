import json
from pathlib import Path

import numpy as np
import pytest

from experiments.scalar_lift import hadamard_frame_search as frame
from experiments.scalar_lift import identity_multiplier as identity
from experiments.scalar_lift import quadratic_sign_lift as signs


ROOT = Path(__file__).resolve().parents[1]


def configuration(p):
    return json.loads((ROOT/f'certificates/scalar-lift/identity-parity-p{p}.json').read_text())


def test_identity_geometry_is_distinct_and_expands_the_exact_control():
    record = configuration(31)
    data = identity.chart(31, 3, 25)
    assert data['n'] == 11 and len(data['shifts']) == 5
    assert not np.any(data['phase'])
    assert all(data['index'][25*x % 31] == data['index'][x] for x in range(31))
    rows = identity.expanded_states(record['control_states'], data)
    assert np.array_equal(rows, frame.expand(record['control_states'], data))
    proof = identity.verify_states(record['control_states'], data, record['control_target_row'])
    assert proof['all_integer_correlations'] == [31]+[0]*30
    changed = list(record['control_states'])
    changed[1] ^= 1
    with pytest.raises(ValueError, match='correlation'):
        identity.verify_states(changed, data)
    # Independent per-site pair sums agree with every compressed chart row.
    rng = np.random.default_rng(31)
    for _ in range(12):
        axes = rng.integers(0, 4, data['n'])
        actual = axes[data['index']]
        for shift, row in zip(data['shifts'], data['rows'], strict=True):
            compressed = sum(count for (i,j,delta),count in row.items() if axes[i] == axes[j])
            assert compressed == np.sum(actual == np.roll(actual, -shift))
            assert all(delta == 0 for i,j,delta in row)


def test_every_public_map_coefficient_and_expanded_shift_replays():
    for p, count, expanded in ((31,211,6330),(223,11027,2447994)):
        result = identity.validate_map(configuration(p))
        assert result['coefficient_inputs'] == count
        assert result['expanded_correlations_checked'] == expanded
        assert result['all_coefficients_verified'] and result['all_expanded_shifts_verified']
    bad = configuration(31)
    bad['inverse_columns'][0] = bad['inverse_columns'][1]
    with pytest.raises(ValueError, match='rank'):
        identity.validate_map(bad)
    with pytest.raises(ValueError, match='identity-action'):
        identity.validate_map(dict(configuration(31), action='C3'))


def test_four_count_affine_kernel_matches_complete_actual_sign_census():
    record = configuration(31)
    data = identity.chart(31,3,25)
    axes = [s//2 for s in record['control_states']]
    system = identity.prepare(axes,data,record['control_target_row'],fixed_signs={1:1})
    n = data['n']
    masks = np.asarray([x for x in range(1 << n) if not x & 1 and x & 2])
    negative = (masks[:,None] >> np.arange(n)[None,:]) & 1
    expanded_signs = (1-2*negative)[:,data['index']]
    expanded_axes = np.asarray(axes)[data['index']]
    actual = np.stack([
        np.sum(expanded_signs*np.roll(expanded_signs,-s,axis=1)*
               (expanded_axes == np.roll(expanded_axes,-s))[None,:],axis=1)
        for s in range(1,31)],axis=1)
    valid = np.all(actual % 4 == 0,axis=1)
    for group,count in zip(system['groups'],system['negative_counts'],strict=True):
        valid &= negative[:,group].sum(axis=1) % 2 == count % 2
    expected = set(map(int,masks[valid]))
    observed = {signs.sign_mask(system,y) for y in range(1 << len(system['basis']))}
    assert expected == observed
    equations = signs.quadratic_mod8(system)
    for y in range(1 << len(system['basis'])):
        mask = signs.sign_mask(system,y)
        correlations = signs.correlations(system,mask)
        assert [signs.evaluate_quadratic(e,y) for e in equations] == [(c//4)%2 for c in correlations]


def test_complete_native_lift_recovers_actual_T31_and_rejects_bad_counts(tmp_path):
    record = configuration(31)
    data = identity.chart(31,3,25)
    axes = [s//2 for s in record['control_states']]
    result = identity.solve_bits(axes,data,record['control_target_row'],first_negative=1,
                                 seconds=2,max_roots=0,workdir=tmp_path)
    assert result['complete_for_fixed_axes'] and result['solutions']
    assert result['native']['assignments_covered'] == 1 << result['affine_dimension']
    for item in result['solutions']:
        proof = identity.verify_states(item['states'],data,record['control_target_row'])
        assert proof['shift_count'] == 31
    with pytest.raises(ValueError,match='norm p'):
        identity.prepare(axes,data,(-3,3,3,3))
    with pytest.raises(identity.InfeasibleAxes,match='positive origin'):
        identity.prepare(axes,data,record['control_target_row'],fixed_signs={0:1})


def test_all_three_T223_rows_and_both_head_gauges_are_supported():
    record = configuration(223)
    assert identity.target_rows(223,3) == [(-14,3,3,3),(13,6,3,3),(-5,9,9,6)]
    data = identity.chart(223,3,183)
    for example in [record['fixed_axis_control']]+record['additional_fixed_axis_controls']:
        system = identity.prepare(example['axes'],data,example['target_row'],
                                  fixed_signs={1:example['first_negative']})
        assert len(system['basis']) == example['affine_dimension']
        assert example['axes'][1:].count(0) == example['head_count']
        assert system['forced_signs'][1] == example['first_negative']
        normalized = identity.normalize_head(example['axes'],data)
        assert normalized[0] == normalized[1] == 0
    assert identity.target_rows(167,83) == []
    assert identity.target_rows(179,89) == []


def test_known_C3_controls_keep_their_exact_sign_bindings():
    # The two successful older charts remain separate coordinate actions.
    p127 = json.loads((ROOT/'certificates/scalar-lift/hadamard-frame-p127.json').read_text())
    p631 = json.loads((ROOT/'artifacts/hadamard-7572/t631.json').read_text())
    for p,m,mu,states,dimension in (
        (127,9,22,p127['completed_T_matrices']['states'],5),
        (631,21,458,p631['states'],13),
    ):
        data = frame.chart(p,m,mu)
        system = signs.prepare([s//2 for s in states],data)
        assert len(system['basis']) == dimension
        mask = sum((s%2) << i for i,s in enumerate(states))
        assert not any(signs.correlations(system,mask))
        rows = frame.expand(states,data)
        assert all(int(np.sum(rows*np.roll(rows,-s,axis=1))) == (p if s == 0 else 0)
                   for s in range(p))


def test_T491_feasibility_controls_have_the_actual_identity_geometry():
    record = configuration(491)
    assert identity.target_rows(491,7) == [(1,21,7,0)]
    assert record['signed_witness_claimed'] is False
    data = identity.chart(491,7,138)
    assert data['n'] == 71 and record['field_degree'] == 35
    example = record['model_controls'][0]
    system = identity.prepare(example['axes'],data,example['target_row'],fixed_signs={1:0})
    assert len(system['basis']) == 31
    equations = signs.quadratic_mod8(system)
    for y in (0,1,(1 << 31)-1):
        mask = signs.sign_mask(system,y)
        states = [2*a+((mask >> i) & 1) for i,a in enumerate(example['axes'])]
        rows = identity.expanded_states(states,data)
        full = [int(np.sum(rows*np.roll(rows,-s,axis=1))) for s in range(491)]
        compressed = signs.correlations(system,mask)
        assert full[0] == 491 and all(c % 4 == 0 for c in full[1:])
        assert list(compressed) == [full[s] for s in data['shifts']]
        assert [signs.evaluate_quadratic(eq,y) for eq in equations] == [(c//4)%2 for c in compressed]

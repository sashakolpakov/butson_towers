"""Independent signed-word census and wide mod8 controls for four counts."""
from itertools import product
import json
from pathlib import Path
import random

import numpy as np
import pytest

from experiments.scalar_lift import unrestricted_parity_map as maps
from experiments.scalar_lift import unrestricted_t_sign_model as model

ROOT = Path(__file__).resolve().parents[1]
T7 = [0,6,0,4,4,5,2]
TARGET7 = (2,1,1,1)


def direct_correlations(states):
    p = len(states)
    return [sum((1-2*(states[x]%2))*(1-2*(states[(x+s)%p]%2))
                for x in range(p) if states[x]//2 == states[(x+s)%p]//2)
            for s in range(p)]


def test_complete_order7_four_count_pipeline_matches_all_signed_words():
    direct = set()
    for word in product(range(8),repeat=6):
        states = (0,)+word
        sums = tuple(sum(1-2*(s%2) for s in states if s//2 == j) for j in range(4))
        if sums == TARGET7 and direct_correlations(states) == [7]+[0]*6:
            direct.add(states)
    found = set()
    patterns = prepared = assignments = 0
    data = maps.periodic_chart(7)
    for labels in product(range(4),repeat=6):
        axes = (0,)+labels
        if any(sum(axes[x] == axes[(x+s)%7] for x in range(7))%2 for s in range(1,7)):
            continue
        patterns += 1
        try:
            system = model.prepare(axes,data,TARGET7)
        except model.InfeasibleAxes:
            continue
        prepared += 1
        equations = model.quadratic_mod8(system)
        for assignment in range(1<<len(system['basis'])):
            assignments += 1
            mask = model.sign_mask(system,assignment)
            states = model.states_from_assignment(system,assignment)
            actual = direct_correlations(states)
            compressed = model.correlations(system,mask)
            assert compressed == actual[1:4]
            assert all(c%4 == 0 for c in compressed)
            observed = [model.evaluate_quadratic(e,assignment) for e in equations]
            assert observed == [(c//4)%2 for c in compressed]
            if any(observed) or not model.exact_count_mask(system,mask) or any(compressed):
                continue
            proof = model.verify_states(states,data,TARGET7)
            assert proof['all_integer_correlations'] == [7]+[0]*6
            assert proof['sequence_row_sums'] == list(TARGET7)
            found.add(tuple(states))
    assert (patterns,prepared,assignments) == (504,120,192)
    assert found == direct
    assert len(found) == 180
    assert tuple(T7) in found


def test_fixed_signs_extreme_counts_and_invalid_inputs():
    axes = [s//2 for s in T7]
    system = model.prepare(axes,7,TARGET7,fixed_signs={5:1})
    assert system['forced_signs'] == {0:0,5:1,2:0,6:0,3:0,4:0,1:0}
    assert len(system['basis']) == 0
    assert model.states_from_assignment(system,0) == T7
    assert model.exact_count_mask(system,system['particular'])
    with pytest.raises(model.InfeasibleAxes,match='origin'):
        model.prepare(axes,7,TARGET7,fixed_signs={0:1})
    with pytest.raises(model.InfeasibleAxes,match='count'):
        model.prepare(axes,7,(-2,1,1,1))
    with pytest.raises(ValueError,match='fixed'):
        model.prepare(axes,7,TARGET7,fixed_signs={7:0})
    with pytest.raises(ValueError,match='integer'):
        model.prepare([0.0]+axes[1:],7,TARGET7)
    with pytest.raises(ValueError,match='norm'):
        model.prepare(axes,7,(1,1,1,1))
    with pytest.raises(ValueError,match='range'):
        model.exact_count_mask(system,1<<7)


@pytest.mark.parametrize('p,target,dimension', [
    (167,(9,6,5,5),80),
    (179,(9,8,5,3),86),
])
def test_wide_four_count_mod4_and_mod8_against_actual_integer_correlations(p,target,dimension):
    config = json.loads((ROOT/f'certificates/scalar-lift/unrestricted-parity-p{p}.json').read_text())
    sampler = maps.DeterminantSampler(config,seed=p+1000)
    for _ in range(100):
        axes = sampler.next_axes(target_row=target,max_draws=100)
        if axes is None:
            continue
        try:
            system = model.prepare(axes,p,target)
        except model.InfeasibleAxes:
            continue
        break
    else:
        pytest.fail('no compatible deterministic control pattern')
    assert len(system['basis']) == dimension > 64
    assert len(system['groups']) == len(system['negative_counts']) == 4
    equations = model.quadratic_mod8(system)
    rng = random.Random(p)
    for assignment in [0]+[rng.getrandbits(dimension) for _ in range(12)]:
        mask = model.sign_mask(system,assignment)
        values = model.correlations(system,mask)
        assert all(c%4 == 0 for c in values)
        assert [model.evaluate_quadratic(e,assignment) for e in equations] == [(c//4)%2 for c in values]
        assert all(sum((mask>>i)&1 for i in group)%2 == count%2
                   for group,count in zip(system['groups'],system['negative_counts']))
    states = model.states_from_assignment(system,0)
    assert model.correlations(system,system['particular']) == direct_correlations(states)[1:(p+1)//2]


def test_final_verifier_ignores_search_pair_tables_and_rejects_nonzero():
    data = maps.periodic_chart(7)
    data['rows'] = [{} for _ in data['rows']]
    proof = model.verify_states(T7,data,TARGET7)
    assert proof['every_expanded_periodic_correlation_verified']
    assert proof['sequence_support_sizes'] == [2,1,3,1]
    with pytest.raises(ValueError,match='correlation'):
        model.verify_states([0]*7,data)
    with pytest.raises(ValueError,match='sums'):
        model.verify_states(T7,7,(-2,1,1,1))
    with pytest.raises(ValueError,match='positive origin'):
        model.verify_states([1]+T7[1:],7,TARGET7)

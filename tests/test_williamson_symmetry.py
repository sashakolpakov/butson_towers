"""Independent finite source census and literal multiplier checks."""
import hashlib
from itertools import product
import json
from pathlib import Path

import numpy as np
import pytest

from experiments.scalar_lift import williamson_symmetry as symmetry
from experiments.scalar_lift.williamson_scalar_endpoints import verify_t_rows

ROOT = Path(__file__).resolve().parents[1]
CERTIFICATES = ROOT/'certificates/scalar-lift'


def independent_rows(axes, signs, data):
    result = np.ones((4,data['p']), dtype=np.int64)
    for orbit, label, sign in zip(data['orbits'],axes,signs,strict=True):
        for k, point in enumerate(orbit):
            coordinate = label
            if data['action'] == 'C3' and label:
                coordinate = 1+(label-1+k) % 3
            result[:,point] = sign
            result[coordinate,point] = -sign
    return result


def actual_correlations(rows):
    return [int(sum(rows[j,x]*rows[j,(x+s) % rows.shape[1]]
                    for j in range(4) for x in range(rows.shape[1])))
            for s in range(rows.shape[1])]


def test_complete_W19_census_is_independent_of_forced_sign_solver():
    data = symmetry.chart(19,6,'C3',8)
    accepted, parity_count = [], 0
    for axes in product(range(4),repeat=3):
        _,_,_,forced = symmetry.coefficients(axes,data)
        matching = []
        for signs in product((-1,1),repeat=3):
            rows = independent_rows(axes,signs,data)
            assert np.array_equal(rows,symmetry.expand(axes,signs,data))
            correlations = actual_correlations(rows)
            if all(c % 16 == 0 for c in correlations[1:]):
                matching.append(list(signs));parity_count += 1
            if correlations != [76]+[0]*18:
                continue
            proof = symmetry.verify_source(rows,data)
            assert proof['row_sums'] == [7,3,3,3]
            index = (np.arange(19)[None,:]-np.arange(19)[:,None]) % 19
            matrices = rows[:,index]
            assert np.array_equal(sum(a@a.T for a in matrices),76*np.eye(19,dtype=np.int64))
            for i in range(4):
                for j in range(i):
                    assert np.array_equal(matrices[i]@matrices[j].T,matrices[j]@matrices[i].T)
            accepted.append([list(axes),list(signs)])
        assert matching == [forced.tolist()]
    record = json.loads((CERTIFICATES/'williamson-symmetry-controls.json').read_text())
    assert parity_count == 64 and len(accepted) == 18
    digest = hashlib.sha256(json.dumps(accepted,separators=(',',':')).encode()).hexdigest()
    assert digest == record['W19_census']['accepted_state_list_sha256']
    assert record['W19_census']['signed_words_checked'] == 512


@pytest.mark.parametrize('p,m,action',[(3,2,'identity'),(223,6,'identity'),
    (223,6,'C3'),(463,22,'identity'),(491,10,'identity'),(739,18,'C3'),
    (283,6,'identity'),(283,6,'C3'),(751,10,'identity'),(1423,18,'C3')])
def test_all_shift_geometry_and_actual_normalized_sign_equations(p,m,action):
    data = symmetry.chart(p,m,action)
    n = data['N']
    assert data['geometry_audit']['all_nonzero_shift_tables_verified'] == p-1
    rng = np.random.default_rng(p)
    axes = rng.integers(0,4,n)
    diagonal,weights,absolute,forced = symmetry.coefficients(axes,data)
    assert np.all(absolute % 2 == 1) and np.all(weights % 2 == 0)
    baseline = (absolute+1)//2 % 2
    columns = []
    for j in range(n):
        signs = np.ones(n,dtype=np.int64);signs[j] = -1
        normalized = diagonal+signs+np.einsum('sij,i,j->s',weights,signs,signs)
        columns.append(normalized//2 % 2 ^ baseline)
    assert np.array_equal(np.stack(columns,axis=1),np.eye(n,dtype=np.int64))
    for signs in (rng.choice((-1,1),n),forced):
        rows = independent_rows(axes,signs,data)
        actual = np.array([np.sum(rows*np.roll(rows,-s,axis=1),dtype=np.int64)
                           for s in range(p)])
        normalized = diagonal+signs+np.einsum('sij,i,j->s',weights,signs,signs)
        assert np.array_equal(actual[1:],4*normalized[data['index'][1:]])
        assert np.all(actual[1:] % 8 == 0)
        if np.array_equal(signs,forced):assert np.all(actual[1:] % 16 == 0)
    serialized = symmetry.coefficient_record(data)
    assert serialized['normalization']['compressed_Gram_scale'] == 4*p
    assert serialized['normalization']['source_at_zero'] == [1,1,1,1]
    assert len(serialized['rows']) == n


@pytest.mark.parametrize('p,m,action',[(3,2,'identity'),(7,2,'identity'),(19,6,'C3')])
def test_joint_row_count_gate_against_complete_signed_state_census(p,m,action):
    data = symmetry.chart(p,m,action)
    expected = set()
    for states in product(range(8),repeat=data['N']):
        axes = [s//2 for s in states]
        signs = [1-2*(s % 2) for s in states]
        row = independent_rows(axes,signs,data).sum(axis=1).tolist()
        if sum(r*r for r in row) == 4*p:
            expected.add(tuple(sorted(row)) if action == 'identity' else tuple(row))
    rows = (symmetry.identity_row_types(p,m,data['N']) if action == 'identity'
            else symmetry.c3_row_types(p,m,data['N']))
    assert expected == {tuple(r['row_sums']) for r in rows}


def test_current_selection_and_positive_T_inputs_replay():
    record = json.loads((CERTIFICATES/'williamson-symmetry-selection.json').read_text())
    controls = json.loads((CERTIFICATES/'williamson-symmetry-controls.json').read_text())
    frontier = record['frontier']['orders']
    assert len(frontier) == len(set(frontier)) == 81
    assert all(n < 10000 for n in frontier)
    assert 4792 in frontier and 5752 in frontier
    assert not {5164,7516,7572,8332} & set(frontier)
    for item in controls['known_T_inputs']:
        proof = verify_t_rows(np.asarray(item['rows'],dtype=np.int64))
        assert proof['order'] == item['t']
    selected = symmetry.select_charts(frontier,[item['t'] for item in controls['known_T_inputs']])
    assert selected == record['charts']
    assert len({row['p'] for row in selected}) == 7 and len(selected) == 9
    assert {r['target_order'] for chart in selected for r in chart['target_routes']} == {
        2676,2956,3004,4460,5692,5892,7924,9260}
    for chart in selected:
        for row in chart['row_types']:
            assert sum(r*r for r in row['row_sums']) == 4*chart['p']
        assert chart['nonzero_orbits'] == (chart['p']-1)//chart['m']


def test_palette_normalization_and_corruption_rejection():
    assert np.array_equal(symmetry.H.T@symmetry.H,4*np.eye(4,dtype=np.int64))
    assert np.array_equal(symmetry.H.T@np.ones(4,dtype=np.int64),2*np.ones(4,dtype=np.int64))
    allowed = {tuple(sign*symmetry.H[:,axis]) for axis in range(4) for sign in (-1,1)}
    assert allowed == {word for word in product((-1,1),repeat=4) if np.prod(word) == -1}
    record = json.loads((CERTIFICATES/'williamson-symmetry-controls.json').read_text())
    rows = np.asarray(record['positive_source']['rows'],dtype=np.int64)
    assert symmetry.verify_source(rows)['component_order'] == 19
    rows[0,1] *= -1
    with pytest.raises(ValueError,match='symmetry'):symmetry.verify_source(rows)
    with pytest.raises(ValueError,match='integer signs'):symmetry.verify_source(rows.astype(float))
    for args in [(9,2,'identity'),(19,3,'C3'),(19,6,'bad'),(True,2,'identity')]:
        with pytest.raises(ValueError):symmetry.chart(*args)
    with pytest.raises(ValueError,match='exact order'):symmetry.chart(19,6,'C3',7)
    data = symmetry.chart(19,6,'C3')
    with pytest.raises(ValueError,match='one axis'):symmetry.coefficients([0],data)
    with pytest.raises(ValueError,match='one axis'):symmetry.coefficients([0,1,4],data)
    with pytest.raises(ValueError,match='Boolean'):symmetry.coefficients([0,1,True],data)
    with pytest.raises(ValueError,match='one sign'):symmetry.expand([0,1,2],[1,0,1],data)
    with pytest.raises(ValueError):symmetry.identity_row_types(19,0,3)
    with pytest.raises(ValueError):symmetry.c3_row_types(19,2,9)
    with pytest.raises(ValueError):symmetry.select_charts([2676],[True])

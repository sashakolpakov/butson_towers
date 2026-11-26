"""Complete coefficient and finite-census checks of unrestricted parity maps."""
from copy import deepcopy
import json
from pathlib import Path

import pytest

from experiments.scalar_lift import hadamard_frame_search as frame
from experiments.scalar_lift import unrestricted_parity_map as maps
from experiments.scalar_lift.large_parity_map import BinaryField

ROOT = Path(__file__).resolve().parents[1]


def actual_parities(axes):
    p = len(axes)
    return [sum(axes[x] == axes[(x+s)%p] for x in range(p))%2
            for s in range(1,p)]


def test_complete_unrestricted_order7_parity_locus():
    config = maps.recover_map(7)
    assert maps.validate_map(config)['coefficients_checked'] == 79
    assert maps.validate_expanded_geometry(config)['exact_correlations_compared'] == 79*6
    field = BinaryField(3,config['field_modulus'])
    images = set()
    for P in range(8):
        for Q in range(8):
            if not(P or Q):
                continue
            for w in range(8):
                U,V = (w,field.mul(1^field.mul(Q,w),field.inverse(P))) if P else (field.inverse(Q),w)
                image = frame.binary_linear(config['inverse_columns'],P|(Q<<3)|(U<<6)|(V<<9))
                axes = [0]+[(image>>(2*j))&3 for j in range(6)]
                assert not any(actual_parities(axes))
                images.add(image)
    direct = set()
    for color in range(1<<12):
        axes = [0]+[(color>>(2*j))&3 for j in range(6)]
        if not any(actual_parities(axes)):
            direct.add(color)
    assert images == direct
    assert len(images) == 8*(8*8-1) == 504


@pytest.mark.parametrize('p,degree,coefficients,comparisons', [
    (167,83,55279,9176314),
    (179,89,63547,11311366),
])
def test_retained_map_full_coefficients_and_expanded_geometry(p,degree,coefficients,comparisons):
    config = json.loads((ROOT/f'certificates/scalar-lift/unrestricted-parity-p{p}.json').read_text())
    assert config['field_degree'] == degree
    proof = maps.validate_map(config)
    assert proof['all_quadratic_coefficients_verified']
    assert proof['coefficients_checked'] == coefficients
    assert (proof['color_map_rank'],proof['equation_map_rank']) == (4*degree,degree)
    geometry = maps.validate_expanded_geometry(config)
    assert geometry['exact_correlations_compared'] == comparisons
    assert geometry['all_expanded_shift_parities_match_orbit_coefficients']
    target = (9,6,5,5) if p == 167 else (9,8,5,3)
    sampler = maps.DeterminantSampler(config,seed=p)
    for _ in range(4):
        axes = sampler.next_axes(target_row=target,max_draws=100)
        assert axes is not None
        assert not any(actual_parities(axes))
        assert all(axes.count(j) >= abs(target[j]) and (axes.count(j)-target[j])%2 == 0
                   for j in range(4))


def test_configuration_checks_and_coefficient_checks_reject_corruption():
    config = maps.recover_map(7)
    bad = deepcopy(config)
    bad['inverse_columns'][0] = bad['inverse_columns'][1]
    with pytest.raises(ValueError,match='singular'):
        maps.validate_map(bad)
    bad = deepcopy(config)
    bad['field_modulus'] = 9
    with pytest.raises(ValueError,match='reducible'):
        maps.validate_map(bad)
    bad = deepcopy(config)
    bad['inverse_columns'][0] ^= bad['inverse_columns'][1]
    assert frame.binary_rank(bad['inverse_columns']) == 12
    with pytest.raises(ValueError,match='coefficient|linear|equation'):
        maps.validate_map(bad)
    with pytest.raises(ValueError,match='unrestricted'):
        maps.validate_map(dict(config,chart_kind='twisted'))
    for p in (True,7.0,9):
        with pytest.raises(ValueError):
            maps.periodic_chart(p)


def test_complete_row_types_and_sampler_draw_limits():
    assert maps.row_sum_types(167) == [
        (9,6,5,5),(9,7,6,1),(9,9,2,1),(10,7,3,3),(11,6,3,1)]
    assert maps.row_sum_types(179) == [
        (9,7,7,0),(9,8,5,3),(9,9,4,1),(11,7,3,0),(12,5,3,1),(13,3,1,0)]
    sampler = maps.DeterminantSampler(maps.recover_map(7),seed=13)
    assert sampler.next_axes(max_draws=0) is None
    assert sampler.draws == 0
    with pytest.raises(ValueError,match='row sums'):
        sampler.next_axes(target_row=(1,1,1,1))
    assert sampler.next_axes(target_row=(2,1,1,1),max_draws=100) is not None

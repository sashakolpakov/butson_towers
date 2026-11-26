"""Independent field, SL2 census, stream replay and expanded-parity checks."""
import json
from pathlib import Path
import random
import shutil

import numpy as np
import pytest

from experiments.scalar_lift import hadamard_frame_search as frame
from experiments.scalar_lift.exact_parity_sampler import NativeDeterminantSampler,MASK64
from experiments.scalar_lift.large_parity_map import recover_map


@pytest.fixture(scope='module')
def library_folder(tmp_path_factory):
    if not(shutil.which('c++') or shutil.which('clang++')):pytest.skip('C++17 compiler required')
    return tmp_path_factory.mktemp('native-parity-sampler')


def make(config,folder,seed=53):
    return NativeDeterminantSampler(config,seed=seed,workdir=folder)


def test_all_504_SL2_GF8_matrices_match_python_coordinate_map(library_folder):
    config=recover_map(p=19,m=3,mu=7,head_net=1,tail_net=1)
    sampler=make(config,library_folder)
    field=sampler.field;d=field.degree
    seen=set()
    for p in range(8):
        for q in range(8):
            if not(p or q):continue
            for w in range(8):
                b0,b1=(w,field.mul(1^field.mul(q,w),field.inverse(p))) if p else (field.inverse(q),w)
                assert field.mul(p,b1)^field.mul(q,b0)==1
                packed=p|(q<<d)|(b0<<(2*d))|(b1<<(3*d))
                want=frame.binary_linear(config['inverse_columns'],packed)
                assert sampler.color_parameters(p,q,w)==want
                seen.add(want)
    assert len(seen)==504


def test_degree37_native_field_and_all_148_bit_coordinates_match_python(library_folder):
    config=json.loads(Path('certificates/scalar-lift/t223-parity-map.json').read_text())
    sampler=make(config,library_folder);rng=random.Random(99)
    for a,b in [(0,0),(1,1)]+[(rng.getrandbits(37),rng.getrandbits(37)) for _ in range(300)]:
        product,reciprocal=sampler.field_operation(a,b)
        assert product==sampler.field.mul(a,b)
        assert reciprocal==(sampler.field.inverse(a) if a else 0)
    for p,q,w in [(0,1,0),(0,(1<<37)-1,(1<<37)-1)]+[(rng.getrandbits(37),rng.getrandbits(37),rng.getrandbits(37)) for _ in range(100)]:
        b0,b1=(w,sampler.field.mul(1^sampler.field.mul(q,w),sampler.field.inverse(p))) if p else (sampler.field.inverse(q),w)
        packed=p|(q<<37)|(b0<<74)|(b1<<111)
        assert sampler.color_parameters(p,q,w)==frame.binary_linear(config['inverse_columns'],packed)


def splitmix(state):
    state=(state+0x9e3779b97f4a7c15)&MASK64
    z=state;z=((z^(z>>30))*0xbf58476d1ce4e5b9)&MASK64
    z=((z^(z>>27))*0x94d049bb133111eb)&MASK64
    return state,z^(z>>31)


def xoshiro(state):
    def rotate(x,k):return ((x<<k)|(x>>(64-k)))&MASK64
    state=list(state);result=(rotate((state[1]*5)&MASK64,7)*9)&MASK64;t=(state[1]<<17)&MASK64
    state[2]^=state[0];state[3]^=state[1];state[1]^=state[2];state[0]^=state[3]
    state[2]^=t;state[3]=rotate(state[3],45)
    return state,result


def test_stream_and_draw_counts_match_independent_Python_generator_and_resume(library_folder):
    config=recover_map(p=19,m=3,mu=7,head_net=1,tail_net=1)
    sampler=make(config,library_folder,seed=0);seed=0;state=[];mask=7
    for _ in range(4):seed,word=splitmix(seed);state.append(word)
    assert tuple(state)==sampler.getstate()[1]
    expected=[]
    for _ in range(50):
        while True:
            state,p=xoshiro(state);state,q=xoshiro(state);p&=mask;q&=mask
            if p or q:break
        state,w=xoshiro(state);w&=mask
        expected.append(sampler.color_parameters(p,q,w))
    got=[sampler.draw_packed() for _ in range(50)]
    assert got==expected and sampler.draws==50
    assert sampler.getstate()==('exact-parity-xoshiro256ss-v1',tuple(state))
    record=json.loads(json.dumps({'state':sampler.getstate(),'draws':sampler.draws}))
    resumed=make(config,library_folder,seed=999)
    resumed.setstate(record['state']);resumed.draws=record['draws']
    assert [sampler.draw_packed() for _ in range(50)]==[resumed.draw_packed() for _ in range(50)]
    assert sampler.draws==resumed.draws==100
    with pytest.raises(ValueError,match='resume state'):resumed.setstate(random.Random(1).getstate())


def test_filtered_T223_axes_have_exact_expanded_parities_and_resume(library_folder):
    config=json.loads(Path('certificates/scalar-lift/t223-parity-map.json').read_text())
    sampler=make(config,library_folder);data=frame.chart(223,3,183)
    for h in (5,7,9,19):
        axes=sampler.next_axes(heads={h},normalize_first_head=True,max_draws=1000000)
        assert axes[0]==axes[1]==0
        assert sum(a==0 for a in axes[1:])==h
        rows=frame.expand([2*a for a in axes],data)
        assert all(int(np.sum(rows*np.roll(rows,-s,axis=1)))%2==0 for s in range(1,223))
    before=sampler.draws
    assert sampler.next_axes(heads=set(),max_draws=7) is None
    assert sampler.draws==before
    assert sampler.next_axes(heads={5},max_draws=0) is None
    assert sampler.draws==before

    # All-head axes would have odd unsigned correlations, so this exact
    # determinant sampler cannot accept them. A rejected batch must still
    # advance and serialize every draw before the next accepted sample.
    assert sampler.next_axes(heads={74},max_draws=7) is None
    assert sampler.draws==before+7
    resumed=make(config,library_folder,seed=999)
    resumed.setstate(json.loads(json.dumps(sampler.getstate())))
    resumed.draws=sampler.draws
    for h in (5,9):
        assert sampler.next_axes(heads={h},normalize_first_head=True,max_draws=1000000)==resumed.next_axes(heads={h},normalize_first_head=True,max_draws=1000000)
        assert sampler.draws==resumed.draws

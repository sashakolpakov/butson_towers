"""Crash/resume boundaries and a real small exact endpoint control."""
import json
import shutil
import sqlite3

import pytest

from experiments.scalar_lift import exact_t_campaign as campaign
from experiments.scalar_lift.large_parity_map import recover_map


def test_resume_continues_after_committed_head_sign_branch(monkeypatch,tmp_path):
    config=recover_map(p=19,m=3,mu=7,head_net=1,tail_net=1)
    axes=[0,0,0,0,1,2,3]
    monkeypatch.setattr(campaign.DeterminantSampler,'next_axes',lambda *a,**k:axes)
    calls=[]
    def checked_branch(*args,**kwargs):
        calls.append(kwargs['fixed_signs'])
        if len(calls)==1:campaign.STOP=True
        return {'complete_for_fixed_axes':True,'solutions':[],
                'native':{'assignments_covered':64},'mod8_roots':0,'exact_count_survivors':0}
    monkeypatch.setattr(campaign,'solve_bits',checked_branch)
    result=campaign.run(config,tmp_path,seconds=3,heads=(3,),branch_policy='all',nice=0)
    assert result['phase']=='stopped'
    assert len(calls)==1
    result=campaign.run(config,tmp_path,seconds=3,heads=(3,),branch_policy='all',
                        nice=0,max_patterns=1)
    assert result['phase']=='pattern_limit'
    assert len(calls)==2 and calls[0]!=calls[1]
    assert result['totals']['patterns_started']==1
    assert result['totals']['patterns_completed']==1
    assert result['totals']['native_affine_assignments_tested']==128
    with sqlite3.connect(tmp_path/'work.sqlite3') as db:
        assert db.execute('SELECT count(*) FROM attempts').fetchone()[0]==2
    assert not (tmp_path/'verified-witness.json').exists()
    with pytest.raises(ValueError,match='identical configuration'):
        campaign.run(config,tmp_path,seconds=1,heads=(3,),seed=99,nice=0)


def test_real_exact_witness_triggers_independent_endpoint_proof(monkeypatch,tmp_path):
    if shutil.which('c++') is None:pytest.skip('C++17 compiler unavailable')
    config=recover_map(p=31,m=3,mu=5,head_net=-1,tail_net=3)
    config['fixed_first_head_negative']=True
    planted=[0,1,2,2,7,6,2,6,2,5,3]
    axes=[s//2 for s in planted]
    monkeypatch.setattr(campaign.DeterminantSampler,'next_axes',lambda *a,**k:axes)
    result=campaign.run(config,tmp_path,seconds=8,heads=(1,),branch_policy='affine',nice=0)
    assert result['phase']=='verified_witness'
    witness=json.loads((tmp_path/'verified-witness.json').read_text())
    assert witness['expanded_verification']['every_expanded_periodic_correlation_verified']
    endpoints=json.loads((tmp_path/'endpoints.json').read_text())
    assert endpoints['orthogonal_design']['full_polynomial_Gram_verified']
    assert [r['order'] for r in endpoints['verified_real_orders']]==[372]
    assert endpoints['verified_real_orders'][0]['entry_alphabet']==[-1,1]


def test_committed_witness_is_recovered_after_output_crash(monkeypatch,tmp_path):
    if shutil.which('c++') is None:pytest.skip('C++17 compiler unavailable')
    config=recover_map(p=31,m=3,mu=5,head_net=-1,tail_net=3)
    config['fixed_first_head_negative']=True
    axes=[s//2 for s in [0,1,2,2,7,6,2,6,2,5,3]]
    monkeypatch.setattr(campaign.DeterminantSampler,'next_axes',lambda *a,**k:axes)
    write=campaign.atomic_json
    def crash(path,value):
        if path.name=='verified-witness.json':raise RuntimeError('simulated output crash')
        write(path,value)
    monkeypatch.setattr(campaign,'atomic_json',crash)
    with pytest.raises(RuntimeError,match='simulated output crash'):
        campaign.run(config,tmp_path,seconds=8,heads=(1,),branch_policy='affine',nice=0)
    assert not (tmp_path/'verified-witness.json').exists()
    monkeypatch.setattr(campaign,'atomic_json',write)
    monkeypatch.setattr(campaign,'solve_bits',lambda *a,**k:pytest.fail('committed witness was resampled'))
    result=campaign.run(config,tmp_path,seconds=8,heads=(1,),branch_policy='affine',nice=0)
    assert result['phase']=='verified_witness'
    assert (tmp_path/'endpoints.json').exists()


def test_detached_resume_rejects_different_configuration(tmp_path):
    from argparse import Namespace
    snapshot=tmp_path/'frozen';snapshot.mkdir()
    (snapshot/'parity-map.json').write_text('{"p":223}')
    requested=tmp_path/'other.json';requested.write_text('{"p":167}')
    with pytest.raises(ValueError,match='different parity configuration'):
        campaign.launch(Namespace(output=tmp_path,config=requested))


def test_native_sampler_state_resumes_in_the_same_work_ledger(monkeypatch,tmp_path):
    if shutil.which('c++') is None:pytest.skip('C++17 compiler unavailable')
    config=recover_map(p=19,m=3,mu=7,head_net=1,tail_net=1)
    monkeypatch.setattr(campaign,'solve_bits',lambda *a,**k:{
        'complete_for_fixed_axes':True,'solutions':[],
        'native':{'assignments_covered':64},'mod8_roots':0,'exact_count_survivors':0})
    for limit in (1,2):
        result=campaign.run(config,tmp_path,seconds=5,heads=(3,),branch_policy='affine',
                            nice=0,max_patterns=limit,sampler_engine='native')
        assert result['totals']['patterns_completed']==limit
        assert result['sampler']=='native'
    with sqlite3.connect(tmp_path/'work.sqlite3') as db:
        state=json.loads(db.execute("SELECT value FROM meta WHERE key='sampler_state'").fetchone()[0])
        assert state[0]=='exact-parity-xoshiro256ss-v1'
        assert db.execute('SELECT count(*) FROM patterns').fetchone()[0]==2

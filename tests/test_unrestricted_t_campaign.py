"""Real positive construction and transactional recovery for large-sign campaigns."""
import json
import sqlite3

import pytest

from experiments.scalar_lift import unrestricted_t_campaign as campaign
from experiments.scalar_lift import unrestricted_parity_map as parity


def small_config():
    return parity.recover_map(7)


def plant(monkeypatch):
    states=[0,6,0,4,4,5,2]
    monkeypatch.setattr(campaign.parity.DeterminantSampler,'next_axes',
                        lambda *a,**k:[s//2 for s in states])


def test_real_small_campaign_constructs_a_verified_endpoint(monkeypatch,tmp_path):
    pytest.importorskip('ortools')
    plant(monkeypatch)
    result=campaign.run(small_config(),tmp_path,seconds=5,pattern_seconds=2,
                        backends=('cp',),nice=0)
    assert result['phase']=='verified_witness'
    proof=json.loads((tmp_path/'verified-witness.json').read_text())
    assert proof['expanded_verification']['every_expanded_periodic_correlation_verified']
    endpoints=json.loads((tmp_path/'endpoints.json').read_text())
    assert endpoints['orthogonal_design']['full_polynomial_Gram_verified']
    assert [r['order'] for r in endpoints['verified_real_orders']]==[84]


def test_committed_positive_result_survives_a_file_output_crash(monkeypatch,tmp_path):
    pytest.importorskip('ortools')
    plant(monkeypatch);write=campaign.common.atomic_json;config=small_config()
    def crash(path,value):
        if path.name=='verified-witness.json':raise RuntimeError('simulated output interruption')
        write(path,value)
    monkeypatch.setattr(campaign.common,'atomic_json',crash)
    with pytest.raises(RuntimeError,match='simulated output interruption'):
        campaign.run(config,tmp_path,seconds=5,backends=('cp',),nice=0)
    monkeypatch.setattr(campaign.common,'atomic_json',write)
    monkeypatch.setattr(campaign,'solver',lambda _backend:lambda *a,**k:pytest.fail('stored positive resampled'))
    result=campaign.run(config,tmp_path,seconds=5,backends=('cp',),nice=0)
    assert result['phase']=='verified_witness'
    assert (tmp_path/'endpoints.json').exists()


def test_timeout_and_resume_do_not_claim_exclusion(monkeypatch,tmp_path):
    plant(monkeypatch);config=small_config()
    def unknown(*args,**kwargs):
        return {'solver_status':'UNKNOWN','complete_for_fixed_axes':False,'solutions':[]}
    monkeypatch.setattr(campaign,'solver',lambda _backend:unknown)
    result=campaign.run(config,tmp_path,seconds=3,backends=('cp',),nice=0,max_patterns=1)
    assert result['phase']=='pattern_limit'
    assert result['totals']['unknown_patterns']==1
    assert result['totals']['proved_infeasible_patterns']==0
    assert not (tmp_path/'verified-witness.json').exists()
    result=campaign.run(config,tmp_path,seconds=3,backends=('cp',),nice=0,max_patterns=1)
    assert result['totals']['patterns_checked']==1
    with sqlite3.connect(tmp_path/'work.sqlite3') as db:
        assert db.execute('SELECT state FROM jobs').fetchone()[0]=='unknown'


def test_mod4_infeasibility_is_recorded_instead_of_killing_worker(monkeypatch,tmp_path):
    plant(monkeypatch)
    def impossible(*args,**kwargs):
        raise campaign.exact.InfeasibleAxes('contradictory modulo-four equations')
    monkeypatch.setattr(campaign.exact,'prepare',impossible)
    monkeypatch.setattr(campaign,'solver',lambda _backend:lambda *a,**k:pytest.fail('infeasible signs reached backend'))
    result=campaign.run(small_config(),tmp_path,seconds=3,backends=('cp',),nice=0,max_patterns=1)
    assert result['totals']['proved_infeasible_patterns']==1
    assert result['totals']['unknown_patterns']==0
    assert result['phase']=='pattern_limit'

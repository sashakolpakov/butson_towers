"""Transactional resume, independent witness gating and frozen-source controls."""
from argparse import Namespace
import json
from pathlib import Path
import shutil
import sqlite3
from types import SimpleNamespace

import pytest

from experiments.scalar_lift import identity_t_campaign as campaign

ROOT = Path(__file__).resolve().parents[1]


def small_config():
    return json.loads((ROOT/'certificates/scalar-lift/identity-parity-p31.json').read_text())


def native_result(complete=True, states=None):
    return {'solver_status': 'COMPLETE' if complete else 'UNFINISHED',
            'complete_for_fixed_axes': complete,
            'solutions': [] if states is None else [{'states': states}],
            'native': {'assignments_covered': 64}, 'mod8_roots': 0, 'exact_count_survivors': 0}


def install_axes(monkeypatch, axes):
    class Sampler:
        def __init__(self):
            import random
            self.rng = random.Random(3);self.draws = 0
        def next_axes(self, **kwargs):
            assert kwargs['normalize_first_head'] is False
            assert kwargs['heads'] == {axes[1:].count(0)}
            self.rng.random();self.draws += 1
            return list(axes)
    monkeypatch.setattr(campaign, 'make_sampler', lambda *a: Sampler())


def test_schedule_interleaves_all_three_identity_rows_and_gauges():
    schedule = campaign.strata(223, 3)
    assert [r['head_count'] for r in schedule] == [5, 4, 2, 7, 6, 4, 9, 8, 6]
    assert [r['target_row'] for r in schedule[:3]] == [
        [-14, 3, 3, 3], [13, 6, 3, 3], [-5, 9, 9, 6]]
    assert [r['first_negative'] for r in schedule[:3]] == [1, 0, 1]
    for offsets in [(1,), (0, 0), (-2,)]:
        with pytest.raises(ValueError):campaign.strata(223, 3, head_offsets=offsets)


def test_resume_starts_after_the_committed_head_sign_branch(monkeypatch, tmp_path):
    config = small_config()
    axes = [0, 0, 0, 0, 1, 2, 2, 2, 3, 3, 3]
    install_axes(monkeypatch, axes)
    calls = []
    def solve(*args, **kwargs):
        calls.append(kwargs['fixed_signs'])
        if len(calls) == 1:campaign.STOP = True
        return native_result()
    monkeypatch.setattr(campaign.identity, 'solve_bits', solve)
    first = campaign.run(config, tmp_path, seconds=5, head_offsets=(2,), nice=0, sampler_engine='python')
    assert first['phase'] == 'stopped' and first['totals']['patterns_completed'] == 0
    second = campaign.run(config, tmp_path, seconds=5, head_offsets=(2,), nice=0,
                          sampler_engine='python', max_patterns=1)
    assert second['phase'] == 'pattern_limit'
    assert len(calls) == 2 and calls[0] != calls[1]
    assert second['totals']['patterns_started'] == second['totals']['patterns_completed'] == 1
    assert second['totals']['branch_attempts'] == second['totals']['branches_completed'] == 2
    assert second['totals']['native_affine_assignments_tested'] == 128
    with pytest.raises(ValueError, match='identical configuration'):
        campaign.run(config, tmp_path, seconds=5, head_offsets=(2,), nice=0, seed=999)


def test_timeout_remains_pending_and_is_replayed_before_new_branches(monkeypatch, tmp_path):
    config = small_config();axes = [0, 0, 0, 0, 1, 2, 2, 2, 3, 3, 3]
    install_axes(monkeypatch, axes);calls = []
    def solve(*args, **kwargs):
        calls.append(kwargs['fixed_signs'])
        if len(calls) == 1:
            campaign.STOP = True
            return native_result(False)
        return native_result()
    monkeypatch.setattr(campaign.identity, 'solve_bits', solve)
    first = campaign.run(config, tmp_path, seconds=5, head_offsets=(2,), nice=0, sampler_engine='python')
    assert first['pending_branches'] == 1 and first['totals']['branches_completed'] == 0
    assert first['totals']['patterns_completed'] == 0
    second = campaign.run(config, tmp_path, seconds=5, head_offsets=(2,), nice=0,
                          sampler_engine='python', max_patterns=1)
    assert calls[0] == calls[1] and calls[1] != calls[2]
    assert second['pending_branches'] == 0
    assert second['totals']['branch_attempts'] == 3 and second['totals']['branches_completed'] == 2
    assert second['totals']['patterns_completed'] == 1
    assert second['totals']['native_affine_assignments_tested'] == 192


def test_real_native_identity_witness_gets_a_separate_endpoint_proof(monkeypatch, tmp_path):
    if shutil.which('c++') is None:pytest.skip('C++17 compiler unavailable')
    config = small_config();states = config['control_states']
    install_axes(monkeypatch, [s//2 for s in states])
    result = campaign.run(config, tmp_path, seconds=10, head_offsets=(0,), nice=0, sampler_engine='python')
    assert result['phase'] == 'verified_witness'
    witness = json.loads((tmp_path/'verified-witness.json').read_text())
    assert witness['target_row'] == [-2, 3, 3, 3]
    assert witness['expanded_verification']['every_expanded_periodic_correlation_verified']
    endpoints = json.loads((tmp_path/'endpoints.json').read_text())
    assert endpoints['orthogonal_design']['full_polynomial_Gram_verified']
    assert [r['order'] for r in endpoints['verified_real_orders']] == [372]


def test_false_positive_backend_result_cannot_be_written(monkeypatch, tmp_path):
    config = small_config();states = config['control_states'].copy();states[2] ^= 1
    install_axes(monkeypatch, [s//2 for s in states])
    monkeypatch.setattr(campaign.identity, 'solve_bits', lambda *a, **k: native_result(states=states))
    with pytest.raises(ValueError):
        campaign.run(config, tmp_path, seconds=5, head_offsets=(0,), nice=0, sampler_engine='python')
    assert not (tmp_path/'verified-witness.json').exists() and not (tmp_path/'endpoints.json').exists()
    with sqlite3.connect(tmp_path/'work.sqlite3') as db:
        assert db.execute("SELECT value FROM meta WHERE key='pending_witness'").fetchone() is None


def test_durable_witness_recovers_after_output_crash_without_resampling(monkeypatch, tmp_path):
    config = small_config();states = config['control_states']
    install_axes(monkeypatch, [s//2 for s in states])
    monkeypatch.setattr(campaign.identity, 'solve_bits', lambda *a, **k: native_result(states=states))
    write = campaign.atomic_json
    def crash(path, value):
        if path.name == 'verified-witness.json':raise RuntimeError('simulated output crash')
        write(path, value)
    monkeypatch.setattr(campaign, 'atomic_json', crash)
    with pytest.raises(RuntimeError, match='simulated output crash'):
        campaign.run(config, tmp_path, seconds=5, head_offsets=(0,), nice=0, sampler_engine='python')
    monkeypatch.setattr(campaign, 'atomic_json', write)
    monkeypatch.setattr(campaign.identity, 'solve_bits', lambda *a, **k: pytest.fail('durable witness was resampled'))
    result = campaign.run(config, tmp_path, seconds=5, head_offsets=(0,), nice=0, sampler_engine='python')
    assert result['phase'] == 'verified_witness' and (tmp_path/'endpoints.json').exists()


def test_detached_worker_freezes_exact_sources_and_detects_tampering(monkeypatch, tmp_path):
    calls = []
    real_popen = campaign.subprocess.Popen
    def popen(command, **kwargs):
        calls.append((command, kwargs));return SimpleNamespace(pid=12345)
    monkeypatch.setattr(campaign.subprocess, 'Popen', popen)
    args = Namespace(output=tmp_path, config=ROOT/'certificates/scalar-lift/identity-parity-p31.json',
                     seconds=10, branch_seconds=1, head_offsets='0,2,4', seed=123,
                     branch_policy='all', sampler='native', nice=0, max_patterns=1,
                     row_targets=[(-2, 3, 3, 3)])
    result = campaign.launch(args)
    assert result['pid'] == 12345
    snapshot = tmp_path/'frozen';campaign.check_snapshot(snapshot)
    manifest = json.loads((snapshot/'source-manifest.json').read_text())
    assert all('.local' not in name and '__pycache__' not in name for name in manifest)
    assert calls[0][1]['env']['PYTHONPATH'] == str(snapshot)
    assert calls[0][1]['env']['PYTHONDONTWRITEBYTECODE'] == '1'
    assert '-B' in calls[0][0]
    assert '--target-row=-2,3,3,3' in calls[0][0]
    launch = json.loads((tmp_path/'launch.json').read_text())
    assert launch['cwd'] == str(snapshot) and launch['environment']['PYTHONPATH'] == str(snapshot)
    if shutil.which('c++') is not None:
        # Run the frozen command synchronously as a bounded control, rather
        # than trusting an editable checkout to supply a missing dependency.
        monkeypatch.setattr(campaign.subprocess, 'Popen', real_popen)
        child = campaign.subprocess.run(calls[0][0], cwd=snapshot, env=calls[0][1]['env'],
                                         text=True, capture_output=True, timeout=20)
        assert child.returncode == 0, child.stdout+child.stderr
        status = json.loads((tmp_path/'status.json').read_text())
        assert status['phase'] in {'pattern_limit', 'verified_witness'}
        assert status['totals']['branch_attempts'] >= 1
        campaign.check_snapshot(snapshot)
    changed = snapshot/'experiments/scalar_lift/identity_multiplier.py'
    changed.chmod(0o644);changed.write_text(changed.read_text()+'\n')
    with pytest.raises(ValueError, match='source hash mismatch'):
        campaign.launch(args)


def test_native_sampler_stream_resumes_in_same_ledger(monkeypatch, tmp_path):
    if shutil.which('c++') is None:pytest.skip('C++17 compiler unavailable')
    config = small_config()
    monkeypatch.setattr(campaign.identity, 'solve_bits', lambda *a, **k: native_result())
    for limit in (1, 2):
        result = campaign.run(config, tmp_path, seconds=10, head_offsets=(0,), nice=0,
                              sampler_engine='native', max_patterns=limit)
        assert result['totals']['patterns_completed'] == limit
    with sqlite3.connect(tmp_path/'work.sqlite3') as db:
        state = json.loads(db.execute("SELECT value FROM meta WHERE key='sampler_state'").fetchone()[0])
        assert state[0] == 'exact-parity-xoshiro256ss-v1'
        assert db.execute('SELECT count(*) FROM patterns').fetchone()[0] == 2


def test_accepted_sample_and_pattern_insert_are_one_transaction(monkeypatch, tmp_path):
    config = small_config();axes = [0, 0, 0, 0, 1, 2, 2, 2, 3, 3, 3]
    install_axes(monkeypatch, axes)
    real_connect = sqlite3.connect
    class CrashInsert(sqlite3.Connection):
        def execute(self, sql, *args, **kwargs):
            if sql.startswith('INSERT INTO patterns('):
                raise RuntimeError('simulated pattern insert crash')
            return super().execute(sql, *args, **kwargs)
    monkeypatch.setattr(campaign.sqlite3, 'connect', lambda *a, **k: real_connect(*a, factory=CrashInsert, **k))
    with pytest.raises(RuntimeError, match='simulated pattern insert crash'):
        campaign.run(config, tmp_path, seconds=5, head_offsets=(2,), nice=0, sampler_engine='python')
    with real_connect(tmp_path/'work.sqlite3') as db:
        assert db.execute('SELECT count(*) FROM patterns').fetchone()[0] == 0
        assert db.execute("SELECT value FROM meta WHERE key='sampler_state'").fetchone() is None
        assert db.execute("SELECT value FROM meta WHERE key='sampler_draws'").fetchone() is None
    monkeypatch.setattr(campaign.sqlite3, 'connect', real_connect)
    monkeypatch.setattr(campaign.identity, 'solve_bits', lambda *a, **k: native_result())
    resumed = campaign.run(config, tmp_path, seconds=5, head_offsets=(2,), nice=0,
                           sampler_engine='python', max_patterns=1)
    assert resumed['totals']['patterns_started'] == resumed['totals']['patterns_completed'] == 1


@pytest.mark.parametrize('flag', [['--target-row=-14,3,3,3'], ['--target-row', '-14,3,3,3']])
def test_explicit_negative_target_row_parses_and_restricts_schedule(tmp_path, flag):
    args = campaign.parse_args(['--config', 'map.json', '--output', str(tmp_path), *flag])
    assert args.row_targets == [(-14, 3, 3, 3)]
    schedule = campaign.strata(223, 3, args.row_targets)
    assert [entry['head_count'] for entry in schedule] == [5, 7, 9]

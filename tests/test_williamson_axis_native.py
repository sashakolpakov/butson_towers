"""Exact native censuses, standalone resume and frozen-input boundaries."""
from argparse import Namespace
import copy
import fcntl
import hashlib
from itertools import product
import json
import os
from pathlib import Path
import shutil
import subprocess

import numpy as np
import pytest

from experiments.scalar_lift import williamson_axis_native as native
from experiments.scalar_lift import williamson_symmetry as symmetry


pytestmark=pytest.mark.skipif(shutil.which('c++') is None,reason='C++17 compiler unavailable')


def options(coefficients,output,*,mode='enumerate',patterns=64,seconds=10):
    return Namespace(coefficients=coefficients,output=output,mode=mode,seconds=seconds,
        seed=123456789,max_patterns=patterns,max_witnesses=0,temperature=1000,
        minimum_temperature=1,epoch_steps=11,nice=0)


def argv(args):
    return ['--coefficients',str(args.coefficients),'--output',str(args.output),'--mode',args.mode,
        '--seconds',str(args.seconds),'--seed',str(args.seed),'--max-patterns',str(args.max_patterns),
        '--max-witnesses',str(args.max_witnesses),'--temperature',str(args.temperature),
        '--minimum-temperature',str(args.minimum_temperature),'--epoch-steps',str(args.epoch_steps),
        '--nice',str(args.nice)]


@pytest.fixture(scope='module')
def coefficients(tmp_path_factory):
    directory=tmp_path_factory.mktemp('native-williamson-inputs')
    result={}
    for p,m,action in ((19,6,'C3'),(3,2,'identity')):
        data=symmetry.chart(p,m,action)
        record=symmetry.coefficient_record(data)
        path=directory/f'p{p}.json';path.write_text(json.dumps(record))
        result[p]=(path,record)
    return result


@pytest.fixture(scope='module')
def p19_census(tmp_path_factory,coefficients):
    output=tmp_path_factory.mktemp('native-williamson-census')
    args=options(coefficients[19][0],output)
    native.main(argv(args))
    return output,args


def independent_rows(axes,signs,p,m,mu,action):
    """Literal residue expansion, independent of coefficient/sign preparation."""
    rows=np.ones((4,p),dtype=np.int64);seen={0};orbit=0
    for representative in range(1,p):
        if representative in seen:continue
        x=representative
        for k in range(m):
            assert x not in seen
            seen.add(x);a=axes[orbit]
            if action=='C3' and a:a=1+(a-1+k)%3
            rows[:,x]=signs[orbit]*(1-2*(np.arange(4)==a))
            x=x*mu%p
        assert x==representative
        orbit+=1
    assert len(seen)==p and orbit==len(axes)
    return rows


def correlations(rows):
    return [int(np.sum(rows*np.roll(rows,-s,axis=1))) for s in range(rows.shape[1])]


def test_native_p19_census_matches_all_512_literal_signed_words(p19_census,coefficients):
    output,_=p19_census;record=coefficients[19][1];expected={}
    for axes in product(range(4),repeat=3):
        for signs in product((-1,1),repeat=3):
            rows=independent_rows(axes,signs,19,6,record['mu'],'C3')
            if correlations(rows)==[76]+[0]*18:expected[(axes,signs)]=rows.tolist()
    witnesses=[json.loads(line) for line in (output/'witnesses.jsonl').read_text().splitlines()]
    assert len(expected)==len(witnesses)==18
    assert {(tuple(w['axes']),tuple(w['signs'])) for w in witnesses}==set(expected)
    for witness in witnesses:
        rows=np.asarray(witness['rows'],dtype=np.int64)
        assert witness['rows']==expected[(tuple(witness['axes']),tuple(witness['signs']))]
        matrices=np.asarray([[np.roll(row,i) for i in range(19)] for row in rows])
        assert np.array_equal(sum(a@a.T for a in matrices),76*np.eye(19,dtype=np.int64))
        assert all(np.array_equal(a,a.T) for a in matrices)
        assert all(np.array_equal(a@b.T,b@a.T) for a in matrices for b in matrices)
    status=json.loads((output/'status.json').read_text())
    assert status['phase']=='exhausted' and status['totals']['axis_patterns_evaluated']==64
    assert status['totals']['verified_witnesses']==18


def test_native_identity_control_exhausts_all_four_W3_words(tmp_path,coefficients):
    args=options(coefficients[3][0],tmp_path,patterns=4)
    native.main(argv(args))
    status=json.loads((tmp_path/'status.json').read_text())
    assert status['phase']=='exhausted' and status['totals']['verified_witnesses']==4
    for line in (tmp_path/'witnesses.jsonl').read_text().splitlines():
        witness=json.loads(line)
        assert witness['signs']==[-1]
        rows=independent_rows(witness['axes'],witness['signs'],3,2,2,'identity')
        assert correlations(rows)==[12,0,0] and witness['rows']==rows.tolist()


def checkpoint_without_elapsed(output):
    lines=(output/'checkpoint.txt').read_text().splitlines()
    counts=lines[2].split();counts[5]='ACTIVE_ELAPSED';lines[2]=' '.join(counts)
    return lines


@pytest.mark.parametrize('mode',['random','anneal'])
def test_frozen_renamed_launcher_resumes_exact_native_stream(tmp_path,coefficients,mode):
    resumed=tmp_path/'resumed';args=options(coefficients[19][0],resumed,mode=mode,patterns=16)
    native.main(argv(args))
    launch=json.loads((resumed/'launch.json').read_text());command=launch['resume_command']
    assert Path(command[2]).parent==resumed/'frozen'
    assert Path(command[2]).with_suffix('.cpp').exists()
    command[command.index('--max-patterns')+1]='32'
    subprocess.run(command,cwd=tmp_path,env=dict(os.environ,PYTHONPATH='',PYTHONNOUSERSITE='1'),
        check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    fresh=tmp_path/'fresh';fresh_args=options(coefficients[19][0],fresh,mode=mode,patterns=32)
    native.main(argv(fresh_args))
    assert checkpoint_without_elapsed(resumed)==checkpoint_without_elapsed(fresh)
    assert (resumed/'witnesses.jsonl').read_text()==(fresh/'witnesses.jsonl').read_text()
    assert json.loads((resumed/'status.json').read_text())['totals']==json.loads((fresh/'status.json').read_text())['totals']


def test_parameters_bind_both_frozen_inputs_and_manifest(p19_census):
    output,args=p19_census;path=output/'frozen/input.txt';manifest=output/'frozen/source-manifest.json'
    original=path.read_bytes();original_manifest=manifest.read_bytes()
    try:
        path.chmod(0o644);path.write_bytes(original+b'0\n')
        with pytest.raises(ValueError,match='frozen source hash mismatch'):native.prepare(args)
        values=json.loads(original_manifest);values['input.txt']=hashlib.sha256(path.read_bytes()).hexdigest()
        manifest.chmod(0o644);manifest.write_text(json.dumps(values))
        with pytest.raises(ValueError,match='manifest differs'):native.prepare(args)
    finally:
        path.write_bytes(original);path.chmod(0o444)
        manifest.chmod(0o644);manifest.write_bytes(original_manifest);manifest.chmod(0o444)
    native.prepare(args)


def test_native_geometry_guard_and_output_lock_do_not_clobber_an_owner(tmp_path,p19_census,coefficients):
    output,_=p19_census;command=json.loads((output/'launch.json').read_text())['command']
    bad=copy.deepcopy(coefficients[19][1]);bad['rows'][0][0][3]+=2
    bad_input=tmp_path/'bad-input.txt';bad_input.write_bytes(native.serialize(bad,native.row_targets(bad)))
    bad_output=tmp_path/'bad';attempt=list(command);attempt[1]=str(bad_input);attempt[2]=str(bad_output)
    result=subprocess.run(attempt,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    assert result.returncode!=0 and b'coefficient audit failed' in result.stderr
    assert json.loads((bad_output/'status.json').read_text())['phase']=='failed'
    locked=tmp_path/'locked';locked.mkdir();status=locked/'status.json';status.write_text('{"phase":"owner-running"}\n')
    with (locked/'worker.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        attempt=list(command);attempt[2]=str(locked)
        result=subprocess.run(attempt,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    assert result.returncode!=0 and b'another native worker owns' in result.stderr
    assert json.loads(status.read_text())=={'phase':'owner-running'}
    assert not (locked/'failure.json').exists()


def test_native_budget_is_cumulative_on_resume(tmp_path,coefficients):
    args=options(coefficients[19][0],tmp_path,mode='random',patterns=0,seconds=0.02)
    native.main(argv(args));first=json.loads((tmp_path/'status.json').read_text())
    native.main(argv(args));second=json.loads((tmp_path/'status.json').read_text())
    assert first['phase']==second['phase']=='budget_exhausted'
    assert first['totals']==second['totals']
    assert second['elapsed_seconds']>=first['elapsed_seconds']>=args.seconds

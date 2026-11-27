"""Independent source census and durable exact CP campaign controls."""
from itertools import product
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np
import pytest

pytest.importorskip('ortools.sat.python.cp_model')
from ortools.sat.python import cp_model
from experiments.scalar_lift import williamson_axis_cp as cp
from experiments.scalar_lift import williamson_symmetry as model
from experiments.scalar_lift import williamson_campaign as campaign


@pytest.mark.parametrize('p,m,action,target',[
    (3,2,'identity',[3,-1,-1,-1]),
    (7,6,'C3',[-5,-1,-1,-1]),
    (19,6,'C3',[7,3,3,3]),
])
def test_positive_controls_are_complete_literal_sources(p,m,action,target):
    result=cp.solve(p,m,action,seconds=3,target=target)
    assert result['source_claimed']
    rows=np.asarray(result['rows'],dtype=np.int64)
    assert set(np.unique(rows)) <= {-1,1}
    assert rows.sum(axis=1).tolist()==target
    assert [int(np.sum(rows*np.roll(rows,-s,axis=1))) for s in range(p)] == [4*p]+[0]*(p-1)
    blocks=np.stack([np.asarray([np.roll(row,i) for i in range(p)]) for row in rows])
    assert np.array_equal(sum(a@a.T for a in blocks),4*p*np.eye(p,dtype=np.int64))


def test_full_CP_census_matches_independent_signed_word_enumeration():
    data=model.chart(19,6,'C3');expected=set()
    # The residue geometry is rebuilt directly here, outside the compressed
    # coefficients or forced-sign equation used by either solver.
    transform=np.ones((4,4),dtype=np.int64)-2*np.eye(4,dtype=np.int64)
    for states in product(range(8),repeat=3):
        rows=np.ones((4,19),dtype=np.int64)
        seen={0};orbit=0
        for representative in range(1,19):
            if representative in seen:continue
            axis,negative=divmod(states[orbit],2);point=representative
            for k in range(6):
                seen.add(point)
                coordinate=0 if axis==0 else 1+(axis-1+k)%3
                rows[:,point]=(1-2*negative)*transform[:,coordinate]
                point=point*8%19
            orbit+=1
        if all(int(np.sum(rows*np.roll(rows,-s,axis=1)))==0 for s in range(1,19)):
            expected.add(states)
    assert len(expected)==18
    instance,axes,negative,_=cp.build(data,gauge=False)
    class Census(cp_model.CpSolverSolutionCallback):
        def __init__(self):super().__init__();self.states=set()
        def on_solution_callback(self):
            self.states.add(tuple(2*self.value(a)+self.value(s) for a,s in zip(axes,negative,strict=True)))
    solver=cp_model.CpSolver();solver.parameters.enumerate_all_solutions=True
    solver.parameters.num_search_workers=1
    census=Census();status=solver.solve(instance,census)
    assert status==cp_model.OPTIMAL and census.states==expected


def test_incompatible_chart_is_infeasible_and_malformed_norm_rejected():
    result=cp.solve(13,6,'C3',seconds=2,gauge=False)
    assert result['status']=='INFEASIBLE' and not result['source_claimed']
    with pytest.raises(ValueError,match='squared norm'):
        cp.build(model.chart(19,6,'C3'),[7,3,3,4])
    with pytest.raises(ValueError,match='exact integer'):
        cp.build(model.chart(19,6,'C3'),[7.0,3,3,3])


def test_committed_positive_source_recovers_without_resampling(tmp_path):
    config={'p':19,'m':6,'mu':8,'action':'C3','target_orders':[228],
            'target_row':[7,3,3,3]}
    for _ in range(2):
        campaign.run(config,tmp_path,budget=3,per_call=2,base_seed=1901,nice=0)
        status=json.loads((tmp_path/'status.json').read_text())
        assert status['phase']=='verified_witness'
    assert len(list(tmp_path.glob('attempt-*.json')))==1
    witness=json.loads((tmp_path/'verified-witness.json').read_text())
    witness['rows'][0][1]*=-1
    (tmp_path/'verified-witness.json').write_text(json.dumps(witness))
    with pytest.raises(ValueError,match='disagree'):
        campaign.run(config,tmp_path,budget=3,per_call=2,base_seed=1901,nice=0)


def frozen_command(tmp_path,config):
    root=tmp_path/'frozen';module=root/'experiments/scalar_lift'
    module.mkdir(parents=True)
    (root/'experiments/__init__.py').touch();(module/'__init__.py').touch()
    names=['williamson_symmetry.py','williamson_axis_cp.py','williamson_campaign.py']
    source=Path(campaign.__file__).parent
    for name in names:shutil.copy2(source/name,module/name)
    (module/'config.json').write_text(json.dumps(config))
    names.append('config.json')
    manifest={'files':{name:hashlib.sha256((module/name).read_bytes()).hexdigest() for name in names}}
    (module/'source-manifest.json').write_text(json.dumps(manifest))
    output=tmp_path/'work'
    command=[sys.executable,'-B','-m','experiments.scalar_lift.williamson_campaign',
             '--config',str(module/'config.json'),'--output',str(output),
             '--seconds','20','--call-seconds','15','--seed','1001','--nice','0']
    env=dict(os.environ,PYTHONPATH=str(root),OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    return root,module,output,command,env


def test_actual_frozen_source_executes_without_checkout_fallback(tmp_path):
    config={'p':19,'m':6,'mu':8,'action':'C3','target_orders':[228],'target_row':[7,3,3,3]}
    root,module,output,command,env=frozen_command(tmp_path,config)
    subprocess.run(command,cwd=root,env=env,check=True,capture_output=True,text=True,timeout=10)
    assert json.loads((output/'status.json').read_text())['phase']=='verified_witness'
    (module/'williamson_symmetry.py').write_text((module/'williamson_symmetry.py').read_text()+'\n# altered\n')
    result=subprocess.run(command,cwd=root,env=env,capture_output=True,text=True,timeout=10)
    assert result.returncode and 'Frozen source hash mismatch' in result.stderr


def test_stop_request_interrupts_a_live_CP_call(tmp_path):
    config={'p':463,'m':22,'mu':38,'action':'identity','target_orders':[9260],
            'target_row':[23,-21,-21,-21]}
    root,module,output,command,env=frozen_command(tmp_path,config)
    process=subprocess.Popen(command,cwd=root,env=env,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    try:
        deadline=time.monotonic()+8
        while not (output/'status.json').exists() and time.monotonic()<deadline:time.sleep(.02)
        assert (output/'status.json').exists()
        (output/'STOP').touch()
        stdout,stderr=process.communicate(timeout=5)
        assert process.returncode==0,(stdout,stderr)
        assert json.loads((output/'status.json').read_text())['phase']=='stopped'
        assert (output/'checkpoint.json').exists()
    finally:
        if process.poll() is None:process.kill();process.wait()

"""Single-core, checkpointed CP-SAT calls on one selected Williamson chart."""
import argparse
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import signal
import threading
import time

import numpy as np
from ortools.sat.python import cp_model
from experiments.scalar_lift.williamson_symmetry import chart, coefficients, expand, verify_source
from experiments.scalar_lift.williamson_axis_cp import build

STOP = False
ACTIVE_SOLVER = None


def utc():
    return datetime.now(timezone.utc).isoformat()


def atomic(path, data):
    temporary = path.with_name(path.name+'.tmp')
    temporary.write_text(json.dumps(data,indent=2)+'\n')
    temporary.replace(path)


def stop(_signal,_frame):
    global STOP
    STOP = True
    if ACTIVE_SOLVER is not None:
        ACTIVE_SOLVER.stop_search()


def run(config, output, budget, per_call, base_seed, nice):
    global ACTIVE_SOLVER
    output = Path(output).resolve();output.mkdir(parents=True,exist_ok=True)
    lock = (output/'worker.lock').open('a')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    if nice:os.nice(nice)
    for sig in (signal.SIGTERM,signal.SIGINT):signal.signal(sig,stop)
    sources = {p.name:hashlib.sha256(p.read_bytes()).hexdigest()
               for p in (Path(__file__),Path(__file__).with_name('williamson_symmetry.py'),
                         Path(__file__).with_name('williamson_axis_cp.py'))}
    fingerprint = hashlib.sha256(json.dumps({'config':config,'sources':sources,
                                'base_seed':base_seed,'per_call':per_call},sort_keys=True).encode()).hexdigest()
    checkpoint_path = output/'checkpoint.json'
    checkpoint = json.loads(checkpoint_path.read_text()) if checkpoint_path.exists() else {
        'fingerprint':fingerprint,'attempts_completed':0,'elapsed_seconds':0.,
        'branches':0,'conflicts':0,'unknown_calls':0}
    if checkpoint['fingerprint'] != fingerprint:
        raise ValueError('Resume requires identical sources, chart and solver settings')
    data = chart(config['p'],config['m'],config['action'])
    if data['mu'] != config['mu']:
        raise ValueError('The selected multiplier differs from the frozen chart')
    before = checkpoint['elapsed_seconds'];started = time.monotonic()
    status = {'schema':'selected-williamson-cp-campaign-v1','p':data['p'],'m':data['m'],
              'mu':data['mu'],'action':data['action'],'source_kind':'Williamson',
              'nonzero_orbits':data['N'],'backend':'CP-SAT','pid':os.getpid(),
              'phase':'running','start_utc':utc(),'fingerprint':fingerprint,
              'target_orders':config['target_orders'],'target_row':config['target_row'],
              'budget_seconds':budget,'budget_seconds_this_session':max(0,budget-before),
              'solver_threads':1,'model_scope':'Selected multiplier chart, exact target row sums and stated residue-scaling gauge.',
              'totals':{'verified_witnesses':0}}
    status_lock = threading.Lock();event = threading.Event()
    def write_status():
        with status_lock:
            status.update(updated_utc=utc(),elapsed_seconds_this_session=time.monotonic()-started,
                          elapsed_seconds=before+time.monotonic()-started)
            atomic(output/'status.json',status)
    def heartbeat():
        while not event.wait(5):write_status()
    write_status();thread=threading.Thread(target=heartbeat,daemon=True);thread.start()
    try:
        if (output/'verified-witness.json').exists():
            saved=json.loads((output/'verified-witness.json').read_text())
            if saved['fingerprint'] != fingerprint:
                raise ValueError('Saved witness has a different source fingerprint')
            saved_rows=expand(saved['axes'],saved['signs'],data)
            if not np.array_equal(saved_rows,np.asarray(saved['rows'],dtype=np.int64)):
                raise ValueError('Saved witness rows disagree with its orbit states')
            verify_source(saved_rows,data)
            status['phase']='verified_witness';status['totals']['verified_witnesses']=1
            return
        while not STOP and not (output/'STOP').exists():
            remaining=budget-before-(time.monotonic()-started)
            if remaining <= 0:
                status['phase']='budget_exhausted';break
            attempt=checkpoint['attempts_completed']
            seed=(base_seed+attempt)%2147483647
            model,axes,negative,description=build(data,config['target_row'],True)
            if attempt == 0:
                (output/'model.pbtxt').write_text(str(model.proto))
            solver=cp_model.CpSolver()
            solver.parameters.num_search_workers=1
            solver.parameters.max_time_in_seconds=min(per_call,remaining)
            solver.parameters.random_seed=seed
            solver.parameters.randomize_search=True
            ACTIVE_SOLVER=solver
            with status_lock:
                status['current']={'attempt':attempt,'seed':seed,'model':description,
                                   'limit_seconds':solver.parameters.max_time_in_seconds}
            write_status()
            finished=threading.Event();answer={}
            def invoke():
                try:answer['code']=solver.solve(model)
                except BaseException as exc:answer['error']=exc
                finally:finished.set()
            solve_thread=threading.Thread(target=invoke,daemon=True);solve_thread.start()
            while not finished.wait(.25):
                if STOP or (output/'STOP').exists():solver.stop_search()
            solve_thread.join();ACTIVE_SOLVER=None
            if 'error' in answer:raise answer['error']
            code=answer['code']
            result={'attempt':attempt,'seed':seed,'status':solver.status_name(code),
                    'wall_seconds':solver.wall_time,'branches':solver.num_branches,
                    'conflicts':solver.num_conflicts,'completed_utc':utc()}
            if code in (cp_model.OPTIMAL,cp_model.FEASIBLE):
                color=[solver.value(v) for v in axes]
                signs=np.asarray([1-2*solver.value(v) for v in negative],dtype=np.int64)
                if not np.array_equal(coefficients(color,data)[3],signs):
                    raise ArithmeticError('Solver signs disagree with independent formula')
                rows=expand(color,signs,data);proof=verify_source(rows,data)
                if proof['row_sums'] != config['target_row']:
                    raise ArithmeticError('Expanded source row sums differ')
                witness={'p':data['p'],'m':data['m'],'mu':data['mu'],'action':data['action'],
                         'source_kind':'Williamson','axes':color,'signs':signs.tolist(),
                         'rows':rows.tolist(),'expanded_verification':proof,
                         'fingerprint':fingerprint,'found_utc':utc(),
                         'target_orders':config['target_orders'],'new_order_claimed':False}
                atomic(output/'verified-witness.json',witness)
                status['phase']='verified_witness';status['totals']['verified_witnesses']=1
            atomic(output/f'attempt-{attempt:06d}.json',result)
            checkpoint.update(attempts_completed=attempt+1,
                              elapsed_seconds=before+time.monotonic()-started,
                              branches=checkpoint['branches']+solver.num_branches,
                              conflicts=checkpoint['conflicts']+solver.num_conflicts,
                              unknown_calls=checkpoint['unknown_calls']+int(code==cp_model.UNKNOWN))
            atomic(checkpoint_path,checkpoint)
            with status_lock:
                status['totals'].update(solver_calls=checkpoint['attempts_completed'],
                                       unknown_calls=checkpoint['unknown_calls'],
                                       branches=checkpoint['branches'],conflicts=checkpoint['conflicts'])
            write_status()
            if status['phase']=='verified_witness':break
            if code == cp_model.INFEASIBLE:
                status['phase']='model_infeasible';break
            if code == cp_model.MODEL_INVALID:
                raise ValueError('CP-SAT rejected its validated model')
        else:
            status['phase']='stopped'
        if STOP or (output/'STOP').exists():status['phase']='stopped'
    except BaseException as exc:
        status.update(phase='failed',error=f'{type(exc).__name__}: {exc}')
        raise
    finally:
        ACTIVE_SOLVER=None
        checkpoint['elapsed_seconds']=before+time.monotonic()-started
        atomic(checkpoint_path,checkpoint)
        event.set();thread.join(timeout=6);write_status();lock.close()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--seconds',type=float,default=259200)
    parser.add_argument('--call-seconds',type=float,default=900)
    parser.add_argument('--seed',type=int,required=True)
    parser.add_argument('--nice',type=int,default=10)
    args=parser.parse_args()
    if args.seconds<=0 or args.call_seconds<=0:parser.error('Positive budgets required')
    snapshot=Path(__file__).resolve().parent
    manifest=snapshot/'source-manifest.json'
    if manifest.exists():
        for name,digest in json.loads(manifest.read_text())['files'].items():
            path=snapshot/name
            if not path.resolve().is_relative_to(snapshot) or hashlib.sha256(path.read_bytes()).hexdigest()!=digest:
                raise ValueError('Frozen source hash mismatch: '+name)
    run(json.loads(args.config.read_text()),args.output,args.seconds,args.call_seconds,args.seed,args.nice)


if __name__=='__main__':main()

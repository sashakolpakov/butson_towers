"""Durable exact sign-solver campaigns for unrestricted T167 and T179 inputs.

The verified determinant map generates parity-compatible coordinate labels.
SAT/CP-SAT then impose mod-four affine links, mod-eight equations, all four
exact counts and every integer correlation. UNKNOWN is never an exclusion.
"""
from __future__ import annotations
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import signal
import shutil
import sqlite3
import subprocess
import sys
import threading
import time

from experiments.scalar_lift import exact_t_campaign as common
from experiments.scalar_lift import unrestricted_parity_map as parity
from experiments.scalar_lift import unrestricted_t_sign_model as exact

SOURCE_NAMES=tuple(dict.fromkeys(common.SOURCE_NAMES+(
    'unrestricted_t_campaign.py','unrestricted_parity_map.py','unrestricted_t_sign_model.py',
    'unrestricted_t_sign_sat.py','unrestricted_t_sign_cp.py')))
STOP=False


def row_types(p):
    """All four-square row sums up to coordinate signs and permutations."""
    from math import isqrt
    result=[]
    for a in range(isqrt(p)+1):
        for b in range(a+1):
            for c in range(b+1):
                rest=p-a*a-b*b-c*c
                if rest<0:continue
                d=isqrt(rest)
                if d<=c and d*d==rest:result.append((a,b,c,d))
    return result


def solver(backend):
    if backend=='sat':
        from experiments.scalar_lift.unrestricted_t_sign_sat import solve_rows
        if shutil.which('cryptominisat5') is None:raise RuntimeError('CryptoMiniSat5 is required for the SAT backend')
        return solve_rows
    if backend=='cp':
        from experiments.scalar_lift.unrestricted_t_sign_cp import solve
        return solve
    raise ValueError('backend must be sat or cp')


def run(config,output,seconds=259200,pattern_seconds=900,seed=20261126010,
        backends=('sat','cp'),nice=10,max_patterns=None):
    global STOP
    STOP=False
    if seconds<=0 or pattern_seconds<=0 or not backends:raise ValueError('positive limits and a backend are required')
    for backend in backends:solver(backend)
    output=Path(output).resolve();output.mkdir(parents=True,exist_ok=True)
    lock=(output/'worker.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    if nice:os.nice(nice)
    data=parity.periodic_chart(config['p']);audit=parity.validate_map(config)
    targets=row_types(data['p'])
    hashes={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in SOURCE_NAMES}
    specification={'config':config,'seed':seed,'backends':list(backends),'targets':targets,'source_sha256':hashes}
    fingerprint=hashlib.sha256(json.dumps(specification,sort_keys=True).encode()).hexdigest()
    db=sqlite3.connect(output/'work.sqlite3')
    db.execute('PRAGMA journal_mode=WAL');db.execute('PRAGMA synchronous=FULL')
    db.executescript('''CREATE TABLE IF NOT EXISTS meta(key TEXT PRIMARY KEY,value TEXT NOT NULL);
    CREATE TABLE IF NOT EXISTS jobs(id INTEGER PRIMARY KEY,digest TEXT UNIQUE NOT NULL,
      axes TEXT NOT NULL,target TEXT NOT NULL,backend TEXT NOT NULL,seed INTEGER NOT NULL,
      state TEXT NOT NULL DEFAULT 'active',result TEXT);''')
    def get(key,default=None):
        row=db.execute('SELECT value FROM meta WHERE key=?',(key,)).fetchone()
        return json.loads(row[0]) if row else default
    def put(key,value):db.execute('INSERT OR REPLACE INTO meta VALUES (?,?)',(key,json.dumps(value)))
    if get('fingerprint',fingerprint)!=fingerprint:
        raise ValueError('resume requires the identical frozen sources and configuration')
    totals=get('totals',{'patterns_started':0,'patterns_checked':0,'proved_infeasible_patterns':0,
                         'unknown_patterns':0,'verified_witnesses':0})
    sampler=parity.DeterminantSampler(config,seed=seed)
    if get('sampler_state') is not None:
        sampler.rng.setstate(common.tuplify(get('sampler_state')));sampler.draws=get('sampler_draws',0)
    with db:
        put('fingerprint',fingerprint);put('specification',specification);put('parity_audit',audit)
    started=time.monotonic();guard=threading.RLock();finished=threading.Event()
    status={'schema':'exact-unrestricted-sign-campaign-v1','pid':os.getpid(),'p':data['p'],
            'phase':'running','start_utc':common.utc(),'fingerprint':fingerprint,
            'budget_seconds_this_session':seconds,'pattern_limit_seconds':pattern_seconds,
            'backends':list(backends),'target_row_types':targets,'totals':totals,
            'method':'SL2 parity -> exact mod4/mod8 sign constraints -> four exact counts and all integer correlations',
            'sampler':'arbitrary-size exact determinant sampler','current':None}
    def checkpoint(phase=None):
        with guard:
            if phase is not None:status['phase']=phase
            status.update(updated_utc=common.utc(),elapsed_seconds_this_session=time.monotonic()-started,
                          sampled_determinant_matrices=sampler.draws)
            if status['current']:
                status['current']['elapsed_seconds']=time.monotonic()-status['current']['monotonic_start']
            common.atomic_json(output/'status.json',status)
    def heartbeat():
        while not finished.wait(15):checkpoint()
    def stop(_signum,_frame):
        global STOP
        STOP=True
    previous={}
    if threading.current_thread() is threading.main_thread():
        for sig in (signal.SIGTERM,signal.SIGINT):previous[sig]=signal.signal(sig,stop)
    checkpoint();thread=threading.Thread(target=heartbeat,daemon=True);thread.start()
    phase='budget_exhausted'
    try:
        pending=get('pending_witness')
        if (output/'verified-witness.json').exists():
            pending=json.loads((output/'verified-witness.json').read_text())
        if pending:
            exact.verify_states(pending['states'],data,pending['target_row'])
            common.atomic_json(output/'verified-witness.json',pending)
            common.verified_endpoints(pending['states'],data,output)
            phase='verified_witness'
        else:
            while time.monotonic()-started<seconds:
                if STOP or (output/'STOP').exists():phase='stopped';break
                if max_patterns is not None and totals['patterns_checked']>=max_patterns:
                    phase='pattern_limit';break
                active=db.execute("SELECT id,axes,target,backend,seed FROM jobs WHERE state='active' ORDER BY id LIMIT 1").fetchone()
                if active is None:
                    index=totals['patterns_started'];target=targets[index%len(targets)]
                    backend=backends[(index//len(targets))%len(backends)]
                    axes=sampler.next_axes(target_row=target,max_draws=10000)
                    if axes is None:checkpoint();continue
                    packed=json.dumps(axes,separators=(',',':'));target_json=json.dumps(target)
                    digest=hashlib.sha256((packed+target_json+backend).encode()).hexdigest()
                    if db.execute('SELECT 1 FROM jobs WHERE digest=?',(digest,)).fetchone():continue
                    with db:
                        cursor=db.execute('INSERT INTO jobs(digest,axes,target,backend,seed) VALUES (?,?,?,?,?)',
                                          (digest,packed,target_json,backend,seed+index))
                        with guard:totals['patterns_started']+=1
                        put('sampler_state',sampler.rng.getstate());put('sampler_draws',sampler.draws);put('totals',totals)
                    active=(cursor.lastrowid,packed,target_json,backend,seed+index)
                job,packed,target_json,backend,job_seed=active
                axes=json.loads(packed);target=tuple(json.loads(target_json))
                result=None
                try:
                    system=exact.prepare(axes,data,target)
                    dimension=len(system['basis']);del system
                except exact.InfeasibleAxes as exc:
                    dimension=None
                    result={'solver_status':'INFEASIBLE_ARITHMETIC','reason':str(exc),
                            'complete_for_fixed_axes':True,'solutions':[]}
                remaining=seconds-(time.monotonic()-started)
                if remaining<=0:break
                with guard:
                    status['current']={'job_id':job,'backend':backend,'target_row':target,
                                       'affine_dimension':dimension,'monotonic_start':time.monotonic(),
                                       'limit_seconds':min(pattern_seconds,remaining)}
                checkpoint()
                if result is None:
                    result=solver(backend)(axes,data,target,seconds=min(pattern_seconds,remaining),
                                           seed=job_seed,max_solutions=1,workdir=output/'solver')
                solutions=result.get('solutions',[])
                for item in solutions:exact.verify_states(item['states'],data,target)
                complete=result.get('complete_for_fixed_axes',False)
                with db:
                    db.execute('UPDATE jobs SET state=?,result=? WHERE id=?',
                               ('witness' if solutions else 'excluded' if complete else 'unknown',json.dumps(result),job))
                    with guard:
                        totals['patterns_checked']+=1
                        totals['proved_infeasible_patterns']+=int(complete and not solutions)
                        totals['unknown_patterns']+=int(not complete and not solutions)
                        totals['verified_witnesses']+=len(solutions)
                    if solutions:
                        witness={'p':data['p'],'m':1,'mu':1,'target_row':target,**solutions[0],
                                 'fingerprint':fingerprint,'found_utc':common.utc()}
                        put('pending_witness',witness)
                    put('totals',totals)
                if solutions:
                    common.atomic_json(output/'verified-witness.json',witness)
                    common.verified_endpoints(witness['states'],data,output)
                    phase='verified_witness';break
                with guard:status['current']=None
                checkpoint()
        checkpoint(phase);print(json.dumps(status),flush=True)
        return status
    except BaseException as exc:
        with guard:status['error']=f'{type(exc).__name__}: {exc}'
        checkpoint('failed');raise
    finally:
        finished.set();thread.join(timeout=2)
        for sig,handler in previous.items():signal.signal(sig,handler)
        db.close();lock.close()


def launch(args):
    output=args.output.resolve();snapshot=output/'frozen';module=snapshot/'experiments'/'scalar_lift'
    if snapshot.exists():
        if json.loads((snapshot/'parity-map.json').read_text())!=json.loads(args.config.read_text()):
            raise ValueError('frozen parity map differs; use the saved launch command or a fresh directory')
    else:
        module.mkdir(parents=True);(snapshot/'experiments'/'__init__.py').touch();(module/'__init__.py').touch()
        for name in SOURCE_NAMES:shutil.copy2(Path(__file__).with_name(name),module/name)
        shutil.copy2(args.config,snapshot/'parity-map.json')
        common.atomic_json(snapshot/'source-manifest.json',{
            str(p.relative_to(snapshot)):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(snapshot.rglob('*')) if p.is_file()})
    command=[sys.executable,'-u','-m','experiments.scalar_lift.unrestricted_t_campaign',
             '--config',str(snapshot/'parity-map.json'),'--output',str(output),'--seconds',str(args.seconds),
             '--pattern-seconds',str(args.pattern_seconds),'--seed',str(args.seed),
             '--backends',args.backends,'--nice',str(args.nice)]
    if args.max_patterns is not None:command+=['--max-patterns',str(args.max_patterns)]
    with (output/'worker.log').open('ab',buffering=0) as log:
        child=subprocess.Popen(command,cwd=snapshot,env=dict(os.environ,PYTHONPATH=str(snapshot)),
                               stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
    common.atomic_json(output/'launch.json',{'pid':child.pid,'command':command,'launched_utc':common.utc()})
    print(json.dumps({'pid':child.pid,'output':str(output)}))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--seconds',type=float,default=259200);parser.add_argument('--pattern-seconds',type=float,default=900)
    parser.add_argument('--backends',default='sat,cp');parser.add_argument('--seed',type=int,default=20261126010)
    parser.add_argument('--nice',type=int,default=10);parser.add_argument('--max-patterns',type=int);parser.add_argument('--detach',action='store_true')
    args=parser.parse_args()
    if args.detach:launch(args)
    else:run(json.loads(args.config.read_text()),args.output,args.seconds,args.pattern_seconds,args.seed,
             tuple(args.backends.split(',')),args.nice,args.max_patterns)


if __name__=='__main__':main()

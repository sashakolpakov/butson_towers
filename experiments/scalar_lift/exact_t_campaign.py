"""Resumable exact T-design searches using determinant parity coordinates.

This is the successful SL2 -> mod4 -> mod8 -> integer correlation pipeline.
SQLite commits completed branch checks and sampler state together. No numerical
objective is used. A partial branch is never reported as a completed pattern.
"""
from __future__ import annotations

import argparse
import datetime
import fcntl
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import shutil
import signal
import sqlite3
import subprocess
import sys
import time

from experiments.scalar_lift import hadamard_frame_search as frame
from experiments.scalar_lift.large_parity_map import DeterminantSampler, validate_map
from experiments.scalar_lift.quadratic_sign_lift import iter_head_sign_branches, solve_bits, verify_states

SOURCE_NAMES = ('exact_t_campaign.py', 'hadamard_frame_search.py',
                'large_parity_map.py', 'quadratic_sign_lift.py',
                'quadratic_sign_lift_bits.cpp', 'williamson_scalar_endpoints.py',
                'williamson_inputs.json', 'exact_parity_sampler.py', 'exact_parity_sampler.cpp')
STOP = False


def atomic_json(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value, indent=2)+'\n')
    temporary.replace(path)


def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def tuplify(value):
    return tuple(map(tuplify, value)) if isinstance(value, list) else value


def stop_handler(_signum, _frame):
    global STOP
    STOP = True


def branch_count(axes, head_net, first_negative, policy):
    if policy == 'affine':
        return 1
    h = sum(a == 0 for a in axes[1:])
    positives = (h+head_net)//2
    if first_negative == 0:
        positives -= 1
    return math.comb(h-1, positives) if 0 <= positives <= h-1 else 0


def branch_at(axes, head_net, first_negative, policy, index):
    if policy == 'affine':
        return None
    return next(itertools.islice(iter_head_sign_branches(axes, head_net, first_negative),
                                 index, None))


def verified_endpoints(states, data, output):
    """Independent OD/source identities and literal endpoint recipes.

    These are construction certificates. Literature status and new Lean
    literals require a subsequent audit and are deliberately not inferred.
    """
    from experiments.scalar_lift.williamson_scalar_endpoints import (
        verify_orthogonal_design, verified_williamson, ScalarEndpoint, row_controls)
    rows = frame.expand(states, data)
    matrices, od_proof = verify_orthogonal_design(rows)
    widths = {223:(3,5,13), 167:(7,13), 179:(11,13), 739:(1,)}.get(data['p'], (3,))
    endpoints = []
    for width in widths:
        source, proof = verified_williamson(width)
        matrix = ScalarEndpoint(matrices, source)
        endpoints.append({'order':matrix.order, 'Williamson_order':width,
                          'source':proof, 'row_controls':row_controls(matrix),
                          'entry_recipe':'sum_k X[k,i,j] * A[k,a,b], row=(i,a), column=(j,b)',
                          'entry_alphabet':[-1,1], 'Gram_scale':matrix.order,
                          'Gram_proof':'Verified OD polynomial identity and commuting symmetric Williamson identities.',
                          'novelty':'Not assessed by this worker.',
                          'Lean_status':'New literal witness not yet added.'})
    result = {'p':data['p'], 'states':states, 'orthogonal_design':od_proof,
              'verified_real_orders':endpoints}
    atomic_json(output/'endpoints.json', result)
    return result


def run(config, output, seconds=259200, heads=(5,7,9), seed=20261126001,
        branch_policy='all', branch_seconds=120, max_patterns=None, nice=10,
        sampler_engine='python'):
    global STOP
    STOP = False
    if seconds <= 0 or branch_seconds <= 0 or not heads:
        raise ValueError('positive limits and at least one coordinate count are required')
    if branch_policy not in ('all','affine'):
        raise ValueError('branch policy must be all or affine')
    if sampler_engine not in ('python','native'):
        raise ValueError('sampler must be python or native')
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    lock = (output/'worker.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    if nice:
        os.nice(nice)
    data = frame.chart(config['p'], config['m'], config['mu'])
    audit = validate_map(config)
    first_negative = int(config.get('fixed_first_head_negative',False))
    source_hashes = {name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                     for name in SOURCE_NAMES}
    specification = {'config':config, 'heads':list(heads), 'seed':seed,'sampler':sampler_engine,
                     'branch_policy':branch_policy, 'source_sha256':source_hashes}
    fingerprint = hashlib.sha256(json.dumps(specification,sort_keys=True).encode()).hexdigest()
    db = sqlite3.connect(output/'work.sqlite3')
    db.execute('PRAGMA journal_mode=WAL')
    db.execute('PRAGMA synchronous=FULL')
    db.executescript('''
        CREATE TABLE IF NOT EXISTS meta(key TEXT PRIMARY KEY,value TEXT NOT NULL);
        CREATE TABLE IF NOT EXISTS patterns(
          id INTEGER PRIMARY KEY, digest TEXT UNIQUE NOT NULL, axes TEXT NOT NULL,
          total_branches INTEGER NOT NULL,next_branch INTEGER NOT NULL DEFAULT 0,
          skipped_branches INTEGER NOT NULL DEFAULT 0,state TEXT NOT NULL DEFAULT 'active');
        CREATE TABLE IF NOT EXISTS attempts(
          id INTEGER PRIMARY KEY,pattern_id INTEGER NOT NULL,branch_index INTEGER NOT NULL,
          result TEXT NOT NULL,utc TEXT NOT NULL);
    ''')
    def get(key, default=None):
        row = db.execute('SELECT value FROM meta WHERE key=?',(key,)).fetchone()
        return json.loads(row[0]) if row else default
    def put(key, value):
        db.execute('INSERT OR REPLACE INTO meta VALUES (?,?)',(key,json.dumps(value)))
    old = get('fingerprint')
    if old is not None and old != fingerprint:
        raise ValueError('resume requires the identical configuration and frozen sources; use a fresh output directory')
    totals = get('totals', {'patterns_started':0,'patterns_completed':0,'branches_checked':0,
                           'unfinished_attempts':0,'native_affine_assignments_tested':0,
                           'mod8_roots':0,'exact_count_survivors':0,'verified_witnesses':0})
    if sampler_engine == 'native':
        from experiments.scalar_lift.exact_parity_sampler import NativeDeterminantSampler
        sampler = NativeDeterminantSampler(config,seed=seed,workdir=output/'native')
    else:
        sampler = DeterminantSampler(config,seed=seed)
    state = get('sampler_state')
    if state is not None:
        sampler.rng.setstate(tuplify(state))
        sampler.draws = get('sampler_draws',0)
    with db:
        put('fingerprint',fingerprint)
        put('specification',specification)
        put('map_audit',audit)
    previous_handlers = {}
    import threading
    if threading.current_thread() is threading.main_thread():
        for sig in (signal.SIGINT,signal.SIGTERM):
            previous_handlers[sig] = signal.signal(sig,stop_handler)
    started = time.monotonic()
    status = {'schema':'resumable-exact-t-search-v1','pid':os.getpid(),'start_utc':utc(),
              'p':data['p'],'m':data['m'],'mu':data['mu'],'heads':list(heads),
              'method':'SL2 parity coordinates -> exact mod4 elimination -> complete native mod8 filter -> exact counts and integer correlations',
              'branch_policy':branch_policy,'head_schedule':'round_robin','fingerprint':fingerprint,
              'sampler':sampler_engine,
              'budget_seconds_this_session':seconds,'totals':totals,
              'work_database':str(output/'work.sqlite3'),
              'assignment_counter_scope':'Executed assignments, including any replay of interrupted branches; not a distinct-assignment census.'}
    def checkpoint(phase, current=None):
        status.update(phase=phase,updated_utc=utc(),elapsed_seconds_this_session=time.monotonic()-started,
                      sampled_determinant_matrices=sampler.draws,current=current)
        atomic_json(output/'status.json',status)
    checkpoint('running')
    phase = 'budget_exhausted'
    try:
        pending = get('pending_witness')
        if (output/'verified-witness.json').exists():
            pending = json.loads((output/'verified-witness.json').read_text())
        if pending is not None:
            verify_states(pending['states'],data,config['head_net'],config['tail_net'])
            atomic_json(output/'verified-witness.json',pending)
            verified_endpoints(pending['states'],data,output)
            checkpoint('verified_witness')
            return status
        while time.monotonic()-started < seconds:
            if STOP or (output/'STOP').exists():
                phase = 'stopped'
                break
            if max_patterns is not None and totals['patterns_completed'] >= max_patterns:
                phase = 'pattern_limit'
                break
            active = db.execute("SELECT id,axes,total_branches,next_branch,skipped_branches FROM patterns WHERE state='active' ORDER BY id LIMIT 1").fetchone()
            if active is None:
                selected_head = heads[totals['patterns_started'] % len(heads)]
                axes = sampler.next_axes(heads={selected_head},normalize_first_head=True,max_draws=100000)
                if axes is None:
                    with db:
                        put('sampler_state',sampler.rng.getstate())
                        put('sampler_draws',sampler.draws)
                    checkpoint('running')
                    continue
                packed = json.dumps(axes,separators=(',',':'))
                digest = hashlib.sha256(packed.encode()).hexdigest()
                if db.execute('SELECT 1 FROM patterns WHERE digest=?',(digest,)).fetchone():
                    continue
                count = branch_count(axes,config['head_net'],first_negative,branch_policy)
                with db:
                    cursor = db.execute('INSERT INTO patterns(digest,axes,total_branches) VALUES (?,?,?)',
                                        (digest,packed,count))
                    totals['patterns_started'] += 1
                    put('sampler_state',sampler.rng.getstate())
                    put('sampler_draws',sampler.draws)
                    put('totals',totals)
                active = (cursor.lastrowid,packed,count,0,0)
            pattern_id, packed, count, index, skipped = active
            axes = json.loads(packed)
            if index >= count:
                with db:
                    db.execute('UPDATE patterns SET state=? WHERE id=?',
                               ('complete' if not skipped else 'partial',pattern_id))
                    if not skipped:
                        totals['patterns_completed'] += 1
                    put('totals',totals)
                checkpoint('running')
                continue
            remaining = seconds-(time.monotonic()-started)
            if remaining <= 0:
                break
            current = {'pattern_id':pattern_id,'head_count':sum(a==0 for a in axes[1:]),
                       'branch_index':index,'total_branches':count}
            checkpoint('running',current)
            signs = branch_at(axes,config['head_net'],first_negative,branch_policy,index)
            result = solve_bits(axes,data,config['head_net'],config['tail_net'],first_negative,
                                seconds=min(branch_seconds,remaining),workdir=output/'native',fixed_signs=signs)
            complete = result['complete_for_fixed_axes']
            time_left = seconds-(time.monotonic()-started)
            # Keep a deadline-interrupted branch active for exact replay on resume.
            advance = complete or time_left > 0
            with db:
                db.execute('INSERT INTO attempts(pattern_id,branch_index,result,utc) VALUES (?,?,?,?)',
                           (pattern_id,index,json.dumps(result),utc()))
                totals['branches_checked'] += 1
                totals['unfinished_attempts'] += int(not complete)
                totals['native_affine_assignments_tested'] += result.get('native',{}).get('assignments_covered',0)
                totals['mod8_roots'] += result.get('mod8_roots',0)
                totals['exact_count_survivors'] += result.get('exact_count_survivors',0)
                if result['solutions']:
                    totals['verified_witnesses'] += len(result['solutions'])
                    witness = {'p':data['p'],'m':data['m'],'mu':data['mu'],
                               **result['solutions'][0],'fingerprint':fingerprint,'found_utc':utc()}
                    put('pending_witness',witness)
                if advance:
                    db.execute('UPDATE patterns SET next_branch=?,skipped_branches=skipped_branches+? WHERE id=?',
                               (index+1,int(not complete),pattern_id))
                put('totals',totals)
            checkpoint('running',current)
            if result['solutions']:
                atomic_json(output/'verified-witness.json',witness)
                checkpoint('verified_witness',current)
                verified_endpoints(witness['states'],data,output)
                phase = 'verified_witness'
                break
        checkpoint(phase)
        print(json.dumps(status),flush=True)
        return status
    except BaseException as exc:
        status['error'] = f'{type(exc).__name__}: {exc}'
        checkpoint('failed')
        raise
    finally:
        for sig, handler in previous_handlers.items():
            signal.signal(sig,handler)
        db.close()
        lock.close()


def launch(args):
    """Freeze a small self-contained source set before detaching a worker."""
    output = args.output.resolve()
    snapshot = output/'frozen'
    module = snapshot/'experiments'/'scalar_lift'
    if snapshot.exists():
        if json.loads((snapshot/'parity-map.json').read_text()) != json.loads(args.config.read_text()):
            raise ValueError('existing frozen worker has a different parity configuration')
        if not (module/'exact_parity_sampler.py').exists() and args.sampler != 'python':
            raise ValueError('existing frozen worker uses Python sampling; resume its saved launch command')
    else:
        module.mkdir(parents=True)
        (snapshot/'experiments'/'__init__.py').touch()
        (module/'__init__.py').touch()
        for name in SOURCE_NAMES:
            shutil.copy2(Path(__file__).with_name(name),module/name)
        shutil.copy2(args.config,snapshot/'parity-map.json')
        atomic_json(snapshot/'source-manifest.json',{
            str(p.relative_to(snapshot)):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(snapshot.rglob('*')) if p.is_file()})
    command = [sys.executable,'-u','-m','experiments.scalar_lift.exact_t_campaign',
               '--config',str(snapshot/'parity-map.json'),'--output',str(output),
               '--seconds',str(args.seconds),'--branch-seconds',str(args.branch_seconds),
               '--heads',args.heads,'--seed',str(args.seed),'--branch-policy',args.branch_policy,
               '--nice',str(args.nice)]
    if (module/'exact_parity_sampler.py').exists():
        command += ['--sampler',args.sampler]
    if args.max_patterns is not None:
        command += ['--max-patterns',str(args.max_patterns)]
    env = dict(os.environ,PYTHONPATH=str(snapshot))
    with (output/'worker.log').open('ab',buffering=0) as log:
        child = subprocess.Popen(command,cwd=snapshot,env=env,stdin=subprocess.DEVNULL,
                                 stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
    atomic_json(output/'launch.json',{'pid':child.pid,'command':command,'launched_utc':utc()})
    print(json.dumps({'pid':child.pid,'output':str(output),'launch':str(output/'launch.json')}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--heads',default='5,7,9')
    parser.add_argument('--branch-policy',choices=('all','affine'),default='all')
    parser.add_argument('--sampler',choices=('native','python'),default='native')
    parser.add_argument('--seed',type=int,default=20261126001)
    parser.add_argument('--seconds',type=float,default=259200)
    parser.add_argument('--branch-seconds',type=float,default=120)
    parser.add_argument('--nice',type=int,default=10)
    parser.add_argument('--max-patterns',type=int)
    parser.add_argument('--detach',action='store_true')
    args = parser.parse_args()
    if args.detach:
        launch(args)
    else:
        config = json.loads(args.config.read_text())
        run(config,args.output,args.seconds,tuple(map(int,args.heads.split(','))),args.seed,
            args.branch_policy,args.branch_seconds,args.max_patterns,args.nice,args.sampler)


if __name__ == '__main__':
    main()

"""Checkpointed exact searches with identity multiplier action on coordinates.

Only completed fixed-pattern, row-sum and head-sign branches are excluded.
Partial native scans remain pending for replay on resume. A positive result
is reconstructed at every actual cyclic shift before its durable checkpoint.
"""
from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import sqlite3
import subprocess
import sys
import threading
import time

from experiments.scalar_lift import identity_multiplier as identity
from experiments.scalar_lift.exact_t_campaign import (
    atomic_json, branch_at, branch_count, tuplify, utc, verified_endpoints,
)
from experiments.scalar_lift.large_parity_map import DeterminantSampler

SOURCE_NAMES = (
    'identity_t_campaign.py', 'identity_multiplier.py', 'exact_t_campaign.py',
    'hadamard_frame_search.py', 'large_parity_map.py', 'quadratic_sign_lift.py',
    'quadratic_sign_lift_bits.cpp', 'exact_parity_sampler.py', 'exact_parity_sampler.cpp',
    'williamson_scalar_endpoints.py', 'williamson_inputs.json',
)
STOP = False


def stop_handler(_signum, _frame):
    global STOP
    STOP = True


def strata(p, m, row_targets=None, head_offsets=(0, 2, 4)):
    """Interleave row types at each increase in the minimum head support."""
    rows = identity.target_rows(p, m) if row_targets is None else row_targets
    rows = [tuple(row) for row in rows]
    if not rows or len(rows) != len(set(rows)):
        raise ValueError('distinct compatible target rows are required')
    if not head_offsets or len(set(head_offsets)) != len(head_offsets) or any(
            type(x) is not int or x < 0 or x % 2 for x in head_offsets):
        raise ValueError('distinct nonnegative even head offsets are required')
    for row in rows:
        if (len(row) != 4 or any(type(x) is not int for x in row)
                or sum(x*x for x in row) != p or (row[0]-1) % m
                or any(x % m for x in row[1:])):
            raise ValueError('target row does not have the required norm or orbit congruences')
    count = (p-1)//m
    result = []
    for offset in head_offsets:
        for row in rows:
            head_net = (row[0]-1)//m
            heads = abs(head_net)+offset
            if heads > count or heads+sum(abs(x//m) for x in row[1:]) > count:
                continue
            result.append({'target_row': list(row), 'head_count': heads,
                           'head_net': head_net,
                           'first_negative': (int(head_net < 0) if heads else None)})
    if not result:
        raise ValueError('no head stratum can realize the requested target rows')
    return result


def make_sampler(config, seed, engine, directory):
    if engine == 'native':
        from experiments.scalar_lift.exact_parity_sampler import NativeDeterminantSampler
        return NativeDeterminantSampler(config, seed=seed, workdir=directory)
    return DeterminantSampler(config, seed=seed)


def run(config, output, seconds=259200, row_targets=None, head_offsets=(0, 2, 4),
        seed=20261127223001, branch_policy='all', branch_seconds=120,
        max_patterns=None, nice=10, sampler_engine='native', heartbeat_seconds=5,
        sampler_batch_draws=2048):
    """Run one locked worker; seconds is the cumulative active-session budget."""
    global STOP
    STOP = False
    if (seconds <= 0 or branch_seconds <= 0 or heartbeat_seconds <= 0
            or type(sampler_batch_draws) is not int or sampler_batch_draws <= 0
            or (max_patterns is not None and (type(max_patterns) is not int or max_patterns < 0))):
        raise ValueError('positive resource bounds and a nonnegative pattern limit are required')
    if branch_policy not in ('all', 'affine') or sampler_engine not in ('native', 'python'):
        raise ValueError('invalid branch policy or sampler engine')
    if config.get('action') != 'identity':
        raise ValueError('an explicitly validated identity-action map is required')
    data = identity.chart(config['p'], config['m'], config['mu'])
    schedule = strata(data['p'], data['m'], row_targets, head_offsets)
    audit = identity.validate_map(config)
    hashes = {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
              for name in SOURCE_NAMES}
    specification = {'config': config, 'schedule': schedule, 'seed': seed,
                     'sampler': sampler_engine, 'branch_policy': branch_policy,
                     'source_sha256': hashes}
    fingerprint = hashlib.sha256(json.dumps(specification, sort_keys=True).encode()).hexdigest()
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    lock = (output/'worker.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    if nice:
        os.nice(nice)
    db = sqlite3.connect(output/'work.sqlite3')
    db.execute('PRAGMA journal_mode=WAL')
    db.execute('PRAGMA synchronous=FULL')
    db.executescript('''
        CREATE TABLE IF NOT EXISTS meta(key TEXT PRIMARY KEY,value TEXT NOT NULL);
        CREATE TABLE IF NOT EXISTS patterns(
          id INTEGER PRIMARY KEY,digest TEXT UNIQUE NOT NULL,axes TEXT NOT NULL,
          stratum INTEGER NOT NULL,total_branches INTEGER NOT NULL,
          next_branch INTEGER NOT NULL DEFAULT 0,state TEXT NOT NULL DEFAULT 'active');
        CREATE TABLE IF NOT EXISTS branches(
          pattern_id INTEGER NOT NULL,branch_index INTEGER NOT NULL,complete INTEGER NOT NULL,
          PRIMARY KEY(pattern_id,branch_index));
        CREATE TABLE IF NOT EXISTS attempts(
          id INTEGER PRIMARY KEY,pattern_id INTEGER NOT NULL,branch_index INTEGER NOT NULL,
          result TEXT NOT NULL,utc TEXT NOT NULL);
    ''')
    def get(key, default=None):
        record = db.execute('SELECT value FROM meta WHERE key=?', (key,)).fetchone()
        return json.loads(record[0]) if record else default
    def put(key, value):
        db.execute('INSERT OR REPLACE INTO meta VALUES (?,?)', (key, json.dumps(value)))
    old = get('fingerprint')
    if old is not None and old != fingerprint:
        db.close();lock.close()
        raise ValueError('resume requires identical configuration, schedule and frozen sources')
    totals = get('totals', dict.fromkeys((
        'patterns_started', 'patterns_completed',
        'branch_attempts', 'branches_completed', 'unfinished_attempts',
        'native_affine_assignments_tested', 'mod8_roots', 'exact_count_survivors',
        'verified_witnesses'), 0))
    sampler = make_sampler(config, seed, sampler_engine, output/'native')
    if get('sampler_state') is not None:
        sampler.rng.setstate(tuplify(get('sampler_state')))
        sampler.draws = get('sampler_draws', 0)
    elapsed_before = get('elapsed_seconds', 0.)
    with db:
        put('fingerprint', fingerprint);put('specification', specification);put('map_audit', audit)
    handlers = {}
    if threading.current_thread() is threading.main_thread():
        for sig in (signal.SIGINT, signal.SIGTERM):
            handlers[sig] = signal.signal(sig, stop_handler)
    started = time.monotonic()
    status = {'schema': 'resumable-identity-t-search-v1', 'pid': os.getpid(), 'start_utc': utc(),
              'p': data['p'], 'm': data['m'], 'mu': data['mu'], 'action': 'identity',
              'schedule': schedule, 'sampler': sampler_engine, 'branch_policy': branch_policy,
              'fingerprint': fingerprint, 'budget_seconds': seconds,
              'budget_seconds_this_session': max(0, seconds-elapsed_before),
              'work_database': str(output/'work.sqlite3'),
              'method': 'Identity-orbit SL2 parity -> exact mod4 elimination -> native mod8 enumeration -> four exact counts -> all integer shifts',
              'assignment_counter_scope': 'Returned native assignments; repeated partial branches count again, uncommitted killed calls do not count.',
              'exclusion_scope': 'Only completed fixed axes, target row and stated first-head/sign-branch constraints; no global nonexistence claim.'}
    status_lock, heartbeat_stop = threading.Lock(), threading.Event()
    def checkpoint(phase, current=None):
        elapsed = elapsed_before+time.monotonic()-started
        with db:
            put('elapsed_seconds', elapsed)
        with status_lock:
            status.update(phase=phase, current=current, updated_utc=utc(),
                          elapsed_seconds=elapsed, elapsed_seconds_this_session=time.monotonic()-started,
                          sampled_determinant_matrices=sampler.draws, totals=dict(totals),
                          pending_branches=db.execute('SELECT count(*) FROM branches WHERE complete=0').fetchone()[0],
                          pending_patterns=db.execute("SELECT count(*) FROM patterns WHERE state='pending'").fetchone()[0])
            atomic_json(output/'status.json', status)
    def heartbeat():
        while not heartbeat_stop.wait(heartbeat_seconds):
            with status_lock:
                if status.get('phase') == 'running':
                    status.update(updated_utc=utc(), elapsed_seconds=elapsed_before+time.monotonic()-started,
                                  elapsed_seconds_this_session=time.monotonic()-started)
                    atomic_json(output/'status.json', status)
    checkpoint('running')
    thread = threading.Thread(target=heartbeat, daemon=True)
    thread.start()
    attempted = set()
    phase = 'budget_exhausted'
    def save_sampler():
        put('sampler_state', sampler.rng.getstate());put('sampler_draws', sampler.draws)
    def accept_witness(witness):
        witness = dict(witness)
        witness['expanded_verification'] = identity.verify_states(witness['states'], data, witness['target_row'])
        atomic_json(output/'verified-witness.json', witness)
        checkpoint('source_verified_endpoint_pending')
        verified_endpoints(witness['states'], data, output)
        checkpoint('verified_witness')
    try:
        pending = get('pending_witness')
        if pending is not None:
            accept_witness(pending)
            return status
        while elapsed_before+time.monotonic()-started < seconds:
            if STOP or (output/'STOP').exists():
                phase = 'stopped';break
            if max_patterns is not None and totals['patterns_completed'] >= max_patterns:
                phase = 'pattern_limit';break
            # Retry each old partial branch at most once this session. Further
            # timeouts remain pending; they never become complete exclusions.
            retry = next((row for row in db.execute('''
                SELECT p.id,p.axes,p.stratum,p.total_branches,b.branch_index FROM branches b
                JOIN patterns p ON p.id=b.pattern_id WHERE b.complete=0 ORDER BY p.id,b.branch_index''')
                          if (row[0], row[4]) not in attempted), None)
            active = retry or db.execute('''
                SELECT id,axes,stratum,total_branches,next_branch FROM patterns
                WHERE state='active' AND next_branch<total_branches ORDER BY id LIMIT 1''').fetchone()
            if active is None:
                stratum_index = totals['patterns_started'] % len(schedule)
                selected = schedule[stratum_index]
                axes = sampler.next_axes(heads={selected['head_count']}, normalize_first_head=False,
                                         max_draws=sampler_batch_draws)
                if axes is None:
                    with db:
                        save_sampler()
                    checkpoint('running');continue
                if selected['head_count']:
                    axes = identity.normalize_head(axes, data)
                if len(axes) != data['n'] or axes[0] != 0 or axes[1:].count(0) != selected['head_count']:
                    raise ArithmeticError('sampler violated the identity-orbit support filter')
                packed = json.dumps(axes, separators=(',', ':'))
                key = hashlib.sha256(json.dumps([selected['target_row'], axes], separators=(',', ':')).encode()).hexdigest()
                if db.execute('SELECT 1 FROM patterns WHERE digest=?', (key,)).fetchone():
                    with db:
                        save_sampler()
                    checkpoint('running');continue
                count = (branch_count(axes, selected['head_net'], selected['first_negative'], branch_policy)
                         if selected['head_count'] else 1)
                with db:
                    cursor = db.execute('INSERT INTO patterns(digest,axes,stratum,total_branches) VALUES (?,?,?,?)',
                                        (key, packed, stratum_index, count))
                    totals['patterns_started'] += 1;put('totals', totals);save_sampler()
                active = (cursor.lastrowid, packed, stratum_index, count, 0)
            pattern, packed, stratum_index, count, index = active
            selected = schedule[stratum_index]
            axes = json.loads(packed)
            fixed = ({0: 0} if not selected['head_count'] else
                     branch_at(axes, selected['head_net'], selected['first_negative'], branch_policy, index))
            remaining = seconds-elapsed_before-(time.monotonic()-started)
            if remaining <= 0:
                break
            current = {'pattern_id': pattern, 'branch_index': index, 'total_branches': count,
                       'head_count': selected['head_count'], 'target_row': selected['target_row']}
            checkpoint('running', current)
            result = identity.solve_bits(axes, data, selected['target_row'],
                                         first_negative=selected['first_negative'], fixed_signs=fixed,
                                         seconds=min(branch_seconds, remaining), workdir=output/'native')
            complete = result['complete_for_fixed_axes'] is True
            solutions = []
            for candidate in result['solutions']:
                proof = identity.verify_states(candidate['states'], data, selected['target_row'])
                solutions.append({**candidate, 'expanded_verification': proof})
            witness = None
            if solutions:
                witness = {'p': data['p'], 'm': data['m'], 'mu': data['mu'], 'action': 'identity',
                           **solutions[0], 'target_row': selected['target_row'], 'fingerprint': fingerprint,
                           'found_utc': utc(), 'pattern_id': pattern, 'branch_index': index}
            attempted.add((pattern, index))
            with db:
                db.execute('INSERT INTO attempts(pattern_id,branch_index,result,utc) VALUES (?,?,?,?)',
                           (pattern, index, json.dumps(result), utc()))
                db.execute('INSERT OR REPLACE INTO branches VALUES (?,?,?)', (pattern, index, int(complete)))
                db.execute('UPDATE patterns SET next_branch=MAX(next_branch,?) WHERE id=?', (index+1, pattern))
                totals['branch_attempts'] += 1;totals['branches_completed'] += int(complete)
                totals['unfinished_attempts'] += int(not complete)
                totals['native_affine_assignments_tested'] += result.get('native', {}).get('assignments_covered', 0)
                totals['mod8_roots'] += result.get('mod8_roots', 0)
                totals['exact_count_survivors'] += result.get('exact_count_survivors', 0)
                if witness is not None:
                    totals['verified_witnesses'] += len(solutions);put('pending_witness', witness)
                next_branch = db.execute('SELECT next_branch FROM patterns WHERE id=?', (pattern,)).fetchone()[0]
                if next_branch >= count:
                    unfinished = db.execute('SELECT count(*) FROM branches WHERE pattern_id=? AND complete=0', (pattern,)).fetchone()[0]
                    state = 'pending' if unfinished else 'complete'
                    db.execute('UPDATE patterns SET state=? WHERE id=?', (state, pattern))
                    if state == 'complete':
                        totals['patterns_completed'] += 1
                put('totals', totals)
            checkpoint('running', current)
            if witness is not None:
                accept_witness(witness)
                phase = 'verified_witness';break
        checkpoint(phase)
        return status
    except BaseException as exc:
        status['error'] = f'{type(exc).__name__}: {exc}'
        checkpoint('failed')
        raise
    finally:
        heartbeat_stop.set();thread.join(timeout=heartbeat_seconds+1)
        for sig, handler in handlers.items():
            signal.signal(sig, handler)
        db.close();lock.close()


def check_snapshot(snapshot):
    snapshot = Path(snapshot).resolve()
    expected = json.loads((snapshot/'source-manifest.json').read_text())
    actual = {str(p.relative_to(snapshot)) for p in snapshot.rglob('*')
              if p.is_file() and p.name != 'source-manifest.json'}
    if actual != set(expected):
        raise ValueError('frozen source file set differs from its manifest')
    for name, digest in expected.items():
        path = snapshot/name
        if path.is_symlink() or not path.resolve().is_relative_to(snapshot) or hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError('frozen source hash mismatch: '+name)


def launch(args):
    """Freeze explicit public sources; parent process performs no search."""
    output = args.output.resolve();output.mkdir(parents=True, exist_ok=True)
    snapshot = output/'frozen';module = snapshot/'experiments/scalar_lift'
    config = json.loads(args.config.read_text())
    if config.get('action') != 'identity':
        raise ValueError('an identity-action parity map is required')
    if snapshot.exists():
        check_snapshot(snapshot)
        if json.loads((snapshot/'parity-map.json').read_text()) != config:
            raise ValueError('existing frozen worker has a different identity configuration')
    else:
        module.mkdir(parents=True)
        (snapshot/'experiments/__init__.py').touch();(module/'__init__.py').touch()
        for name in SOURCE_NAMES:
            shutil.copy2(Path(__file__).with_name(name), module/name)
        atomic_json(snapshot/'parity-map.json', config)
        atomic_json(snapshot/'source-manifest.json', {
            str(p.relative_to(snapshot)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(snapshot.rglob('*')) if p.is_file()})
        for path in snapshot.rglob('*'):
            if path.is_file():path.chmod(0o444)
    command = [sys.executable, '-B', '-u', '-m', 'experiments.scalar_lift.identity_t_campaign',
               '--config', str(snapshot/'parity-map.json'), '--output', str(output),
               '--seconds', str(args.seconds), '--branch-seconds', str(args.branch_seconds),
               '--head-offsets', args.head_offsets, '--seed', str(args.seed),
               '--branch-policy', args.branch_policy, '--sampler', args.sampler, '--nice', str(args.nice)]
    for row in getattr(args, 'row_targets', None) or []:
        command.append('--target-row='+','.join(map(str, row)))
    if args.max_patterns is not None:command += ['--max-patterns', str(args.max_patterns)]
    env = dict(os.environ, PYTHONPATH=str(snapshot), PYTHONDONTWRITEBYTECODE='1',
               OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1')
    with (output/'worker.log').open('ab', buffering=0) as log:
        child = subprocess.Popen(command, cwd=snapshot, env=env, stdin=subprocess.DEVNULL,
                                 stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
    public_env = {key: env[key] for key in ('PYTHONPATH', 'PYTHONDONTWRITEBYTECODE',
                  'OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS')}
    atomic_json(output/'launch.json', {'pid': child.pid, 'command': command, 'cwd': str(snapshot),
                                      'environment': public_env, 'launched_utc': utc()})
    return {'pid': child.pid, 'output': str(output), 'launch': str(output/'launch.json')}


def parse_row(value):
    try:
        row = tuple(map(int, value.split(',')))
    except ValueError as exc:
        raise argparse.ArgumentTypeError('target row must contain four comma-separated integers') from exc
    if len(row) != 4:
        raise argparse.ArgumentTypeError('target row must contain four comma-separated integers')
    return row


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--head-offsets', default='0,2,4')
    parser.add_argument('--target-row', type=parse_row, action='append', dest='row_targets',
                        help='four comma-separated integer sums; repeat for several rows (default: all compatible rows)')
    parser.add_argument('--branch-policy', choices=('all', 'affine'), default='all')
    parser.add_argument('--sampler', choices=('native', 'python'), default='native')
    parser.add_argument('--seed', type=int, default=20261127223001)
    parser.add_argument('--seconds', type=float, default=259200)
    parser.add_argument('--branch-seconds', type=float, default=120)
    parser.add_argument('--nice', type=int, default=10)
    parser.add_argument('--max-patterns', type=int)
    parser.add_argument('--detach', action='store_true')
    # argparse treats a negative comma-separated tuple as an option unless it
    # is attached with '='. Accept both spellings and always freeze the latter.
    raw = list(sys.argv[1:] if argv is None else argv)
    normalized, index = [], 0
    while index < len(raw):
        if raw[index] == '--target-row' and index+1 < len(raw) and ',' in raw[index+1]:
            normalized.append('--target-row='+raw[index+1]);index += 2
        else:
            normalized.append(raw[index]);index += 1
    return parser.parse_args(normalized)


def main():
    args = parse_args()
    if args.detach:
        result = launch(args)
    else:
        frozen = Path(__file__).resolve().parents[2]
        if (frozen/'source-manifest.json').exists():check_snapshot(frozen)
        result = run(json.loads(args.config.read_text()), args.output, seconds=args.seconds,
                     row_targets=args.row_targets,
                     head_offsets=tuple(map(int, args.head_offsets.split(','))), seed=args.seed,
                     branch_policy=args.branch_policy, branch_seconds=args.branch_seconds,
                     max_patterns=args.max_patterns, nice=args.nice, sampler_engine=args.sampler)
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()

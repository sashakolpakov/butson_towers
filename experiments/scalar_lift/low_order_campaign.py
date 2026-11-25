"""Bounded, checkpointed low-order T-matrix searches from exact parity maps.

All candidate matrices are checked over every expanded shift. A completed
fixed-coordinate search is not a nonexistence result for its target order.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time

from experiments.scalar_lift.hadamard_frame_search import chart
from experiments.scalar_lift.large_parity_map import DeterminantSampler, validate_map
from experiments.scalar_lift.quadratic_sign_lift import solve_bits, iter_head_sign_branches


def write_record(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value, indent=2)+'\n')
    temporary.replace(path)


def run(config, output, seconds=600, pattern_seconds=15, heads=(5,), seed=20261125223,
        max_patterns=None):
    if seconds <= 0 or pattern_seconds <= 0 or not heads:
        raise ValueError('positive time limits and a nonempty head-count selection are required')
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    data = chart(config['p'], config['m'], config['mu'])
    map_check = validate_map(config)
    sampler = DeterminantSampler(config, seed=seed)
    source_paths = [Path(__file__), Path(__file__).with_name('quadratic_sign_lift.py'),
                    Path(__file__).with_name('quadratic_sign_lift_bits.cpp'),
                    Path(__file__).with_name('large_parity_map.py')]
    report = {
        'schema':'bounded-low-order-t-matrix-search-v1',
        'target_order':12*data['p'], 'p':data['p'], 'm':data['m'], 'mu':data['mu'],
        'head_net':config['head_net'], 'tail_net':config['tail_net'],
        'first_negative':int(config.get('fixed_first_head_negative',False)),
        'seed':seed, 'heads':list(heads), 'wall_limit_seconds':seconds,
        'per_pattern_native_limit_seconds':pattern_seconds,
        'map_validation':map_check, 'source_sha256':{
            p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths},
        'patterns':[], 'solutions':[], 'finished':False,
        'scope':'The sampled coordinate patterns and stated sign/count gauges only.',
    }
    seen = set()
    while time.monotonic()-started < seconds and (max_patterns is None or len(report['patterns']) < max_patterns):
        axes = sampler.next_axes(heads=set(heads), normalize_first_head=True, max_draws=100000)
        if axes is None or tuple(axes) in seen:
            continue
        seen.add(tuple(axes))
        index = len(report['patterns'])
        remaining = seconds-(time.monotonic()-started)
        if remaining <= 0:
            break
        branches = list(iter_head_sign_branches(axes, config['head_net'], report['first_negative']))
        branch_results = []
        for fixed_signs in branches:
            remaining = seconds-(time.monotonic()-started)
            if remaining <= 0:
                break
            branch = solve_bits(axes, data, config['head_net'], config['tail_net'],
                                report['first_negative'], min(pattern_seconds, remaining),
                                workdir=output/'native', fixed_signs=fixed_signs)
            branch_results.append(branch)
            if branch['solutions']:
                break
        complete = len(branch_results) == len(branches) and all(r['complete_for_fixed_axes'] for r in branch_results)
        result = {
            'solver_status':'COMPLETE' if complete else 'UNFINISHED',
            'complete_for_fixed_axes':complete,
            'head_sign_branches':len(branches), 'completed_head_branches':len(branch_results),
            'branches':branch_results,
            'mod8_roots':sum(r.get('mod8_roots',0) for r in branch_results),
            'exact_count_survivors':sum(r.get('exact_count_survivors',0) for r in branch_results),
            'affine_assignments_covered':sum(r.get('native',{}).get('assignments_covered',0) for r in branch_results),
            'solutions':[solution for r in branch_results for solution in r['solutions']],
        }
        record = {'index':index, 'head_count':sum(a==0 for a in axes[1:]),
                  'axes':axes, 'result':result}
        report['patterns'].append(record)
        report['solutions'].extend(result['solutions'])
        report['elapsed_seconds'] = time.monotonic()-started
        report['sampled_determinant_matrices'] = sampler.draws
        report['completed_fixed_patterns'] = sum(p['result']['complete_for_fixed_axes'] for p in report['patterns'])
        report['mod8_roots'] = sum(p['result'].get('mod8_roots',0) for p in report['patterns'])
        report['count_survivors'] = sum(p['result'].get('exact_count_survivors',0) for p in report['patterns'])
        report['affine_assignments_covered'] = sum(p['result'].get('affine_assignments_covered',0) for p in report['patterns'])
        write_record(output/'campaign.json',report)
        if result['solutions'] or index%10==0:
            print(json.dumps({k:report.get(k) for k in ('target_order','elapsed_seconds','completed_fixed_patterns','mod8_roots','count_survivors','affine_assignments_covered')}|{'solutions':len(report['solutions'])}),flush=True)
        if result['solutions']:
            write_record(output/'verified-witness.json',{
                'p':data['p'],'m':data['m'],'mu':data['mu'],
                'states':result['solutions'][0]['states'],
                'expanded_verification':result['solutions'][0]['expanded_verification'],
            })
            break
    report['finished'] = True
    report['elapsed_seconds'] = time.monotonic()-started
    report['sampled_determinant_matrices'] = sampler.draws
    write_record(output/'campaign.json',report)
    return report


def run_projection(config, output, seconds=3600, trial_seconds=60, seed=22320261125,
                   initial=None):
    """Cycle through spectral searches; retain and exactly verify the best word.

    Trials alternate short and long restart intervals. A retained nonzero word
    is search state only. The separate witness file requires every exact
    expanded correlation to vanish.
    """
    from experiments.scalar_lift import t_sequence_projection as projection
    if seconds <= 0 or trial_seconds <= 0:
        raise ValueError('positive campaign and trial durations are required')
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    assignment, library = projection.build_assignment(output/'native')
    controls = projection.controls(assignment)
    data = chart(config['p'], config['m'], config['mu'])
    target = projection.Projections(data, None, config['head_net'], config['tail_net'],
                                    assignment).target
    best = None
    if initial is not None:
        states = json.loads(Path(initial).read_text())['states']
        best = projection.verify(states, data, target)
        write_record(output/'best.json', best)
    source_paths = [Path(__file__), Path(projection.__file__),
                    Path(projection.__file__).with_suffix('.cpp'),
                    Path(__file__).with_name('hadamard_frame_search.py')]
    report = {
        'schema':'bounded-low-order-projection-search-v1',
        'target_order':12*data['p'], 'p':data['p'], 'm':data['m'], 'mu':data['mu'],
        'head_net':config['head_net'], 'tail_net':config['tail_net'],
        'seed':seed, 'wall_limit_seconds':seconds, 'trial_seconds':trial_seconds,
        'controls':controls, 'source_sha256':{
            p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths},
        'trials':[], 'finished':False,
        'scope':'Heuristic search only; every retained best is verified with integers.',
    }
    profiles = [(None, beta, orientation)
                for orientation in ('discrete-first', 'spectral-first')
                for beta in (.2, .35, .5, .65, .8)]
    profiles += [(5, .35, 'discrete-first'), (5, .5, 'spectral-first')]
    started = time.monotonic()
    while time.monotonic()-started < seconds:
        index = len(report['trials'])
        heads, beta, orientation = profiles[index % len(profiles)]
        restart_after = 2000 if (index//len(profiles)) % 2 == 0 else 100000
        warm = None
        if index % 5 == 0 and best is not None and (heads is None or
                sum(s//2 == 0 for s in best['states'][1:]) == heads):
            warm = best['states']
        remaining = seconds-(time.monotonic()-started)
        if remaining <= 0:
            break
        trial = projection.search(
            data, heads, config['head_net'], config['tail_net'], assignment,
            min(trial_seconds, remaining), seed+index, beta,
            output/f'trial-{index:05d}.json', restart_after, warm, orientation)
        if best is None or trial['best_energy'] < best['best_energy']:
            best = trial
            write_record(output/'best.json', best)
        report['trials'].append({
            'index':index, 'heads':heads, 'beta':beta, 'orientation':orientation,
            'restart_after':restart_after, 'seed':seed+index,
            'initialization':'previous_best' if warm is not None else 'fresh_random',
            'best_energy':trial['best_energy'], 'exact_T_design':trial['exact_T_design'],
            'iterations':trial['total_iterations'], 'seconds':trial['total_seconds'],
        })
        report['elapsed_seconds'] = time.monotonic()-started
        report['best_energy'] = best['best_energy']
        report['exact_T_design'] = best['exact_T_design']
        write_record(output/'campaign.json', report)
        print(json.dumps({'trial':index, 'target_order':report['target_order'],
                          'best_energy':best['best_energy'],
                          'exact_T_design':best['exact_T_design']}), flush=True)
        if best['exact_T_design']:
            write_record(output/'verified-witness.json', best)
            break
    report['finished'] = True
    report['elapsed_seconds'] = time.monotonic()-started
    write_record(output/'campaign.json', report)
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--seconds',type=float,default=600)
    parser.add_argument('--pattern-seconds',type=float,default=15)
    parser.add_argument('--heads',default='5')
    parser.add_argument('--seed',type=int,default=20261125223)
    parser.add_argument('--max-patterns',type=int)
    parser.add_argument('--engine',choices=('signs','projection'),default='signs')
    parser.add_argument('--trial-seconds',type=float,default=60)
    parser.add_argument('--initial',type=Path)
    args=parser.parse_args()
    config=json.loads(args.config.read_text())
    if args.engine == 'projection':
        report=run_projection(config,args.output,args.seconds,args.trial_seconds,
                              args.seed,args.initial)
        print(json.dumps({k:report.get(k) for k in
            ('finished','best_energy','exact_T_design','elapsed_seconds')}))
    else:
        report=run(config,args.output,args.seconds,args.pattern_seconds,
                   tuple(map(int,args.heads.split(','))),args.seed,args.max_patterns)
        print(json.dumps({'finished':report['finished'],'solutions':len(report['solutions']),
                          'patterns':len(report['patterns']),'elapsed_seconds':report['elapsed_seconds']}))


if __name__=='__main__':main()

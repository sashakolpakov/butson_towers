"""Exact CP-SAT for a selected symmetric Williamson multiplier chart."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from experiments.scalar_lift.williamson_symmetry import chart, rotate, coefficients, expand, verify_source, H


def _solver_module():
    try:
        from ortools.sat.python import cp_model
    except ImportError as exc:
        raise RuntimeError("CP-SAT requires the optional search dependencies: pip install '.[search]'") from exc
    return cp_model


def integer_nets(target, data):
    """Convert four exact source row sums into signed orbit counts."""
    target = np.asarray(target)
    if target.shape != (4,) or not np.issubdtype(target.dtype,np.integer):
        raise ValueError('Four exact integer row sums are required')
    target = target.astype(np.int64)
    if int(target @ target) != 4*data['p']:
        raise ValueError('Four requested row sums must have squared norm 4p')
    numerator = H.T @ target - 2
    denominator = 4*data['m']
    if data['action'] == 'C3':
        if not np.all(numerator[1:] == numerator[1]):
            raise ValueError('C3 row sums require equal transformed tails')
        numerator[1:] *= 3
    if np.any(numerator % denominator):
        raise ValueError('Row sums contradict orbit congruences')
    nets = (numerator // denominator).tolist()
    return nets


def build(data, target=None, gauge=True):
    """Encode every integer correlation and the forced-sign XOR identities.

    The returned axis and negative-sign variables describe complete orbit
    states. Optional row sums and the residue-scaling gauge remain explicit
    in the returned model description.
    """
    cp_model = _solver_module()
    started = time.monotonic()
    n = data['N']
    model = cp_model.CpModel()
    axes = [model.new_int_var(0, 3, f'axis{i}') for i in range(n)]
    negative = [model.new_bool_var(f'neg{i}') for i in range(n)]
    head = [model.new_bool_var(f'head{i}') for i in range(n)]
    for i in range(n):
        model.add(axes[i] == 0).only_enforce_if(head[i])
        model.add(axes[i] != 0).only_enforce_if(head[i].Not())
    rotated = {(i, 0): axes[i] for i in range(n)}
    if data['action'] == 'C3':
        for i in range(n):
            for phase in (1, 2):
                value = model.new_int_var(0, 3, f'rot{i}_{phase}')
                model.add_element(axes[i], [rotate(j, phase) for j in range(4)], value)
                rotated[i, phase] = value
    pairs = sorted({(i,j,d) for row in data['rows'] for i,j,d in row if i != j})
    differences, equal, both = {}, {}, {}
    for i,j,d in pairs:
        if (i,j) not in differences:
            value = model.new_bool_var(f'diff{i}_{j}')
            model.add_bool_xor([negative[i], negative[j], value.Not()])
            differences[i,j] = value
        match = model.new_bool_var(f'equal{i}_{j}_{d}')
        model.add(axes[i] == rotated[j,d]).only_enforce_if(match)
        model.add(axes[i] != rotated[j,d]).only_enforce_if(match.Not())
        cut = model.new_bool_var(f'cut{i}_{j}_{d}')
        model.add(cut <= match)
        model.add(cut <= differences[i,j])
        model.add(cut >= match + differences[i,j] - 1)
        equal[i,j,d], both[i,j,d] = match, cut
    for ell,row in enumerate(data['rows']):
        constant = 0
        diagonal = [0]*n
        linear = []
        parity = [negative[ell]]
        for (i,j,d),count in row.items():
            if i == j:
                if d == 0:
                    constant += count
                else:
                    diagonal[i] += count
            else:
                if count % 2:
                    raise ValueError('Inversion-paired offdiagonal coefficient must be even')
                linear.append(count*(equal[i,j,d] - 2*both[i,j,d]))
                if (count//2) % 2:
                    parity.append(equal[i,j,d])
        if constant % 2 != 1 or any(value % 2 for value in diagonal):
            raise ValueError('The sign-elimination coefficient parity failed')
        for i,count in enumerate(diagonal):
            if count:
                linear.append(count*head[i])
                if (count//2) % 2:
                    parity.append(head[i])
        model.add(constant + sum(linear) + 1 - 2*negative[ell] == 0)
        # xor(parity) = (constant+1)/2 mod 2. CP-SAT requires odd parity.
        wanted = ((constant+1)//2) % 2
        model.add_bool_xor(parity + [bool(1-wanted)])
    nets = integer_nets(target, data) if target is not None else None
    if nets is not None:
        coordinates = []
        for i in range(n):
            state = model.new_int_var(0, 7, f'state{i}')
            model.add(state == 2*axes[i] + negative[i])
            values = []
            for j in range(4):
                value = model.new_int_var(-1, 1, f'value{i}_{j}')
                model.add_element(state, [int(k//2 == j)*(1-2*(k%2)) for k in range(8)], value)
                values.append(value)
            coordinates.append(values)
        model.add(sum(row[0] for row in coordinates) == nets[0])
        if data['action'] == 'identity':
            for j in (1,2,3):
                model.add(sum(row[j] for row in coordinates) == nets[j])
        else:
            model.add(sum(sum(row[1:]) for row in coordinates) == nets[1])
    if gauge:
        # Every source with p>1 has a negative orbit sign: all-positive
        # signs give normalized correlation 1+A>0. Residue scaling sends
        # one negative orbit to residue1 and preserves multiplier symmetry.
        model.add(negative[0] == 1)
        if nets is not None and nets[0] < 0:
            model.add(axes[0] == 0)
        elif data['action'] == 'C3':
            if nets is not None and nets[1] < 0:
                model.add(axes[0] == 1)
            else:
                model.add(axes[0] <= 1)
        elif nets is None:
            model.add(axes[0] == 0)  # Full S4 coordinate permutation.
    failure = model.validate()
    if failure:
        raise ValueError(failure)
    return model, axes, negative, {'model_seconds':time.monotonic()-started,
                                   'orbit_count':n,'pair_variables':len(pairs),
                                   'target_orbit_nets':nets,'gauges_enabled':gauge}


def solve(p, m, action, *, seconds=10, seed=1, target=None, gauge=True,
          output=None, hint=None):
    """Run one single-thread call and independently verify any literal source.

    Solver status applies to the supplied chart, row sums and gauges.
    A positive result contains all four sign rows and expanded verification.
    """
    cp_model = _solver_module()
    data = chart(p,m,action)
    model, axes, negative, description = build(data,target,gauge)
    if hint is not None:
        for variable,value in zip(axes,hint['axes'],strict=True):
            model.add_hint(variable,int(value))
        for variable,value in zip(negative,hint['signs'],strict=True):
            model.add_hint(variable,int(value < 0))
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = 1
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.random_seed = seed
    solver.parameters.randomize_search = True
    status = solver.solve(model)
    result = {'p':p,'m':m,'mu':data['mu'],'action':action,'backend':'CP-SAT',
              'status':solver.status_name(status),'seed':seed,'seconds_limit':seconds,
              'wall_seconds':solver.wall_time,'branches':solver.num_branches,
              'conflicts':solver.num_conflicts,'model':description,
              'target_row':target,'scope':'Selected dense symmetric Williamson multiplier chart only.',
              'source_claimed':False}
    if status in (cp_model.OPTIMAL,cp_model.FEASIBLE):
        color = [solver.value(v) for v in axes]
        signs = np.asarray([1-2*solver.value(v) for v in negative],dtype=np.int64)
        _,_,_,forced = coefficients(color,data)
        if not np.array_equal(signs,forced):
            raise ArithmeticError('CP signs disagree with independent forced-sign formula')
        rows = expand(color,signs,data)
        proof = verify_source(rows,data)
        if target is not None and proof['row_sums'] != target:
            raise ArithmeticError('CP source fails requested exact row sums')
        result.update(axes=color, signs=signs.tolist(), rows=rows.tolist(),
                      expanded_verification=proof,source_claimed=True)
    if output is not None:
        output=Path(output);output.parent.mkdir(parents=True,exist_ok=True)
        output.write_text(json.dumps(result,indent=2)+'\n')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--p',type=int,required=True)
    parser.add_argument('--m',type=int,required=True)
    parser.add_argument('--action',choices=('identity','C3'),required=True)
    parser.add_argument('--seconds',type=float,default=10)
    parser.add_argument('--seed',type=int,default=1)
    parser.add_argument('--target',help='Four comma-separated integer row sums')
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--no-gauge',action='store_true')
    args=parser.parse_args()
    target=None if args.target is None else list(map(int,args.target.split(',')))
    result=solve(args.p,args.m,args.action,seconds=args.seconds,seed=args.seed,
                 target=target,gauge=not args.no_gauge,output=args.output)
    print(json.dumps({k:v for k,v in result.items() if k not in
                      ('axes','signs','rows','expanded_verification')}))


if __name__=='__main__':main()

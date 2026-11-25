"""Exact CP-SAT repair of a neighborhood in a twisted periodic T-design.

All states outside the selected multiplication orbits are fixed. Free
states retain their full signed-coordinate meaning. The model imposes the
actual Pi-twisted integer correlations and exact head/tail sign counts.
A state 2*a+epsilon denotes (-1)^epsilon times standard vector e_a.
Coordinate 0 orbits are called heads; the other three coordinates are tails.
The coordinate permutation fixes 0 and sends 1 to 2 to 3 to 1.
Every returned candidate is independently expanded at all p residues.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
from math import isqrt
from pathlib import Path
import random
import time

import numpy as np

from experiments.scalar_lift import hadamard_frame_search as frame


def dot_state(left, right, delta):
    a, sa = divmod(left, 2)
    b, sb = divmod(right, 2)
    rotated = 0 if b == 0 else 1+(b-1+delta) % 3
    return (1 if sa == sb else -1) if a == rotated else 0


def exact_state_report(states, data, head_net, tail_net):
    states = list(map(int, states))
    if len(states) != data['n'] or states[0] != 0 or any(not 0 <= s < 8 for s in states):
        raise ValueError('complete signed-coordinate states with b(0)=e0 are required')
    b = frame.expand(states, data)
    sums = b.sum(axis=1).tolist()
    expected = [1+data['m']*head_net]+[data['m']//3*tail_net]*3
    if sum(x*x for x in expected) != data['p']:
        raise ValueError('requested coordinate sums fail the exact norm condition')
    if sums != expected:
        raise ValueError('candidate does not have the declared coordinate sums')
    if not np.all(np.sum(b*b, axis=0) == 1):
        raise ValueError('candidate is not a signed-coordinate function')
    correlations = [int(np.sum(b*np.roll(b, -s, axis=1))) for s in range(data['p'])]
    compressed = [sum(count*dot_state(states[i],states[j],delta)
                      for (i,j,delta),count in row.items()) for row in data['rows']]
    if compressed != [correlations[s] for s in data['shifts']]:
        raise ValueError('twisted model and expanded correlations disagree')
    return dict(states=states, row_sums=sums,
                head_orbits=sum(s//2 == 0 for s in states[1:]),
                correlations=compressed,
                energy=sum(c*c for c in compressed),
                l1=sum(abs(c) for c in compressed),
                max_expanded_nonzero_shift_residual=max(map(abs,correlations[1:]),default=0),
                all_expanded_shift_correlations_checked=True,
                exact_T_design=all(c == 0 for c in correlations[1:]))


def choose_free(states, count, seed=0, all_heads=True):
    if not 1 <= count < len(states):
        raise ValueError('neighborhood size must be between 1 and n-1')
    rng = random.Random(seed)
    heads = [i for i in range(1,len(states)) if states[i]//2 == 0]
    if all_heads:
        selected = rng.sample(heads,min(count,len(heads)))
    else:
        selected = []
    remaining = [i for i in range(1,len(states)) if i not in selected]
    selected += rng.sample(remaining,count-len(selected))
    return sorted(selected)


def build_model(states, data, free_orbits, head_net=-5, tail_net=3,
                objective='exact', improve=True, encoding='table'):
    from ortools.sat.python import cp_model
    initial = exact_state_report(states,data,head_net,tail_net)
    states = initial['states']
    free = sorted(set(map(int,free_orbits)))
    if not free or free[0] < 1 or free[-1] >= len(states):
        raise ValueError('free orbits must be nonzero valid indices')
    heads = initial['head_orbits']
    tails = len(states)-1-heads
    if (heads-head_net)%2 or (tails-tail_net)%2:
        raise ValueError('exact counts are not integral')
    wanted_negative_head = (heads-head_net)//2
    wanted_negative_tail = (tails-tail_net)//2
    allowed = list(range(8))
    if wanted_negative_head == heads:
        allowed.remove(0)
    elif wanted_negative_head == 0:
        allowed.remove(1)
    model = cp_model.CpModel()
    variables, indicators = {}, {}
    for i in free:
        variables[i] = model.new_int_var_from_domain(cp_model.Domain.from_values(allowed),f's{i}')
        row = {}
        for state in allowed:
            literal = model.new_bool_var(f's{i}_is{state}')
            model.add(variables[i] == state).only_enforce_if(literal)
            model.add(variables[i] != state).only_enforce_if(literal.Not())
            row[state] = literal
            model.add_hint(literal,int(states[i] == state))
        model.add_exactly_one(list(row.values()))
        indicators[i] = row
        model.add_hint(variables[i],states[i])
    fixed = set(range(1,len(states)))-set(free)
    requirements = [
        (lambda s:s//2 == 0, heads),
        (lambda s:s == 1, wanted_negative_head),
        (lambda s:s//2 != 0 and s%2 == 1,wanted_negative_tail),
    ]
    for predicate,target in requirements:
        fixed_count = sum(predicate(states[i]) for i in fixed)
        terms = [literal for i in free for state,literal in indicators[i].items() if predicate(state)]
        model.add(sum(terms) == target-fixed_count)

    constant = [0]*len(data['rows'])
    boundary = [defaultdict(lambda:[0]*8) for _ in data['rows']]
    interacting = defaultdict(set)
    for ell,row in enumerate(data['rows']):
        for (i,j,delta),count in row.items():
            if i in variables and j in variables and i != j:
                interacting[i,j].add(delta)
            elif i in variables or j in variables:
                free_i = i if i in variables else j
                values = boundary[ell][free_i]
                for state in allowed:
                    a = state if i == free_i else states[i]
                    b = state if j == free_i else states[j]
                    values[state] += count*dot_state(a,b,delta)
            else:
                constant[ell] += count*dot_state(states[i],states[j],delta)
    dots = {}
    if encoding not in ('table', 'bits'):
        raise ValueError('pair encoding must be table or bits')
    bit_states, rotated_bits = {}, {}
    def xor_bit(a, b, name, hint):
        out = model.new_bool_var(name)
        model.add_bool_xor([a,b,out.Not()])
        model.add_hint(out,int(hint))
        return out
    def and_bit(a, b, name, hint):
        out = model.new_bool_var(name)
        model.add_implication(out,a)
        model.add_implication(out,b)
        model.add_bool_or([out,a.Not(),b.Not()])
        model.add_hint(out,int(hint))
        return out
    if encoding == 'bits':
        for i in free:
            bit_states[i] = []
            for bit in range(3):
                out = model.new_bool_var(f'bit{bit}_{i}')
                model.add(out == sum(literal for state,literal in indicators[i].items() if (state >> bit) & 1))
                model.add_hint(out,(states[i] >> bit) & 1)
                bit_states[i].append(out)
            a,b = bit_states[i][1:]
            c = xor_bit(a,b,f'axisxor{i}',((states[i] >> 1) ^ (states[i] >> 2)) & 1)
            rotated_bits[i] = [(a,b),(b,c),(c,a)]
    for (i,j),deltas in interacting.items():
        deltas = sorted(deltas)
        if encoding == 'table':
            products = [model.new_int_var(-1,1,f'dot{i}_{j}_{delta}') for delta in deltas]
            tuples = [[a,b,*[dot_state(a,b,delta) for delta in deltas]] for a in allowed for b in allowed]
            model.add_allowed_assignments([variables[i],variables[j],*products],tuples)
            for delta,variable in zip(deltas,products,strict=True):
                dots[i,j,delta] = variable
                model.add_hint(variable,dot_state(states[i],states[j],delta))
        else:
            different = xor_bit(bit_states[i][0],bit_states[j][0],f'signdiff{i}_{j}',(states[i] ^ states[j]) & 1)
            a0,a1 = bit_states[i][1:]
            for delta in deltas:
                b0,b1 = rotated_bits[j][delta]
                rotated = 0 if states[j]//2 == 0 else 1+(states[j]//2-1+delta) % 3
                difference = (states[i]//2) ^ rotated
                d0 = xor_bit(a0,b0,f'ax0diff{i}_{j}_{delta}',difference & 1)
                d1 = xor_bit(a1,b1,f'ax1diff{i}_{j}_{delta}',difference >> 1)
                match = and_bit(d0.Not(),d1.Not(),f'match{i}_{j}_{delta}',difference == 0)
                opposite = and_bit(match,different,f'opp{i}_{j}_{delta}',difference == 0 and bool((states[i] ^ states[j]) & 1))
                dots[i,j,delta] = match-2*opposite
    residuals = []
    initial_cost = initial['energy'] if objective == 'energy' else initial['l1']
    bound = max(0,initial_cost-int(improve))
    for ell,row in enumerate(data['rows']):
        expr = constant[ell]+sum(coefficient*indicators[i][state]
               for i,values in boundary[ell].items() for state in allowed
               if (coefficient := values[state]))
        expr += sum(count*dots[i,j,delta] for (i,j,delta),count in row.items()
                    if (i,j,delta) in dots)
        if objective == 'exact':
            model.add(expr == 0)
        elif objective == 'energy':
            radius = min(data['p'],isqrt(bound))
            residual = model.new_int_var(-radius,radius,f'c{ell}')
            square = model.new_int_var(0,radius*radius,f'energy{ell}')
            model.add(residual == expr)
            model.add_allowed_assignments([residual,square],[[v,v*v] for v in range(-radius,radius+1)])
            residuals.append(square)
            if abs(initial['correlations'][ell]) <= radius:
                model.add_hint(residual,initial['correlations'][ell])
                model.add_hint(square,initial['correlations'][ell]**2)
        elif objective == 'l1':
            absolute = model.new_int_var(0,min(data['p'],bound),f'abs{ell}')
            model.add_abs_equality(absolute,expr)
            residuals.append(absolute)
            if abs(initial['correlations'][ell]) <= bound:
                model.add_hint(absolute,abs(initial['correlations'][ell]))
        else:
            raise ValueError('objective must be exact, energy or l1')
    if residuals:
        model.add(sum(residuals) <= bound)
        model.minimize(sum(residuals))
    error = model.validate()
    if error:raise ValueError(error)
    info = dict(free_orbits=free, neighborhood_size=len(free),
                total_model_variables=len(model.proto.variables),
                total_model_constraints=len(model.proto.constraints),
                coupled_state_pairs=len(interacting), pair_encoding=encoding,
                coupled_state_pair_tables=len(interacting) if encoding == 'table' else 0,
                exact_head_orbits=heads, exact_negative_head_orbits=wanted_negative_head,
                exact_negative_tail_orbits=wanted_negative_tail,
                objective=objective, strict_improvement=bool(improve and residuals),
                initial_energy=initial['energy'], initial_l1=initial['l1'])
    return model,variables,info


def repair(states,data,free_orbits,head_net=-5,tail_net=3,seconds=5,
           objective='exact',improve=True,seed=0,encoding='table'):
    from ortools.sat.python import cp_model
    if seconds <= 0:raise ValueError('a positive time bound is required')
    start = time.monotonic()
    model,variables,info = build_model(states,data,free_orbits,head_net,tail_net,objective,improve,encoding)
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_search_workers = 1
    solver.parameters.random_seed = seed
    solver.parameters.max_presolve_iterations = 1
    solver.parameters.cp_model_presolve = False
    solver.parameters.cp_model_probing_level = 0
    solver.parameters.linearization_level = 0
    solver.parameters.repair_hint = False
    solver.parameters.hint_conflict_limit = 100
    status = solver.solve(model)
    result = dict(info, solver_status=solver.status_name(status),
                  build_and_solve_seconds=time.monotonic()-start,
                  solver_seconds=solver.wall_time,
                  branches=solver.num_branches,conflicts=solver.num_conflicts,
                  completed_neighborhood_decision=status in (cp_model.OPTIMAL,cp_model.INFEASIBLE))
    if status in (cp_model.OPTIMAL,cp_model.FEASIBLE):
        candidate = list(states)
        for i,variable in variables.items():candidate[i] = solver.value(variable)
        verification = exact_state_report(candidate,data,head_net,tail_net)
        if objective == 'exact' and not verification['exact_T_design']:
            raise ValueError('solver claimed an invalid exact completion')
        if improve:
            if objective == 'energy' and verification['energy'] >= info['initial_energy']:
                raise ValueError('energy did not strictly improve')
            if objective == 'l1' and verification['l1'] >= info['initial_l1']:
                raise ValueError('L1 did not strictly improve')
        result.update(verification)
        result['unchanged_fixed_orbits_verified'] = all(candidate[i] == states[i] for i in range(len(states)) if i not in variables)
    return result


def parity_projection_model(states, data, changes):
    """Keep the head set fixed and change at most `changes` tail colors.

    This encodes every unsigned-correlation parity exactly. Signs are not
    search variables at this stage; they are retained only for reporting
    the projected candidate's actual integer correlations.
    """
    from experiments.scalar_lift.boolean_quadratic_sat import XorCNF
    from experiments.scalar_lift.large_parity_map import CorrelationForms
    if not isinstance(changes,int) or changes < 0:
        raise ValueError('a nonnegative integer color-change bound is required')
    axes = [state//2 for state in states]
    if len(axes) != data['n'] or axes[0] != 0:
        raise ValueError('a complete signed-coordinate seed is required')
    forms = CorrelationForms(data)
    cnf = XorCNF()
    color_bits, changed = [], []
    for i in range(1,len(axes)):
        if axes[i] == 0:
            color_bits.extend([-cnf.true,-cnf.true])
            continue
        a,b = cnf.var(),cnf.var()
        color_bits.extend([a,b])
        cnf.clause([a,b])
        da = -a if axes[i] & 1 else a
        db = -b if axes[i] & 2 else b
        changed.append(-cnf.and_literal(-da,-db))
    slack = [cnf.var() for _ in range(changes.bit_length())]
    cnf.weighted_exact([(literal,1) for literal in changed] +
                       [(literal,1 << bit) for bit,literal in enumerate(slack)],changes)
    terms = [[] for _ in data['rows']]
    for i,coefficient in enumerate(forms.linear):
        for ell in _set_bits(coefficient):terms[ell].append(color_bits[i])
    for i,row in enumerate(forms.cross):
        for j in range(i):
            coefficient = row[j]
            if coefficient:
                product = cnf.and_literal(color_bits[i],color_bits[j])
                for ell in _set_bits(coefficient):terms[ell].append(product)
    for ell,equation in enumerate(terms):
        cnf.xor(equation,(forms.constant >> ell) & 1)
    return cnf,color_bits,dict(color_change_bound=changes,head_set_fixed=True,
                               head_orbits=sum(a == 0 for a in axes[1:]),
                               tail_orbits=len(changed),
                               exact_unsigned_parity_equations=len(data['rows']))


def _set_bits(value):
    while value:
        low = value & -value
        yield low.bit_length()-1
        value ^= low


def project_parity(states, data, changes=8, head_net=-5, tail_net=3,
                   seconds=10, seed=0, workdir=None):
    before = exact_state_report(states,data,head_net,tail_net)
    cnf,color_bits,info = parity_projection_model(states,data,changes)
    start = time.monotonic()
    result = cnf.solve(seconds,seed,workdir)
    values = result.pop('values')
    result.update(info,projection_seconds=time.monotonic()-start,
                  initial_energy=before['energy'],initial_l1=before['l1'],
                  initial_odd_correlation_classes=sum(c % 2 for c in before['correlations']))
    if result['status'] != 'SAT':return result
    axes = [0]
    def literal_value(literal):
        value = values[abs(literal)]
        return value if literal > 0 else 1-value
    for i in range(0,len(color_bits),2):
        axes.append(literal_value(color_bits[i]) | (literal_value(color_bits[i+1]) << 1))
    original_axes = [s//2 for s in states]
    changed = [i for i in range(1,len(axes)) if axes[i] != original_axes[i]]
    if len(changed) > changes or any((a == 0) != (b == 0) for a,b in zip(axes,original_axes,strict=True)):
        raise ValueError('parity projection changed the fixed head set or exceeded its bound')
    for row in data['rows']:
        if sum(count for (i,j,delta),count in row.items()
               if axes[i] == (0 if axes[j] == 0 else 1+(axes[j]-1+delta) % 3)) % 2:
            raise ValueError('parity projection failed an actual unsigned correlation')
    candidate = [2*a+(s & 1) for a,s in zip(axes,states,strict=True)]
    preferred_sign = int(head_net < 0)
    chosen = next((i for i in range(1,len(candidate)) if candidate[i]//2 == 0 and candidate[i] % 2 == preferred_sign),None)
    if chosen is None:
        raise ValueError('there is no head of the requested sign to normalize')
    scale = data['orbits'][chosen][0][0]
    normalized = []
    normalized_seed_axes = []
    for orbit in data['orbits']:
        point = scale*orbit[0][0] % data['p']
        source = int(data['index'][point]);phase = int(data['phase'][point])
        axis,negative = divmod(candidate[source],2)
        normalized.append(2*(0 if axis == 0 else 1+(axis-1+phase) % 3)+negative)
        old_axis = states[source]//2
        normalized_seed_axes.append(0 if old_axis == 0 else 1+(old_axis-1+phase) % 3)
    expanded = frame.expand(candidate,data)
    direct = frame.expand(normalized,data)
    if not np.array_equal(direct,expanded[:,(scale*np.arange(data['p'])) % data['p']]):
        raise ValueError('normalization is not the declared position permutation')
    after = exact_state_report(normalized,data,head_net,tail_net)
    if any(c % 2 for c in after['correlations']):
        raise ValueError('normalized projection lost parity feasibility')
    result.update(axes=[s//2 for s in normalized],states=normalized,
                  normalized_seed_axes=normalized_seed_axes,
                  head_set_fixed_before_declared_gauge=True,
                  changed_tail_orbits=changed,changed_tail_orbit_count=len(changed),
                  gauge_scale=scale,first_head_negative=bool(preferred_sign),
                  retained_original_signs=True,original_signs_energy=after['energy'],
                  original_signs_l1=after['l1'],original_signs_correlations=after['correlations'],
                  all_expanded_original_sign_correlations_checked=True,
                  row_sums=after['row_sums'],all_unsigned_parities_verified=True,
                  exact_T_design_with_original_signs=after['exact_T_design'])
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--free',type=int,default=20)
    parser.add_argument('--seconds',type=float,default=5)
    parser.add_argument('--trials',type=int,default=1)
    parser.add_argument('--seed',type=int,default=0)
    parser.add_argument('--objective',choices=['exact','energy','l1'],default='exact')
    parser.add_argument('--allow-equal',action='store_true')
    parser.add_argument('--encoding',choices=['table','bits'],default='table')
    parser.add_argument('--parity-projection',action='store_true')
    parser.add_argument('--changes',type=int,default=8)
    args=parser.parse_args()
    record=json.loads(args.input.read_text())
    data=frame.chart(record.get('p',223),record.get('m',3),record.get('mu',183))
    states=record['states']
    if args.parity_projection:
        result=project_parity(states,data,args.changes,record.get('head_net',-5),
                              record.get('tail_net',3),args.seconds,args.seed)
        result.update(p=data['p'],m=data['m'],mu=data['mu'],
                      head_net=record.get('head_net',-5),tail_net=record.get('tail_net',3))
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps({k:v for k,v in result.items() if k not in ('states','axes','original_signs_correlations','changed_tail_orbits')}),flush=True)
        return
    outputs=[]
    for trial in range(args.trials):
        free=choose_free(states,args.free,args.seed+trial)
        result=repair(states,data,free,record.get('head_net',-5),record.get('tail_net',3),
                      args.seconds,args.objective,not args.allow_equal,args.seed+trial,args.encoding)
        outputs.append(result)
        if 'states' in result:states=result['states']
        print(json.dumps({k:v for k,v in result.items() if k not in ('states','free_orbits','correlations')}),flush=True)
        if result.get('exact_T_design'):break
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(dict(p=data['p'],m=data['m'],mu=data['mu'],
                                           head_net=record.get('head_net',-5),tail_net=record.get('tail_net',3),
                                           states=states,trials=outputs),indent=2)+'\n')


if __name__=='__main__':main()

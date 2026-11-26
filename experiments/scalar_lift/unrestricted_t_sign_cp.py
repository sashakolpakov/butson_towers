"""Bounded exact CP-SAT sign lifts for unrestricted periodic T-matrix axes.

The coordinate label at each cyclic site is fixed. Signs satisfy all four
exact row sums and every integer periodic correlation. Affine mod-four links
and the quadratic mod-eight equations are optional redundant constraints;
they use Python integer masks and impose no 64-variable limit.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import operator
from pathlib import Path
import time

from experiments.scalar_lift import unrestricted_t_sign_model as exact


def bits(value):
    while value:
        low=value&-value
        yield low.bit_length()-1
        value^=low


def xor_equal(model,literals,rhs):
    if not literals:
        if rhs:model.add(False)
    else:
        model.add_bool_xor(literals if rhs else [*literals,True])


def build_model(system,*,add_mod8=True,affine_links=True,fixed_free=None):
    """Encode exact counts and correlations, plus the supplied affine/MQ cuts."""
    from ortools.sat.python import cp_model
    if add_mod8 and not affine_links:
        raise ValueError('mod-eight equations require affine links')
    model=cp_model.CpModel()
    n=len(system['axes']);dimension=len(system['basis'])
    signs=[model.new_bool_var(f'x{i}') for i in range(n)]
    free=[model.new_bool_var(f'y{j}') for j in range(dimension)] if affine_links else []
    if affine_links:
        for i,form in enumerate(system['forms']):
            xor_equal(model,[signs[i]]+[free[j] for j in bits(form)],(system['particular']>>i)&1)
    for i,value in system['forced_signs'].items():model.add(signs[i]==value)
    fixed_free=dict(fixed_free or {})
    for index,value in fixed_free.items():
        if isinstance(index,bool) or isinstance(value,bool):raise ValueError('fixed free signs require integer indices and bits')
        index,value=operator.index(index),operator.index(value)
        if not affine_links or not 0<=index<dimension or value not in (0,1):
            raise ValueError('fixed free signs require in-range indices and binary values')
        model.add(free[index]==value)
    for group,target in zip(system['groups'],system['negative_counts'],strict=True):
        model.add(sum(signs[i] for i in group)==target)
    cuts={}
    for diagonal,weight in zip(system['diagonal'],system['weights'],strict=True):
        absolute=int(diagonal);terms=[]
        for i in range(n):
            for j in range(i+1,n):
                w=int(weight[i,j])
                if not w:continue
                absolute+=w
                if (i,j) not in cuts:
                    cut=model.new_bool_var(f'cut{i}_{j}')
                    xor_equal(model,[signs[i],signs[j],cut],0)
                    cuts[i,j]=cut
                terms.append(w*cuts[i,j])
        if absolute%2:raise exact.InfeasibleAxes('unsigned correlation has odd parity')
        model.add(sum(terms)==absolute//2)
    equations=exact.quadratic_mod8(system) if add_mod8 else []
    products={}
    for equation in equations:
        literals=[free[j] for j in bits(equation['linear'])]
        for u,v in equation['quadratic']:
            pair=tuple(sorted((u,v)))
            if pair not in products:
                z=model.new_bool_var(f'product{pair[0]}_{pair[1]}')
                model.add_multiplication_equality(z,[free[pair[0]],free[pair[1]]])
                products[pair]=z
            literals.append(products[pair])
        xor_equal(model,literals,equation['constant'])
    error=model.validate()
    if error:raise ValueError(error)
    return model,signs,free,{
        'model_variables':len(model.proto.variables),
        'model_constraints':len(model.proto.constraints),
        'cyclic_sites':n,'affine_dimension':dimension,
        'affine_mod4_links_in_model':bool(affine_links),
        'quadratic_mod8_equations':len(equations),
        'cut_variables':len(cuts),'quadratic_product_variables':len(products),
        'exact_row_sum_constraints':4,
        'integer_correlation_constraints':len(system['weights']),
        'exact_negative_counts':list(system['negative_counts']),
        'fixed_free_variables':{int(k):int(v) for k,v in fixed_free.items()},
    }


def solve(axes,data_or_p,target_row,seconds=30,seed=0,max_solutions=1,workdir=None,*,
          all_solutions=False,add_mod8=True,affine_links=True,fixed_signs=None,
          fixed_free=None,log_search=False):
    """Run one CP-SAT worker and independently verify every positive word."""
    from ortools.sat.python import cp_model
    if not math.isfinite(seconds) or seconds<=0 or max_solutions<0:raise ValueError('finite positive duration and nonnegative solution limit required')
    started=time.monotonic()
    common={'method':'CP-SAT exact unrestricted periodic T sign lift',
            'requested_wall_seconds':seconds,'requested_seed':int(seed),
            'solver_workers':1,'solutions':[],'stats':{},
            'scope':'Fixed cyclic coordinate labels, four target sums, and recorded fixed signs/free variables only.'}
    try:
        system=exact.prepare(axes,data_or_p,target_row,fixed_signs=fixed_signs)
    except exact.InfeasibleAxes as exc:
        return dict(common,solver_status='INFEASIBLE_ARITHMETIC',reason=str(exc),
                    complete_for_fixed_axes=True,complete_for_fixed_constraints=True,elapsed_seconds=time.monotonic()-started)
    common.update(p=system['data']['p'],axes=list(system['axes']),target_row=list(system['target_row']),
                  fixed_signs={int(k):int(v) for k,v in (fixed_signs or {}).items()})
    model,signs,free,metadata=build_model(system,add_mod8=add_mod8,
                                         affine_links=affine_links,fixed_free=fixed_free)
    build_seconds=time.monotonic()-started
    remaining=seconds-build_seconds
    if remaining<=0:
        return dict(common,**metadata,solver_status='UNKNOWN',
                    termination='budget exhausted during exact model construction',
                    complete_for_fixed_axes=False,complete_for_fixed_constraints=False,build_seconds=build_seconds,
                    solver_seconds=0,elapsed_seconds=time.monotonic()-started)
    solver=cp_model.CpSolver()
    solver.parameters.max_time_in_seconds=remaining
    solver.parameters.num_search_workers=1
    solver.parameters.random_seed=int(seed)%2147483647
    solver.parameters.enumerate_all_solutions=bool(all_solutions)
    solver.parameters.log_search_progress=bool(log_search)
    solutions=[];limited=False
    class Collector(cp_model.CpSolverSolutionCallback):
        def on_solution_callback(self):
            nonlocal limited
            mask=sum(int(self.value(x))<<i for i,x in enumerate(signs))
            if any(sum((mask>>i)&1 for i in group)!=target
                   for group,target in zip(system['groups'],system['negative_counts'],strict=True)):
                raise ArithmeticError('CP-SAT candidate violated an exact negative count')
            if any(exact.correlations(system,mask)):
                raise ArithmeticError('CP-SAT candidate violated an integer correlation')
            states=[2*a+((mask>>i)&1) for i,a in enumerate(system['axes'])]
            proof=exact.verify_states(states,system['data'],system['target_row'])
            solutions.append({'states':states,'negative_mask':mask,
                              'expanded_verification':proof})
            if not all_solutions or (max_solutions and len(solutions)>=max_solutions):
                limited=True;self.stop_search()
    status=solver.solve(model,Collector())
    complete=not limited and status in (cp_model.OPTIMAL,cp_model.INFEASIBLE)
    common.pop('stats',None)
    return dict(common,**metadata,solver_status=solver.status_name(status),solutions=solutions,
                complete_for_fixed_axes=complete,complete_for_fixed_constraints=complete,solution_limit_reached=limited,
                stats={'branches':solver.num_branches,'conflicts':solver.num_conflicts,
                       'solver_seconds':solver.wall_time,'affine_dimension':metadata['affine_dimension']},
                build_seconds=build_seconds,solver_seconds=solver.wall_time,
                solver_branches=solver.num_branches,solver_conflicts=solver.num_conflicts,
                elapsed_seconds=time.monotonic()-started)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--axes',type=Path,required=True)
    parser.add_argument('--p',type=int,required=True)
    parser.add_argument('--target',help='four comma-separated integer row sums')
    parser.add_argument('--seconds',type=float,default=30)
    parser.add_argument('--seed',type=int,default=0)
    parser.add_argument('--all-solutions',action='store_true')
    parser.add_argument('--max-solutions',type=int,default=1)
    parser.add_argument('--without-mod8',action='store_true')
    parser.add_argument('--report',type=Path,required=True)
    args=parser.parse_args();record=json.loads(args.axes.read_text())
    axes=record['axes'] if isinstance(record,dict) else record
    target=tuple(map(int,args.target.split(','))) if args.target else record.get('target_row') if isinstance(record,dict) else None
    if target is None:raise ValueError('provide --target or target_row in the axes record')
    report=solve(axes,args.p,target,seconds=args.seconds,seed=args.seed,
                 all_solutions=args.all_solutions,max_solutions=args.max_solutions,
                 add_mod8=not args.without_mod8)
    report['source_sha256']={Path(p).name:hashlib.sha256(Path(p).read_bytes()).hexdigest()
                             for p in (__file__,exact.__file__)}
    args.report.parent.mkdir(parents=True,exist_ok=True)
    temporary=args.report.with_suffix(args.report.suffix+'.tmp')
    temporary.write_text(json.dumps(report,indent=2)+'\n');temporary.replace(args.report)
    print(json.dumps({key:report.get(key) for key in
          ('p','solver_status','complete_for_fixed_constraints','elapsed_seconds')}
          |{'exact_solutions':len(report['solutions'])}))


if __name__=='__main__':main()

"""Compare the exact twisted-neighborhood model with independent finite checks."""
from itertools import product

import pytest

from ortools.sat.python import cp_model

from experiments.scalar_lift import hadamard_frame_search as frame
from experiments.scalar_lift.t_sequence_neighborhood import (
    build_model, choose_free, exact_state_report, repair,
)

T19=[0,0,1,0,3,2,6]


@pytest.mark.parametrize('encoding',['table','bits'])
def test_fixed_neighborhood_complete_solution_set_matches_expanded_arithmetic(encoding):
    data=frame.chart(19,3,7)
    free=[1,2,4,5]
    expected=set()
    for values in product(range(8),repeat=len(free)):
        states=T19.copy()
        for i,state in zip(free,values):states[i]=state
        try:report=exact_state_report(states,data,1,1)
        except ValueError:continue
        if report['exact_T_design'] and report['head_orbits']==3:expected.add(tuple(values))
    model,variables,info=build_model(T19,data,free,1,1,encoding=encoding)
    class Collector(cp_model.CpSolverSolutionCallback):
        def __init__(self):super().__init__();self.values=set()
        def on_solution_callback(self):self.values.add(tuple(self.value(variables[i]) for i in free))
    collector=Collector();solver=cp_model.CpSolver()
    solver.parameters.enumerate_all_solutions=True
    solver.parameters.num_search_workers=1
    solver.parameters.max_time_in_seconds=5
    status=solver.solve(model,collector)
    assert status==cp_model.OPTIMAL
    assert collector.values==expected
    assert expected
    assert info['exact_head_orbits']==3


def test_neighborhood_repair_returns_actual_all_shift_witness():
    data=frame.chart(19,3,7)
    damaged=T19.copy();damaged[3],damaged[6]=damaged[6],damaged[3]
    before=exact_state_report(damaged,data,1,1)
    assert before['energy']>0
    result=repair(damaged,data,[3,6],1,1,seconds=5)
    assert result['solver_status'] in ('FEASIBLE','OPTIMAL')
    assert result['exact_T_design']
    assert result['unchanged_fixed_orbits_verified']
    assert result['all_expanded_shift_correlations_checked']
    assert result['row_sums']==[4,1,1,1]


def test_neighborhood_energy_objective_is_a_strict_verified_improvement():
    data=frame.chart(19,3,7)
    damaged=T19.copy();damaged[3],damaged[6]=damaged[6],damaged[3]
    result=repair(damaged,data,[3,6],1,1,seconds=5,objective='energy')
    assert result['solver_status'] in ('FEASIBLE','OPTIMAL')
    assert result['energy']<result['initial_energy']
    assert result['unchanged_fixed_orbits_verified']
    assert result['energy']==0
    free=choose_free(T19,4,seed=0)
    assert set([1,2,3])<=set(free)
    assert 0 not in free


def test_parity_projection_zero_radius_and_nearby_feasible_axes():
    from experiments.scalar_lift.t_sequence_neighborhood import project_parity
    data=frame.chart(19,3,7)
    fixed=project_parity(T19,data,changes=0,head_net=1,tail_net=1,seconds=5)
    assert fixed['status']=='SAT'
    assert fixed['changed_tail_orbit_count']==0
    assert fixed['all_unsigned_parities_verified']
    assert fixed['exact_T_design_with_original_signs']
    damaged=T19.copy();damaged[4]=5
    # A tail color changes, retaining its sign and every exact coordinate sum.
    assert exact_state_report(damaged,data,1,1)['row_sums']==[4,1,1,1]
    excluded=project_parity(damaged,data,changes=0,head_net=1,tail_net=1,seconds=5)
    assert excluded['status']=='UNSAT'
    repaired=project_parity(damaged,data,changes=1,head_net=1,tail_net=1,seconds=5)
    assert repaired['status']=='SAT'
    assert repaired['changed_tail_orbit_count']==1
    assert all((a==0)==(b==0) for a,b in zip(repaired['axes'],repaired['normalized_seed_axes']))
    assert repaired['all_unsigned_parities_verified']
    assert repaired['row_sums']==[4,1,1,1]
    assert all(c%2==0 for c in repaired['original_signs_correlations'])
    negative=project_parity([0,1,2],frame.chart(7,3,2),changes=0,head_net=-1,tail_net=1,seconds=5)
    assert negative['status']=='SAT' and negative['states'][1]==1
    assert negative['first_head_negative']

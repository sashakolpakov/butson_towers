"""Compare exact CP-SAT sign solutions with independent finite cyclic sums."""
from itertools import combinations,product

import numpy as np
import pytest

pytest.importorskip('ortools.sat.python.cp_model')

from experiments.scalar_lift import hadamard_frame_search as frame
from experiments.scalar_lift import unrestricted_t_sign_cp as backend


def exact_masks(axes,target):
    """Enumerate only exact-count sign words and compute literal cyclic sums."""
    n=len(axes);groups=[[i for i,a in enumerate(axes) if a==j] for j in range(4)]
    choices=[]
    for group,wanted in zip(groups,target,strict=True):
        negative=(len(group)-wanted)//2
        choices.append(list(combinations([i for i in group if i],negative)))
    answer=set()
    for blocks in product(*choices):
        negative={i for block in blocks for i in block}
        signs=[-1 if i in negative else 1 for i in range(n)]
        if all(sum(signs[x]*signs[(x+s)%n] for x in range(n)
                   if axes[x]==axes[(x+s)%n])==0 for s in range(1,n)):
            answer.add(sum(1<<i for i in negative))
    return answer


@pytest.mark.parametrize('affine,mod8',[(True,True),(True,False),(False,False)])
def test_all_T7_solutions_equal_independent_complete_enumeration(affine,mod8):
    planted=[0,6,0,4,4,5,2];axes=[s//2 for s in planted];target=(2,1,1,1)
    expected=exact_masks(axes,target)
    result=backend.solve(axes,7,target,seconds=5,all_solutions=True,max_solutions=0,
                         affine_links=affine,add_mod8=mod8)
    assert result['solver_status']=='OPTIMAL'
    assert result['complete_for_fixed_axes']
    assert {s['negative_mask'] for s in result['solutions']}==expected
    assert len(expected)==2
    for solution in result['solutions']:
        proof=solution['expanded_verification']
        assert proof['every_expanded_periodic_correlation_verified']
        assert proof['all_integer_correlations']==[7]+[0]*6
        assert proof['sequence_row_sums']==list(target)


def test_unrestricted_T19_solution_set_matches_every_exact_count_sign_word():
    data=frame.chart(19,3,7)
    b=frame.expand([0,0,1,0,3,2,6],data)
    axes=np.argmax(abs(b),axis=0).tolist();target=(4,1,1,1)
    expected=exact_masks(axes,target)
    result=backend.solve(axes,19,target,seconds=10,all_solutions=True,max_solutions=0)
    assert result['complete_for_fixed_axes']
    assert {s['negative_mask'] for s in result['solutions']}==expected
    assert expected
    assert result['exact_row_sum_constraints']==4
    assert result['integer_correlation_constraints']==9
    assert result['solver_workers']==1


def test_first_solution_limit_and_exhausted_wall_budget_do_not_claim_completion():
    axes=[0,3,0,2,2,2,1];target=(2,1,1,1)
    result=backend.solve(axes,7,target,seconds=5,seed=20261126179)
    assert result['solutions'] and not result['complete_for_fixed_axes']
    assert result['solution_limit_reached']
    unfinished=backend.solve(axes,7,target,seconds=1e-12)
    assert unfinished['solver_status']=='UNKNOWN'
    assert unfinished['solutions']==[]
    assert not unfinished['complete_for_fixed_axes']
    with pytest.raises(ValueError,match='positive duration'):
        backend.solve(axes,7,target,seconds=0)


def test_count_obstruction_is_distinguished_from_solver_timeout():
    # A valid norm-seven target cannot be carried by seven coordinate-zero sites.
    result=backend.solve([0]*7,7,(2,1,1,1),seconds=1)
    assert result['solver_status']=='INFEASIBLE_ARITHMETIC'
    assert result['complete_for_fixed_axes']
    assert result['solutions']==[]


@pytest.mark.parametrize('p,target,dimension,word',[
    (167,(9,6,5,5),80,
     '00323223210132001320010230020221022331203121332023111021303313310021300020231110120000132210201211230302200102121203102310023002222132032000202322020232032221213113032'),
    (179,(9,8,5,3),86,
     '01120231132012333210321012132331203011330103310203033001123003313131321201020133323303220232013030313031030103300201300130232233330303100030121013123110332301031002111001100113033'),
])
def test_actual_frontier_models_retain_more_than_64_affine_bits(p,target,dimension,word):
    # These are parity-compatible coordinate words, not completed T-designs.
    axes=list(map(int,word))
    assert len(axes)==p
    assert all(sum(axes[x]==axes[(x+s)%p] for x in range(p))%2==0 for s in range(1,p))
    system=backend.exact.prepare(axes,p,target)
    assert len(system['basis'])==dimension
    assert any(form>>64 for form in system['forms'])
    model,signs,free,metadata=backend.build_model(system)
    assert model.validate()==''
    assert len(signs)==p and len(free)==dimension
    assert metadata['quadratic_mod8_equations']==(p-1)//2
    assert metadata['integer_correlation_constraints']==(p-1)//2
    assert metadata['exact_row_sum_constraints']==4

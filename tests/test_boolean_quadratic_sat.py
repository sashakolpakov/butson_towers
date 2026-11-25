"""Exhaustive arithmetic encodings and independent signed-matrix controls."""
import itertools
import json
from pathlib import Path
import shutil

import pytest

from experiments.scalar_lift.boolean_quadratic_sat import XorCNF, solve
from experiments.scalar_lift import hadamard_frame_search as frame
from experiments.scalar_lift import quadratic_sign_lift as lift

ROOT = Path(__file__).resolve().parents[1]
pytestmark = pytest.mark.skipif(shutil.which('cryptominisat5') is None,
                                reason='optional CryptoMiniSat executable is unavailable')


def test_weighted_binary_adders_against_exhaustive_integer_sums():
    for weights in ((1,1,1), (1,2,3), (2,2,4)):
        for values in itertools.product((0,1), repeat=3):
            expected = sum(w*v for w,v in zip(weights,values))
            for target in (expected, expected+1):
                problem = XorCNF()
                variables = [problem.var() for _ in values]
                for variable, value in zip(variables,values):
                    problem.clause([variable if value else -variable])
                problem.weighted_exact(list(zip(variables,weights)),target)
                result = problem.solve(2)
                assert (result['status']=='SAT') == (target==expected)


def test_weighted_complementary_literals_and_native_xor_semantics():
    for value in (0,1):
        for target in (2,3,4):
            problem = XorCNF()
            x = problem.var()
            problem.clause([x if value else -x])
            # 2x + 3(1-x) + x = 3 identically.
            problem.weighted_exact([(x,2),(-x,3),(x,1)],target)
            result = problem.solve(2)
            assert (result['status']=='SAT') == (target==3)
    problem = XorCNF()
    x,y = problem.var(),problem.var()
    # x XOR x XOR NOT y XOR TRUE = y.
    problem.xor([x,x,-y,problem.true],1)
    result = problem.solve(2)
    assert result['status']=='SAT' and result['values'][y]==1


def test_negative_head_gauge_recovers_a_full_t31_control():
    states = [0,1,2,2,7,6,2,6,2,5,3]
    data = frame.chart(31,3,5)
    system = lift.prepare([s//2 for s in states],data,-1,3,1)
    result = solve(system,seconds=3)
    assert result['status']=='SAT' and result['exact_zero']
    assert result['expanded_verification']['shift_count']==31
    assert result['expanded_verification']['sequence_row_sums']==[-2,3,3,3]


def test_complete_planted_t631_signs_recovered():
    record = json.loads((ROOT/'artifacts/hadamard-7572/t631.json').read_text())
    data = frame.chart(631,21,458)
    system = lift.prepare([s//2 for s in record['states']],data)
    result = solve(system,seconds=10)
    assert result['status']=='SAT' and result['exact_zero']
    assert result['states']==record['states']
    assert result['expanded_verification']['shift_count']==631

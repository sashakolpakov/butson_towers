"""Independent integer controls for the bounded unrestricted SAT backend."""
import itertools
import json
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
import pytest

from experiments.scalar_lift import unrestricted_t_sign_sat as sat

pytestmark = pytest.mark.skipif(shutil.which('cryptominisat5') is None,
                                reason='optional CryptoMiniSat executable is unavailable')


def test_weighted_adders_against_exhaustive_integer_sums(tmp_path):
    for weights in ((1, 1, 1), (1, 2, 3), (2, 2, 4)):
        for values in itertools.product((0, 1), repeat=3):
            expected = sum(w*v for w, v in zip(weights, values))
            for target in (expected, expected+1):
                problem = sat.XorCNF()
                variables = [problem.var() for _ in values]
                for variable, value in zip(variables, values):
                    problem.clause([variable if value else -variable])
                problem.weighted_exact(list(zip(variables, weights)), target)
                result = problem.solve(2, workdir=tmp_path)
                assert (result['status'] == 'SAT') == (target == expected)


def test_complementary_literals_and_native_xor_semantics(tmp_path):
    for value in (0, 1):
        for target in (2, 3, 4):
            problem = sat.XorCNF()
            x = problem.var()
            problem.clause([x if value else -x])
            problem.weighted_exact([(x, 2), (-x, 3), (x, 1)], target)
            result = problem.solve(2, workdir=tmp_path)
            assert (result['status'] == 'SAT') == (target == 3)
    problem = sat.XorCNF()
    x, y = problem.var(), problem.var()
    problem.xor([x, x, -y, problem.true], 1)
    result = problem.solve(2, workdir=tmp_path)
    assert result['status'] == 'SAT' and result['values'][y] == 1


def _brute_states(axes, target):
    """All signs with the origin positive, checked directly at every shift."""
    p = len(axes)
    found = set()
    for bits in itertools.product((0, 1), repeat=p-1):
        signs = [0, *bits]
        b = np.zeros((4, p), dtype=np.int64)
        b[axes, np.arange(p)] = [1-2*bit for bit in signs]
        if b.sum(axis=1).tolist() != list(target):
            continue
        if all(int(np.sum(b*np.roll(b, -s, axis=1))) == 0 for s in range(1, p)):
            found.add(tuple(2*a+bit for a, bit in zip(axes, signs)))
    return found


@pytest.mark.parametrize('lazy', (False, True))
def test_complete_unrestricted_t7_matches_all_signs(tmp_path, lazy):
    from experiments.scalar_lift import unrestricted_t_sign_model as lift
    states = [0, 6, 0, 4, 4, 5, 2]
    axes, target = [s//2 for s in states], (2, 1, 1, 1)
    expected = _brute_states(axes, target)
    assert tuple(states) in expected
    result = sat.solve_rows(axes, 7, target, seconds=5, max_solutions=0,
                            all_solutions=True, exact_correlations=not lazy,
                            workdir=tmp_path)
    assert result['status'] == 'SAT' and result['complete_for_fixed_axes']
    assert {tuple(item['states']) for item in result['solutions']} == expected
    for item in result['solutions']:
        assert lift.verify_states(item['states'], lift.periodic_chart(7), target)


def test_native_sat_preserves_86_independent_affine_bits(tmp_path):
    from experiments.scalar_lift import unrestricted_t_sign_model as lift
    n = 86
    groups = [list(range(i, n, 4)) for i in range(4)]
    targets = [3, 4, 5, 6]
    system = dict(axes=[i % 4 for i in range(n)], particular=0,
                  basis=[1 << i for i in range(n)], forms=[1 << i for i in range(n)],
                  groups=groups, negative_counts=targets,
                  diagonal=np.zeros(1, dtype=np.int64),
                  weights=np.zeros((1, n, n), dtype=np.int64))
    problem, free, _ = sat.encode_sign_system(system)
    result = problem.solve(3, workdir=tmp_path)
    assert len(free) == 86 and result['status'] == 'SAT'
    assignment = sum(result['values'][variable] << i for i, variable in enumerate(free))
    mask = lift.sign_mask(system, assignment)
    assert [sum((mask >> i) & 1 for i in group) for group in groups] == targets


@pytest.mark.parametrize('p,m,mu,folder', [
    (127, 9, 22, 't-matrices-127'),
    (631, 21, 458, 'hadamard-7572'),
])
def test_existing_planted_integer_systems_recover_signs(p, m, mu, folder, tmp_path):
    from experiments.scalar_lift import hadamard_frame_search as frame
    from experiments.scalar_lift import quadratic_sign_lift as old
    from experiments.scalar_lift import unrestricted_t_sign_model as lift
    root = Path(__file__).resolve().parents[1]
    states = json.loads((root/'artifacts'/folder/f't{p}.json').read_text())['states']
    data = frame.chart(p, m, mu)
    system = old.prepare([state//2 for state in states], data)
    system['groups'] = [system['heads'], system['tails'], [], []]
    system['negative_counts'] = [*system['negative_counts'], 0, 0]
    problem, free, _ = sat.encode_sign_system(system)
    result = problem.solve(5, workdir=tmp_path)
    assert result['status'] == 'SAT'
    assignment = sum(result['values'][variable] << i for i, variable in enumerate(free))
    mask = old.sign_mask(system, assignment)
    recovered = [2*a+((mask >> i) & 1) for i, a in enumerate(system['axes'])]
    assert recovered == states
    b = frame.expand(recovered, data)
    full = (2*np.argmax(abs(b), axis=0)+(b.sum(axis=0) < 0)).tolist()
    proof = lift.verify_states(full, p, b.sum(axis=1).tolist())
    assert proof['every_expanded_periodic_correlation_verified']


def test_arithmetic_unsat_and_bounded_unknown_are_distinct(tmp_path):
    # Seven positions cannot have coordinate sums of this absolute size.
    result = sat.solve_rows([0]*7, 7, (2, 1, 1, 1), seconds=1, workdir=tmp_path)
    assert result['status'] == 'UNSAT' and result['complete_for_fixed_axes']
    states = [0, 6, 0, 4, 4, 5, 2]
    result = sat.solve_rows([s//2 for s in states], 7, (2, 1, 1, 1),
                            seconds=1e-12, workdir=tmp_path)
    assert result['status'] == 'UNKNOWN' and not result['complete_for_fixed_axes']
    assert result['solutions'] == []


def test_prepared_small_pattern_has_a_complete_unsat_result(tmp_path):
    from experiments.scalar_lift import unrestricted_t_sign_model as lift
    axes, target = [0, 0, 0, 1, 0, 2, 3], (2, 1, 1, 1)
    system = lift.prepare(axes, 7, target)
    assert _brute_states(axes, target) == set()
    result = sat.solve(system, seconds=3, workdir=tmp_path)
    assert result['status'] == 'UNSAT' and result['complete_for_fixed_axes']
    assert result['attempts'][-1]['status'] == 'UNSAT'


def test_lazy_rejection_adds_a_full_integer_equation(tmp_path):
    # Abstract weighted-pair control: C=8 is zero modulo8 but not as an
    # integer. It is deliberately not a periodic construction candidate.
    weights = np.zeros((1, 3, 3), dtype=np.int64)
    weights[0, 1, 2] = weights[0, 2, 1] = 8
    system = dict(axes=[0, 1, 1], data={'p': 3}, particular=0, basis=[],
                  forms=[0, 0, 0], groups=[[0], [1, 2], [], []],
                  negative_counts=[0, 0, 0, 0], target_row=(1, 2, 0, 0),
                  diagonal=np.zeros(1, dtype=np.int64), weights=weights)
    result = sat.solve(system, seconds=3, workdir=tmp_path, exact_correlations=False)
    assert result['status'] == 'UNSAT' and result['complete_for_fixed_axes']
    assert result['nonzero_correlations_rejected'] == 1
    assert result['exact_correlations_encoded_finally'] == 1
    assert [attempt['status'] for attempt in result['attempts']] == ['SAT', 'UNSAT']


def test_real_cli_accepts_a_large_production_seed(tmp_path):
    root = Path(__file__).resolve().parents[1]
    axes_path, report_path = tmp_path/'axes.json', tmp_path/'report.json'
    axes_path.write_text(json.dumps({'axes': [0, 3, 0, 2, 2, 2, 1]}))
    process = subprocess.run([
        sys.executable, '-m', 'experiments.scalar_lift.unrestricted_t_sign_sat',
        '--axes', str(axes_path), '--p', '7', '--target-row', '2', '1', '1', '1',
        '--seed', '20261126167', '--seconds', '3', '--workdir', str(tmp_path/'native'),
        '--report', str(report_path),
    ], cwd=root, check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    assert json.loads(process.stdout)['status'] == 'SAT'
    result = json.loads(report_path.read_text())
    assert result['status'] == 'SAT' and result['solutions']
    assert result['solutions'][0]['expanded_verification'][
        'every_expanded_periodic_correlation_verified']

"""Independent finite controls for the unrestricted-width exact sign lift."""
import itertools
import json
from pathlib import Path
import random

import numpy as np
import pytest

from experiments.scalar_lift import hadamard_frame_search as frame
from experiments.scalar_lift import quadratic_sign_lift as lift

ROOT = Path(__file__).resolve().parents[1]


def test_complete_t19_sign_sets_and_quadratic_reduction():
    data = frame.chart(19, 3)
    patterns = total = affine_assignments = 0
    for additional_heads in itertools.combinations(range(2, 7), 2):
        tails = [i for i in range(2, 7) if i not in additional_heads]
        for colors in itertools.product((1, 2, 3), repeat=3):
            axes = [0]*7
            for i, color in zip(tails, colors):
                axes[i] = color
            unsigned = frame.expand([2*a for a in axes], data)
            if any(int(np.sum(unsigned*np.roll(unsigned, -s, axis=1))) % 2
                   for s in range(1, 19)):
                continue
            patterns += 1
            expected = set()
            for mask in range(1 << 5):
                states = [2*a+((mask >> (i-2)) & 1 if i >= 2 else 0)
                          for i, a in enumerate(axes)]
                rows = frame.expand(states, data)
                if rows.sum(axis=1).tolist() != [4, 1, 1, 1]:
                    continue
                if all(int(np.sum(rows*np.roll(rows, -s, axis=1))) == 0
                       for s in range(1, 19)):
                    expected.add(tuple(states))
            try:
                system = lift.prepare(axes, data)
            except lift.InfeasibleAxes:
                assert not expected
            else:
                equations = lift.quadratic_mod8(system)
                for assignment in range(1 << len(system['basis'])):
                    values = lift.correlations(system, lift.sign_mask(system, assignment))
                    assert all(v % 4 == 0 for v in values)
                    assert [lift.evaluate_quadratic(e, assignment) for e in equations] == [
                        (v//4) % 2 for v in values]
                    affine_assignments += 1
            result = lift.solve(axes, data, seconds=3, all_solutions=True, max_solutions=0)
            assert result['complete_for_fixed_axes']
            found = {tuple(item['states']) for item in result['solutions']}
            assert found == expected
            total += len(found)
    assert (patterns, total, affine_assignments) == (90, 150, 210)


@pytest.mark.parametrize('p,m,mu,folder', [
    (127, 9, 22, 't-matrices-127'),
    (631, 21, 458, 'hadamard-7572'),
])
def test_planted_signs_recovered_and_mod8_matches_integer_correlations(p, m, mu, folder, tmp_path):
    states = json.loads((ROOT/'artifacts'/folder/f't{p}.json').read_text())['states']
    data = frame.chart(p, m, mu)
    axes = [state//2 for state in states]
    system = lift.prepare(axes, data)
    equations = lift.quadratic_mod8(system)
    rng = random.Random(p)
    dimension = len(system['basis'])
    assignments = [0]+[1 << u for u in range(dimension)]+[
        rng.randrange(1 << dimension) for _ in range(80)]
    for assignment in assignments:
        values = lift.correlations(system, lift.sign_mask(system, assignment))
        assert all(v % 4 == 0 for v in values)
        assert [lift.evaluate_quadratic(e, assignment) for e in equations] == [
            (v//4) % 2 for v in values]
    result = lift.solve(axes, data, seconds=10, all_solutions=True, max_solutions=0)
    assert result['complete_for_fixed_axes']
    assert [item['states'] for item in result['solutions']] == [states]
    assert result['solutions'][0]['expanded_verification'][
        'every_expanded_periodic_correlation_verified']
    for block_variables in (6, 11):
        native = lift.solve_bits(axes, data, seconds=5, workdir=tmp_path,
                                 block_variables=block_variables)
        assert native['complete_for_fixed_axes']
        assert [item['states'] for item in native['solutions']] == [states]
        assert native['native']['assignments_covered'] == 1 << dimension


def test_binary_elimination_above_machine_word_width():
    variables = 131
    rows = [
        (1 << 129) | (1 << 65) | (1 << variables),
        (1 << 130) | (1 << 64),
        (1 << 100) | (1 << 2) | (1 << variables),
    ]
    particular, basis = lift._eliminate(rows, variables)
    mask = (1 << variables)-1
    assert len(basis) == 128
    for row in rows:
        assert ((particular & row & mask).bit_count() & 1) == ((row >> variables) & 1)
        assert all((vector & row & mask).bit_count() % 2 == 0 for vector in basis)
    with pytest.raises(lift.InfeasibleAxes):
        lift._eliminate([1, 1 | (1 << variables)], variables)


def test_negative_head_control_and_mq_only_complete_solution_set():
    states = [0, 1, 2, 2, 7, 6, 2, 6, 2, 5, 3]
    data = frame.chart(31, 3, 5)
    result = lift.solve([s//2 for s in states], data, -1, 3, 1,
                        seconds=5, all_solutions=True, max_solutions=0,
                        exact_correlations=False)
    assert result['complete_for_fixed_axes']
    found = {tuple(item['states']) for item in result['solutions']}
    assert tuple(states) in found
    for item in result['solutions']:
        assert item['expanded_verification']['sequence_row_sums'] == [-2, 3, 3, 3]


def test_native_streaming_mq_against_exhaustive_polynomials(tmp_path):
    import shutil
    import subprocess
    if shutil.which('c++') is None:
        pytest.skip('C++17 compiler unavailable')
    binary = lift.build_bits_binary(tmp_path)
    rng = random.Random(223)
    for d, count in [(0, 1), (1, 2), (5, 9), (6, 10), (7, 12),
                     (13, 14), (13, 70), (13, 131)]:
        planted = rng.randrange(1 << d)
        equations = []
        for _ in range(count):
            linear = rng.randrange(1 << d)
            rows = [rng.randrange(1 << u) for u in range(d)]
            constant = (linear & planted).bit_count() % 2
            for u in range(d):
                if (planted >> u) & 1:
                    constant ^= (rows[u] & planted).bit_count() % 2
            equations.append((constant, linear, rows))
        expected = set()
        for assignment in range(1 << d):
            satisfied = True
            for constant, linear, rows in equations:
                value = constant ^ ((linear & assignment).bit_count() % 2)
                for u in range(d):
                    if (assignment >> u) & 1:
                        value ^= (rows[u] & assignment).bit_count() % 2
                if value:
                    satisfied = False
                    break
            if satisfied:
                expected.add(assignment)
        lines = [f'{d} {count}']
        for constant, linear, rows in equations:
            lines.extend([f'{constant} {linear}', ' '.join(map(str, rows))])
        source, output = tmp_path/'input.txt', tmp_path/'output.json'
        source.write_text('\n'.join(lines)+'\n')
        # Exercise each block width and both sides of packed-equation word
        # boundaries. Scalar enumeration supplies the complete expected set.
        widths = range(6, 13) if d >= 13 else (11,)
        for block_variables in widths:
            subprocess.run([str(binary), str(source), str(output), '5', '0',
                            str(block_variables)], check=True,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            result = json.loads(output.read_text())
            assert result['complete']
            assert result['assignments_covered'] == 1 << d
            assert set(result['roots']) == expected
        assert planted in expected


def test_exact_head_branches_are_complete_and_preserve_solutions():
    for count, expected_count in [(5, 1), (7, 6), (9, 28)]:
        axes = [0]*(count+1)+[1]*(74-count)
        branches = list(lift.iter_head_sign_branches(axes, -5, 1))
        expected = {bits for bits in itertools.product((0, 1), repeat=count)
                    if bits[0] == 1 and count-2*sum(bits) == -5}
        assert {tuple(branch[i] for i in range(1, count+1)) for branch in branches} == expected
        assert len(branches) == expected_count
    axes = [0, 0, 0, 0, 1, 2, 3]
    data = frame.chart(19, 3)
    whole = lift.solve(axes, data, seconds=3, all_solutions=True, max_solutions=0)
    assert whole['complete_for_fixed_axes'] and len(whole['solutions']) == 6
    by_branches = set()
    for fixed in lift.iter_head_sign_branches(axes, 1, 0):
        result = lift.solve(axes, data, seconds=3, all_solutions=True, max_solutions=0,
                            fixed_signs=fixed)
        assert result['complete_for_fixed_axes']
        by_branches.update(tuple(item['states']) for item in result['solutions'])
    assert by_branches == {tuple(item['states']) for item in whole['solutions']}


def test_native_streaming_covers_all_32_bit_assignments(tmp_path):
    import shutil
    import subprocess
    if shutil.which('c++') is None:
        pytest.skip('C++17 compiler unavailable')
    binary = lift.build_bits_binary(tmp_path)
    dimension = 32
    planted = 0xD56A4C93
    lines = [f'{dimension} {dimension}']
    for index in range(dimension):
        lines.extend([f'{(planted >> index) & 1} {1 << index}',
                      ' '.join(['0']*dimension)])
    source, output = tmp_path/'wide.txt', tmp_path/'wide.json'
    source.write_text('\n'.join(lines)+'\n')
    subprocess.run([str(binary), str(source), str(output), '10', '0'], check=True,
                   stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    result = json.loads(output.read_text())
    assert result['complete'] and result['assignments_covered'] == 1 << dimension
    assert result['roots'] == [planted]


def test_native_root_limit_is_reported_as_unfinished(tmp_path):
    import subprocess
    binary = lift.build_bits_binary(tmp_path)
    source, output = tmp_path/'all-roots.txt', tmp_path/'limited.json'
    source.write_text('15 1\n0 0\n'+' '.join(['0']*15)+'\n')
    subprocess.run([str(binary), str(source), str(output), '5', '3', '11'],
                   check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    result = json.loads(output.read_text())
    assert not result['complete']
    assert result['root_limit_reached'] and not result['time_limit_reached']
    assert result['roots'] == [0, 1, 2]
    assert result['assignments_covered'] == 2048

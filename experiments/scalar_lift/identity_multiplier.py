"""Exact periodic T-designs constant on odd multiplier orbits.

The multiplier fixes all four coordinate labels. This geometry differs
from the three-cycle action in hadamard_frame_search. A validated binary
map parametrizes unsigned correlation parity by SL(2,2**d); signs still
require four exact counts and every integer periodic correlation.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
from math import isqrt
import operator
from pathlib import Path
import subprocess
import tempfile
import time

import numpy as np

from experiments.scalar_lift import hadamard_frame_search as base
from experiments.scalar_lift.large_parity_map import BinaryField
from experiments.scalar_lift.quadratic_sign_lift import (
    InfeasibleAxes, _eliminate, build_bits_binary, correlations,
    evaluate_quadratic, quadratic_mod8, sign_mask,
)


def _integer(value, name):
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f'{name} must be an integer, not a Boolean')
    try:
        return operator.index(value)
    except TypeError as exc:
        raise ValueError(f'{name} must be an integer') from exc


def chart(p, m, mu=None):
    """Build actual site orbits and all unordered correlation-pair weights."""
    p, m = _integer(p, 'prime length'), _integer(m, 'multiplier order')
    if (p < 3 or any(p % k == 0 for k in range(2, isqrt(p)+1)) or
            m < 1 or m % 2 != 1 or (p-1) % m):
        raise ValueError('an odd prime and an odd divisor of p-1 are required')
    factors = [q for q in range(2, m+1) if m % q == 0 and
               all(q % k for k in range(2, isqrt(q)+1))]
    def has_order(value):
        return pow(value, m, p) == 1 and all(pow(value, m//q, p) != 1 for q in factors)
    if mu is None:
        mu = next(value for value in range(1, p) if has_order(value))
    mu = _integer(mu, 'multiplier')
    if not 1 <= mu < p or not has_order(mu):
        raise ValueError('the multiplier must have the exact stated order')
    index = np.zeros(p, dtype=np.int64)
    seen, orbits = {0}, [[(0, 0)]]
    for representative in range(1, p):
        if representative in seen:
            continue
        orbit, value = [], representative
        for _ in range(m):
            if value in seen:
                raise ArithmeticError('multiplier orbits overlap')
            seen.add(value)
            index[value] = len(orbits)
            orbit.append((value, 0))
            value = value*mu % p
        if value != representative:
            raise ArithmeticError('multiplier orbit did not close')
        orbits.append(orbit)
    shifts, covered = [], {0}
    for shift in range(1, p):
        if shift in covered:
            continue
        shifts.append(shift)
        value = shift
        for _ in range(m):
            covered.update((value, -value % p))
            value = value*mu % p
    rows = []
    for shift in shifts:
        counts = defaultdict(int)
        for x in range(p):
            i, j = sorted((int(index[x]), int(index[(x+shift) % p])))
            counts[i, j, 0] += 1
        rows.append(dict(counts))
    return dict(chart_kind='identity-multiplier', p=p, m=m, mu=mu,
                n=len(orbits), index=index, phase=np.zeros(p, dtype=np.int64),
                orbits=orbits, shifts=shifts, rows=rows)


identity_chart = chart


def target_rows(p, m):
    """All compatible row sums up to signs/permutations of coordinates1..3.

    The origin is +e0, so row0=1 modulo m and the other rows are multiples
    of m. Orbit support counts must realize the corresponding signed nets.
    This is necessary arithmetic, not a design existence assertion.
    """
    data = chart(p, m)
    result = []
    for a in range(-isqrt(p), isqrt(p)+1):
        if (a-1) % m:
            continue
        for b in range(0, isqrt(p-a*a)+1, m):
            for c in range(0, min(b, isqrt(p-a*a-b*b))+1, m):
                rest = p-a*a-b*b-c*c
                d = isqrt(rest)
                if d > c or d*d != rest or d % m:
                    continue
                nets = ((a-1)//m, b//m, c//m, d//m)
                if sum(abs(x) for x in nets) <= data['n']-1 and (data['n']-1-sum(nets)) % 2 == 0:
                    result.append((a, b, c, d))
    return sorted(result, key=lambda row: (-abs(row[0]), row))


def normalize_head(axes, data):
    """Send a nonzero coordinate-zero orbit to residue1 by x -> r*x."""
    axes = list(map(lambda a: _integer(a, 'coordinate label'), axes))
    if len(axes) != data['n'] or axes[0] != 0 or any(a not in range(4) for a in axes):
        raise ValueError('one four-color label per orbit, with zero first, is required')
    selected = next((i for i in range(1, len(axes)) if axes[i] == 0), None)
    if selected is None:
        raise InfeasibleAxes('no nonzero head to normalize')
    scale = data['orbits'][selected][0][0]
    return [axes[int(data['index'][scale*orbit[0][0] % data['p']])]
            for orbit in data['orbits']]


def validate_map(config):
    """Check every quadratic coefficient and every expanded nonzero shift.

    The asserted identity is F(T(P,Q,U,V))=L(PV+QU+1), with invertible
    binary maps T,L. Degree at most two makes these coefficient checks a
    proof for all binary inputs, rather than random testing of the map.
    """
    started = time.monotonic()
    if config.get('action') != 'identity':
        raise ValueError('an explicitly identity-action map is required')
    data = chart(config['p'], config['m'], config['mu'])
    degree = _integer(config['field_degree'], 'field degree')
    field = BinaryField(degree, _integer(config['field_modulus'], 'field modulus'))
    columns = config['inverse_columns']
    n = 4*degree
    if (degree < 1 or not field.irreducible() or data['n'] != 2*degree+1 or
            len(data['rows']) != degree or len(columns) != n or
            any(type(c) is not int or not 0 <= c < 1 << n for c in columns) or
            base.binary_rank(columns) != n):
        raise ValueError('field, orbit dimensions or full binary coordinate rank failed')
    pairs = [(u, v) for u in range(n) for v in range(u)]
    inputs = [0]+[1 << u for u in range(n)]+[(1 << u) | (1 << v) for u, v in pairs]
    colors = [base.binary_linear(columns, value) for value in inputs]
    axes = np.asarray([[0]+[(c >> (2*j)) & 3 for j in range(2*degree)] for c in colors], dtype=np.uint8)
    values = np.zeros((len(inputs), degree), dtype=np.uint8)
    for ell, row in enumerate(data['rows']):
        for (i, j, delta), count in row.items():
            if delta != 0:
                raise ArithmeticError('identity geometry cannot rotate coordinates')
            if count & 1:
                values[:, ell] ^= axes[:, i] == axes[:, j]
    parities = [sum(int(b) << ell for ell, b in enumerate(row)) for row in values]
    if any(value != parities[0] for value in parities[1:1+n]):
        raise ValueError('a linear coefficient of the determinant identity failed')
    lookup = {pair:k for k, pair in enumerate(pairs)}
    def polar(u, v):
        u, v = max(u, v), min(u, v)
        return parities[1+n+lookup[u, v]] ^ parities[1+u] ^ parities[1+v] ^ parities[0]
    output = [polar(3*degree+j, 0) for j in range(degree)]
    if base.binary_rank(output) != degree or output[0] != parities[0]:
        raise ValueError('constant coefficient or output-map rank failed')
    for u, v in pairs:
        expected = 0
        if ((3*degree <= u < 4*degree and v < degree) or
                (2*degree <= u < 3*degree and degree <= v < 2*degree)):
            expected = field.mul(1 << (u % degree), 1 << (v % degree))
        if polar(u, v) != base.binary_linear(output, expected):
            raise ValueError(f'quadratic coefficient {(u,v)} failed')
    shift_class = {}
    for ell, shift in enumerate(data['shifts']):
        value = shift
        for _ in range(data['m']):
            shift_class[value] = shift_class[-value % data['p']] = ell
            value = value*data['mu'] % data['p']
    actual = axes[:, data['index']]
    for shift in range(1, data['p']):
        observed = np.sum(actual == np.roll(actual, -shift, axis=1), axis=1) % 2
        if not np.array_equal(observed, values[:, shift_class[shift]]):
            raise ValueError(f'expanded coefficient geometry failed at shift{shift}')
    return dict(coefficient_inputs=len(inputs), color_map_rank=n, equation_map_rank=degree,
                all_coefficients_verified=True, all_expanded_shifts_verified=True,
                expanded_correlations_checked=len(inputs)*(data['p']-1),
                seconds=time.monotonic()-started)


def prepare(axes, data, target_row, *, fixed_signs=None, force_extremes=True):
    """Exact mod4 signs with all four count parities and extreme counts."""
    axes = tuple(_integer(a, 'coordinate label') for a in axes)
    n, length = data['n'], len(data['rows'])
    if data.get('chart_kind') != 'identity-multiplier':
        raise ValueError('explicit identity-action orbit data are required')
    if len(axes) != n or axes[0] != 0 or any(a not in range(4) for a in axes):
        raise ValueError('one coordinate label per orbit and positive origin are required')
    target = tuple(_integer(r, 'row sum') for r in target_row)
    m = data['m']
    if (len(target) != 4 or sum(r*r for r in target) != data['p'] or
            (target[0]-1) % m or any(r % m for r in target[1:])):
        raise ValueError('row sums must have norm p and the multiplier congruences')
    orbit_targets = (1+(target[0]-1)//m,)+tuple(r//m for r in target[1:])
    weights = np.zeros((length, n, n), dtype=np.int32)
    diagonal = np.zeros(length, dtype=np.int32)
    for ell, row in enumerate(data['rows']):
        for (i, j, delta), count in row.items():
            if delta != 0:
                raise ValueError('an identity-action pair cannot rotate coordinates')
            if axes[i] != axes[j]:
                continue
            if i == j:
                diagonal[ell] += count
            else:
                weights[ell, i, j] += count
                weights[ell, j, i] += count
    absolute = diagonal+weights.sum(axis=(1, 2))//2
    if np.any(absolute % 2):
        raise InfeasibleAxes('unsigned correlation parity failed')
    fixed = {0:0}
    for i, bit in (fixed_signs or {}).items():
        i, bit = _integer(i, 'fixed orbit'), _integer(bit, 'fixed sign')
        if not 0 <= i < n or bit not in (0, 1):
            raise ValueError('fixed signs need valid orbit indices and binary values')
        if i in fixed and fixed[i] != bit:
            raise InfeasibleAxes('a fixed sign contradicts the positive origin')
        fixed[i] = bit
    groups = [[i for i, a in enumerate(axes) if a == j] for j in range(4)]
    wanted = []
    for group, net in zip(groups, orbit_targets, strict=True):
        if len(group) < abs(net) or (len(group)-net) % 2:
            raise InfeasibleAxes('support sizes cannot realize the requested sums')
        count = (len(group)-net)//2
        wanted.append(count)
        unknown = [i for i in group if i not in fixed]
        needed = count-sum(fixed[i] for i in group if i in fixed)
        if not 0 <= needed <= len(unknown):
            raise InfeasibleAxes('fixed signs contradict an exact coordinate count')
        if force_extremes and needed in (0, len(unknown)):
            fixed.update((i, int(needed != 0)) for i in unknown)
    equations = []
    for ell in range(length):
        degrees = weights[ell].sum(axis=1)
        equations.append(sum(int(d % 2) << i for i, d in enumerate(degrees)) |
                         (int(absolute[ell]//2 % 2) << n))
    equations.extend((1 << i) | (bit << n) for i, bit in fixed.items())
    equations.extend(sum(1 << i for i in group) | ((count % 2) << n)
                     for group, count in zip(groups, wanted, strict=True))
    particular, basis = _eliminate(equations, n)
    forms = [sum(((vector >> i) & 1) << j for j, vector in enumerate(basis)) for i in range(n)]
    return dict(axes=axes, data=data, weights=weights, diagonal=diagonal,
                absolute_correlations=absolute, particular=particular, basis=basis,
                forms=forms, groups=groups, negative_counts=wanted, target_row=target,
                forced_signs=fixed)


def exact_count_mask(system, mask):
    mask = _integer(mask, 'negative-sign mask')
    if not 0 <= mask < 1 << len(system['axes']):
        raise ValueError('negative-sign mask is outside the orbit range')
    return all(sum((mask >> i) & 1 for i in group) == target
               for group, target in zip(system['groups'], system['negative_counts'], strict=True))


def expanded_states(states, data):
    """Expand each signed coordinate orbit without any coordinate rotation."""
    states = tuple(_integer(s, 'signed coordinate state') for s in states)
    if len(states) != data['n'] or states[0] != 0 or any(s not in range(8) for s in states):
        raise ValueError('one state in0..7 per orbit and positive origin are required')
    values = np.asarray(states, dtype=np.int64)[data['index']]
    rows = np.zeros((4, data['p']), dtype=np.int64)
    rows[values//2, np.arange(data['p'])] = 1-2*(values % 2)
    return rows


def verify_states(states, data, target_row=None):
    """Rebuild the site geometry independently and check all integer shifts."""
    fresh = chart(data['p'], data['m'], data['mu'])
    rows = expanded_states(states, fresh)
    values = [int(np.sum(rows*np.roll(rows, -s, axis=1), dtype=np.int64)) for s in range(fresh['p'])]
    if values != [fresh['p']]+[0]*(fresh['p']-1):
        raise ValueError('candidate failed an expanded integer periodic correlation')
    sums = list(map(int, rows.sum(axis=1)))
    if target_row is not None:
        target = [_integer(r, 'row sum') for r in target_row]
        if len(target) != 4 or sums != target:
            raise ValueError('candidate failed the four requested coordinate sums')
    return dict(every_expanded_periodic_correlation_verified=True, shift_count=fresh['p'],
                sequence_row_sums=sums, all_integer_correlations=values,
                sequence_support_sizes=list(map(int, np.count_nonzero(rows, axis=1))))


def solve_bits(axes, data, target_row, *, first_negative=None, fixed_signs=None,
               seconds=10, max_roots=100000, workdir=None, block_variables=11):
    """Completely filter bounded MQ signs, then test exact counts and norms."""
    max_roots = _integer(max_roots, 'root limit')
    block_variables = _integer(block_variables, 'block variables')
    if seconds <= 0 or max_roots < 0 or not 6 <= block_variables <= 12:
        raise ValueError('positive time, nonnegative root limit and block size6..12 required')
    started = time.monotonic()
    fixed = dict(fixed_signs or {})
    if first_negative is not None:
        first_negative = _integer(first_negative, 'first-head sign')
        if first_negative not in (0, 1) or len(axes) < 2 or axes[1] != 0:
            raise ValueError('the optional first-head gauge needs a head and a binary sign')
        if 1 in fixed and fixed[1] != first_negative:
            raise InfeasibleAxes('additional fixed signs contradict the first-head gauge')
        fixed[1] = first_negative
    try:
        system = prepare(axes, data, target_row, fixed_signs=fixed)
    except InfeasibleAxes as exc:
        return dict(solver_status='INFEASIBLE_ARITHMETIC', reason=str(exc),
                    complete_for_fixed_axes=True, solutions=[], elapsed_seconds=time.monotonic()-started,
                    scope='Fixed identity-action axes, target row and supplied fixed signs only.')
    dimension = len(system['basis'])
    if dimension > 60:
        raise ValueError('the streaming native filter supports at most60 free variables')
    equations = quadratic_mod8(system)
    with tempfile.TemporaryDirectory(prefix='identity-sign-bits-') as temporary:
        directory = Path(temporary) if workdir is None else Path(workdir)
        directory.mkdir(parents=True, exist_ok=True)
        binary = build_bits_binary(directory)
        source, target = directory/'mq-input.txt', directory/'mq-result.json'
        lines = [f'{dimension} {len(equations)}']
        for equation in equations:
            rows = [0]*dimension
            for u, v in equation['quadratic']:
                rows[v] |= 1 << u
            lines.extend((f"{equation['constant']} {equation['linear']}", ' '.join(map(str, rows))))
        source.write_text('\n'.join(lines)+'\n')
        subprocess.run([str(binary), str(source), str(target), str(seconds), str(max_roots), str(block_variables)],
                       check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        native = json.loads(target.read_text())
    solutions, passed_counts = [], 0
    for assignment in native['roots']:
        if any(evaluate_quadratic(eq, assignment) for eq in equations):
            raise ArithmeticError('native MQ root failed independent evaluation')
        mask = sign_mask(system, assignment)
        if not exact_count_mask(system, mask):
            continue
        passed_counts += 1
        if any(correlations(system, mask)):
            continue
        states = [2*a+((mask >> i) & 1) for i, a in enumerate(axes)]
        proof = verify_states(states, data, target_row)
        solutions.append(dict(states=states, expanded_verification=proof))
    return dict(solver_status='COMPLETE' if native['complete'] else 'UNFINISHED',
                complete_for_fixed_axes=native['complete'], affine_dimension=dimension,
                forced_sign_count=len(system['forced_signs']), quadratic_mod8_equations=len(equations),
                exact_negative_counts=system['negative_counts'], mod8_roots=len(native['roots']),
                exact_count_survivors=passed_counts, solutions=solutions, native=native,
                elapsed_seconds=time.monotonic()-started,
                scope='Fixed identity-action coordinate pattern, four exact sums and stated gauges only.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(validate_map(json.loads(args.config.read_text()))))


if __name__ == '__main__':
    main()

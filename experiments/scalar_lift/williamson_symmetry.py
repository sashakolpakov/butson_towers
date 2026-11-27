"""Symmetric Williamson sources with residue-multiplier covariance.

At residue zero the four sign functions have value (1,1,1,1). Off zero,
their product is negative, so their vectors are signed columns of H=J-2I.
An axis assignment on the multiplier orbits determines all remaining signs
through an exact congruence. Only the full integer source verifier certifies
a completed source. Neither a chart nor its row-sum screen asserts existence.
"""
from __future__ import annotations

from collections import Counter
import argparse
import hashlib
from itertools import combinations_with_replacement
import json
from math import isqrt
import operator
from pathlib import Path

import numpy as np

H = np.ones((4, 4), dtype=np.int64) - 2*np.eye(4, dtype=np.int64)
H.setflags(write=False)
KIND = 'symmetric-Williamson-multiplier'


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _integer(value, name):
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(name+' must be an integer, not a Boolean')
    try:
        return operator.index(value)
    except TypeError as exc:
        raise ValueError(name+' must be an integer') from exc


def _prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n)+1))


def _prime_divisors(n):
    return [d for d in range(2, n+1) if n % d == 0 and _prime(d)]


def rotate(axis, phase):
    """Apply the stated three-cycle, fixing coordinate zero."""
    axis, phase = _integer(axis, 'axis'), _integer(phase, 'phase')
    _require(0 <= axis < 4, 'axis must be in 0..3')
    return 0 if axis == 0 else 1+(axis-1+phase) % 3


def chart(p, m, action, mu=None):
    """Build literal orbits and verify every expanded coefficient table.

    The identity action permits any even divisor m of p-1. The three-cycle
    action requires 6|m, so the element -1 acts trivially in both cases.
    The zero residue is separate; N counts only the nonzero orbits.
    """
    p, m = _integer(p, 'prime length'), _integer(m, 'multiplier order')
    _require(p > 2 and _prime(p), 'an odd prime length is required')
    _require(m >= 2 and m % 2 == 0 and (p-1) % m == 0,
             'an even divisor of p-1 is required')
    _require(action in ('identity', 'C3'), 'coordinate action must be identity or C3')
    _require(action != 'C3' or m % 6 == 0,
             'the three-cycle must be trivial on the inversion element')
    factors = _prime_divisors(m)

    def correct_order(value):
        return (pow(value, m, p) == 1 and
                all(pow(value, m//q, p) != 1 for q in factors))

    if mu is None:
        mu = next(value for value in range(2, p) if correct_order(value))
    mu = _integer(mu, 'multiplier')
    _require(1 <= mu < p and correct_order(mu), 'multiplier has the wrong exact order')
    _require(pow(mu, m//2, p) == p-1, 'inversion is missing from the multiplier group')
    index = np.full(p, -1, dtype=np.int64)
    phase = np.zeros(p, dtype=np.int64)
    seen, orbits = {0}, []
    for representative in range(1, p):
        if representative in seen:
            continue
        orbit, value = [], representative
        for k in range(m):
            _require(value not in seen, 'multiplier orbits overlap')
            seen.add(value)
            index[value] = len(orbits)
            phase[value] = k % 3 if action == 'C3' else 0
            orbit.append(value)
            value = value*mu % p
        _require(value == representative, 'multiplier orbit did not close')
        orbits.append(orbit)
    _require(len(seen) == p, 'multiplier orbits do not cover every residue')

    def counts(shift):
        row = Counter()
        for x in range(1, p):
            y = (x+shift) % p
            if y == 0:
                continue
            i, j = int(index[x]), int(index[y])
            delta = int(phase[y]-phase[x]) % 3
            if i > j:
                i, j, delta = j, i, -delta % 3
            row[i, j, delta] += 1
        return row

    rows = [counts(orbit[0]) for orbit in orbits]
    for shift in range(1, p):
        _require(counts(shift) == rows[int(index[shift])],
                 'expanded shift coefficients disagree with their orbit representative')
    for row in rows:
        _require(all(count % 2 == 0 for (i,j,_),count in row.items() if i != j),
                 'an off-diagonal pair coefficient is odd')
        _require(sum(count for (i,j,d),count in row.items() if i == j and d == 0) % 2 == 1,
                 'the fixed-point diagonal coefficient is not odd')
        for i in range(len(orbits)):
            _require((row[i,i,1]+row[i,i,2]) % 2 == 0,
                     'the rotated diagonal coefficients have incorrect parity')
    index.setflags(write=False)
    phase.setflags(write=False)
    return dict(chart_kind=KIND, p=p, m=m, mu=mu, action=action, N=len(orbits),
                index=index, phase=phase, orbits=orbits, rows=rows,
                geometry_audit=dict(all_nonzero_shift_tables_verified=p-1,
                                    all_offdiagonal_coefficients_even=True,
                                    diagonal_parity_verified=True),
                normalization=dict(value_at_zero=[1,1,1,1],
                                   off_zero_vectors='signed columns of J4-2I4',
                                   weighted_value_at_zero=['1/2']*4,
                                   weighted_coordinate_norm_scale=p,
                                   sign_source_summed_norm_scale=4*p))


def _axes(axes, data):
    _require(data.get('chart_kind') == KIND, 'a verified Williamson chart is required')
    axes = tuple(_integer(a, 'axis') for a in axes)
    _require(len(axes) == data['N'] and all(0 <= a < 4 for a in axes),
             'one axis in 0..3 is required for every nonzero orbit')
    return axes


def coefficients(axes, data):
    """Return D,W,A and the unique signs satisfying the normalized mod4 cut.

    For shift orbit ell, C_F/4 = sigma_ell + D_ell +
    sum_(i<j) W_ell,i,j sigma_i sigma_j. The off-diagonal W are even and
    A=D+sum W is odd, giving sigma_ell=(-1)^((A_ell+1)/2).
    """
    axes = _axes(axes, data)
    n = data['N']
    diagonal = np.zeros(n, dtype=np.int64)
    weights = np.zeros((n,n,n), dtype=np.int64)
    for ell, row in enumerate(data['rows']):
        for (i,j,phase),count in row.items():
            if axes[i] != rotate(axes[j], phase):
                continue
            if i == j:
                diagonal[ell] += count
            else:
                weights[ell,i,j] += count
    absolute = diagonal+weights.sum(axis=(1,2))
    _require(np.all(absolute % 2 == 1) and np.all(weights % 2 == 0),
             'forced-sign coefficient parity failed')
    signs = 1-2*((absolute+1)//2 % 2)
    return diagonal, weights, absolute, signs


def expand(axes, signs, data):
    """Expand actual sign functions, retaining the special all-ones value at zero."""
    axes = _axes(axes, data)
    signs = tuple(_integer(sign, 'sign') for sign in signs)
    _require(len(signs) == data['N'] and all(s in (-1,1) for s in signs),
             'one sign is required for every nonzero orbit')
    rows = np.ones((4,data['p']), dtype=np.int64)
    for x in range(1, data['p']):
        i = int(data['index'][x])
        rows[:,x] = signs[i]*H[:,rotate(axes[i], int(data['phase'][x]))]
    return rows


def verify_source(rows, data=None):
    """Verify every literal integer source identity independently of the solver.

    Complete cyclic coefficients prove the full circulant matrix Grams.
    Dense matrix products are unnecessary; small independent tests also
    materialize those products. A supplied chart additionally checks the
    actual multiplier covariance, without trusting compressed Gram values.
    """
    rows = np.asarray(rows)
    _require(rows.ndim == 2 and rows.shape[0] == 4 and rows.shape[1] > 0,
             'four nonempty sign rows are required')
    _require(np.issubdtype(rows.dtype, np.integer) and np.all(abs(rows) == 1),
             'source entries must be integer signs')
    rows = rows.astype(np.int64)
    p = rows.shape[1]
    _require(p % 2 == 1 and np.all(rows[:,0] == 1),
             'an odd length and normalized all-ones value at zero are required')
    _require(np.array_equal(rows[:,(-np.arange(p)) % p], rows),
             'literal Williamson symmetry failed')
    _require(np.all(np.prod(rows[:,1:], axis=0) == -1),
             'literal Williamson product condition failed')
    if data is not None:
        _require(data.get('chart_kind') == KIND and data['p'] == p,
                 'source and multiplier chart do not agree')
        permutation = [0,3,1,2] if data['action'] == 'C3' else [0,1,2,3]
        _require(np.array_equal(rows[:,(data['mu']*np.arange(p)) % p], rows[permutation]),
                 'literal multiplier covariance failed')
    correlations = [int(np.sum(rows*np.roll(rows,-s,axis=1), dtype=np.int64))
                    for s in range(p)]
    _require(correlations == [4*p]+[0]*(p-1), 'complete integer Williamson norm failed')
    for i in range(4):
        for j in range(i):
            forward = [int(np.dot(rows[i], np.roll(rows[j],-s))) for s in range(p)]
            reverse = [int(np.dot(rows[j], np.roll(rows[i],-s))) for s in range(p)]
            _require(forward == reverse, 'complete amicability coefficients failed')
    return dict(row_sums=rows.sum(axis=1).tolist(), entry_alphabet=[-1,1],
                component_order=p, vector_dimension=4,
                all_integer_correlations=correlations,
                all_expanded_periodic_correlations_verified=True,
                sign_entries=True, symmetry=True, product_theorem_condition=True,
                pairwise_amicability=True, amicability_coefficients_checked=6*p,
                full_sum_of_matrix_Grams=True, dense_matrix_products_performed=False,
                source_Gram_identity=f'sum_j W_j W_j^T={4*p}I_{p}',
                source_summed_Gram_scale=4*p,
                weighted_coordinate_correlation_scale=p,
                coordinate_transform='F=(J4-2I4)b; b(0)=ones/2 and b(x)=signed unit axis off zero',
                transformation_Gram_scale=4,
                actual_off_zero_vectors=sorted({tuple(map(int,rows[:,x])) for x in range(1,p)}))


def coefficient_record(data):
    """Serialize the verified literal geometry for a frozen solver input."""
    _require(data.get('chart_kind') == KIND, 'a verified Williamson chart is required')
    result = {key:data[key] for key in ('p','m','mu','action','N','orbits')}
    result.update(index=data['index'].tolist(), phase=data['phase'].tolist(),
                  rows=[[[i,j,d,count] for (i,j,d),count in sorted(row.items())]
                        for row in data['rows']],
                  normalization=dict(source_at_zero=[1,1,1,1],
                                     source_off_zero='sigma*(ones-2*e_axis)',
                                     compressed_Gram_scale=4*data['p'],
                                     row_equation_scale=4, entry_alphabet=[-1,1]),
                  geometry_audit=data['geometry_audit'],
                  sign_elimination='A_ell=sum(count*equal_rotated_axes); sigma_ell=1-2*(((A_ell+1)//2)%2)')
    return result


def _row_parameters(p, m, n, action):
    p, m, n = (_integer(p, 'prime length'), _integer(m, 'multiplier order'),
               _integer(n, 'orbit count'))
    _require(p > 2 and _prime(p), 'an odd prime length is required')
    _require(m >= 2 and m % 2 == 0 and (p-1) % m == 0 and n == (p-1)//m,
             'row-count parameters do not describe the multiplier orbits')
    _require(action != 'C3' or m % 6 == 0, 'three-cycle row counts require 6|m')
    return p, m, n


def identity_row_types(p, m, n):
    """Complete row-norm and joint sign-state counts, up to coordinate permutations."""
    p, m, n = _row_parameters(p, m, n, 'identity')
    bound = isqrt(4*p)
    values = [r for r in range(-bound,bound+1) if (r-1) % m == 0
              and abs((r-1)//m) <= n and ((r-1)//m-n) % 2 == 0]
    result = []
    for row in combinations_with_replacement(values,4):
        if sum(r*r for r in row) != 4*p:
            continue
        nets = [(r-1)//m for r in row]
        if sum(nets) % 2:
            continue
        total = sum(nets)//2
        if abs(total) > n or (total-n) % 2 or any((total-k) % 2 for k in nets):
            continue
        groups = [(total-k)//2 for k in nets]
        minimum = sum(abs(g) for g in groups)
        if minimum <= n and (n-minimum) % 2 == 0:
            result.append(dict(row_sums=list(row), signed_axis_sums=groups,
                               total_orbit_sign_sum=total))
    return result


def c3_row_types(p, m, n):
    """Complete three-cycle row norms and joint signed-axis orbit counts."""
    p, m, n = _row_parameters(p, m, n, 'C3')
    result = []
    unit = m//3
    for head in range(-isqrt(4*p),isqrt(4*p)+1):
        if (head-1) % m:
            continue
        a = (head-1)//m
        if abs(a) > n or (a-n) % 2:
            continue
        remaining = 4*p-head*head
        if remaining % 3:
            continue
        tail = isqrt(remaining//3)
        if tail*tail != remaining//3:
            continue
        for value in sorted({tail,-tail}):
            if (value-1) % unit:
                continue
            b = (value-1)//unit
            if abs(b) > 3*n or (b-n) % 2 or (b-a) % 4 or (3*a+b) % 4:
                continue
            h, t = (b-a)//4, (3*a+b)//4
            if abs(h)+abs(t) <= n and (n-abs(h)-abs(t)) % 2 == 0:
                result.append(dict(row_sums=[head,value,value,value],
                                   signed_head_axis_sum=h, signed_tail_axis_sum=t,
                                   total_orbit_sign_sum=h+t))
    return result


def select_charts(orders, known_t_lengths, *, max_orbits=80, order_ceiling=10000):
    """Screen prime sources, then retain the strongest admissible multiplier per prime.

    This screens only the stated identity and three-cycle actions, row norms
    and finite sign-state counts. For the retained 81-order frontier every
    prime route is3 modulo4, so its quotient by inversion has odd order and
    these are the possible coordinate actions. No correlations are solved.
    """
    max_orbits = _integer(max_orbits, 'orbit limit')
    order_ceiling = _integer(order_ceiling, 'order ceiling')
    _require(max_orbits > 0 and order_ceiling > 0, 'positive selection limits are required')
    lengths = sorted({_integer(t, 'known T length') for t in known_t_lengths})
    _require(lengths and all(t > 0 for t in lengths), 'positive known T lengths are required')
    routes = {}
    for value in orders:
        order = _integer(value, 'target order')
        _require(order > 0 and order % 4 == 0, 'target orders must be positive multiples of four')
        if order >= order_ceiling:
            continue
        for t in lengths:
            if order % (4*t):
                continue
            p = order//(4*t)
            if p > 2 and _prime(p):
                route = dict(target_order=order, known_T_length=t)
                if route not in routes.setdefault(p,[]):
                    routes[p].append(route)
    candidates = []
    for p, outputs in routes.items():
        for m in range(2,p,2):
            if (p-1) % m or (p-1)//m > max_orbits:
                continue
            n = (p-1)//m
            actions = [('identity', identity_row_types(p,m,n))]
            if m % 6 == 0:
                actions.append(('C3', c3_row_types(p,m,n)))
            for action, row_types in actions:
                if row_types:
                    candidates.append(dict(p=p,m=m,nonzero_orbits=n,action=action,
                                           row_types=row_types,target_routes=sorted(outputs,key=lambda r:r['target_order'])))
    strongest = []
    for p in sorted({r['p'] for r in candidates}):
        minimum = min(r['nonzero_orbits'] for r in candidates if r['p'] == p)
        strongest.extend(r for r in candidates if r['p'] == p and r['nonzero_orbits'] == minimum)
    return sorted(strongest,key=lambda r:(r['nonzero_orbits'],r['p'],r['action']))


def main(argv=None):
    """Generate frozen coefficients or replay the compact selection certificate."""
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--p',type=int)
    parser.add_argument('--m',type=int)
    parser.add_argument('--mu',type=int)
    parser.add_argument('--action',choices=('identity','C3'))
    parser.add_argument('--output',type=Path)
    parser.add_argument('--check-selection',type=Path,nargs='?',
                        const=root/'certificates/scalar-lift/williamson-symmetry-selection.json')
    args = parser.parse_args(argv)
    try:
        if args.check_selection is not None:
            from experiments.scalar_lift.williamson_scalar_endpoints import verify_t_rows
            record = json.loads(args.check_selection.read_text())
            controls_path = root/record['known_T_inputs_certificate']
            controls_raw = controls_path.read_bytes()
            _require(hashlib.sha256(controls_raw).hexdigest() == record['known_T_inputs_sha256'],
                     'known T-input certificate hash differs')
            controls = json.loads(controls_raw)
            for item in controls['known_T_inputs']:
                verify_t_rows(np.asarray(item['rows'],dtype=np.int64))
            orders = record['frontier']['orders']
            _require(len(orders) == len(set(orders)) == 81 and orders == sorted(orders),
                     'the retained 81-order frontier is malformed')
            order_hash = hashlib.sha256(json.dumps(orders,separators=(',',':')).encode()).hexdigest()
            _require(order_hash == record['frontier']['ordered_list_sha256'],
                     'frontier order-list hash differs')
            selected = select_charts(orders,[r['t'] for r in controls['known_T_inputs']],
                                     max_orbits=record['limits']['max_nonzero_orbits'],
                                     order_ceiling=record['limits']['target_order_exclusive'])
            _require(selected == record['charts'],'selection certificate differs from exact row-count replay')
            print(json.dumps(dict(ok=True,prime_inputs=len({r['p'] for r in selected}),
                                  charts=len(selected),
                                  target_orders=sorted({x['target_order'] for r in selected for x in r['target_routes']}))))
            return 0
        if None in (args.p,args.m,args.action) or args.output is None:
            parser.error('coefficient generation requires --p, --m, --action and --output')
        data = chart(args.p,args.m,args.action,args.mu)
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(coefficient_record(data),separators=(',',':'))+'\n')
        print(json.dumps(dict(ok=True,p=data['p'],m=data['m'],mu=data['mu'],
                              nonzero_orbits=data['N'],output=str(args.output))))
        return 0
    except (ValueError,OSError,KeyError,TypeError) as exc:
        parser.exit(1,'Williamson symmetry audit failed: '+str(exc)+'\n')


if __name__ == '__main__':
    raise SystemExit(main())

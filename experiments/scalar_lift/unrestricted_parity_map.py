"""Exact SL2 coordinates for unrestricted periodic correlation parities.

There is one coordinate label at each residue, with label zero at the
origin. Labels 0,1,2,3 represent 0,1,w,1+w in GF(4), w^2+w+1=0.
The four field-coordinate blocks are P,Q,U,V, in that order. An exported
binary isomorphism T satisfies F(T(P,Q,U,V))=L(P*V+Q*U+1), where F lists
the unsigned correlation parities and L is another binary isomorphism.

Recovery uses the exact polar-form and hyperbolic-basis algorithm of
large_parity_map. Every exported coefficient is checked against actual
site pairs; a second check expands every coefficient probe at every shift.
The separate sign solver must still satisfy exact counts and correlations.
"""
from __future__ import annotations
import argparse
from collections import Counter
import json
from math import isqrt
import operator
from pathlib import Path
import random
import time
import numpy as np
from experiments.scalar_lift import hadamard_frame_search as base
from experiments.scalar_lift.large_parity_map import (
    BinaryField,BinaryCoordinates,CorrelationForms,bits,inverse_matrix,
    matrix_apply,DeterminantSampler as _DeterminantSampler,
)


def _integer(value, name):
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f'{name} must be an integer, not a Boolean')
    try:
        return operator.index(value)
    except TypeError as exc:
        raise ValueError(f'{name} must be an integer') from exc


def periodic_chart(p):
    """Return actual site-pair correlation rows, with no covariance imposed."""
    p = _integer(p, 'p')
    if p<3 or p%2==0 or any(p%d==0 for d in range(3,isqrt(p)+1,2)):
        raise ValueError('p must be an odd prime')
    shifts=list(range(1,(p+1)//2))
    rows=[]
    for shift in shifts:
        row=Counter()
        for x in range(p):
            y=(x+shift)%p
            row[min(x,y),max(x,y),0]+=1
        if sum(row.values())!=p:raise ArithmeticError('chart count failed')
        rows.append(dict(row))
    return dict(p=p,m=1,mu=1,n=p,index=np.arange(p,dtype=np.int64),
                phase=np.zeros(p,dtype=np.int64),orbits=[[(x,0)] for x in range(p)],
                shifts=shifts,rows=rows,chart_kind='unrestricted-periodic')


def row_sum_types(p):
    """All nonnegative decreasing integer rows of squared norm p."""
    p = periodic_chart(p)['p']
    result = []
    for a in range(isqrt(p)+1):
        for b in range(a+1):
            for c in range(b+1):
                remaining = p-a*a-b*b-c*c
                if remaining < 0:
                    continue
                d = isqrt(remaining)
                if d <= c and d*d == remaining:
                    result.append((a,b,c,d))
    return result


def recover_map(p, seed=None):
    """Construct one determinant factor, or fail rather than omit factors."""
    start = time.monotonic()
    data = periodic_chart(p)
    p = data['p']
    seed = p if seed is None else seed
    forms = CorrelationForms(data)
    n, degree = forms.variables, len(data['rows'])
    if n != 4*degree or degree < 3 or degree % 2 != 1:
        raise ValueError('this construction currently requires a single odd-degree factor of degree at least 3')
    inv0 = inverse_matrix(forms.polar_rows[0])
    transform = None
    field = None
    for which in range(1, degree):
        # B0*T=Bwhich; all polar matrices are symmetric.
        columns = [matrix_apply(inv0, row) for row in forms.polar_rows[which]]
        value, krylov = 1, []
        for _ in range(degree):
            krylov.append(value)
            value = base.binary_linear(columns, value)
        try:
            coordinate = BinaryCoordinates(krylov)
            polynomial = (1 << degree) ^ coordinate.coordinates(value)
        except ValueError:
            continue
        candidate = BinaryField(degree, polynomial)
        if not candidate.irreducible():
            continue
        # Check the operator relation on every binary coordinate vector.
        for j in range(n):
            total, power = 0, 1 << j
            for k in range(degree+1):
                if polynomial >> k & 1:
                    total ^= power
                power = base.binary_linear(columns, power)
            if total:
                break
        else:
            transform, field = columns, candidate
            break
    if transform is None:
        raise ValueError('could not recover a full-degree scalar field')

    def times_t(v): return base.binary_linear(transform, v)
    def scale(a, v):
        result = 0
        while a:
            if a & 1:
                result ^= v
            a >>= 1
            if a:
                v = times_t(v)
        return result

    u = 1
    v = 1 << next(bits(forms.polar_rows[0][0]))
    output_columns, power = [], u
    for _ in range(degree):
        output_columns.append(forms.polar(power, v))
        power = times_t(power)
    output_coordinates = BinaryCoordinates(output_columns)
    constant = output_coordinates.coordinates(forms.constant)
    inv_constant = field.inverse(constant)

    def q(v): return field.mul(output_coordinates.coordinates(forms.value(v)), inv_constant)
    def pairing(u, v):
        return field.mul(output_coordinates.coordinates(forms.polar(u, v)), inv_constant)

    def normalize_pair(a, b):
        cross = pairing(a, b)
        if not cross:
            raise ValueError('degenerate candidate plane')
        b = scale(field.inverse(cross), b)
        qa, qb = q(a), q(b)
        if qb:
            y = field.artin_schreier(field.mul(qa, qb))
            if y is None:
                return None
            t = field.mul(y, field.inverse(qb))
        else:
            t = qa
        e = a ^ scale(t, b)
        if not e or q(e):
            raise ValueError('isotropic-vector construction failed')
        f = b
        f = scale(field.inverse(pairing(e, f)), f)
        f ^= scale(q(f), e)
        if q(f) or pairing(e, f) != 1:
            raise ValueError('hyperbolic pair normalization failed')
        return e, f

    rng = random.Random(seed)
    first = None
    for _ in range(100):
        a, b = rng.getrandbits(n), rng.getrandbits(n)
        if pairing(a, b):
            first = normalize_pair(a, b)
            if first:
                break
    if first is None:
        raise ValueError('failed to find first hyperbolic plane')
    e1, f1 = first

    def complement(w): return w ^ scale(pairing(w, f1), e1) ^ scale(pairing(w, e1), f1)
    second = None
    for _ in range(100):
        a, b = complement(rng.getrandbits(n)), complement(rng.getrandbits(n))
        if pairing(a, b):
            second = normalize_pair(a, b)
            if second:
                break
    if second is None:
        raise ValueError('failed to find complementary hyperbolic plane')
    e2, f2 = second
    columns = []
    for vector in (e1, e2, f2, f1):
        power = vector
        for _ in range(degree):
            columns.append(power)
            power = times_t(power)
    if base.binary_rank(columns) != n:
        raise ValueError('constructed hyperbolic coordinate map is singular')
    config = dict(schema='unrestricted-periodic-parity-map-v1',
                  p=p, m=1, mu=1, field_degree=degree, field_modulus=field.modulus,
                  inverse_columns=columns,
                  construction='ratios of polar forms and exact hyperbolic splitting',
                  construction_seed=seed,
                  construction_seconds=round(time.monotonic()-start, 6))
    config['chart_kind']='unrestricted-periodic'
    config['field_coordinate_order'] = ['P','Q','U','V']
    return config


def _configuration(config):
    if config.get('chart_kind') != 'unrestricted-periodic' or config.get('m',1) != 1 or config.get('mu',1) != 1:
        raise ValueError('an unrestricted periodic parity map is required')
    if config.get('field_coordinate_order', ['P','Q','U','V']) != ['P','Q','U','V']:
        raise ValueError('field coordinates must be ordered P,Q,U,V')
    data = periodic_chart(config['p'])
    degree = _integer(config['field_degree'], 'field degree')
    if degree < 1 or data['n'] != 2*degree+1:
        raise ValueError('incorrect field degree for a single determinant factor')
    field = BinaryField(degree, _integer(config['field_modulus'], 'field modulus'))
    if not field.irreducible():
        raise ValueError('field modulus is reducible')
    variables = 4*degree
    columns = config['inverse_columns']
    if len(columns) != variables or any(type(c) is not int or not 0 <= c < (1 << variables) for c in columns) or base.binary_rank(columns) != variables:
        raise ValueError('singular or incorrectly sized map')
    return data, degree, field, columns


def validate_map(config):
    """Prove the full quadratic identity by every Boolean coefficient."""
    start = time.monotonic()
    data, degree, field, columns = _configuration(config)
    variables = 4*degree
    if data['n'] != 2*degree+1 or len(data['rows']) != degree:
        raise ValueError('chart is not a single determinant factor')
    pairs = [(u, v) for u in range(variables) for v in range(u)]
    inputs = [0]+[1 << u for u in range(variables)]+[(1 << u)|(1 << v) for u,v in pairs]
    colors = [base.binary_linear(columns, value) for value in inputs]
    axes = np.array([[0]+[(value >> (2*j)) & 3 for j in range(2*degree)] for value in colors], dtype=np.uint8)
    values = np.zeros((len(inputs), degree), dtype=np.uint8)
    for ell, row in enumerate(data['rows']):
        for (i,j,delta), count in row.items():
            if count & 1:
                rotated = np.where(axes[:,j] == 0, 0, 1+(axes[:,j]-1+delta)%3)
                values[:,ell] ^= axes[:,i] == rotated
    parity = [sum(int(bit) << ell for ell,bit in enumerate(row)) for row in values]
    if any(value != parity[0] for value in parity[1:1+variables]):
        raise ValueError('nonzero unexpected linear coefficient')
    pair_index = {pair:i for i,pair in enumerate(pairs)}
    def polar(u,v):
        if u < v:u,v=v,u
        return parity[1+variables+pair_index[(u,v)]] ^ parity[1+u] ^ parity[1+v] ^ parity[0]
    determinant_columns = [polar(3*degree+j,0) for j in range(degree)]
    if base.binary_rank(determinant_columns) != degree or determinant_columns[0] != parity[0]:
        raise ValueError('determinant equation map has wrong rank/constant')
    for i,(u,v) in enumerate(pairs):
        term = 0
        if (3*degree <= u < 4*degree and 0 <= v < degree) or (2*degree <= u < 3*degree and degree <= v < 2*degree):
            term = field.mul(1 << (u % degree), 1 << (v % degree))
        observed = parity[1+variables+i] ^ parity[1+u] ^ parity[1+v] ^ parity[0]
        if observed != base.binary_linear(determinant_columns, term):
            raise ValueError(f'quadratic coefficient failed at {u},{v}')
    return dict(all_quadratic_coefficients_verified=True,
                coefficients_checked=len(inputs), color_map_rank=variables,
                equation_map_rank=degree, field_modulus_irreducible=True,
                exact_identity='F(T(P,Q,U,V)) = L(P*V+Q*U+1), T and L invertible',
                validation_seconds=round(time.monotonic()-start, 6))


def validate_expanded_geometry(config):
    """Compare all coefficient probes with every directly expanded shift.

    This independent geometry check uses actual length-p color words and
    equality counts rather than the compressed orbit coefficient formulas.
    """
    start = time.monotonic()
    data, _degree, _field, _columns = _configuration(config)
    n = 4*config['field_degree']
    inputs = [0]+[1 << i for i in range(n)]+[(1 << i)|(1 << j)
                                               for i in range(n) for j in range(i)]
    colors = [base.binary_linear(config['inverse_columns'], v) for v in inputs]
    axes = np.array([[0]+[(v >> (2*j)) & 3 for j in range(data['n']-1)]
                     for v in colors], dtype=np.uint8)
    expanded = axes[:, data['index']]
    expanded = np.where(expanded == 0, 0,
                        1+(expanded-1+data['phase'][None,:]) % 3).astype(np.uint8)
    expected = np.zeros((len(inputs), len(data['rows'])), dtype=np.uint8)
    for ell,row in enumerate(data['rows']):
        for (i,j,delta),count in row.items():
            if count & 1:
                rotated = np.where(axes[:,j] == 0, 0, 1+(axes[:,j]-1+delta) % 3)
                expected[:,ell] ^= axes[:,i] == rotated
    shift_class = {}
    for ell, representative in enumerate(data['shifts']):
        for k in range(data['m']):
            shift = representative*pow(data['mu'], k, data['p']) % data['p']
            shift_class[shift] = shift_class[-shift % data['p']] = ell
    if len(shift_class) != data['p']-1:
        raise ValueError('shift representatives do not cover nonzero shifts')
    for shift in range(1, data['p']):
        actual = np.count_nonzero(expanded == np.roll(expanded, -shift, axis=1), axis=1) % 2
        if not np.array_equal(actual, expected[:,shift_class[shift]]):
            raise ValueError(f'expanded correlation geometry failed at shift {shift}')
    return dict(all_expanded_shift_parities_match_orbit_coefficients=True,
                coefficient_inputs=len(inputs), nonzero_shifts=data['p']-1,
                exact_correlations_compared=len(inputs)*(data['p']-1),
                validation_seconds=round(time.monotonic()-start, 6))


class DeterminantSampler(_DeterminantSampler):
    """Sample the exact parity locus, optionally filtering four support sizes.

    Field operations use arbitrary-size Python integers. No first-nonzero
    site sign or coordinate is fixed: the sole imposed label is at residue0.
    """
    def __init__(self, config, seed=None, rng=None):
        _configuration(config)
        super().__init__(config, seed=seed, rng=rng)

    def next_axes(self, *, target_row=None, max_draws=None):
        if max_draws is not None:
            max_draws = _integer(max_draws, 'draw limit')
            if max_draws < 0:
                raise ValueError('draw limit must be nonnegative')
        target = None
        if target_row is not None:
            target = tuple(_integer(x, 'row sum') for x in target_row)
            if len(target) != 4 or sum(r*r for r in target) != self.config['p']:
                raise ValueError('four row sums of squared norm p are required')
        stop = None if max_draws is None else self.draws+max_draws
        while stop is None or self.draws < stop:
            color = self.draw_packed()
            axes = [0]+[(color>>(2*j))&3 for j in range(self.nonorigin_orbits)]
            if target is not None:
                support = [axes.count(j) for j in range(4)]
                if any(w < abs(r) or (w-r)%2 for w,r in zip(support,target)):
                    continue
            return axes
        return None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--p', type=int, required=True)
    parser.add_argument('--seed', type=int)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    config = recover_map(args.p, seed=args.seed)
    config['validation'] = validate_map(config)
    config['expanded_geometry_validation'] = validate_expanded_geometry(config)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(config, indent=2)+'\n')
    print(json.dumps({key:value for key,value in config.items() if key != 'inverse_columns'}))


if __name__ == '__main__':
    main()

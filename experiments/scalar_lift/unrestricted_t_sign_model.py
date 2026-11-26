"""Exact four-count modular sign systems for unrestricted periodic T-designs.

At site i, the entry is (-1)^x_i e_{a_i}. Coordinate labels a_i are fixed;
the binary variables x_i remain unknown. All site pairs are retained. The
mod4 equations give x=particular+sum_j basis[j]*y_j over F2. The imported
quadratic_mod8 then computes each integer correlation divided by4 modulo2.
Four exact negative counts and the full integer correlations remain required.
"""
from __future__ import annotations

from collections.abc import Mapping
import operator

import numpy as np

from experiments.scalar_lift.quadratic_sign_lift import (
    InfeasibleAxes, _eliminate, correlations, evaluate_quadratic,
    quadratic_mod8, sign_mask,
)
from experiments.scalar_lift.unrestricted_parity_map import periodic_chart


def _integer(value, name):
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f'{name} must be an integer, not a Boolean')
    try:
        return operator.index(value)
    except TypeError as exc:
        raise ValueError(f'{name} must be an integer') from exc


def _chart(data_or_p):
    if not isinstance(data_or_p, Mapping):
        return periodic_chart(data_or_p)
    data = data_or_p
    p = _integer(data['p'], 'p')
    if (data.get('chart_kind') != 'unrestricted-periodic' or
            data.get('m') != 1 or data.get('mu') != 1 or data.get('n') != p or
            not np.array_equal(data.get('index'), np.arange(p)) or
            not np.array_equal(data.get('phase'), np.zeros(p)) or
            list(data.get('shifts', ())) != list(range(1,(p+1)//2))):
        raise ValueError('actual-site unrestricted periodic data are required')
    return data


def prepare(axes, data_or_p, target_row, *, fixed_signs=None,
            force_extremes=True):
    """Build mod4 signs with four count parities and forced exact-count signs.

    ``weights[ell,i,j]`` is the symmetric off-diagonal pair multiplicity;
    ``diagonal[ell]`` is sign independent. ``groups`` contains four lists of
    actual sites, including the origin in group0. ``negative_counts[j]`` is
    the exact requested number of negative signs on that coordinate. Only
    its parity is used by binary elimination; SAT backends impose it exactly.

    Origin sign x_0=0 is imposed. No sign at another site is gauged unless
    supplied explicitly. Invalid inputs raise ValueError; inconsistent
    parity/count arithmetic raises InfeasibleAxes.
    """
    data = _chart(data_or_p)
    n, length = data['n'], len(data['rows'])
    axes = tuple(_integer(a, 'coordinate label') for a in axes)
    if len(axes) != n or axes[0] != 0 or any(a < 0 or a > 3 for a in axes):
        raise ValueError('one coordinate label in0..3 per site, with zero first, is required')
    target = tuple(_integer(r, 'row sum') for r in target_row)
    if len(target) != 4 or sum(r*r for r in target) != data['p']:
        raise ValueError('four integer row sums of squared norm p are required')

    weights = np.zeros((length,n,n), dtype=np.int32)
    diagonal = np.zeros(length, dtype=np.int32)
    for ell, row in enumerate(data['rows']):
        for (i,j,delta), count in row.items():
            if delta != 0:
                raise ValueError('an unrestricted site pair cannot rotate coordinates')
            if axes[i] != axes[j]:
                continue
            if i == j:
                diagonal[ell] += count
            else:
                weights[ell,i,j] += count
                weights[ell,j,i] += count
    absolute = diagonal+weights.sum(axis=(1,2))//2
    if np.any(absolute%2):
        raise InfeasibleAxes('unsigned correlation parity failed')

    fixed = {0:0}
    for index, bit in (fixed_signs or {}).items():
        index = _integer(index, 'fixed site')
        bit = _integer(bit, 'fixed sign bit')
        if not 0 <= index < n or bit not in (0,1):
            raise ValueError('fixed signs need valid sites and binary sign bits')
        if index in fixed and fixed[index] != bit:
            raise InfeasibleAxes('a fixed sign contradicts the positive origin')
        fixed[index] = bit

    groups = [[i for i,a in enumerate(axes) if a == j] for j in range(4)]
    negative_counts = []
    for group, row_sum in zip(groups,target, strict=True):
        if len(group) < abs(row_sum) or (len(group)-row_sum)%2:
            raise InfeasibleAxes('support sizes do not permit the requested row')
        count = (len(group)-row_sum)//2
        negative_counts.append(count)
        known = sum(fixed[i] for i in group if i in fixed)
        unknown = [i for i in group if i not in fixed]
        needed = count-known
        if not 0 <= needed <= len(unknown):
            raise InfeasibleAxes('fixed signs contradict an exact coordinate count')
        if force_extremes and needed in (0,len(unknown)):
            bit = int(needed != 0)
            fixed.update((i,bit) for i in unknown)

    rows = []
    for ell in range(length):
        degrees = weights[ell].sum(axis=1)
        rows.append(sum(int(d%2)<<i for i,d in enumerate(degrees)) |
                    (int(absolute[ell]//2%2)<<n))
    rows.extend((1<<i)|(bit<<n) for i,bit in fixed.items())
    rows.extend(sum(1<<i for i in group)|((count%2)<<n)
                for group,count in zip(groups,negative_counts, strict=True))
    particular, basis = _eliminate(rows,n)
    forms = [sum(((vector>>i)&1)<<j for j,vector in enumerate(basis))
             for i in range(n)]
    return {
        'axes':axes, 'data':data, 'weights':weights, 'diagonal':diagonal,
        'absolute_correlations':absolute, 'particular':particular, 'basis':basis,
        'forms':forms, 'groups':groups, 'negative_counts':negative_counts,
        'target_row':target, 'forced_signs':fixed,
        'scope':'Fixed coordinate labels, four requested sums and positive origin.',
    }


prepare_rows = prepare


def exact_count_mask(system, mask):
    """Test all four exact coordinate counts on a packed negative-sign mask."""
    mask = _integer(mask, 'negative-sign mask')
    if not 0 <= mask < (1<<len(system['axes'])):
        raise ValueError('negative-sign mask is outside the site range')
    return all(sum((mask>>i)&1 for i in group) == count
               for group,count in zip(system['groups'],system['negative_counts'], strict=True))


def states_from_assignment(system, assignment):
    """Recover actual-site states2*a_i+x_i from free affine variables."""
    assignment = _integer(assignment, 'affine assignment')
    mask = sign_mask(system,assignment)
    return [2*a+((mask>>i)&1) for i,a in enumerate(system['axes'])]


def verify_states(states, data_or_p, target_row=None):
    """Independently reconstruct every integer entry and check every shift."""
    data = _chart(data_or_p)
    p = data['p']
    states = tuple(_integer(s, 'signed coordinate state') for s in states)
    if len(states) != p or states[0] != 0 or any(s < 0 or s > 7 for s in states):
        raise ValueError('one state in0..7 per site, with positive origin, is required')
    rows = np.zeros((4,p), dtype=np.int64)
    for site,state in enumerate(states):
        rows[state//2,site] = 1-2*(state%2)
    values = [int(np.sum(rows*np.roll(rows,-shift,axis=1), dtype=np.int64))
              for shift in range(p)]
    if values != [p]+[0]*(p-1):
        raise ValueError('candidate failed an expanded integer periodic correlation')
    sums = list(map(int,rows.sum(axis=1)))
    if target_row is not None:
        target = [_integer(r, 'row sum') for r in target_row]
        if len(target) != 4 or sums != target:
            raise ValueError('candidate failed the four requested coordinate sums')
    return {
        'every_expanded_periodic_correlation_verified':True,
        'shift_count':p, 'sequence_row_sums':sums,
        'sequence_support_sizes':list(map(int,np.count_nonzero(rows,axis=1))),
        'all_integer_correlations':values,
        'target_row_verified':target_row is not None,
    }

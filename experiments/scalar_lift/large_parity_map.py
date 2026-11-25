"""Exact determinant coordinates for larger twisted periodic T-design searches.

The color of an orbit is an element of F4, encoded by its two polynomial
basis bits (0,1,w,w+1). The permutation of the three nonzero colors is
multiplication by w. The unsigned-correlation parities are quadratic forms.

Rather than choosing an auxiliary Fourier field, this module recovers the
common scalar field from ratios of their polar forms, then constructs a
hyperbolic basis. Every exported map is independently checked against the
actual expanded correlation coefficients, on all constant, linear and
quadratic inputs. No table of all elements of the scalar field is created.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import random
import time

import numpy as np

from experiments.scalar_lift import hadamard_frame_search as base


def poly_mod(a: int, modulus: int) -> int:
    degree = modulus.bit_length() - 1
    while a.bit_length() > degree:
        a ^= modulus << (a.bit_length() - degree - 1)
    return a


def poly_gcd(a: int, b: int) -> int:
    while b:
        a, b = b, poly_mod(a, b)
    return a


@dataclass(frozen=True)
class BinaryField:
    degree: int
    modulus: int

    def mul(self, a: int, b: int) -> int:
        result = 0
        while b:
            if b & 1:
                result ^= a
            a <<= 1
            if a >> self.degree:
                a ^= self.modulus
            b >>= 1
        return result

    def square(self, a: int) -> int:
        return self.mul(a, a)

    def inverse(self, a: int) -> int:
        if not a:
            raise ZeroDivisionError('zero has no field inverse')
        u, v, g1, g2 = a, self.modulus, 1, 0
        while u != 1:
            j = u.bit_length() - v.bit_length()
            if j < 0:
                u, v, g1, g2, j = v, u, g2, g1, -j
            u ^= v << j
            g1 ^= g2 << j
        return poly_mod(g1, self.modulus)

    def trace(self, a: int) -> int:
        value, power = 0, a
        for _ in range(self.degree):
            value ^= power
            power = self.square(power)
        if value not in (0, 1):
            raise ValueError('trace did not land in F2')
        return value

    def artin_schreier(self, c: int) -> int | None:
        """Solve y^2+y=c by the half trace, for odd field degree."""
        if self.degree % 2 != 1:
            raise ValueError('half trace requires odd extension degree')
        if self.trace(c):
            return None
        value, power = 0, c
        for _ in range((self.degree + 1) // 2):
            value ^= power
            power = self.square(self.square(power))
        if self.square(value) ^ value != c:
            raise ValueError('half-trace identity failed')
        return value

    def irreducible(self) -> bool:
        if self.modulus.bit_length() != self.degree + 1 or not self.modulus & 1:
            return False
        if self.degree == 1:
            return self.modulus == 3
        factors = []
        remaining = self.degree
        for q in range(2, self.degree + 1):
            if remaining % q == 0:
                factors.append(q)
                while remaining % q == 0:
                    remaining //= q
        checks = {self.degree // q for q in factors}
        power = 2
        for i in range(1, self.degree + 1):
            power = self.square(power)
            if i in checks and poly_gcd(power ^ 2, self.modulus) != 1:
                return False
        return power == 2


class BinaryCoordinates:
    """Coordinates in an explicitly supplied independent binary vector list."""
    def __init__(self, columns):
        self.columns = tuple(columns)
        self.pivots = {}
        for i, column in enumerate(columns):
            value, coordinates = column, 1 << i
            while value:
                bit = value.bit_length() - 1
                if bit not in self.pivots:
                    self.pivots[bit] = (value, coordinates)
                    break
                pivot, coeff = self.pivots[bit]
                value ^= pivot
                coordinates ^= coeff
            if not value:
                raise ValueError('dependent coordinate columns')

    def coordinates(self, value):
        answer = 0
        while value:
            bit = value.bit_length() - 1
            if bit not in self.pivots:
                raise ValueError('vector outside coordinate span')
            pivot, coeff = self.pivots[bit]
            value ^= pivot
            answer ^= coeff
        return answer


def inverse_matrix(rows):
    n = len(rows)
    work = [row | (1 << (n+i)) for i, row in enumerate(rows)]
    for col in range(n):
        at = next((i for i in range(col, n) if work[i] >> col & 1), None)
        if at is None:
            raise ValueError('singular binary matrix')
        work[col], work[at] = work[at], work[col]
        for i in range(n):
            if i != col and work[i] >> col & 1:
                work[i] ^= work[col]
    return [row >> n for row in work]


def matrix_apply(rows, value):
    return sum(((row & value).bit_count() & 1) << i for i, row in enumerate(rows))


class CorrelationForms:
    def __init__(self, data):
        self.data = data
        self.variables = 2*(data['n']-1)
        n = self.variables
        self.constant = 0
        self.linear = [0]*n
        self.cross = [[0]*n for _ in range(n)]
        for ell, row in enumerate(data['rows']):
            bit = 1 << ell
            for (i, j, delta), count in row.items():
                if not count & 1:
                    continue
                self.constant ^= bit
                a = 0 if i == 0 else 1 << (2*(i-1))
                b = 0 if i == 0 else 1 << (2*(i-1)+1)
                c = 0 if j == 0 else 1 << (2*(j-1))
                d = 0 if j == 0 else 1 << (2*(j-1)+1)
                if delta == 1:
                    c, d = d, c ^ d
                elif delta == 2:
                    c, d = c ^ d, c
                x, y = a ^ c, b ^ d
                # 1 + norm(x+y*w), norm=x+y+x*y over F2.
                for variable in bits(x ^ y):
                    self.linear[variable] ^= bit
                for u in bits(x):
                    for v in bits(y):
                        if u == v:
                            self.linear[u] ^= bit
                        else:
                            self.cross[max(u, v)][min(u, v)] ^= bit
        self.polar_rows = []
        for ell in range(len(data['rows'])):
            rows = [0]*n
            for i in range(n):
                for j in range(i):
                    if self.cross[i][j] >> ell & 1:
                        rows[i] |= 1 << j
                        rows[j] |= 1 << i
            self.polar_rows.append(rows)

    def value(self, v):
        result = 0
        for i in bits(v):
            result ^= self.linear[i]
            lower = v & ((1 << i)-1)
            for j in bits(lower):
                result ^= self.cross[i][j]
        return result

    def polar(self, u, v):
        return sum((((matrix_apply(rows, u) & v).bit_count() & 1) << ell)
                   for ell, rows in enumerate(self.polar_rows))


def bits(value):
    while value:
        lowest = value & -value
        yield lowest.bit_length()-1
        value ^= lowest


def recover_map(p=223, m=3, mu=183, head_net=-5, tail_net=3, seed=223):
    start = time.monotonic()
    data = base.chart(p, m, mu)
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
    config = dict(p=p, m=m, mu=mu, field_degree=degree, field_modulus=field.modulus,
                  inverse_columns=columns, head_net=head_net, tail_net=tail_net,
                  fixed_first_head_negative=head_net < 0,
                  construction='ratios of polar forms and exact hyperbolic splitting',
                  construction_seed=seed,
                  construction_seconds=round(time.monotonic()-start, 6))
    return config


def validate_map(config):
    """Check every coefficient against actual orbit correlation parities."""
    start = time.monotonic()
    data = base.chart(config['p'], config['m'], config['mu'])
    degree = config['field_degree']
    field = BinaryField(degree, config['field_modulus'])
    if not field.irreducible():
        raise ValueError('field modulus is reducible')
    variables = 4*degree
    columns = config['inverse_columns']
    if len(columns) != variables or any(type(c) is not int or not 0 <= c < (1 << variables) for c in columns) or base.binary_rank(columns) != variables:
        raise ValueError('singular or incorrectly sized map')
    if data['n'] != 2*degree+1 or len(data['rows']) != degree:
        raise ValueError('chart is not a single determinant factor')
    bh, bt = 1+data['m']*config['head_net'], data['m']//3*config['tail_net']
    if bh*bh+3*bt*bt != data['p']:
        raise ValueError('DC norm condition failed')
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
                exact_identity='unsigned correlation parity = L(p*b1+q*b0+1), L invertible',
                validation_seconds=round(time.monotonic()-start, 6))


def validate_expanded_geometry(config):
    """Compare all coefficient probes with every directly expanded shift.

    This independent geometry check uses actual length-p color words and
    equality counts rather than the compressed orbit coefficient formulas.
    """
    start = time.monotonic()
    data = base.chart(config['p'], config['m'], config['mu'])
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


class DeterminantSampler:
    """Rejection-free determinant-one draws, with optional color filtering.

    Before color conditioning, each determinant-one field matrix has the
    same sampling probability. Multiplicative normalization is an optional
    search symmetry reduction, not a claim of independent uniform draws in
    a quotient.
    """
    def __init__(self, config, seed=None, rng=None):
        self.config = config
        self.degree = config['field_degree']
        self.field = BinaryField(self.degree, config['field_modulus'])
        self.columns = config['inverse_columns']
        self.rng = rng or random.Random(seed)
        self.draws = 0
        self.nonorigin_orbits = 2*self.degree
        self.even_mask = ((1 << (2*self.nonorigin_orbits))-1)//3
        self._chart = None

    def draw_packed(self):
        degree, rng, field = self.degree, self.rng, self.field
        while True:
            p, q = rng.getrandbits(degree), rng.getrandbits(degree)
            if p or q:
                break
        w = rng.getrandbits(degree)
        if p:
            b0, b1 = w, field.mul(1 ^ field.mul(q,w), field.inverse(p))
        else:
            b0, b1 = field.inverse(q), w
        packed = p | (q << degree) | (b0 << (2*degree)) | (b1 << (3*degree))
        self.draws += 1
        return base.binary_linear(self.columns, packed)

    def head_count(self, color):
        occupied = (color | (color >> 1)) & self.even_mask
        return self.nonorigin_orbits - occupied.bit_count()

    def next_axes(self, heads=None, first_head=False, normalize_first_head=False,
                  max_draws=None):
        stop = None if max_draws is None else self.draws + max_draws
        while stop is None or self.draws < stop:
            color = self.draw_packed()
            if heads is not None and self.head_count(color) not in heads:
                continue
            if first_head and not normalize_first_head and color & 3:
                continue
            axes = [0]+[(color >> (2*j)) & 3 for j in range(self.nonorigin_orbits)]
            if normalize_first_head:
                if self._chart is None:
                    c = self.config
                    self._chart = base.chart(c['p'], c['m'], c['mu'])
                axes, _scale = normalize_head(axes, self._chart)
            return axes
        return None


def normalize_head(axes, data):
    """Translate a nonzero head representative to1 by x->r*x.

    This permutes positions, preserves all correlations and commutes with
    the declared multiplier covariance. The function also transports each
    representative's stored three-cycle phase; it is not merely an orbit
    list permutation.
    """
    if len(axes) != data['n'] or axes[0] != 0 or any(type(a) is not int or not 0 <= a < 4 for a in axes):
        raise ValueError('a complete four-color orbit assignment is required')
    selected = next((i for i in range(1,len(axes)) if axes[i] == 0), None)
    if selected is None:
        raise ValueError('there is no nonorigin head to normalize')
    scale = data['orbits'][selected][0][0]
    result = []
    for orbit in data['orbits']:
        point = (scale*orbit[0][0]) % data['p']
        axis = axes[int(data['index'][point])]
        phase = int(data['phase'][point])
        result.append(0 if axis == 0 else 1+(axis-1+phase)%3)
    if result[0] != 0 or result[1] != 0:
        raise ValueError('head normalization failed')
    if sum(x == 0 for x in result) != sum(x == 0 for x in axes):
        raise ValueError('head normalization changed support counts')
    return result, scale


def sample(config, rng=None, count=1, first_head=False, heads=None,
           normalize_first_head=False):
    sampler = DeterminantSampler(config,rng=rng)
    for _ in range(count):
        yield sampler.next_axes(heads,first_head,normalize_first_head)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--p',type=int,default=223)
    parser.add_argument('--m',type=int,default=3)
    parser.add_argument('--mu',type=int,default=183)
    parser.add_argument('--head-net',type=int,default=-5)
    parser.add_argument('--tail-net',type=int,default=3)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    config=recover_map(args.p,args.m,args.mu,args.head_net,args.tail_net)
    config['validation']=validate_map(config)
    config['expanded_geometry_validation']=validate_expanded_geometry(config)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(config,indent=2)+'\n')
    print(json.dumps({k:v for k,v in config.items() if k!='inverse_columns'}))


if __name__=='__main__':main()

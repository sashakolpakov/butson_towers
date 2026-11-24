"""Faithful new period-field scalar restrictions, with explicitly nonflat outputs.

The additional period generator is retained.  Source entries still lie in
Q(zeta_6), so their numerical action repeats; the field certificate does
not assert existence of a flat frame.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from dataclasses import dataclass
from functools import cache
from pathlib import Path

import numpy as np
from flint import fmpz_mat
from sympy import isprime, n_order, primitive_root, totient

from experiments.scalar_lift.cyclotomic_towers import (
    cyclotomic_orders_with_degree_dividing,
    representation_obstructions,
)
from experiments.scalar_lift.forward_source_registry import load_source as registry_source

ROOT = Path(__file__).resolve().parents[2]
FRONTIER = ROOT / 'certificates/scalar-lift/forward-novelty-frontier-2026-11-22.json'
PRESETS = {167: (294923, (3850, 233, 9)),
           179: (268501, (721, 186, 5)),
           191: (1431737, (27861, 788, 25)),
           223: (189997, (5369, 342, 17)),
           283: (748819, (5384, 165, 7)),
           347: (3191707, (143224, 3241, 87)),
           359: (3856379, (17582, 95, 9)),
           491: (2524723, (1472, 27, 1))}


@cache
def load_ambient_source(source_id):
    if source_id.startswith('catalogue-'):
        from experiments.scalar_lift.catalogue_source_acquisition import load_source
        return load_source(source_id)
    if source_id == 'saved-unreal9-q6':
        from experiments.scalar_lift.nine_source_acquisition import load_source
        return load_source(source_id)
    return registry_source(source_id)


def _hash(array):
    data = np.asarray(array, dtype='<i8')
    return hashlib.sha256(data.tobytes()).hexdigest()


def _exact(array):
    return fmpz_mat(np.asarray(array).tolist())


def first_split_conductors(r, count=2, maximum_h=50000):
    """A bounded conductor enumeration, not an assertion of exhaustive fields."""
    result = []
    for h in range(2, maximum_h + 1, 2):
        conductor = r*h + 1
        if isprime(conductor) and pow(2, h, conductor) == 1:
            result.append({'conductor': conductor, 'subgroup_order': h,
                           'order_of_two': int(n_order(2, conductor))})
            if len(result) == count:
                break
    return result


def target_preflight():
    encoded = FRONTIER.read_bytes()
    frontier = json.loads(encoded)
    targets = frontier['retained_source_listed_unresolved_orders']
    result = []
    for n in (3, 7, 13, 15, 35):
        # The fixed common120 sources generate the full entry field, also
        # after a global phase shift; do not silently replace them by BH(n,6).
        phase = 6 if n in (3, 7, 13) else 120
        result.append({'source_order': n, 'entry_field_phase': phase,
                       'entry_field_degree': int(totient(phase)),
                       'order_divisible_targets': [N for N in targets if N % n == 0],
                       'field_degree_compatible_targets': [N for N in targets
                           if N % n == 0 and (N//n) % int(totient(phase)) == 0]})
    return {'frontier_sha256': hashlib.sha256(encoded).hexdigest(),
            'attributed_snapshot_only': True, 'source_routes': result,
            'fixed_source_routes_above_20000': [(n, N) for n in (3, 7, 13, 15, 35)
                                               for N in targets if N > 20000 and N % n == 0]}


@dataclass(frozen=True)
class PeriodField:
    degree: int
    conductor: int
    subgroup_order: int
    primitive_root: int
    subgroup: np.ndarray
    coset_representatives: np.ndarray
    multiplication: np.ndarray
    trace_gram: np.ndarray
    basis_numerator: np.ndarray
    basis_denominator: int
    norm_witness: tuple[int, int, int]

    def proof(self):
        r, ell, h = self.degree, self.conductor, self.subgroup_order
        full_fields = [q for q in cyclotomic_orders_with_degree_dividing(4*r)
                       if int(totient(q)) == 4*r and q % 3 == 0]
        return {'real_period_degree': r, 'prime_conductor': ell, 'subgroup_order': h,
                'primitive_root': self.primitive_root,
                'two_splits_completely': pow(2, h, ell) == 1,
                'order_of_two': int(n_order(2, ell)),
                'norm_witness': list(self.norm_witness),
                'trace_gram_formula': f'{ell} I - {h} J',
                'multiplication_definition': 'M[k,j]=count(x in H: g^k-x in g^j H)-h*[j=0]',
                'multiplication_sha256': _hash(self.multiplication),
                'trace_gram_sha256': _hash(self.trace_gram),
                'basis_numerator_sha256': _hash(self.basis_numerator),
                'basis_denominator': self.basis_denominator,
                'trace_basis_exactly_orthonormal': True,
                'trace_basis_2_adically_integral_and_identity_modulo_two': True,
                'period_generator_self_adjoint_for_trace_metric': True,
                'period_generator_is_primitive': True,
                'primitive_reason': 'prime field degree; period variance r*ell-h*r-1 is nonzero',
                'compositum': 'Q(zeta_12) L', 'compositum_degree': 4*r,
                'faithful_auxiliary_action_retained': True,
                'full_cyclotomic_degree': 4*(ell-1),
                'full_cyclotomic_action_asserted': False,
                'full_cyclotomic_same_degree_sextic_fields': [
                    {'phase': q, 'standard_metric_obstructions': representation_obstructions(4*r, q)}
                    for q in full_fields]}


@cache
def period_field(r, conductor=None, norm_witness=None):
    r = int(r)
    if conductor is None:
        conductor, norm_witness = PRESETS[r]
    ell = int(conductor)
    a, b, t = map(int, norm_witness)
    if not (isprime(r) and r % 4 == 3 and isprime(ell) and (ell-1) % (2*r) == 0):
        raise ValueError('prime real-period parameters required')
    h = (ell-1)//r
    if pow(2, h, ell) != 1 or a*a+r*b*b != ell*t*t or t % 2 != 1:
        raise ValueError('split-at-two and odd-denominator norm witnesses required')
    g = int(primitive_root(ell))
    subgroup = np.array([pow(g, r*j, ell) for j in range(h)], dtype=np.int64)
    representatives = np.array([pow(g, k, ell) for k in range(r)], dtype=np.int64)
    lookup = np.full(ell, -1, dtype=np.int64)
    for j, rep in enumerate(representatives):
        indices = subgroup*rep % ell
        if np.any(lookup[indices] != -1):
            raise ValueError('period cosets overlap')
        lookup[indices] = j
    if np.any(lookup[1:] < 0) or lookup[ell-1] != 0 or lookup[2] != 0:
        raise ValueError('period subgroup is not real and split at two')
    multiplication = np.empty((r, r), dtype=np.int64)
    for k, rep in enumerate(representatives):
        labels = lookup[(rep-subgroup) % ell]
        multiplication[k] = np.bincount(labels[labels >= 0], minlength=r)
    multiplication[:, 0] -= h
    unit, ones = np.eye(r, dtype=np.int64), np.ones((r, r), dtype=np.int64)
    gram = ell*unit-h*ones
    character = np.zeros(r, dtype=np.int64)
    character[1:] = [1 if pow(j, (r-1)//2, r) == 1 else -1 for j in range(1, r)]
    jacobsthal = character[(np.arange(r)[None, :]-np.arange(r)[:, None]) % r]
    numerator = r*a*unit-r*b*jacobsthal+(t-a+h*r*t)*ones
    denominator = r*t*ell
    # FLINT arbitrary integers: no int64 overflow is possible in proof products.
    Q, G, S, identity, J, M = map(_exact, (jacobsthal, gram, numerator, unit, ones, multiplication))
    if Q.transpose() != -Q or Q*Q.transpose() != r*identity-J or Q*J != 0*J:
        raise ValueError('Jacobi trace-frame identities failed')
    if S.transpose()*G*S != denominator**2*identity:
        raise ValueError('the exact trace normalization failed')
    if denominator % 2 != 1 or not np.array_equal(numerator % 2, unit):
        raise ValueError('the trace basis is not identity modulo two')
    if M.transpose()*G != G*M or sum(int(M[i, i]) for i in range(r)) != -1:
        raise ValueError('the period multiplication trace identities failed')
    # Conjugation eta_j -> eta_(j+1) supplies another multiplication operator.
    shifted = np.roll(np.roll(multiplication, 1, axis=0), 1, axis=1)
    if M*_exact(shifted) != _exact(shifted)*M:
        raise ValueError('conjugate period multiplication does not commute')
    for array in (subgroup, representatives, multiplication, gram, numerator):
        array.setflags(write=False)
    return PeriodField(r, ell, h, g, subgroup, representatives, multiplication,
                       gram, numerator, denominator, (a, b, t))


def _quaternion(a, b, c, d):
    return np.array([[a, -b, -c, -d], [b, a, -d, c],
                     [c, d, a, -b], [d, -c, b, a]], dtype=np.int64)


def _four_squares(n):
    for a in range(math.isqrt(n), -1, -1):
        for b in range(math.isqrt(n-a*a), -1, -1):
            for c in range(math.isqrt(n-a*a-b*b), -1, -1):
                d = math.isqrt(n-a*a-b*b-c*c)
                if a*a+b*b+c*c+d*d == n:
                    return a, b, c, d
    raise AssertionError('four-square enumeration failed')


def source_coordinate_view(field, source_id='unreal13-q6'):
    """A compact literal integer matrix N=N0 tensor I_r; not a sign frame."""
    source = load_ambient_source(source_id)
    if source.phase not in (6, 12):
        raise ValueError('this coordinate view needs a verified sixth/twelfth root source')
    r = field.degree
    unit = np.eye(4, dtype=np.int64)
    D = _quaternion(0, 1, 1, 1)
    J = np.array([[0, -1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, -1, 0]], dtype=np.int64)
    V = J@(unit-D)  # twice a primitive twelfth root, with V²/4=(I+D)/2.
    if not (np.array_equal(D@D, -3*unit) and np.array_equal(J@J, -unit)
            and np.array_equal(J@D, D@J)
            and np.array_equal(V@V@V@V-4*V@V+16*unit, np.zeros((4, 4), dtype=np.int64))):
        raise ValueError('the faithful quartic field action failed')
    twelfth_powers = [2*unit]
    for _ in range(11):
        next_numerator = twelfth_powers[-1]@V
        if np.any(next_numerator % 2):
            raise ValueError('twelfth root action left the half-integral table')
        twelfth_powers.append(next_numerator//2)
    twice_powers = np.array(twelfth_powers)[::12//source.phase]
    norm_coefficients = _four_squares(r)
    frame = 2*_quaternion(*norm_coefficients)
    if not np.array_equal(frame@frame.T, 4*r*unit):
        raise ValueError('the rational normalization frame failed')
    blocks = twice_powers@frame
    if np.any(blocks % 2):
        raise ValueError('the displayed normalized source is not integral')
    blocks //= 2
    compact = blocks[source.exponents].transpose(0, 2, 1, 3).reshape(4*source.order, 4*source.order)
    target = source.order*4*r
    if _exact(compact)*_exact(compact).transpose() != target*_exact(np.eye(len(compact), dtype=np.int64)):
        raise ValueError('the complete actual source-derived Gram failed')
    metadata = {'source_id': source.id, 'source_hash': source.source_hash,
                'source_order': source.order, 'source_phase': source.phase,
                'scope': 'Supplementary norm-frame normalization, not a pure restriction stage',
                'source_complete_Gram_verified': source.proof['full_Gram_verified'],
                'rational_endpoint_order': target, 'entry_formula': 'N[i,j]=N0[i//r,j//r]*[i%r=j%r]',
                'compact_matrix_sha256': _hash(compact), 'compact_shape': list(compact.shape),
                'integer_entry_alphabet': sorted(set(map(int, compact.flat)) | {0}),
                'complete_exact_Gram_scale': target, 'complete_exact_Gram_verified': True,
                'normalization_frame_quaternion_coefficients': [2*x for x in norm_coefficients],
                'normalization_frame_Gram_scale': 4*r,
                'faithful_ambient_generators': ['J tensor I_r', 'D tensor I_r', 'I_4 tensor M_eta0'],
                'source_action_ignores_auxiliary_generator': True,
                'right_norm_frame_inserted': True,
                'same_as_repeated_four_dimensional_source_blocks': True,
                'flat_entries': False, 'new_real_Hadamard_constructed': False}
    return compact, metadata, {'source_phase_action_twice': twice_powers, 'sqrt_minus_three': D,
                               'imaginary_unit': J, 'norm_frame': frame}


def family(r=179, source_id='unreal13-q6', bundle=None):
    field = period_field(r)
    compact, source, scalar_arrays = source_coordinate_view(field, source_id)
    frontier = json.loads(FRONTIER.read_text())['retained_source_listed_unresolved_orders']
    report = {'schema': 'new-period-ambient-restriction-v1', 'field': field.proof(),
              'source_coordinate_view': source,
              'target_listed_unresolved_in_used_snapshot': source['rational_endpoint_order'] in frontier,
              'missing_step': 'A joint frame whose images of every used source entry are signs',
              'new_Hadamard_order_claimed': False}
    from experiments.scalar_lift.pure_period_restriction import tower_report
    report['pure_restriction_towers'] = tower_report(field, load_ambient_source(source_id))
    if bundle is not None:
        path = Path(bundle)
        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(path, source_exponents=load_ambient_source(source_id).exponents, subgroup=field.subgroup, representatives=field.coset_representatives,
                            period_multiplication=field.multiplication, trace_gram=field.trace_gram,
                            trace_basis_numerator=field.basis_numerator,
                            trace_basis_denominator=np.array(field.basis_denominator, dtype=np.int64),
                            compact_integer_matrix=compact, **scalar_arrays,
                            metadata_json=np.array(json.dumps(report, sort_keys=True)))
    return report



RECIPES = [(179, 'unreal13-q6'), (223, 'unreal13-q6'), (167, 'unreal13-q6'),
           (191, 'saved-unreal9-q6')] + [
               (r, source) for r in (167, 283, 347, 359)
               for source in ('catalogue-p7-q6', 'catalogue-p7-q12')]


def build_report(bundle_directory=None):
    reports = []
    paths = []
    for r, source in RECIPES:
        filename = f'period-{r}-{source}.npz'
        path = None if bundle_directory is None else Path(bundle_directory)/filename
        reports.append(family(r, source, path))
        paths.append(filename)
    states = []
    for report in reports:
        tree = report['pure_restriction_towers']
        states.append(tree['initial_source']['literal_action_metric_sha256'])
        states.extend(stage['literal_action_metric_sha256'] for path in tree['paths'] for stage in path)
    return {'preflight': target_preflight(), 'new_field_families': reports,
            'scalar_tree_execution': {'status': 'complete_for_declared_finite_recipes',
                'recipe_count': len(RECIPES), 'source_to_rational_paths': 2*len(RECIPES),
                'literal_node_visits': len(states), 'literal_distinct_nodes': len(set(states)),
                'literal_duplicate_visits': len(states)-len(set(states)),
                'coordinate_permutation_only_nodes_deduplicated': False,
                'pending_declared_recipes': 0, 'unbounded_ambient_field_search_completed': False,
                'actual_bundle_filenames': paths,
                'retained_target_orders': sorted({r['source_coordinate_view']['rational_endpoint_order']
                                                 for r in reports}),
                'pure_restriction_flat_endpoint_found': False,
                'product_basis_family_excluded_by_structural_zeros': True,
                'no_further_conductor_search_pending_in_this_family': True,
                'nonproduct_joint_frame_not_constructed': True}}


class SavedPeriodBundle:
    """Replay one complete saved definition and access any pure stage entry."""
    def __init__(self, path):
        from experiments.scalar_lift.pure_period_restriction import period_towers
        with np.load(path, allow_pickle=False) as data:
            saved = {key: np.array(data[key]) for key in data.files}
        self.report = json.loads(str(saved.pop('metadata_json')))
        field_record = self.report['field']
        self.field = period_field(field_record['real_period_degree'], field_record['prime_conductor'],
                                  tuple(field_record['norm_witness']))
        source_id = self.report['source_coordinate_view']['source_id']
        self.source = load_ambient_source(source_id)
        expected = family(self.field.degree, source_id)
        if self.report != expected:
            raise ValueError('saved period metadata does not replay')
        compact, _, scalar = source_coordinate_view(self.field, source_id)
        arrays = {'source_exponents': self.source.exponents,
                  'subgroup': self.field.subgroup,
                  'representatives': self.field.coset_representatives,
                  'period_multiplication': self.field.multiplication,
                  'trace_gram': self.field.trace_gram,
                  'trace_basis_numerator': self.field.basis_numerator,
                  'trace_basis_denominator': np.array(self.field.basis_denominator),
                  'compact_integer_matrix': compact, **scalar}
        if set(saved) != set(arrays) or any(not np.array_equal(saved[key], value) for key, value in arrays.items()):
            raise ValueError('saved period arrays do not replay')
        self.paths = period_towers(self.field, self.source)
        self.order = 4*self.field.degree*self.source.order
        self.compact_supplement = compact

    def entry(self, row, column, path=0, stage=-1):
        return self.paths[path][stage].entry(row, column)

    def supplementary_normalized_entry(self, row, column):
        if int(row) != row or int(column) != column or not 0 <= row < self.order or not 0 <= column < self.order:
            raise ValueError('matrix indices outside the normalized supplement')
        i, a = divmod(int(row), self.field.degree)
        j, b = divmod(int(column), self.field.degree)
        return int(self.compact_supplement[i, j]) if a == b else 0

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--bundle-directory', type=Path)
    parser.add_argument('--check', type=Path)
    args = parser.parse_args()
    report = build_report(args.bundle_directory)
    reports = report['new_field_families']
    if args.check and json.loads(args.check.read_text()) != report:
        raise SystemExit('ambient-field report differs')
    if args.output:
        args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'endpoint_orders': [r['source_coordinate_view']['rational_endpoint_order'] for r in reports],
                      'flat': False, 'new_Hadamard': False}))


if __name__ == '__main__':
    main()

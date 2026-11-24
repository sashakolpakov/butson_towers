"""Literal partial third/fifth-root restrictions, interleaved with dyadic descent.

For q=5c, 4|c and 5 not dividing c, the degree-four representation over
Q(zeta_c) has an explicit two-basis frame P_s=2 rho(zeta_q^s). It is flat
exactly on exponents e with e+s nonzero modulo5. Every child below is an
actual Butson matrix, saved with all its entries; no source search occurs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter, deque
from functools import cache
from pathlib import Path

import numpy as np
from sympy import I, Matrix, Poly, Symbol, cyclotomic_poly, eye, totient

from experiments.scalar_lift.fifth_dyadic_tower import fifth_numerators
from experiments.scalar_lift.fifth_source_assembly import root_vectors
from experiments.scalar_lift.forward_butson_descent import (
    available_steps,
    exponent_hash,
    scalar_blocks,
    verify_source,
)
from experiments.scalar_lift.forward_source_registry import VerifiedSource, load_source
from experiments.scalar_lift.large_butson_descendants import fast_scalar_certificate
from experiments.scalar_lift.scalar_restriction_obstructions import alphabet_obstruction
from hadamard.inverse_lift import is_sign_hadamard

ROOT = Path(__file__).resolve().parents[2]
SEEDS = ('common120-15', 'common120-35', 'fifth-5-q20', 'fifth-40-q40',
         'fifth-40-q80', 'fifth-80-q40', 'fifth-160-q20')


def _json_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


@cache
def gaussian_powers():
    """Return integer Gaussian numerators N_k=2 M^k, including N_0."""
    realified = [2*np.eye(8, dtype=np.int64), *fifth_numerators()]
    result = []
    for array in realified:
        re, im = array[::2, ::2], array[1::2, ::2]
        reconstructed = np.empty((8, 8), dtype=np.int64)
        reconstructed[::2, ::2], reconstructed[1::2, 1::2] = re, re
        reconstructed[1::2, ::2], reconstructed[::2, 1::2] = im, -im
        if not np.array_equal(array, reconstructed):
            raise ArithmeticError('the stored action is not Gaussian-linear')
        result.append(Matrix(re.tolist())+I*Matrix(im.tolist()))
    for k in range(5):
        for ell in range(5):
            if result[k]*result[ell] != 2*result[(k+ell) % 5]:
                raise ArithmeticError('the Gaussian fifth powers do not multiply')
            if result[k]*result[ell].H != 2*result[(k-ell) % 5]:
                raise ArithmeticError('the Gaussian fifth powers do not preserve adjoints')
    if sum(result, Matrix.zeros(4)) != Matrix.zeros(4):
        raise ArithmeticError('the Gaussian action fails Phi5')
    return tuple(result)


@cache
def scalar_certificate(upper_phase):
    q = int(upper_phase)
    c, remainder = divmod(q, 5)
    if remainder or c % 4 or c % 5 == 0:
        raise ValueError('need q=5c with4|c and5 not dividing c')
    if int(totient(q)) != 4*int(totient(c)):
        raise ArithmeticError('the relative field degree is not four')
    a, b = pow(5, -1, c), pow(c, -1, 5)
    numerators = gaussian_powers()
    roots = root_vectors(c)
    coefficients = np.zeros((q, 4, 4, roots.shape[1]), dtype=np.int64)
    for e in range(q):
        matrix = numerators[b*e % 5]
        for i in range(4):
            for j in range(4):
                re, im = matrix[i, j].as_real_imag()
                coefficients[e, i, j] = (int(re)*roots[a*e % c]
                    + int(im)*roots[(a*e+c//4) % c])
    x = Symbol('x')
    phi = Poly(cyclotomic_poly(q, x), x)
    total = np.zeros_like(coefficients[0])
    for (power,), value in phi.terms():
        total += int(value)*coefficients[power]
    if np.any(total):
        raise ArithmeticError('the represented primitive root fails Phi_q')
    expected = np.einsum('ij,k->ijk', 2*np.eye(4, dtype=np.int64), roots[1])
    if not np.array_equal(coefficients[5], expected):
        raise ArithmeticError('the fifth power does not recover zeta_c')
    # The cyclic vector identifies this action with ordinary field multiplication.
    v = eye(4)[:, 0]
    cyclic = Matrix.hstack(v, *(numerators[k]*v/2 for k in (1, 2, 3)))
    companion = Matrix([[0, 0, 0, -1], [1, 0, 0, -1], [0, 1, 0, -1], [0, 0, 1, -1]])
    if cyclic.det() == 0 or numerators[1]*cyclic/2 != cyclic*companion:
        raise ArithmeticError('the explicit restriction basis is not faithful')
    flat = []
    units = {1, -1, I, -I}
    for e in range(q):
        if all(value in units for value in numerators[b*e % 5]):
            flat.append(e)
    if flat != [e for e in range(q) if e % 5]:
        raise ArithmeticError('the partial flat alphabet is not exactly the nonzero fifth cosets')
    return {'upper_phase': q, 'lower_phase': c, 'relative_field_degree': 4,
            'CRT_exponents_a_b': [a, b], 'generator_action': 'rho(zeta_q)=zeta_c^a M^b',
            'Gaussian_fifth_numerators': [[[str(v) for v in row] for row in matrix.tolist()]
                                        for matrix in numerators],
            'output_basis_in_fifth_power_basis': [[str(v) for v in row] for row in cyclic.inv().tolist()],
            'output_basis_definition': 'alpha_i=sum_j zeta5^j (cyclic^-1)[j,i]',
            'input_basis_at_shift_s': '2 zeta_q^s alpha_i',
            'input_output_frame': 'P_s=2 rho(zeta_q^s), P_s P_s*=4I',
            'Phi_upper_of_represented_primitive_root_exact_zero': True,
            'generator_fifth_power_is_zeta_lower': True,
            'faithful_field_basis_exactly_identified': True,
            'Gaussian_multiplication_and_adjoint_pairs_checked': 25,
            'all_scalar_cross_Grams_proved': q*q,
            'cross_Gram_identity': 'B_e B_f*=4 rho(zeta_q^(e-f)), independent of frame shift s',
            'unshifted_flat_exponents': flat,
            'shifted_flat_domain': '(e+s) mod5 !=0',
            'forbidden_domain_witness': 'when(e+s) mod5=0, B_e has offdiagonal zeros and diagonal twice a root',
            'no_dense_pair_enumeration_needed': True,
            'operation': 'scalar restriction with explicitly declared unequal input/output field bases'}


@cache
def block_table(upper_phase, frame_shift=0):
    certificate = scalar_certificate(upper_phase)
    q, c = certificate['upper_phase'], certificate['lower_phase']
    a, b = certificate['CRT_exponents_a_b']
    lookup = {1: 0, I: 1, -1: 2, -I: 3}
    table = np.full((q, 4, 4), -1, dtype=np.int64)
    for e in range(q):
        t = e+frame_shift
        k = b*t % 5
        if not k:
            continue
        matrix = gaussian_powers()[k]
        for i in range(4):
            for j in range(4):
                table[e, i, j] = (a*t+(c//4)*lookup[matrix[i, j]]) % c
    return table


@cache
def third_certificate(upper_phase):
    q = int(upper_phase)
    c, remainder = divmod(q, 3)
    if remainder or c % 4 or c % 3 == 0:
        raise ValueError('need q=3c with4|c and3 not dividing c')
    a, b = pow(3, -1, c), pow(c, -1, 3)
    m = Matrix([[-1+I, -1-I], [1-I, -1-I]])/2
    powers = (eye(2), m, (m*m).applyfunc(lambda v: v.expand()))
    frame = Matrix([[1, 1], [1, -1]])
    if (m*m+m+eye(2)).applyfunc(lambda v: v.expand()) != Matrix.zeros(2) or (m*m.H).applyfunc(lambda v: v.expand()) != eye(2):
        raise ArithmeticError('the Gaussian third action failed')
    for k in range(3):
        for ell in range(3):
            if (powers[k]*frame*(powers[ell]*frame).H-2*powers[(k-ell) % 3]).applyfunc(lambda v: v.expand()) != Matrix.zeros(2):
                raise ArithmeticError('a third scalar cross Gram failed')
    cyclic = Matrix.hstack(eye(2)[:, 0], m[:, 0])
    companion = Matrix([[0, -1], [1, -1]])
    if cyclic.det() == 0 or (m*cyclic-cyclic*companion).applyfunc(lambda v: v.expand()) != Matrix.zeros(2):
        raise ArithmeticError('the third restriction basis failed')
    roots = root_vectors(c)
    coefficients = np.zeros((q, 2, 2, roots.shape[1]), dtype=np.int64)
    for e in range(q):
        matrix = 2*powers[b*e % 3]
        for i in range(2):
            for j in range(2):
                re, im = matrix[i, j].expand().as_real_imag()
                coefficients[e, i, j] = (int(re)*roots[a*e % c]
                    +int(im)*roots[(a*e+c//4) % c])
    x = Symbol('x')
    total = np.zeros_like(coefficients[0])
    for (power,), value in Poly(cyclotomic_poly(q, x), x).terms():
        total += int(value)*coefficients[power]
    if np.any(total) or int(totient(q)) != 2*int(totient(c)):
        raise ArithmeticError('the full third cyclotomic representation failed')
    return {'upper_phase': q, 'lower_phase': c, 'relative_field_degree': 2,
            'CRT_exponents_a_b': [a, b], 'generator_action': 'rho(zeta_q)=zeta_c^a M3^b',
            'Gaussian_third_action': [[str(v) for v in row] for row in m.tolist()],
            'Gaussian_frame': [[int(v) for v in row] for row in frame.tolist()],
            'output_basis_in_third_power_basis': [[str(v) for v in row] for row in cyclic.inv().tolist()],
            'input_basis_at_shift_s': 'zeta_q^s times (alpha0+alpha1, alpha0-alpha1)',
            'input_output_frame': 'P_s=rho(zeta_q^s) H2, P_s P_s*=2I',
            'Phi_upper_of_represented_primitive_root_exact_zero': True,
            'faithful_field_basis_exactly_identified': True,
            'Gaussian_cross_Gram_pairs_checked': 9,
            'all_scalar_cross_Grams_proved': q*q,
            'cross_Gram_identity': 'B_e B_f*=2rho(zeta_q^(e-f)), independent of shift s',
            'shifted_flat_domain': '(b*(e+s)) mod3 in{0,1}',
            'forbidden_domain_witness': 'M3^2 H2 has zeros and entries of norm2',
            'operation': 'scalar restriction with explicitly declared unequal field bases'}


@cache
def third_blocks(upper_phase, frame_shift=0):
    certificate = third_certificate(upper_phase)
    q, c = certificate['upper_phase'], certificate['lower_phase']
    a, b = certificate['CRT_exponents_a_b']
    m = Matrix([[-1+I, -1-I], [1-I, -1-I]])/2
    frame = Matrix([[1, 1], [1, -1]])
    powers = (frame, m*frame)
    lookup = {1: 0, I: 1, -1: 2, -I: 3}
    table = np.full((q, 2, 2), -1, dtype=np.int64)
    for e in range(q):
        t = e+frame_shift
        k = b*t % 3
        if k == 2:
            continue
        for i in range(2):
            for j in range(2):
                table[e, i, j] = (a*t+(c//4)*lookup[powers[k][i, j].expand()]) % c
    return table


def restrict_third(source, upper_phase, frame_shift=0):
    array = np.asarray(source, dtype=np.int64)
    table = third_blocks(upper_phase, frame_shift)
    if np.any(table[array] < 0):
        raise ValueError('the source uses a forbidden third coset for this input basis')
    n = len(array)
    return table[array].transpose(0, 2, 1, 3).reshape(2*n, 2*n)


@cache
def input_ids():
    from experiments.scalar_lift.broad_source_acquisition import dataset, root_gauge, sources
    q12 = [row for row in sources() if row['phase'] == 12]
    compatible = [row['id'] for row in q12 if root_gauge(
        dataset(row['order'], 12)[row['catalogue_class_index']], 6)['exists']]
    return (*SEEDS, *(f'cubic13-q24-t{k}' for k in range(3)), 'catalogue-p7-q12',
            *(row['id'] for row in q12), *('gauge3-'+identifier for identifier in compatible))


@cache
def load_input(identifier):
    if identifier.startswith('gauge3-'):
        from experiments.scalar_lift.broad_source_acquisition import dataset, root_gauge
        original = identifier.removeprefix('gauge3-')
        _prefix, n, q, index = original.split('-')
        raw = dataset(int(n), int(q))[int(index)]
        gauge = root_gauge(raw, 6)
        if not gauge['exists']:
            raise ValueError('the declared fixed third-coset gauge does not exist')
        matrix = (raw+np.asarray(gauge['row_shifts'])[:, None]
                  +np.asarray(gauge['column_shifts'])[None, :]) % int(q)
        proof = {**verify_source(matrix, int(q)), 'fixed_input_gauge': gauge,
                 'raw_source_hash': exponent_hash(raw), 'published_source_id': original,
                 'source_construction_search_performed': False}
        array = matrix.astype(np.uint16)
        array.setflags(write=False)
        return VerifiedSource(identifier, len(array), int(q), exponent_hash(array),
                              tuple(map(int, np.unique(array))), proof, array,
                              lambda i, j: int(array[i, j]))
    if identifier.startswith('aalto-'):
        from experiments.scalar_lift.broad_source_acquisition import load_source as broad_source
        return broad_source(identifier)
    if identifier.startswith('catalogue-'):
        from experiments.scalar_lift.catalogue_source_acquisition import (
            load_source as catalogue_source,
        )
        return catalogue_source(identifier)
    return load_source(identifier)


def restrict(source, upper_phase, frame_shift=0):
    array = np.asarray(source, dtype=np.int64)
    table = block_table(upper_phase, frame_shift)
    if np.any(table[array] < 0):
        raise ValueError('the source uses a forbidden fifth coset for this input basis')
    n = len(array)
    return table[array].transpose(0, 2, 1, 3).reshape(4*n, 4*n)


def _matrix_key(array, phase):
    return _json_hash({'phase': phase, 'order': len(array), 'literal_exponent_hash': exponent_hash(array)})


def build_report(bundle_directory, source_ids=None):
    directory = Path(bundle_directory)
    directory.mkdir(parents=True, exist_ok=True)
    nodes, edges, rejections, queue = {}, [], [], deque()
    scalar_proofs = {}
    selected_sources = input_ids() if source_ids is None else tuple(source_ids)

    def admit(array, phase, origin, parent=None, edge=None, registry_proof=None):
        key = _matrix_key(array, phase)
        if key in nodes:
            return key, True
        digest = exponent_hash(array)
        proof = {'order': len(array), 'phase': phase, 'exponent_sha256': digest,
                 'full_Gram_verified': True,
                 'verification': 'verified source and exact scalar cross-Gram induction',
                 'parent_node': parent, 'edge': edge}
        independent = None
        if registry_proof is not None:
            proof['registry_source_proof'] = registry_proof
        elif phase == 2:
            if not is_sign_hadamard((1-2*array).tolist()):
                raise ArithmeticError('a complete terminal sign Gram failed')
            independent = {'all_rows_checked_by_exact_bit_distances': True,
                           'full_sign_Gram_verified': True}
        elif len(array) <= 160:
            independent = verify_source(array, phase)
        record = {'id': key, 'order': len(array), 'phase': phase,
                  'exponent_sha256': digest, 'source_origin': origin, 'parent_node': parent,
                  'edge': edge, 'proof': proof, 'independent_full_Gram_check': independent,
                  'used_alphabet': list(map(int, np.unique(array))),
                  'bundle': str(directory/(key+'.npz')),
                  'new_Hadamard_order_claimed': False}
        record['alphabet_gate'] = alphabet_obstruction(phase, record['used_alphabet'])
        np.savez_compressed(record['bundle'], exponents=np.asarray(array, dtype=np.uint16),
                            metadata_json=np.array(json.dumps(record, sort_keys=True)))
        nodes[key] = record
        can_fifth = phase % 5 == 0 and (phase//5) % 4 == 0 and (phase//5) % 5
        can_third = phase % 3 == 0 and (phase//3) % 4 == 0 and (phase//3) % 3
        if available_steps(phase) or can_fifth or can_third:
            queue.append((key, np.asarray(array, dtype=np.int64)))
        return key, False

    for identifier in selected_sources:
        source = load_input(identifier)
        if source.proof['full_Gram_verified'] is not True or exponent_hash(source.exponents) != source.source_hash:
            raise ArithmeticError('a source proof is not bound to its literal entries')
        admit(source.exponents, source.phase, identifier, registry_proof=source.proof)
    while queue:
        parent_key, array = queue.popleft()
        parent = nodes[parent_key]
        q, origin = parent['phase'], parent['source_origin']
        if q % 3 == 0 and (q//3) % 4 == 0 and (q//3) % 3:
            proof_key = f'partial3-{q}'
            scalar_proofs[proof_key] = third_certificate(q)
            b = pow(q//3, -1, 3)
            for shift in range(3):
                bad = [e for e in parent['used_alphabet'] if b*(e+shift) % 3 == 2]
                recipe = {'kind': 'partial-third', 'upper_phase': q, 'lower_phase': q//3,
                          'frame_shift': shift, 'scalar_proof_key': proof_key}
                if bad:
                    rejections.append({'parent': parent_key, 'source_origin': origin, **recipe,
                                       'excluded_exponents': bad,
                                       'literal_source_occurrence': np.argwhere(array == bad[0])[0].tolist(),
                                       'scalar_block_zero_witness': [0, 0],
                                       'reason': 'actual restricted scalar block has zero entries'})
                    continue
                child_array = restrict_third(array, q, shift)
                child, duplicate = admit(child_array, q//3, origin, parent_key, recipe)
                edges.append({'parent': parent_key, 'child': child, **recipe,
                              'literal_duplicate_pruned': duplicate})
        if q % 5 == 0 and (q//5) % 4 == 0 and (q//5) % 5:
            proof_key = f'partial5-{q}'
            scalar_proofs[proof_key] = scalar_certificate(q)
            for shift in range(5):
                bad = [e for e in parent['used_alphabet'] if (e+shift) % 5 == 0]
                recipe = {'kind': 'partial-fifth', 'upper_phase': q, 'lower_phase': q//5,
                          'frame_shift': shift, 'scalar_proof_key': proof_key}
                if bad:
                    position = np.argwhere(array == bad[0])[0].tolist()
                    rejections.append({'parent': parent_key, 'source_origin': origin, **recipe,
                                       'excluded_exponents': bad,
                                       'literal_source_occurrence': position,
                                       'scalar_block_zero_witness': [0, 1],
                                       'reason': 'actual restricted scalar block has zero offdiagonal entries'})
                    continue
                child_array = restrict(array, q, shift)
                child, duplicate = admit(child_array, q//5, origin, parent_key, recipe)
                edges.append({'parent': parent_key, 'child': child, **recipe,
                              'literal_duplicate_pruned': duplicate})
        for p in available_steps(q):
            for mode in ('primary', 'primitive'):
                proof_key = f'repeated-{q}-{p}-{mode}'
                scalar_proofs[proof_key] = fast_scalar_certificate(q, p, mode)
                recipe = {'kind': 'repeated-prime', 'upper_phase': q, 'lower_phase': q//p,
                          'prime': p, 'generator_mode': mode, 'scalar_proof_key': proof_key}
                old = len(array)
                child_array = scalar_blocks(q, p, mode)[array].transpose(0, 2, 1, 3).reshape(p*old, p*old)
                child, duplicate = admit(child_array, q//p, origin, parent_key, recipe)
                edges.append({'parent': parent_key, 'child': child, **recipe,
                              'literal_duplicate_pruned': duplicate})
    real = [n for n in nodes.values() if n['phase'] == 2]
    return {'schema': 'literal-partial-third-fifth-descent-v1', 'source_ids': list(selected_sources),
            'nodes': list(nodes.values()), 'edges': edges, 'domain_rejections': rejections,
            'scalar_proofs': scalar_proofs, 'source_count': len(selected_sources),
            'literal_node_count': len(nodes), 'successful_scalar_edge_attempts': len(edges),
            'exact_literal_duplicates_pruned': sum(e['literal_duplicate_pruned'] for e in edges),
            'partial_fifth_edges': sum(e['kind'] == 'partial-fifth' for e in edges),
            'partial_third_edges': sum(e['kind'] == 'partial-third' for e in edges),
            'phase_order_inventory': [{'phase': q, 'order': n, 'count': count}
                for (q, n), count in sorted(Counter((x['phase'], x['order']) for x in nodes.values()).items())],
            'real_orders': sorted({n['order'] for n in real}),
            'real_literal_nodes': len(real),
            'all_terminal_sign_Grams_independently_checked': all(
                n['independent_full_Gram_check']['full_sign_Gram_verified'] for n in real),
            'literal_storage_bytes': sum(Path(n['bundle']).stat().st_size for n in nodes.values()),
            'pending_nodes': 0, 'finite_declared_tree_complete': True,
            'only_operations': 'scalar restriction in declared field bases and exact literal duplicate pruning',
            'new_source_construction_or_solver_used': False, 'new_Hadamard_order_found': False,
            'scope': 'the declared source_ids, repeated-prime modes, all five fifth-coset and all three third-coset input-basis positions'}


class SavedNode:
    """A complete literal Butson node, with recursive exact restriction replay."""
    def __init__(self, path):
        self.path = Path(path)
        with np.load(self.path, allow_pickle=False) as data:
            self.exponents = np.asarray(data['exponents'], dtype=np.int64)
            self.metadata = json.loads(str(data['metadata_json']))
        self.order, self.phase = self.metadata['order'], self.metadata['phase']
        if self.exponents.shape != (self.order, self.order) or _matrix_key(self.exponents, self.phase) != self.metadata['id']:
            raise ValueError('the saved literal matrix hash does not match')
        self.exponents.setflags(write=False)

    def replay(self):
        record = self.metadata
        if record['parent_node'] is None:
            source = load_input(record['source_origin'])
            if not np.array_equal(self.exponents, source.exponents) or source.proof != record['proof']['registry_source_proof']:
                raise ValueError('the root source does not replay')
        else:
            parent = SavedNode(self.path.parent/(record['parent_node']+'.npz'))
            parent.replay()
            edge = record['edge']
            if edge['upper_phase'] != parent.phase or edge['lower_phase'] != self.phase:
                raise ValueError('the saved edge does not have the declared fields')
            if edge['kind'] == 'partial-fifth':
                scalar_certificate(parent.phase)
                expected = restrict(parent.exponents, parent.phase, edge['frame_shift'])
            elif edge['kind'] == 'partial-third':
                third_certificate(parent.phase)
                expected = restrict_third(parent.exponents, parent.phase, edge['frame_shift'])
            else:
                p, mode, old = edge['prime'], edge['generator_mode'], parent.order
                fast_scalar_certificate(parent.phase, p, mode)
                expected = scalar_blocks(parent.phase, p, mode)[parent.exponents].transpose(0, 2, 1, 3).reshape(p*old, p*old)
            if not np.array_equal(self.exponents, expected):
                raise ValueError('the saved child fails actual scalar restriction replay')
        return {'literal_matrix_replayed': True, 'all_ancestor_scalar_edges_replayed': True}

    def entry(self, row, column):
        if type(row) is not int or type(column) is not int or not 0 <= row < self.order or not 0 <= column < self.order:
            raise ValueError('matrix entry is outside the saved node')
        return int(self.exponents[row, column])

    def as_verified_source(self):
        self.replay()
        array = self.exponents.astype(np.uint16)
        array.setflags(write=False)
        return VerifiedSource('partial-odd-scalar-'+self.metadata['id'], self.order, self.phase,
                              exponent_hash(array), tuple(self.metadata['used_alphabet']),
                              self.metadata['proof'], array, lambda i, j: int(array[i, j]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle-directory', type=Path, default=ROOT/'artifacts/partial-fifth-descent')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--output', type=Path)
    group.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = build_report(args.bundle_directory)
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    elif json.loads(args.check.read_text()) != result:
        raise ValueError('the partial fifth descent report differs from exact replay')
    print(json.dumps({k: result[k] for k in ('literal_node_count', 'successful_scalar_edge_attempts',
                      'exact_literal_duplicates_pruned', 'partial_fifth_edges', 'partial_third_edges', 'real_orders',
                      'real_literal_nodes', 'pending_nodes', 'new_Hadamard_order_found')}))


if __name__ == '__main__':
    main()

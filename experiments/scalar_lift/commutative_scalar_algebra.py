"""Exact scalar algebras of commuting symmetric Williamson-type sources.

The algebra Q[T] is a product of totally real fields. Its degree may be
smaller than the dimension of the actual module V=Q^w. All field factors,
representation multiplicities, source-polynomial reconstruction recipes
and the transported metric are retained. No cyclic group algebra is assumed.

python-flint is needed only when running these optional algebra audits;
ordinary endpoint construction and verification can remain NumPy-only.
"""
from __future__ import annotations
import argparse
import gc
import hashlib
import json
import math
from pathlib import Path
import time
import numpy as np
try:
    from flint import fmpz_mat, fmpq_mat, fmpz_poly, fmpq_poly
except ImportError:
    fmpz_mat = fmpq_mat = fmpz_poly = fmpq_poly = None
ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / 'experiments/scalar_lift/near_williamson_inputs.json'

def _require_flint():
    if fmpz_mat is None:
        raise RuntimeError('the optional scalar-algebra audit requires python-flint')

def _check(condition, message='exact scalar-algebra identity failed'):
    if not condition:
        raise ArithmeticError(message)

def I(n):
    return fmpz_mat([[int(i == j) for j in range(n)] for i in range(n)])

def Z(n):
    return fmpz_mat(n, n)

def intpoly(poly):
    return [str(int(c)) for c in poly.coeffs()]

def rational_coefficients(poly):
    """Return ascending integer numerators and one positive denominator."""
    coeff = list(poly.coeffs())
    den = math.lcm(*(int(c.denominator) for c in coeff)) if coeff else 1
    return ([int(c.numerator) * (den // int(c.denominator)) for c in coeff], den)

def polynomial_record(poly):
    """Hash canonical ascending rational coefficients without a large dump."""
    coeff = list(poly.coeffs())
    h = hashlib.sha256(f'{len(coeff)}\n'.encode())
    for c in coeff:
        h.update(f'{c.numerator}/{c.denominator}\n'.encode())
    return {'coefficient_count': len(coeff), 'coefficient_sha256': h.hexdigest(), 'hash_encoding': 'ASCII coefficient count newline, then each ascending reduced numerator/denominator newline'}

def matrix_sha(A):
    header = f'{A.nrows()} {A.ncols()}\n'
    h = hashlib.sha256(header.encode())
    for value in A.entries():
        h.update(f'{value.numerator}/{value.denominator}\n'.encode())
    return h.hexdigest()

def eval_poly(poly, powers, n):
    nums, den = rational_coefficients(fmpq_poly(poly))
    value = Z(n)
    for k, c in enumerate(nums):
        if c:
            value += c * powers[k]
    return fmpq_mat(value) / den

def column_matrix(columns, n):
    return fmpq_mat([[column[i, 0] for column in columns] for i in range(n)])

def column(A, j):
    return fmpq_mat([[A[i, j]] for i in range(A.nrows())])

def newton_traces(coeff, count):
    d = len(coeff) - 1
    values = [d]
    for k in range(1, count):
        if k <= d:
            value = -sum((coeff[d - j] * values[k - j] for j in range(1, k))) - k * coeff[d - k]
        else:
            value = -sum((coeff[d - j] * values[k - j] for j in range(1, d + 1)))
        values.append(value)
    return values

def audit(matrices, *, generator_coefficients=(1, 3, 5, 7), with_basis=True, progress=None):
    """Audit an integer symmetric commuting Williamson-type quadruple.

    The compact recipe regenerates every rational source polynomial and CRT
    projector. Every polynomial identity is checked on the full matrix.
    Actual representation multiplicities are retained.
    """
    _require_flint()
    started = time.monotonic()
    array = np.asarray(matrices)
    if array.ndim != 3 or array.shape[0] != 4 or array.shape[1] != array.shape[2] or (array.dtype.kind not in 'iu'):
        raise ValueError('four square integer matrices are required')
    w = array.shape[1]
    if w < 1 or not np.all(np.abs(array) == 1):
        raise ValueError('four sign matrices are required')
    coefficients_tuple = tuple(generator_coefficients)
    if len(coefficients_tuple) != 4 or any((type(c) is not int for c in coefficients_tuple)):
        raise ValueError('four integer generator coefficients are required')
    A = [fmpz_mat(M.tolist()) for M in array]
    identity = I(w)
    if not all((M == M.transpose() for M in A)):
        raise ValueError('source matrices must be symmetric')
    if not all((A[i] * A[j] == A[j] * A[i] for i in range(4) for j in range(i))):
        raise ValueError('source matrices do not commute')
    if sum((M * M.transpose() for M in A), Z(w)) != 4 * w * identity:
        raise ValueError('summed source Gram is not 4wI')
    T = sum((c * M for c, M in zip(coefficients_tuple, A)), Z(w))
    chi = T.charpoly()
    unit, factors = chi.factor()
    _check(int(unit) == 1, 'exact scalar-algebra identity failed')
    factors = sorted(factors, key=lambda fe: (fe[0].degree(), [int(c) for c in fe[0].coeffs()]))
    minimal = fmpz_poly([1])
    rebuilt = fmpz_poly([1])
    for f, multiplicity in factors:
        minimal *= f
        rebuilt *= f ** multiplicity
    _check(chi == rebuilt and minimal.gcd(minimal.derivative()) == 1, 'exact scalar-algebra identity failed')
    _check(T.minpoly() == minimal, 'exact scalar-algebra identity failed')
    D = minimal.degree()
    powers = [identity]
    for _ in range(D):
        powers.append(powers[-1] * T)
    _check(eval_poly(minimal, powers, w) == fmpq_mat(w, w), 'exact scalar-algebra identity failed')
    if progress:
        progress({'stage': 'powers', 'width': w, 'algebra_dimension': D, 'factor_degrees_multiplicities': [[f.degree(), int(m)] for f, m in factors], 'seconds': time.monotonic() - started})
    v_index = None
    for j in range(w):
        K = fmpz_mat([[powers[k][i, j] for k in range(D)] for i in range(w)])
        if K.rank() == D:
            v_index = j
            break
    if v_index is None:
        raise ArithmeticError('no standard-column cyclic vector; choose an explicit generic vector')
    gram = K.transpose() * K
    rhs = fmpz_mat([[A[j][i, v_index] for j in range(4)] for i in range(w)])
    coefficients = fmpq_mat(gram).solve(fmpq_mat(K.transpose() * rhs))
    polynomials = []
    for j in range(4):
        poly = fmpq_poly([coefficients[k, j] for k in range(D)])
        _check(eval_poly(poly, powers, w) == fmpq_mat(A[j]), ('FULL source polynomial identity failed', w, j))
        polynomials.append(poly)
    if progress:
        progress({'stage': 'FULL_polynomials', 'width': w, 'seconds': time.monotonic() - started})
    traces = [sum((int(powers[k][i, i]) for i in range(w))) for k in range(D)]
    factor_traces = [newton_traces([int(c) for c in f.coeffs()], D) for f, _ in factors]
    _check(all((traces[k] == sum((int(mult) * tt[k] for (_, mult), tt in zip(factors, factor_traces))) for k in range(D))), 'exact scalar-algebra identity failed')
    _check(sum((f.degree() * int(mult) for f, mult in factors)) == w, 'exact scalar-algebra identity failed')
    output = {
        'schema': 'exact-WT-operator-module-algebra-v1',
        'width': w,
        'source_matrices_sha256': hashlib.sha256(array.astype(np.int8).tobytes()).hexdigest(),
        'generator_coefficients': list(coefficients_tuple),
        'characteristic_polynomial': intpoly(chi),
        'squarefree_minimal_polynomial': intpoly(minimal),
        'algebra_dimension': D,
        'representation_dimension': w,
        'factors': [{'polynomial': intpoly(f), 'degree': f.degree(), 'representation_multiplicity': int(mult)} for f, mult in factors],
        'source_polynomials': [polynomial_record(p) for p in polynomials],
        'source_polynomial_reconstruction': 'K=(T^k e_j)_(0<=k<degree); solve(K^T K)c=K^T A_i e_j exactly, then verify sum_k c_k T^k=A_i on the FULL matrix; j is cyclic_test_vector_standard_index.',
        'cyclic_test_vector_standard_index': v_index,
        'verification': {
            'full_source_sign_Gram_symmetry_commutation': True,
            'characteristic_factorization_complete': True,
            'minimal_polynomial_squarefree': True,
            'full_minimal_polynomial_matrix_annihilation': True,
            'FLINT_minimal_polynomial_agrees': True,
            'full_four_source_polynomial_matrix_identities': True,
            'source_generated_algebra_equals_Q_T': True,
            'trace_identity_checked_on_all_power_basis_elements': D
        },
        'involution': 'Identity on Q[T], since all source matrices are symmetric and commute.',
        'trace_functional': 'tau(a)=Tr_V(rho(a))/w = sum_i multiplicity_i*Tr_{K_i/Q}(a_i)/w; K_i=Q[t]/factor_i.',
        'source_to_endpoint': {
            'weighted_source_order': 2524,
            'weighted_source_Gram_scale': 2524 * w,
            'operator_module_dimension': w,
            'module_metric': 'I_w',
            'source_entry_palette': 'signed copies of the four source algebra elements p_j(t)',
            'terminal_frame': 'I_w',
            'actual_scalar_image_palette': 'signed copies of the four verified sign matrices A_j',
            'expanded_order': 2524 * w,
            'entry_and_Gram_scales_agree': True,
            'novelty_claim': False,
            'conditional_on_retained_T631_orthogonal_design': True,
            'endpoint_constructed_by_this_algebra_audit': False
        }
    }
    if with_basis:
        M = fmpq_poly(minimal)
        Tq = fmpq_mat(T)
        projectors = []
        all_columns = []
        blocks = []
        for factor_index, (f, mult) in enumerate(factors):
            F = fmpq_poly(f)
            complement = M // F
            g, s, t = complement.xgcd(F)
            _check(g == 1, 'exact scalar-algebra identity failed')
            e = complement * s % M
            _check(all((e % fmpq_poly(other) == (1 if k == factor_index else 0) for k, (other, _) in enumerate(factors))), 'exact scalar-algebra identity failed')
            E = eval_poly(e, powers, w)
            _check(E == E.transpose() and E * E == E, 'exact scalar-algebra identity failed')
            _check(all((E * P == fmpq_mat(w, w) for P in projectors)), 'exact scalar-algebra identity failed')
            projectors.append(E)
            degree = f.degree()
            selected = []
            selected_columns = []
            for j in range(w):
                v = column(E, j)
                orbit = []
                for _ in range(degree):
                    orbit.append(v)
                    v = Tq * v
                candidate = column_matrix(selected_columns + orbit, w)
                rank = candidate.rank()
                if rank == len(selected_columns):
                    continue
                _check(rank == len(selected_columns) + degree, 'exact scalar-algebra identity failed')
                selected.append(j)
                selected_columns += orbit
                if len(selected) == int(mult):
                    break
            _check(len(selected) == int(mult), 'exact scalar-algebra identity failed')
            all_columns += selected_columns
            blocks.append({'factor_index': factor_index, 'degree': degree, 'multiplicity': int(mult), 'projector_polynomial': polynomial_record(e), 'generator_standard_column_indices': selected, 'basis_order': 'for each selected j, e_i(T)e_j,T e_i(T)e_j,...,T^(degree-1)e_i(T)e_j'})
        _check(sum(projectors, fmpq_mat(w, w)) == fmpq_mat(identity), 'exact scalar-algebra identity failed')
        B = column_matrix(all_columns, w)
        _check(B.ncols() == w and B.rank() == w, 'exact scalar-algebra identity failed')
        C = fmpz_mat(w, w)
        offset = 0
        for f, mult in factors:
            degree = f.degree()
            coeff = [int(c) for c in f.coeffs()]
            for _ in range(int(mult)):
                for k in range(degree - 1):
                    C[offset + k + 1, offset + k] = 1
                for k in range(degree):
                    C[offset + k, offset + degree - 1] = -coeff[k]
                offset += degree
        _check(Tq * B == B * fmpq_mat(C), 'exact scalar-algebra identity failed')
        G = B.transpose() * B
        _check(G == G.transpose() and fmpq_mat(C).transpose() * G == G * fmpq_mat(C), 'exact scalar-algebra identity failed')
        output['module_basis_recipe'] = {
            'projector_reconstruction': 'For M=minimal_polynomial and F=factor_i, set L=M/F and e_i=L*(L inverse modF) modM. Use the unique representative of degree<degree(M).',
            'companion_convention': 'Power basis1,t,...,t^(d-1): ones below the diagonal and final column minus the ascending coefficients of the monic factor.',
            'projector_blocks': blocks,
            'B_sha256': matrix_sha(B),
            'G_equals_B_transpose_B_sha256': matrix_sha(G),
            'hash_encoding': 'ASCII rows cols newline, then every row-major rational entry as numerator/denominator newline',
            'B_full_rank': w,
            'full_T_B_equals_B_companions': True,
            'equivalent_B_inverse_T_B_is_repeated_companions': True,
            'full_companion_transpose_G_equals_G_companion': True,
            'G_positive_definite_reason': 'B is invertible over Q and G=B^T B.'
        }
        output['verification']['full_CRT_projector_idempotence_orthogonality_and_sum'] = True
    output['seconds'] = time.monotonic() - started
    if progress:
        progress({'stage': 'COMPLETE', 'width': w, 'seconds': output['seconds'], 'algebra_dimension': D})
    return output

def stable_record(record):
    """Remove timing fields recursively for an exact replay comparison."""
    if isinstance(record, dict):
        return {key: stable_record(value) for key, value in record.items() if key != 'seconds'}
    if isinstance(record, list):
        return [stable_record(value) for value in record]
    return record

def verify(matrices, record):
    """Regenerate every exact object and require the saved compact recipe."""
    actual = audit(matrices, generator_coefficients=tuple(record['generator_coefficients']), with_basis='module_basis_recipe' in record)
    expected = {key: record[key] for key in actual}
    if stable_record(actual) != stable_record(expected):
        raise ArithmeticError('scalar-algebra record does not replay exactly')
    return actual['verification']

def matrices_from_entry(entry):
    """Reconstruct the public literal source, including its first reversal."""
    w = entry['width']
    words = entry['first_rows_pm']
    if len(words) != 4 or any((len(word) != w or not set(word) <= {'+', '-'} for word in words)):
        raise ValueError('four literal sign words of the declared width are required')
    rows = np.array([[1 if sign == '+' else -1 for sign in word] for word in words], dtype=np.int64)
    index = (np.arange(w)[None, :] - np.arange(w)[:, None]) % w
    matrices = np.array([row[index] for row in rows], dtype=np.int64)
    if entry['first_matrix_reversed']:
        matrices[0] = matrices[0][:, ::-1]
    expected = entry['verification']['WT_matrices_int8_row_major_sha256']
    if hashlib.sha256(matrices.astype(np.int8).tobytes()).hexdigest() != expected:
        raise ValueError('literal WT source matrix hash differs')
    return matrices

def main(argv=None):
    """Write or replay one compact certificate for selected literal sources."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, default=INPUT)
    parser.add_argument('--widths', nargs='+', type=int, help='default: every retained source width')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--check', type=Path)
    parser.add_argument('--skip-basis', action='store_true', help='omit optional CRT module-basis reconstruction')
    parser.add_argument('--progress', action='store_true')
    args = parser.parse_args(argv)
    input_bytes = args.inputs.read_bytes()
    inputs = json.loads(input_bytes)
    entries = {entry['width']: entry for entry in inputs['entries']}
    widths = sorted(entries) if args.widths is None else args.widths
    unknown = set(widths) - entries.keys()
    if unknown:
        parser.error(f'no literal source for widths {sorted(unknown)}')
    callback = (lambda result: print(json.dumps(result), flush=True)) if args.progress else None
    records = []
    for width in widths:
        entry = entries[width]
        result = audit(matrices_from_entry(entry), with_basis=not args.skip_basis, progress=callback)
        result['first_matrix_reversed'] = entry['first_matrix_reversed']
        records.append(result)
        gc.collect()
    report = {
        'schema': 'literal-near-Williamson-scalar-algebras-v1',
        'literal_inputs_sha256': hashlib.sha256(input_bytes).hexdigest(),
        'scope': 'Exact factor, multiplicity, source-polynomial and CRT module reconstruction for complete literal sources; no novelty claim.',
        'polynomial_coefficient_convention': 'Ascending integer coefficients; rational polynomials are regenerated and canonically hashed instead of dumped.',
        'records': records,
        'verified_widths': widths,
        'all_source_polynomial_identities_verified_as_FULL_matrices': True
    }
    if args.check:
        expected = json.loads(args.check.read_text())
        if stable_record(report) != stable_record(expected):
            raise ArithmeticError('saved scalar-algebra certificate differs from exact replay')
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + '\n')
    if not args.output and (not args.check):
        print(json.dumps(report, indent=2))
    else:
        print(json.dumps({'verified_widths': widths, 'check_passed': bool(args.check)}))
    return 0
if __name__ == '__main__':
    raise SystemExit(main())

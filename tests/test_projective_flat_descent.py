import itertools
import json
from pathlib import Path

import numpy as np
import pytest
from sympy import Matrix, Poly, eye, kronecker_product, symbols

from experiments.scalar_lift.forward_butson_descent import verify_source
from experiments.scalar_lift.projective_flat_descent import (
    ROOTS,
    descend,
    dita_family,
    endpoint_pair,
    flatness_pullback,
    haagerup_histogram,
    quaternary_source_family,
    real_flatness_pullback,
    restriction_tower,
    sign_row_difference_rank,
    verify_laurent_family,
    verify_quaternary_gram,
    verify_sign_gram,
    verify_sixth_gram,
)
from experiments.scalar_lift.relative_scalar_restriction import (
    RationalArray,
    conjugate_transpose,
    regular_matrix,
)
from hadamard.coefficients import CyclotomicField
from hadamard.hadamard.family import MonomialHadamardFamily

ROOT = Path(__file__).resolve().parents[1]


def regular(data):
    return regular_matrix(RationalArray(np.array(data['numerators'], dtype=np.int64), data['denominator']), 6)


def adjoint_from_json(data):
    value = RationalArray(np.array(data['numerators'], dtype=np.int64), data['denominator'])
    return regular_matrix(conjugate_transpose(value, 6), 6)


@pytest.mark.parametrize('variant', [1, 2])
def test_each_flatness_pullback_is_an_exact_projective_polynomial_identity(variant):
    report = flatness_pullback(variant)
    _, _, blocks = endpoint_pair(variant)
    x = symbols('x0:6')
    h = symbols('h')
    lam = (h-sum(x))/7
    diagonal = sum(v*v-h*v-2*h*lam for v in x)
    for k, weights in enumerate(report['cross_coefficients_in_target_row_major_order']):
        i, j = divmod(k, 6)
        coefficients = ROOTS[blocks[:6, i, j]]
        a = sum(int(coefficients[t, 0])*x[t] for t in range(6))
        b = sum(int(coefficients[t, 1])*x[t] for t in range(6))
        # Norm(a+b*zeta6)=a²+ab+b², computed independently of orbit multiplication.
        norm = a*a+a*b+b*b-h*h
        certificate = diagonal+sum(weight*(x[r]*x[s]-h*lam)
                                   for weight, (r, s) in zip(weights, report['cross_relation_pairs'], strict=True))
        assert Poly(norm-certificate, *x, h).is_zero
    # This is genuinely an endpoint relation, not flatness on all input coordinates.
    assert Poly(a*a+a*b+b*b-h*h, *x, h) != Poly(0, *x, h)


@pytest.mark.parametrize('variant', [1, 2])
def test_nonflat_tower_and_terminal_maps_reconstruct_a_complete_source(variant):
    report = restriction_tower(variant)
    raw = regular(report['raw_zeta7_action'])
    metric = regular(report['transported_metric'])
    c, d = regular(report['left_endpoint_map_C']), regular(report['right_endpoint_frame_D'])
    assert metric == adjoint_from_json(report['left_endpoint_map_C'])*c
    assert d*adjoint_from_json(report['right_endpoint_frame_D']) == 6*metric.inv()
    assert metric != regular(report['canonical_unweighted_trace_metric'])
    assert adjoint_from_json(report['raw_zeta7_action'])*metric*raw == metric
    assert any(value == 0 for value in raw)
    source = np.outer(np.arange(7), np.arange(7)) % 7
    restricted = Matrix.vstack(*(Matrix.hstack(*(raw**int(e) for e in row)) for row in source))
    terminal = kronecker_product(eye(7), c)*restricted*kronecker_product(eye(7), d)
    result = descend(source, 7, variant=variant)
    assert terminal == regular_matrix(RationalArray.exact(ROOTS[result]), 6)
    assert verify_source(result, 6)['full_Gram_verified']


def test_base_extension_uses_the_actual_lower_root_alphabet():
    source = (12*np.outer(np.arange(7), np.arange(7))+7*np.arange(7)[:, None]) % 84
    result = descend(source, 84, 12)
    assert result.shape == (42, 42)
    assert np.any(result % 2)  # Actual twelfth roots survive the lower-field gauge.
    assert verify_source(result, 12)['full_Gram_verified']
    ordinary = descend((12*np.outer(np.arange(7), np.arange(7))) % 84, 84, 12)
    expected = (ordinary+np.repeat(np.arange(7), 6)[:, None]) % 12
    assert np.array_equal(result, expected)


@pytest.mark.parametrize('source,upper,lower', [([[0.5]], 42, 6), ([[True]], 42, 6),
                                                ([[0]], 42, 7), ([[0]], 42, 42),
                                                ([[0]], 30, 6), ([[0, 0], [0, 0]], 42, 6)])
def test_unverified_sources_and_invalid_field_edges_are_rejected(source, upper, lower):
    with pytest.raises(ValueError):
        descend(source, upper, lower)


def test_fast_laurent_certificate_agrees_with_generic_exact_family_engine():
    report = json.loads((ROOT/'certificates/scalar-lift/projective-flat-descent.json').read_text())
    with np.load(ROOT/report['literal_matrix_bundle']['path'], allow_pickle=False) as stored:
        base = stored['family42_base']
        labels = stored['family42_labels']
    assert verify_laurent_family(base, labels)['full_row_and_column_Laurent_Grams_verified']
    field = CyclotomicField(6)
    roots = [field.zeta()**j for j in range(6)]
    family = MonomialHadamardFamily('scalar-reduced order 42', field, ('t',),
        tuple(tuple(roots[e] for e in row) for row in base),
        tuple(tuple((int(k >= 0),) for k in row) for row in labels))
    assert family.verify().valid
    # This corruption is invisible at t=1: only the symbolic coefficient test detects it.
    broken = labels.copy()
    broken[1, 1] = 0 if broken[1, 1] < 0 else -1
    with pytest.raises(ValueError, match='Laurent Gram coefficient'):
        verify_laurent_family(base, broken)


def test_order294_family_has_36_dephased_independent_parameters_and_exact_specializations():
    _, base, labels = dita_family()
    proof = verify_laurent_family(base, labels)
    assert proof['independent_unit_phase_parameters'] == 36
    assert proof['dephased_for_all_parameters']
    assert set(labels.ravel()) == set(range(-1, 36))
    verify_sixth_gram(base)
    # Independent cyclotomic zero-sum implementation checks the complete column Gram.
    assert verify_source(base, 6)['full_Gram_verified']
    damaged = base.copy()
    damaged[1, 1] = (damaged[1, 1]+1) % 6
    with pytest.raises(ValueError, match='Gram'):
        verify_sixth_gram(damaged)


def test_haagerup_counts_match_independent_quadruples_and_ignore_monomial_gauges():
    matrix = 3*np.array([[0, 0, 0, 0], [0, 1, 0, 1], [0, 0, 1, 1], [0, 1, 1, 0]])
    brute = [0]*6
    for i, j, k, column in itertools.product(range(4), repeat=4):
        brute[int((matrix[i, j]+matrix[k, column]-matrix[i, column]-matrix[k, j]) % 6)] += 1
    assert haagerup_histogram(matrix) == brute
    gauged = (matrix[[3, 1, 0, 2]][:, [1, 3, 2, 0]]
              +np.arange(4)[:, None]+2*np.arange(4)) % 6
    assert haagerup_histogram(gauged) == brute


def test_saved_examples_prove_lower_bounds_without_claiming_a_complete_classification():
    report = json.loads((ROOT/'certificates/scalar-lift/projective-flat-descent.json').read_text())
    with np.load(ROOT/report['literal_matrix_bundle']['path'], allow_pickle=False) as stored:
        fingerprints = []
        for row in report['order294_examples']:
            matrix = stored[row['matrix_key']]
            assert haagerup_histogram(matrix) == row['haagerup_histogram']
            verify_sixth_gram(matrix)
            fingerprints.append(tuple(row['haagerup_histogram']))
        assert len(set(fingerprints)) == 5
        assert len({tuple(row['haagerup_histogram']) for row in report['order42_examples']}) == 4
    assert report['order42_monomial_class_lower_bound'] == 4
    assert report['order294_monomial_class_lower_bound'] == 5
    assert not report['scope']['previously_unknown_order_claimed']
    assert report['scope']['real_sign_endpoint_claimed']
    assert not report['scope']['complete_output_classification_claimed']


def test_real_endpoint_flatness_and_norm_factor_through_two_source_quadrics():
    report = real_flatness_pullback()
    x, y, h = symbols('x y h')
    relations = [x*x+y*y-h*h, x*y]
    raw = Matrix([[x, y], [-y, x]])
    metric = Matrix(report['transported_canonical_trace_metric'])
    c, d = Matrix(report['left_endpoint_map_C']), Matrix(report['right_endpoint_map_D'])
    assert raw.T*metric*raw-h*h*metric == relations[0]*metric
    assert c.T*c == metric and d*d.T == 2*metric.inv()
    for value, (a, b) in zip(c*raw*d, report['flatness_ideal_coefficients'], strict=True):
        assert Poly(value*value-h*h-a*relations[0]-b*relations[1], x, y, h).is_zero


def test_real_family_contains_five_fully_verified_inequivalent_H768():
    base, labels, family = quaternary_source_family()
    assert set(labels.ravel()) == set(range(-1, 285))
    assert family['independent_mu4_parameters'] == 285
    verify_quaternary_gram(base)
    report = json.loads((ROOT/'certificates/scalar-lift/projective-flat-descent.json').read_text())
    real = report['real768']
    assert real['binary_code_dimensions'] == [17, 18, 19, 20, 21]
    with np.load(ROOT/report['literal_matrix_bundle']['path'], allow_pickle=False) as stored:
        for i, item in enumerate(real['examples']):
            matrix = stored[item['matrix_key']]
            assert matrix.shape == (768, 768)
            assert verify_sign_gram(matrix)['sign_matrix_sha256'] == item['verification']['sign_matrix_sha256']
            assert sign_row_difference_rank(matrix) == 17+i
            assert verify_quaternary_gram(stored[f'bh384_source{i}']) == item['source_proof']
    assert not real['previously_unknown_real_order_claimed']


def test_real_code_rank_ignores_signed_permutations_and_bad_source_grams_fail():
    matrix = np.array([[1, 1], [1, -1]], dtype=np.int64)
    for _ in range(2):
        matrix = np.kron(matrix, [[1, 1], [1, -1]])
    reference = sign_row_difference_rank(matrix)
    changed = matrix[[3, 6, 0, 2, 1, 5, 7, 4]][:, [5, 3, 0, 6, 7, 1, 2, 4]]
    changed = changed*np.array([1, -1]*4)[:, None]*np.array([-1, -1, 1, 1]*2)
    assert sign_row_difference_rank(changed) == reference
    assert sign_row_difference_rank(changed.T) == reference
    broken = matrix.copy()
    broken[1, 1] *= -1
    with pytest.raises(ValueError, match='Gram'):
        verify_sign_gram(broken)
    bad_source = np.zeros((2, 2), dtype=np.int64)
    with pytest.raises(ValueError, match='source Gram'):
        verify_quaternary_gram(bad_source)

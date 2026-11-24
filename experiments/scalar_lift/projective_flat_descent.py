"""Exact endpoint flatness identities and explicit scalar-reduced families.

The nonmonomial (7,6) pairs are from Egan--Ó Catháin, arXiv:2610.28394,
Theorem 3.1. Their input data are checked here, not trusted as certificates.
Restriction passes through a nonflat quadratic/cubic tower. A weighted trace
form is transported throughout, and a final frame makes the entries flat.
No claim of a previously unknown order or unpublished matrix class is made.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import operator
from functools import cache
from pathlib import Path

import numpy as np
from sympy import Matrix, Rational, eye, kronecker_product

from experiments.scalar_lift.forward_butson_descent import exponent_hash, verify_source
from experiments.scalar_lift.relative_scalar_restriction import (
    RationalArray,
    conjugate_transpose,
    edge,
    regular_matrix,
)
from experiments.scalar_lift.vanishing_collision_search import integral_roots
from hadamard.dense_lift import QUATERNARY_REAL_BLOCKS

ROOT = Path(__file__).resolve().parents[2]
PAPER = "https://arxiv.org/abs/2610.28394"
ROOTS = np.asarray(integral_roots(6), dtype=np.int64)
ROOTS.setflags(write=False)
FRAME_EXPONENTS = 2 * np.array([
    [0, 0, 0, 0, 0, 0], [0, 0, 1, 1, 2, 2], [0, 1, 0, 2, 1, 2],
    [0, 2, 1, 2, 0, 1], [0, 2, 2, 1, 1, 0], [0, 1, 2, 0, 2, 1],
], dtype=np.int64)
FRAME_EXPONENTS.setflags(write=False)


def _regular(coefficients):
    return regular_matrix(RationalArray.exact(coefficients), 6)


def _coefficients(matrix):
    rows, columns = matrix.rows // 2, matrix.cols // 2
    result = np.array([[[matrix[2*i+a, 2*j] for a in range(2)]
                        for j in range(columns)] for i in range(rows)], dtype=object)
    if _regular(result) != matrix:
        raise ValueError("matrix is not regular over the declared Eisenstein field")
    return result


def _adjoint(matrix):
    return regular_matrix(conjugate_transpose(RationalArray.exact(_coefficients(matrix)), 6), 6)


def _json_matrix(matrix):
    return RationalArray.exact(_coefficients(matrix)).to_json()


def _roots_from_regular(matrix):
    coefficients = _coefficients(matrix)
    lookup = {tuple(map(int, root)): exponent for exponent, root in enumerate(ROOTS)}
    try:
        return np.array([[lookup[tuple(value)] for value in row]
                         for row in coefficients], dtype=np.int64)
    except KeyError as error:
        raise ValueError("a final entry is not a sixth root") from error


def _integer(value, name):
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be an exact integer")
    try:
        return operator.index(value)
    except TypeError as error:
        raise ValueError(f"{name} must be an exact integer") from error


def _square(values, phase):
    phase = _integer(phase, "phase")
    if phase < 1:
        raise ValueError("phase must be positive")
    values = np.asarray(values, dtype=object)
    if values.ndim != 2 or not len(values) or values.shape[0] != values.shape[1]:
        raise ValueError("expected a nonempty square exponent matrix")
    return np.array([[_integer(value, "exponent") % phase for value in row]
                     for row in values], dtype=np.int64)


@cache
def endpoint_pair(variant=1):
    """Return an exactly verified pair and its seven flat scalar blocks."""
    variant = _integer(variant, "variant")
    if variant == 1:
        exponents = np.array([
            [-1, -1, -1, -1, -1, 0], [0, 2, 0, 5, -1, -1],
            [-1, 0, 2, 5, 0, -1], [3, 3, -1, 4, 1, -1],
            [4, 2, 2, -1, 4, -1], [4, -1, 0, 1, 0, -1],
        ])
        amplitudes = (exponents >= 0).astype(np.int64)
        amplitudes[0, 5] = 2
        denominator = 2
    elif variant == 2:
        exponents = np.array([
            [0, 0, 0, 4, 1, 4], [2, 2, 5, 0, 0, 0],
            [2, 2, 2, 3, 0, 0], [0, 3, 0, 4, 4, 4],
            [3, 0, 0, 4, 4, 4], [2, 2, 2, 0, 0, 3],
        ])
        amplitudes = np.ones((6, 6), dtype=np.int64)
        for i, j in ((0, 4), (1, 2), (2, 3), (3, 1), (4, 0), (5, 5)):
            amplitudes[i, j] = 2
        denominator = 3
    else:
        raise ValueError("the published pair variant must be 1 or 2")
    numerators = ROOTS[exponents % 6] * amplitudes[:, :, None]
    u = _regular(np.array([Rational(int(value), denominator) for value in numerators.flat],
                         dtype=object).reshape(6, 6, 2))
    z = _regular(ROOTS[FRAME_EXPONENTS])
    if u*_adjoint(u) != eye(12) or z*_adjoint(z) != 6*eye(12):
        raise ValueError("the published action/frame Gram identities failed")
    if sum((u**j for j in range(7)), Matrix.zeros(12)) != Matrix.zeros(12):
        raise ValueError("the scalar action fails Phi_7")
    blocks = np.array([_roots_from_regular(u**j*z) for j in range(7)])
    blocks.setflags(write=False)
    return u, z, blocks


def descend(exponents, upper_phase=42, lower_phase=6, variant=1):
    """BH(n,q) -> BH(6n,s), q | 7s, 6 | s and 7 does not divide s."""
    upper_phase = _integer(upper_phase, "upper phase")
    lower_phase = _integer(lower_phase, "lower phase")
    if upper_phase < 1 or lower_phase < 1 or lower_phase % 6 or lower_phase % 7 == 0:
        raise ValueError("the lower alphabet must contain mu_6 and have no factor seven")
    if (7*lower_phase) % upper_phase:
        raise ValueError("the declared upper alphabet must divide 7 times the lower alphabet")
    source = _square(exponents, upper_phase)
    verify_source(source, upper_phase)
    _, _, blocks = endpoint_pair(variant)
    embedded = ((7*lower_phase)//upper_phase)*source
    seventh = (pow(lower_phase, -1, 7)*embedded) % 7
    lower = (pow(7, -1, lower_phase)*embedded) % lower_phase
    result = (lower[..., None, None] + (lower_phase//6)*blocks[seventh]) % lower_phase
    result = result.transpose(0, 2, 1, 3).reshape(6*len(source), 6*len(source))
    result.setflags(write=False)
    return result


def flatness_pullback(variant=1):
    """Give coefficientwise identities, not numerical flatness residuals.

    Put lambda=(h-sum(x))/7, g_i=x_i²-h*x_i-2*h*lambda and
    g_ij=x_i*x_j-h*lambda. On the seven projective source points
    (h,x)=(1,e_i),(1,-1,...,-1), every g is zero. For EACH final
    coordinate b=sum_i c_i*x_i we verify the polynomial identity

        b*conjugate(b)-h² = sum_i g_i + sum_(i<j) k_ij*g_ij,

    with k_ij=2 Re(c_i conjugate(c_j)). Its only remaining coefficient
    identity is sum k_ij=-5, verified exactly over Q(zeta_6).
    """
    _, _, blocks = endpoint_pair(variant)
    pairs = [(i, j) for i in range(6) for j in range(i+1, 6)]
    weights = []
    for row in range(6):
        for column in range(6):
            differences = [(int(blocks[i, row, column])-int(blocks[j, row, column])) % 6
                           for i, j in pairs]
            # The trace from Q(zeta_6) to Q is 2*a+b in basis (1,zeta_6).
            values = [int(2*ROOTS[e, 0]+ROOTS[e, 1]) for e in differences]
            if sum(values) != -5:
                raise ValueError("a pulled-back projective flatness identity failed")
            weights.append(values)
    return {"source_projective_coordinates": ["h", *[f"x{i}" for i in range(6)]],
            "source_points_on_h_equals_one": [*np.eye(6, dtype=int).tolist(), [-1]*6],
            "lambda": "(h-sum(x))/7", "diagonal_relations": "x_i^2-h*x_i-2*h*lambda",
            "cross_relations": "x_i*x_j-h*lambda", "cross_relation_pairs": [list(pair) for pair in pairs],
            "target_identity": "|sum(c_i*x_i)|^2-h^2=sum(g_i)+sum(k_ij*g_ij)",
            "cross_coefficients_in_target_row_major_order": weights,
            "all_36_polynomial_identities_verified": True,
            "flatness_only_required_at_endpoint": True,
            "lower_scalar_unit_phases_preserve_the_identity": True}


def restriction_tower(variant=1):
    """Transport a chosen positive trace form through a genuine 2-by-3 tower.

    L=K(zeta_7), E=K(eta), K=Q(zeta_6), eta=zeta_7+zeta_7^-1,
    xi=zeta_7-zeta_7^-1. The source starts with form h*I, not I:
    H*H=nI implies H*(hI)H=n(hI) for every central Hermitian h.
    We explicitly derive h and never reset a metric during restriction.
    """
    u, z, _ = endpoint_pair(variant)
    relative = edge(42, 6, (0, 6, 12, 18, 24, 30))
    raw = regular_matrix(relative.action(6), 6)
    vector = _regular(np.eye(6, dtype=int)[:, 0, None, None]*ROOTS[0])
    p = Matrix.hstack(*(u**j*vector for j in range(6)))
    if not p.det() or u*p != p*raw:
        raise ValueError("cyclic basis does not intertwine the actual restriction")
    gram_power = _adjoint(p)*p
    traces = [6, -1, -1, -1, -1, -1, -1]
    trace_system = Matrix([[traces[(a+b) % 7] for b in range(6)] for a in range(6)])
    h = np.array((trace_system.inv()*Matrix(_coefficients(gram_power)[0].tolist())).tolist(),
                 dtype=object)
    h_operator = sum((kronecker_product(eye(6), _regular(h[j][None, None, :]))*raw**j
                      for j in range(6)), Matrix.zeros(12))
    canonical = regular_matrix(relative.metric, 6)
    if canonical*h_operator != gram_power or _adjoint(h_operator)*canonical != canonical*h_operator:
        raise ValueError("the chosen incoming form is not the advertised Hermitian trace weight")

    def multiply(a, b):
        values = np.convolve(np.array(a, dtype=object).ravel(), np.array(b, dtype=object).ravel()).tolist()
        for j in range(len(values)-1, 5, -1):
            leading = values[j]
            for k in range(7):
                values[j-6+k] -= leading
        return Matrix(values[:6])

    one, x = eye(6)[:, 0], eye(6)[:, 1]
    inverse_x = -Matrix.ones(6, 1)
    eta, xi = x+inverse_x, x-inverse_x
    eta2 = multiply(eta, eta)
    basis = Matrix.hstack(one, eta, eta2, xi, multiply(eta, xi), multiply(eta2, xi))
    t = Matrix([[0, 0, 1], [1, 0, 2], [0, 1, -1]])
    a = t*t-4*eye(3)
    r = Matrix.vstack(Matrix.hstack(t, a), Matrix.hstack(eye(3), t))/2
    basis_regular = kronecker_product(basis, eye(2))
    r_regular = kronecker_product(r, eye(2))
    if basis_regular.inv()*raw*basis_regular != r_regular:
        raise ValueError("quadratic/cubic tower differs from direct scalar multiplication")
    c = p*basis_regular
    metric = _adjoint(c)*c
    h_tower = np.array((basis.inv()*Matrix(h.tolist())).tolist(), dtype=object)

    def cubic_operator(poly):
        return sum((kronecker_product(t**j, _regular(poly[j][None, None, :]))
                    for j in range(3)), Matrix.zeros(6))

    h0, h1 = cubic_operator(h_tower[:3]), cubic_operator(h_tower[3:])
    aa = kronecker_product(a, eye(2))
    first_metric = 2*Matrix.vstack(Matrix.hstack(h0, h1*aa), Matrix.hstack(-h1*aa, -h0*aa))
    cubic_trace = Matrix([[int((t**(i+j)).trace()) for j in range(3)] for i in range(3)])
    transferred = kronecker_product(eye(2), cubic_trace, eye(2))*first_metric
    first_adjoint = kronecker_product(
        Matrix.vstack(Matrix.hstack(t, eye(3)), Matrix.hstack(a, t))/2, eye(2))
    if first_adjoint*first_metric*r_regular != first_metric or transferred != metric:
        raise ValueError("incoming-form transport or trace transitivity failed")
    if _adjoint(r_regular)*metric*r_regular != metric:
        raise ValueError("the final weighted norm identity failed")
    d = c.inv()*z
    if d*_adjoint(d) != 6*metric.inv():
        raise ValueError("the terminal frame has the wrong norm scale")
    for j in range(7):
        if c*r_regular**j*d != u**j*z:
            raise ValueError("the terminal map differs from the certified flat endpoint")
    return {"fields": {"K": "Q(zeta_6)", "E": "K(eta)", "L": "E(xi)=K(zeta_7)",
                       "eta_minimal_polynomial": [-1, -2, 1, 1],
                       "xi_squared": "eta^2-4", "relative_degrees": [2, 3]},
            "initial_source_form": "h I; H^*(h I)H=n(h I) follows from the verified source",
            "h_in_zeta7_power_basis": RationalArray.exact(h).to_json(),
            "h_in_tower_basis": RationalArray.exact(h_tower).to_json(),
            "h_positive": "its transferred form is C^*C with C nonsingular at both embeddings of K",
            "tower_basis_in_zeta7_coordinates": np.array(basis.tolist(), dtype=int).tolist(),
            "raw_zeta7_action": _json_matrix(r_regular),
            "first_metric_in_regular_E_coordinates": _json_matrix(first_metric),
            "cubic_trace_metric": [[int(value) for value in row] for row in cubic_trace.tolist()],
            "transported_metric": _json_matrix(metric),
            "canonical_unweighted_trace_metric": _json_matrix(basis_regular.T*canonical*basis_regular),
            "unweighted_and_chosen_weighted_forms_kept_distinct": True,
            "left_endpoint_map_C": _json_matrix(c), "right_endpoint_frame_D": _json_matrix(d),
            "first_weighted_norm_verified": True, "trace_transitivity_verified": True,
            "final_weighted_norm_verified": True, "C_star_C_equals_transported_metric": True,
            "D_D_star_equals_6_inverse_metric": True, "all_terminal_scalar_blocks_reconstructed": True,
            "raw_quadratic_action_has_zero_entries": any(value == 0 for value in r),
            "no_intermediate_flatness_assumed": True}


def _cross_gram(left, right):
    a, b = left[..., 0], left[..., 1]
    c, d = right[..., 0], right[..., 1]
    return np.stack((a@(c+d).T+b@d.T, b@c.T-a@d.T), axis=2)


def verify_sixth_gram(exponents):
    matrix = _square(exponents, 6)
    expected = np.zeros((len(matrix), len(matrix), 2), dtype=np.int64)
    expected[np.arange(len(matrix)), np.arange(len(matrix)), 0] = len(matrix)
    for value in (matrix, matrix.T):
        if not np.array_equal(_cross_gram(ROOTS[value], ROOTS[value]), expected):
            raise ValueError("the complete sixth-root row or column Gram failed")
    return {"order": len(matrix), "phase": 6, "both_full_exact_Grams_verified": True,
            "exponent_sha256": exponent_hash(matrix)}


def verify_laurent_family(base, labels):
    """Check every Laurent Gram coefficient, using sparse coefficient supports.

    Each entry is zeta_6^base[i,j]*t_labels[i,j], with t_-1=1.
    Distinct nonconstant quotients t_k/t_l have distinct exponent vectors.
    Coefficient support intersections avoid expanding a 36-variable ideal.
    """
    base = _square(base, 6)
    labels = np.asarray(labels)
    if labels.shape != base.shape or not np.issubdtype(labels.dtype, np.integer) or np.any(labels < -1):
        raise ValueError("expected integer parameter labels, with -1 for constants")
    n = len(base)
    checked = 0
    for exponents, indices in ((base, labels), (base.T, labels.T)):
        components = []
        total = np.zeros((n, n, 2), dtype=np.int64)
        for label in sorted(map(int, np.unique(indices))):
            support = indices == label
            rows, columns = np.flatnonzero(support.any(axis=1)), np.flatnonzero(support.any(axis=0))
            value = ROOTS[exponents[np.ix_(rows, columns)]]*support[np.ix_(rows, columns)][..., None]
            total[np.ix_(rows, rows)] += _cross_gram(value, value)
            components.append((columns, value))
        expected = np.zeros_like(total)
        expected[np.arange(n), np.arange(n), 0] = n
        if not np.array_equal(total, expected):
            raise ValueError("the constant Laurent Gram coefficient failed")
        for i, (left_columns, left) in enumerate(components):
            for right_columns, right in components[i+1:]:
                _, a, b = np.intersect1d(left_columns, right_columns, return_indices=True)
                if not len(a):
                    continue
                left_value, right_value = left[:, a], right[:, b]
                left_value = left_value[np.any(left_value, axis=(1, 2))]
                right_value = right_value[np.any(right_value, axis=(1, 2))]
                if len(left_value) and len(right_value) and np.any(_cross_gram(left_value, right_value)):
                    raise ValueError("a nonconstant Laurent Gram coefficient failed")
                checked += 1
    return {"order": n, "independent_unit_phase_parameters": len(set(labels.ravel())-{-1}),
            "every_entry_is_a_root_times_a_unit_phase_monomial": True,
            "full_row_and_column_Laurent_Grams_verified": True,
            "nonconstant_coefficient_products_checked": checked,
            "dephased_for_all_parameters": bool(np.all(base[0] == 0) and np.all(base[:, 0] == 0)
                                                 and np.all(labels[0] == -1) and np.all(labels[:, 0] == -1))}


def haagerup_histogram(exponents):
    """Counts of all h_ij*h_kl*conj(h_il*h_kj), including repeated indices.

    Different histograms prove monomial inequivalence; equality is no claim
    of equivalence. Counts are also invariant under transpose/conjugation.
    """
    matrix = _square(exponents, 6)
    histogram = np.zeros(6, dtype=np.int64)
    for row in matrix:
        differences = (row-matrix) % 6
        counts = np.array([(differences == e).sum(axis=1) for e in range(6)]).T
        histogram += np.array([sum(int((counts[:, e]*counts[:, (e-t) % 6]).sum())
                                   for e in range(6)) for t in range(6)], dtype=np.int64)
    if int(histogram.sum()) != len(matrix)**4:
        raise ValueError("the invariant did not count every ordered quadruple")
    return list(map(int, histogram))


def dita_family(variant=1):
    """A dephased order-294 family with 36 independent lower-field phases."""
    a, b = np.repeat(np.arange(7), 7), np.tile(np.arange(7), 7)
    source = (a[:, None]*a+b[:, None]*b) % 7
    base = descend(source, upper_phase=7, variant=variant)
    source_labels = np.where((b[:, None] > 0) & (a > 0), 6*(b[:, None]-1)+a-1, -1)
    labels = np.repeat(np.repeat(source_labels, 6, axis=0), 6, axis=1)
    return source, base, labels


def _specialize(base, labels, parameters, phase=6):
    values = np.array([0, *list(parameters)], dtype=np.int64)
    return ((phase//6)*base+values[labels+1]) % phase


def real_flatness_pullback():
    """Four real endpoint identities modulo x²+y²-h² and xy."""
    action = Matrix([[0, 1], [-1, 0]])  # Basis (1,-i) for Q(i)/Q.
    c, d, metric = Matrix([[1, 1], [1, -1]]), eye(2), 2*eye(2)
    if action.T*metric*action != metric or action**2 != -eye(2):
        raise ValueError("the real restriction fails the transported metric or Phi_4")
    if c.T*c != metric or d*d.T != 2*metric.inv():
        raise ValueError("the real endpoint normalization failed")
    if any(c*action**j*d != Matrix(QUATERNARY_REAL_BLOCKS[j]) for j in range(4)):
        raise ValueError("the real restriction differs from the verified quaternary transfer")
    linear = [[1, -1], [1, 1], [1, 1], [-1, 1]]
    return {"field_extension": "Q(i)/Q", "basis": ["1", "-i"],
            "relative_degree": 2, "incoming_source_form": "I",
            "transported_canonical_trace_metric": [[2, 0], [0, 2]],
            "raw_i_action": [[0, 1], [-1, 0]],
            "left_endpoint_map_C": [[1, 1], [1, -1]], "right_endpoint_map_D": [[1, 0], [0, 1]],
            "source_projective_relations": ["x^2+y^2-h^2", "x*y"],
            "terminal_entry_linear_coefficients": linear,
            "flatness_ideal_coefficients": [[1, 2*a*b] for a, b in linear],
            "identity": "(a*x+b*y)^2-h^2=(x^2+y^2-h^2)+2*a*b*x*y",
            "all_four_polynomial_identities_verified": all(a*a == b*b == 1 for a, b in linear),
            "raw_action_nonflat": True, "endpoint_entries_are_real_signs": True}


@cache
def quaternary_source_family():
    """285 independent mu_4 choices in verified order-384 Diţă sources."""
    from experiments.scalar_lift.inverse_blocks import known_order_twelve_pair

    small = np.asarray(known_order_twelve_pair()[0], dtype=np.int64)
    if not np.array_equal(small@small.T, 12*np.eye(12, dtype=np.int64)):
        raise ValueError("the actual order-twelve input failed its complete Gram")
    seed = small.copy()
    for _ in range(3):
        seed = np.kron(seed, [[1, 1], [1, -1]])
    seed = seed*seed[:, 0, None]
    seed = seed*seed[0]
    if not np.array_equal(seed@seed.T, 96*np.eye(96, dtype=np.int64)):
        raise ValueError("the actual tensor input failed its complete Gram")
    a, b = np.repeat(np.arange(96), 4), np.tile(np.arange(4), 96)
    base = (2*(seed[a[:, None], a] < 0)+b[:, None]*b) % 4
    labels = np.where((b[:, None] > 0) & (a > 0), 95*(b[:, None]-1)+a-1, -1)
    base.setflags(write=False)
    labels.setflags(write=False)
    return base, labels, {"order": 384, "phase": 4, "independent_mu4_parameters": 285,
                         "number_of_literal_parameter_choices": "4^285",
                         "formula": "S[(a,b),(c,d)]=H96[a,c]*i^(b*d+p[b,c])",
                         "H96_construction": "verified H12 tensor Sylvester H8, then dephased",
                         "H12_sha256": exponent_hash(small), "H96_sha256": exponent_hash(seed),
                         "full_source_Gram_for_every_parameter_choice":
                             "sum_d i^((b-b')d)=4 delta_bb'; remaining parameter phases cancel; H96 H96^T=96 I",
                         "source_Gram_scale": 384, "relative_degree": 2}


def verify_quaternary_gram(exponents):
    source = _square(exponents, 4)
    a = np.array([1, 0, -1, 0], dtype=np.int64)[source]
    b = np.array([0, 1, 0, -1], dtype=np.int64)[source]
    expected = len(source)*np.eye(len(source), dtype=np.int64)
    for x, y in ((a, b), (a.T, b.T)):
        if not np.array_equal(x@x.T+y@y.T, expected) or np.any(y@x.T-x@y.T):
            raise ValueError("a complete quaternary source Gram failed")
    return {"order": len(source), "phase": 4, "full_exact_row_and_column_Grams_verified": True,
            "actual_used_exponents": np.unique(source).tolist(), "exponent_sha256": exponent_hash(source)}


def _sign_matrix(values):
    matrix = np.asarray(values)
    if (matrix.ndim != 2 or not len(matrix) or matrix.shape[0] != matrix.shape[1]
            or not np.issubdtype(matrix.dtype, np.integer) or not np.all((matrix == 1) | (matrix == -1))):
        raise ValueError("expected a square exact sign matrix")
    return matrix.astype(np.int64)


def verify_sign_gram(values):
    matrix = _sign_matrix(values)
    expected = len(matrix)*np.eye(len(matrix), dtype=np.int64)
    if not np.array_equal(matrix@matrix.T, expected) or not np.array_equal(matrix.T@matrix, expected):
        raise ValueError("the complete real Hadamard Gram failed")
    return {"order": len(matrix), "phase": 2, "all_entries_signs": True,
            "both_full_exact_Grams_verified": True, "sign_matrix_sha256": exponent_hash(matrix)}


def sign_row_difference_rank(values):
    """Dimension over F2 of span(1, b_i+b_0), b_i the negative-entry rows.

    Column sign changes cancel in row differences; row sign changes add 1.
    Row permutations change the reference within the same span, and column
    permutations permute coordinates. Thus differing ranks prove monomial
    inequivalence. No assumption about a fixed dephasing is needed.
    """
    matrix = _sign_matrix(values)
    bits = matrix < 0
    rows = [int.from_bytes(np.packbits(row ^ bits[0], bitorder='little').tobytes(), 'little')
            for row in bits]
    rows.append((1 << len(matrix))-1)
    basis = {}
    for row in rows:
        while row:
            pivot = row.bit_length()-1
            if pivot in basis:
                row ^= basis[pivot]
            else:
                basis[pivot] = row
                break
    return len(basis)


def real_examples():
    base, labels, family_proof = quaternary_source_family()
    patterns = [[], [[1, 1, 2]], [[1, 1, 1]], [[1, 1, 1], [2, 2, 2]], [[1, 1, 1], [2, 2, 1]]]
    a, b = np.repeat(np.arange(96), 4), np.tile(np.arange(4), 96)
    blocks = np.asarray(QUATERNARY_REAL_BLOCKS, dtype=np.int8)
    arrays = {'real_source384_base': base, 'real_source384_labels': labels}
    outputs = []
    for index, pattern in enumerate(patterns):
        parameter = np.zeros((4, 96), dtype=np.int64)
        for row, column, value in pattern:
            parameter[row, column] = value
        source = (base+parameter[b[:, None], a]) % 4
        target = blocks[source].transpose(0, 2, 1, 3).reshape(768, 768)
        key = f'h768_example{index}'
        arrays[key] = target
        arrays[f'bh384_source{index}'] = source
        outputs.append({"matrix_key": key, "nonzero_source_parameter_exponents": pattern,
                        "source_proof": verify_quaternary_gram(source),
                        "verification": verify_sign_gram(target),
                        "binary_row_difference_code_dimension": sign_row_difference_rank(target)})
    ranks = [row['binary_row_difference_code_dimension'] for row in outputs]
    if len(set(ranks)) != len(outputs):
        raise ValueError("the selected real examples do not have distinct code ranks")
    return {"flatness": real_flatness_pullback(), "source_family": family_proof,
            "examples": outputs, "pairwise_monomially_inequivalent_examples": len(outputs),
            "binary_code_dimensions": ranks, "inequivalence_proof": "dimension of span_F2(1,b_i+b_0)",
            "quaternary_to_real_morphism_is_classical": True,
            "previously_unknown_real_order_claimed": False}, arrays


def build_report():
    pairs, arrays, outputs = [], {}, []
    source_path = ROOT/'certificates/finite-butson/bh7-42.json'
    catalogue = json.loads(source_path.read_text())
    if not catalogue['classification_complete'] or len(catalogue['classes']) != 9:
        raise ValueError("the complete order-seven input catalogue changed")
    for variant in (1, 2):
        u, z, blocks = endpoint_pair(variant)
        pairs.append({"variant": variant, "unitary_scalar_action": _json_matrix(u),
                      "frame": _json_matrix(z), "flat_seventh_root_blocks": blocks.tolist(),
                      "Phi_7_verified": True, "unitarity_verified": True,
                      "frame_Gram_equals_6_I": True,
                      "projective_flatness": flatness_pullback(variant),
                      "nonflat_restriction_tower": restriction_tower(variant)})
        for i, item in enumerate(catalogue['classes']):
            source = _square(item['representative'], 42)
            result = descend(source, variant=variant)
            key = f'bh42_pair{variant}_source{i}'
            arrays[key] = result
            outputs.append({"matrix_key": key, "source_class_index": i,
                            "source_class_digest": item['class_digest'],
                            "source_proof": verify_source(source, 42),
                            "verification": verify_sixth_gram(result),
                            "haagerup_histogram": haagerup_histogram(result)})
    # A known order-seven Laurent source is extracted from two stored points;
    # its entire Laurent Gram is verified before the endpoint is constructed.
    first = np.array(catalogue['classes'][0]['representative'], dtype=np.int64)
    mask = (np.array(catalogue['classes'][3]['representative'])-first) % 42
    if not np.all((mask == 0) | (mask == 1)) or np.any(first % 7):
        raise ValueError("the declared one-parameter source chart changed")
    source_family_proof = verify_laurent_family(first//7, np.where(mask, 0, -1))
    _, _, blocks = endpoint_pair(1)
    family42_base = ((first//7)[..., None, None]+blocks[mask]).transpose(0, 2, 1, 3).reshape(42, 42) % 6
    family42_labels = np.repeat(np.repeat(np.where(mask, 0, -1), 6, axis=0), 6, axis=1)
    family42_proof = verify_laurent_family(family42_base, family42_labels)
    arrays.update(family42_base=family42_base, family42_labels=family42_labels)
    phase12 = _specialize(family42_base, family42_labels, [1], 12)
    source84 = (14*(first//7)+19*mask) % 84  # t=zeta12*zeta7.
    if not np.array_equal(phase12, descend(source84, 84, 12)):
        raise ValueError("lower-field specialization does not commute with restriction")
    phase12_proof = verify_source(phase12, 12)
    arrays['bh42_q12'] = phase12

    source49, family294_base, family294_labels = dita_family()
    family294_proof = verify_laurent_family(family294_base, family294_labels)
    arrays.update(family294_base=family294_base, family294_labels=family294_labels)
    parameters = {"zero": np.zeros((7, 7), dtype=np.int64)}
    value = np.zeros((7, 7), dtype=np.int64)
    value[1, 1] = 1
    parameters['one'] = value
    value = np.zeros((7, 7), dtype=np.int64)
    value[1, 1:3] = 1
    parameters['two'] = value
    value = np.eye(7, dtype=np.int64)
    value[0, 0] = 0
    parameters['diagonal'] = value
    value = np.zeros((7, 7), dtype=np.int64)
    value[1:, 1:] = [[0, 1, 3, 2, 5, 4], [1, 2, 0, 5, 3, 1], [4, 0, 2, 1, 4, 5],
                     [3, 5, 1, 0, 2, 3], [5, 2, 4, 3, 1, 0], [2, 4, 5, 1, 0, 2]]
    parameters['mixed'] = value
    large_outputs = []
    a, b = np.repeat(np.arange(7), 7), np.tile(np.arange(7), 7)
    for name, parameter in parameters.items():
        source = (6*source49+7*parameter[b[:, None], a]) % 42
        result = descend(source)
        if not np.array_equal(result, _specialize(family294_base, family294_labels, parameter[1:, 1:].ravel())):
            raise ValueError("the 36-parameter family disagrees with its actual source reduction")
        key = f'bh294_{name}'
        arrays[key] = result
        large_outputs.append({"matrix_key": key, "parameters": parameter.tolist(),
                              "source_proof": verify_source(source, 42),
                              "verification": verify_sixth_gram(result),
                              "haagerup_histogram": haagerup_histogram(result)})
    lower42 = len({tuple(item['haagerup_histogram']) for item in outputs})
    lower294 = len({tuple(item['haagerup_histogram']) for item in large_outputs})
    real_report, real_arrays = real_examples()
    arrays.update(real_arrays)
    return {"schema": "projective-flat-descent-v1", "published_pair_source": PAPER,
            "published_pair_location": "Theorem 3.1, arXiv v1, 20 October 2026",
            "general_edge": "BH(n,q) -> BH(6n,s) for q dividing 7s, 6 dividing s, and 7 not dividing s",
            "input_catalogue": str(source_path.relative_to(ROOT)),
            "input_catalogue_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
            "endpoint_pairs": pairs, "order42_examples": outputs,
            "order42_monomial_class_lower_bound": lower42,
            "order42_family": {"source": source_family_proof, "endpoint": family42_proof,
                               "phase12_specialization": phase12_proof},
            "order294_family": family294_proof, "order294_examples": large_outputs,
            "order294_monomial_class_lower_bound": lower294,
            "real768": real_report,
            "inequivalence_proof": "different exact ordered Haagerup multiplicity histograms",
            "equal_fingerprints_do_not_assert_equivalence": True,
            "array_definitions": {key: {"shape": list(value.shape), "sha256": exponent_hash(value)}
                                  for key, value in sorted(arrays.items())},
            "scope": {"new_to_this_repository": True, "published_morphism_reproduced": True,
                      "previously_unknown_order_claimed": False, "publication_priority_established": False,
                      "complete_output_classification_claimed": False, "real_sign_endpoint_claimed": True}}, arrays


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument('--output', type=Path)
    action.add_argument('--check', type=Path)
    parser.add_argument('--matrices', type=Path, default=Path('artifacts/projective-flat-descent/examples.npz'))
    args = parser.parse_args()
    report, arrays = build_report()
    if args.output:
        args.matrices.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(args.matrices, **arrays)
        matrix_path = args.matrices
    else:
        old = json.loads(args.check.read_text())
        matrix_path = ROOT/old['literal_matrix_bundle']['path']
        with np.load(matrix_path, allow_pickle=False) as stored:
            if set(stored.files) != set(arrays) or any(not np.array_equal(stored[key], value) for key, value in arrays.items()):
                raise ValueError("the literal matrix bundle differs from the certified construction")
    matrix_path = matrix_path.resolve()
    report['literal_matrix_bundle'] = {"path": str(matrix_path.relative_to(ROOT)) if matrix_path.is_relative_to(ROOT) else str(matrix_path),
                                       "file_sha256": hashlib.sha256(matrix_path.read_bytes()).hexdigest()}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2)+'\n')
    elif json.loads(args.check.read_text()) != report:
        raise ValueError("the flatness/descent certificate failed replay")
    print(json.dumps({"orders": [42, 294, 768], "order42_class_lower_bound": report['order42_monomial_class_lower_bound'],
                      "order294_class_lower_bound": report['order294_monomial_class_lower_bound'],
                      "pairwise_inequivalent_real_H768": report['real768']['pairwise_monomially_inequivalent_examples'],
                      "real_binary_code_dimensions": report['real768']['binary_code_dimensions'],
                      "unit_phase_parameters_at_order294": report['order294_family']['independent_unit_phase_parameters'],
                      "final_flatness_certified_symbolically": True, "new_order_claimed": False}))


if __name__ == '__main__':
    main()

"""Alphabet-preserving quadratic descent of restricted fifth-root sources.

For d=4,8, map restricted BH(n,10d) to restricted BH(2n,5d).
This specializes a known Butson morphism; the additional audit preserves
the fifth-root alphabet needed by the terminal real lift. No new order
or missing large source is claimed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from functools import cache
from numbers import Integral
from pathlib import Path

import numpy as np
from sympy import totient

from experiments.scalar_lift.fifth_source_assembly import (
    complete_tensor_source,
    exact_column_gram,
    root_vectors,
)


def parameters(upper_phase):
    if type(upper_phase) is not int or upper_phase not in (40, 80):
        raise ValueError("the certified descents are phase40 to20 and80 to40")
    lower_phase = upper_phase // 2
    dyadic = lower_phase // 5
    return lower_phase, dyadic


def crt_coordinates(exponent, upper_phase):
    """zeta_(10d)^e = zeta5^c * zeta_(2d)^b, with canonical c,b."""
    _, d = parameters(upper_phase)
    if not isinstance(exponent, Integral) or isinstance(exponent, bool):
        raise ValueError("an exact integer exponent is required")
    exponent = int(exponent) % upper_phase
    color = pow(2 * d, -1, 5) * exponent % 5
    dyadic_exponent = pow(5, -1, 2 * d) * exponent % (2 * d)
    if (2 * d * color + 5 * dyadic_exponent) % upper_phase != exponent:
        raise ValueError("the CRT source reconstruction failed")
    return color, dyadic_exponent


def monomial_action(exponent, upper_phase):
    """Return rho(h) by row permutation and lower-root exponent per row."""
    lower, d = parameters(upper_phase)
    color, dyadic = crt_coordinates(exponent, upper_phase)
    base = (d * color + 5 * (dyadic // 2)) % lower
    if dyadic % 2:
        return (1, 0), ((base + 5) % lower, base)
    return (0, 1), (base, base)


def multiply_actions(first, second, lower_phase):
    p, a = first
    q, b = second
    return tuple(q[p[i]] for i in range(2)), tuple((a[i] + b[p[i]]) % lower_phase for i in range(2))


def adjoint_action(action, lower_phase):
    permutation, phases = action
    inverse = tuple(permutation.index(i) for i in range(2))
    return inverse, tuple(-phases[inverse[i]] % lower_phase for i in range(2))


@cache
def block_exponents(upper_phase):
    """The dense block rho(h)H2; every entry is a lower-field root."""
    lower, _ = parameters(upper_phase)
    h2 = np.array([[0, 0], [0, lower // 2]], dtype=np.int64)
    blocks = []
    for exponent in range(upper_phase):
        permutation, phases = monomial_action(exponent, upper_phase)
        blocks.append((h2[list(permutation)] + np.array(phases)[:, None]) % lower)
    return tuple(blocks)


def validate_source(source, phase):
    rows = [list(row) for row in source]
    if not rows or any(len(row) != len(rows) for row in rows):
        raise ValueError("a nonempty square source is required")
    if any(not isinstance(e, Integral) or isinstance(e, bool) for row in rows for e in row):
        raise ValueError("source exponents must be exact integers")
    matrix = np.array([[int(e) % phase for e in row] for row in rows], dtype=np.int64)
    if np.any(matrix % 5 == 0):
        raise ValueError("the source contains an excluded dyadic root")
    return matrix


def check_source(source, phase):
    matrix = validate_source(source, phase)
    n = len(matrix)
    for candidate in (matrix, matrix.T):
        gram = exact_column_gram(candidate, phase)
        expected = np.zeros_like(gram)
        expected[np.arange(n), np.arange(n), 0] = n
        if not np.array_equal(gram, expected):
            raise ValueError("the source fails its exact cyclotomic Gram")
    return {"order": n, "phase": phase, "all_entries_in_restricted_alphabet": True,
            "complete_row_and_column_Grams_verified": True,
            "source_exponents_sha256": hashlib.sha256(
                json.dumps(matrix.tolist(), separators=(",", ":")).encode()).hexdigest()}


def split_coefficients(source, upper_phase):
    """H=A+eta B; -1 marks a zero, all other integers encode lower roots."""
    lower, d = parameters(upper_phase)
    source = validate_source(source, upper_phase)
    a, b = np.full_like(source, -1), np.full_like(source, -1)
    for i in range(len(source)):
        for j in range(len(source)):
            color, dyadic = crt_coordinates(source[i, j], upper_phase)
            coefficient = (d * color + 5 * (dyadic // 2)) % lower
            (b if dyadic % 2 else a)[i, j] = coefficient
            if coefficient % 5 == 0:
                raise ValueError("a coefficient left the restricted lower alphabet")
            reconstruction = (2 * coefficient + 5 * (dyadic % 2)) % upper_phase
            if reconstruction != source[i, j]:
                raise ValueError("the disjoint coefficient decomposition failed")
    return a, b


def coefficient_product(a, b, phase, shift=0):
    """Exact zeta_phase^shift A B*, with -1 denoting a zero entry."""
    roots = root_vectors(phase)
    n = len(a)
    answer = np.zeros((n, n, len(roots[0])), dtype=np.int64)
    for i in range(n):
        for j in range(n):
            occupied = (a[i] >= 0) & (b[j] >= 0)
            if np.any(occupied):
                answer[i, j] = roots[(a[i, occupied] - b[j, occupied] + shift) % phase].sum(axis=0)
    return answer


def check_relative_equations(source, upper_phase):
    lower, _ = parameters(upper_phase)
    a, b = split_coefficients(source, upper_phase)
    n = len(a)
    if not np.all((a >= 0) ^ (b >= 0)):
        raise ValueError("the two supports do not partition all source entries")
    diagonal = coefficient_product(a, a, lower) + coefficient_product(b, b, lower)
    expected = np.zeros_like(diagonal)
    expected[np.arange(n), np.arange(n), 0] = n
    cross = coefficient_product(b, a, lower) + coefficient_product(a, b, lower, shift=-5)
    if not np.array_equal(diagonal, expected) or np.any(cross):
        raise ValueError("the relative coefficient Gram equations failed")
    return {"order": n, "upper_phase": upper_phase, "lower_phase": lower,
            "eta_upper_exponent": 5, "xi_lower_exponent": 5,
            "disjoint_support_partition_verified": True,
            "A_nonzero_entries": int(np.count_nonzero(a >= 0)),
            "B_nonzero_entries": int(np.count_nonzero(b >= 0)),
            "both_channels_retain_the_same_restricted_alphabet": True,
            "AAstar_plus_BBstar_equals_nI": True,
            "BAstar_plus_xi_inverse_ABstar_equals_zero": True,
            "ABstar_is_nonzero": bool(np.any(coefficient_product(a, b, lower)))}


def descend_source(source, upper_phase):
    lower, _ = parameters(upper_phase)
    source = validate_source(source, upper_phase)
    check_source(source, upper_phase)
    blocks = block_exponents(upper_phase)
    result = np.block([[blocks[int(exponent)] for exponent in row] for row in source])
    check_source(result, lower)
    return result


def scalar_certificate(upper_phase):
    lower, d = parameters(upper_phase)
    actions = [monomial_action(e, upper_phase) for e in range(upper_phase)]
    blocks = block_exponents(upper_phase)
    roots = root_vectors(lower)
    if int(totient(upper_phase)) != 2 * int(totient(lower)):
        raise ValueError("the cyclotomic extension is not quadratic")
    r = actions[5]
    if r != ((1, 0), (5, 0)) or multiply_actions(r, r, lower) != ((0, 1), (5, 5)):
        raise ValueError("the dyadic scalar does not satisfy R^2=xi I")
    for e, action in enumerate(actions):
        color, dyadic = crt_coordinates(e, upper_phase)
        if adjoint_action(action, lower) != actions[-e % upper_phase]:
            raise ValueError("the scalar action does not preserve conjugation")
        if np.any((pow(d, -1, 5) * blocks[e]) % 5 != color):
            raise ValueError("a dense scalar block changed its fifth color")
        for f, other in enumerate(actions):
            if multiply_actions(action, other, lower) != actions[(e + f) % upper_phase]:
                raise ValueError("the scalar action is not multiplicative")
            permutation, phases = actions[(e - f) % upper_phase]
            expected = np.zeros((2, 2, len(roots[0])), dtype=np.int64)
            for i in range(2):
                expected[i, permutation[i]] = 2 * roots[phases[i]]
            gram = np.array([[roots[(blocks[e][i] - blocks[f][j]) % lower].sum(axis=0)
                              for j in range(2)] for i in range(2)])
            if not np.array_equal(gram, expected):
                raise ValueError("a dense scalar cross-Gram identity failed")
    return {"upper_phase": upper_phase, "lower_phase": lower,
            "upper_field_degree": int(totient(upper_phase)),
            "lower_field_degree": int(totient(lower)), "relative_degree": 2,
            "minimal_polynomial": "X^2-zeta_d", "R_squared_equals_xi_I": True,
            "star_compatible_unitary_action_verified": True,
            "CRT_reconstruction_verified_for_all_phases": True,
            "both_coefficient_channels_preserve_fifth_colors": True,
            "all_dense_block_cross_Grams_verified": upper_phase**2,
            "all_multiplication_identities_verified": upper_phase**2,
            "all_scalar_block_exponents": [block.tolist() for block in blocks],
            "admitted_upper_exponents": [e for e in range(upper_phase) if e % 5],
            "new_order_claimed": False}


def certificate():
    original = complete_tensor_source(80)
    middle = descend_source(original, 80)
    terminal = descend_source(middle, 40)
    stages = [check_source(source, phase) for source, phase in ((original, 80), (middle, 40), (terminal, 20))]
    if [(stage["order"], stage["phase"]) for stage in stages] != [(40, 80), (80, 40), (160, 20)]:
        raise ValueError("the iterated source dimensions changed")
    return {"schema": "hadamard.alphabet-preserving-fifth-source-descent.v1",
            "known_morphism_reference": "https://arxiv.org/html/1904.10771#S1.Thmtheorem1.2",
            "attribution": "Specialization of O Cathain-Swartz; no novelty claim for the morphism",
            "scalar_steps": [scalar_certificate(phase) for phase in (40, 80)],
            "complete_source_control": {
                "initial_source": "existing BH40 from the fifth and dyadic order8 controls",
                "stages": stages,
                "relative_coefficient_checks": [check_relative_equations(original, 80),
                                                check_relative_equations(middle, 40)],
                "source_exponents": original.tolist(),
                "descended_phase40_exponents": middle.tolist(),
                "descended_phase20_exponents": terminal.tolist(),
                "terminal_known_real_lift_width": 8,
                "implied_known_real_order": 1280,
                "final_real_matrix_Gram_not_repeated": True,
                "known_orders_only": True},
            "conditional_high_routes": [
                {"initial_source_order": n, "initial_phase": phase,
                 "terminal_source_order": n * (phase // 20), "terminal_phase": 20,
                 "real_order": 8 * n * (phase // 20), "initial_source_constructed": False}
                for n, phase in ((1499, 40), (1499, 80), (2999, 80))],
            "large_Butson_source_or_new_real_order_constructed": False,
            "scope": "Positive iterative descent preserving the admitted alphabet; missing high sources are not constructed"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", type=Path)
    args = parser.parse_args()
    report = certificate()
    if args.check and json.loads(args.check.read_text()) != report:
        raise ValueError("the alphabet-preserving descent certificate did not replay")
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"restricted_source_control_stages": [[s["order"], s["phase"]]
                                                            for s in report["complete_source_control"]["stages"]],
                      "implied_known_real_order": 1280, "new_real_order_constructed": False}))


if __name__ == "__main__":
    main()

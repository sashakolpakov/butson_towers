"""Complete coefficient audit of the order-716 constant tensor template.

There are 22 normalized coefficient matrices and 16 quaternion orientations.
Exact modular moment certificates exclude 324 cases; 28 commuting cases
survive these tests. Surviving a moment test is not a Hadamard construction.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations_with_replacement, permutations, product
from math import isqrt

import numpy as np
from sympy import Matrix

from experiments.scalar_lift.paley_frame_family import modular_inconsistency
from experiments.scalar_lift.paley_sextic_operator import quaternion_orientation_audit

PRIME = 179
MODULI = (3, 5, 7, 11, 13, 17, 19, 23, 1009)
SIGNS = np.array(list(product((-1, 1), repeat=4)), dtype=np.int64)
PAIRS = tuple((i, j) for i in range(4) for j in range(i, 4))


def coefficient_classification_716():
    norm, bound = 3 * PRIME, isqrt(3 * PRIME)
    first_rows = [
        (a, *r)
        for a in range(1, bound + 1, 2)
        for r in combinations_with_replacement(range(0, bound + 1, 2), 3)
        if a * a + sum(x * x for x in r) == norm
    ]
    matrices, trials = [], 0
    for a, u, v, w in first_rows:
        equations = Matrix([
            [u, 0, 0, v, w, 0, -a * u],
            [0, v, 0, u, 0, w, -a * v],
            [0, 0, w, 0, u, v, -a * w],
            [1, 1, 1, 0, 0, 0, -a],
        ])
        reduced, pivots = equations.rref()
        if len(pivots) != 4 or any(j >= 6 for j in pivots):
            raise ValueError("the normalized coefficient equations do not have rank four")
        free = [j for j in range(6) if j not in pivots]
        coefficients = [[Fraction(reduced[i, j]) for j in range(7)] for i in range(4)]
        ranges = [[x for x in range(-bound, bound + 1) if x % 2 == (j < 3)] for j in free]
        for values in product(*ranges):
            trials += 1
            entries = [None] * 6
            for j, value in zip(free, values, strict=True):
                entries[j] = value
            for i, j in enumerate(pivots):
                value = coefficients[i][6] - sum(coefficients[i][k] * entries[k] for k in free)
                if value.denominator != 1 or abs(value) > bound or value.numerator % 2 != (j < 3):
                    break
                entries[j] = int(value)
            else:
                b, c, d, x, y, z = entries
                matrix = [[a, u, v, w], [u, b, x, y], [v, x, c, z], [w, y, z, d]]
                if all(
                    sum(matrix[i][k] * matrix[k][j] for k in range(4))
                    == (norm if i == j else 0)
                    for i in range(4) for j in range(4)
                ):
                    matrices.append(matrix)
    if len(first_rows) != 21 or trials != 11109 or len(matrices) != 22:
        raise ValueError("the complete coefficient classification changed")
    return {
        "norm": norm,
        "entry_bound": bound,
        "normalized_first_rows": [list(row) for row in first_rows],
        "bounded_free_assignments": trials,
        "normalized_matrices": matrices,
        "scope": "All symmetric integral B=I mod2, B^2=537I, up to signed permutation conjugacy",
    }


def required_sums_716(b, d):
    result = []
    for large in range(4):
        for signs in SIGNS:
            k = signs.copy()
            k[large] *= 3
            numerator = b @ k
            if np.any(numerator % 3):
                continue
            sigma = numerator // 3
            if sigma[0] < 0 or np.any((sigma + d @ sigma) % 4 != 2):
                continue
            if int(sigma @ sigma) != 4 * PRIME:
                raise ValueError("an extremal direction has the wrong norm")
            result.append(sigma)
    return result


def moment_arrays_716(b, d, sigma):
    offset = b @ d @ sigma
    if np.any(offset % PRIME):
        return None, None
    offset //= PRIME
    columns = []
    for v in SIGNS:
        for y in SIGNS:
            z = b @ (2 * y - v) - offset
            if np.any(z % 6):
                continue
            z //= 3
            columns.append(
                [1] + v.tolist() + y.tolist() + z.tolist()
                + [int(v[i] * z[j] + v[j] * z[i]) for i, j in PAIRS]
                + [int(z[i] * z[j] - PRIME * v[i] * v[j]) for i, j in PAIRS]
            )
    matrix = np.array(columns, dtype=np.int64).T if columns else np.zeros((33, 0), dtype=np.int64)
    rhs = np.array(
        [PRIME] + sigma.tolist() + ((sigma + d @ sigma) // 2).tolist()
        + [0] * 14 + [-int(sigma[i] * sigma[j]) for i, j in PAIRS], dtype=np.int64
    )
    return matrix, rhs


def first_modular_obstruction(matrix, rhs):
    for modulus in MODULI:
        rank, witness = modular_inconsistency(matrix, rhs, modulus)
        if witness is not None:
            return {"modulus": modulus, "coefficient_rank": rank, **witness}
    return None


def commuting_factor(b, d):
    if not np.array_equal(b @ d, d @ b):
        raise ValueError("the surviving coefficient matrices do not commute")
    numerator = -d @ b
    if np.any(numerator % 3):
        raise ValueError("the surviving quaternion factor is not integral")
    e = numerator // 3
    unit = np.eye(4, dtype=np.int64)
    if not np.array_equal(e.T, -e) or not np.array_equal(e @ e, -PRIME * unit):
        raise ValueError("the quaternion factor has the wrong skew norm")
    if not np.array_equal(d @ e, b) or not np.array_equal(d @ e, e @ d):
        raise ValueError("the commuting factorization B=D E failed")
    return e


def integer_moment_controls(first_survivor):
    """Actual nonnegative integer counts for the first case's four directions.

    These meet all local count moments but do not arrange the states into
    Paley-compatible columns. They are explicitly weaker than a frame.
    """
    sparse_counts = (
        {2: 22, 4: 24, 8: 17, 14: 26, 17: 19, 23: 28, 27: 21, 29: 22},
        {2: 17, 4: 19, 8: 22, 14: 24, 17: 21, 23: 22, 27: 26, 29: 28},
        {2: 26, 4: 17, 8: 21, 14: 22, 17: 22, 23: 24, 27: 28, 29: 19},
        {0: 17, 2: 4, 4: 5, 6: 17, 8: 2, 10: 17, 13: 17, 17: 11,
         18: 17, 21: 17, 23: 9, 25: 17, 27: 5, 29: 7, 31: 17},
    )
    b, d = np.array(first_survivor["B"]), np.array(first_survivor["D"])
    records = []
    for sigma, counts in zip(first_survivor["remaining_directions"], sparse_counts, strict=True):
        matrix, rhs = moment_arrays_716(b, d, np.array(sigma))
        vector = np.array([counts.get(j, 0) for j in range(matrix.shape[1])], dtype=np.int64)
        if np.any(vector < 0) or not np.array_equal(matrix @ vector, rhs):
            raise ValueError("the integer local-moment witness failed")
        records.append({"sigma": sigma, "state_counts": vector.tolist(), "exact_moments_verified": True})
    return records


def survivor_normal_forms(survivors):
    """Signed-permutation orbits, then Paley frame-existence equivalences.

    Global B sign is induced by j->2j on F_179, since two is a nonsquare.
    Simultaneous B,D sign reverses C and replaces R by its transpose.
    These operations permit independent global signs for existence tests.
    """
    if pow(2, (PRIME - 1) // 2, PRIME) != PRIME - 1:
        raise ValueError("doubling does not negate the Paley core")
    signed_groups, broad_groups = {}, {}
    all_permutations = tuple(permutations(range(4)))
    all_signs = tuple(product((-1, 1), repeat=4))
    for index, case in enumerate(survivors):
        b, d = np.array(case["B"]), np.array(case["D"])
        narrow_keys, broad_keys = [], []
        for permutation in all_permutations:
            bp, dp = b[np.ix_(permutation, permutation)], d[np.ix_(permutation, permutation)]
            for signs in all_signs:
                factors = np.outer(signs, signs)
                bb, dd = bp * factors, dp * factors
                narrow_keys.append(tuple(bb.flat) + tuple(dd.flat))
                for b_sign, d_sign in product((-1, 1), repeat=2):
                    broad_keys.append(tuple((b_sign * bb).flat) + tuple((d_sign * dd).flat))
        signed_groups.setdefault(min(narrow_keys), []).append(index)
        broad_groups.setdefault(min(broad_keys), []).append(index)
    classes = []
    for members in broad_groups.values():
        representative = members[0]
        target_b = np.array(survivors[representative]["B"])
        target_d = np.array(survivors[representative]["D"])
        witnesses = []
        for index in members:
            b, d = np.array(survivors[index]["B"]), np.array(survivors[index]["D"])
            found = None
            for permutation in all_permutations:
                bp, dp = b[np.ix_(permutation, permutation)], d[np.ix_(permutation, permutation)]
                for signs in all_signs:
                    factors = np.outer(signs, signs)
                    for b_sign, d_sign in product((-1, 1), repeat=2):
                        if np.array_equal(b_sign * bp * factors, target_b) and np.array_equal(
                            d_sign * dp * factors, target_d
                        ):
                            found = {
                                "survivor_index": index, "permutation": list(permutation),
                                "coordinate_signs": list(signs),
                                "global_B_sign": b_sign, "global_D_sign": d_sign,
                            }
                            break
                    if found is not None:
                        break
                if found is not None:
                    break
            if found is None:
                raise ValueError("a claimed frame-existence equivalence has no exact witness")
            witnesses.append(found)
        e = survivors[representative]["integral_E_with_B_equals_D_E"]
        classes.append({
            "representative_survivor_index": representative,
            "absolute_quaternion_coefficients": sorted(abs(int(x)) for x in e[0][1:]),
            "members_and_equivalence_witnesses": witnesses,
        })
    if len(signed_groups) != 10 or len(classes) != 3:
        raise ValueError("the surviving normal-form counts changed")
    return {
        "simultaneous_signed_permutation_orbit_count": len(signed_groups),
        "simultaneous_signed_permutation_orbits": list(signed_groups.values()),
        "Paley_frame_existence_class_count": len(classes),
        "additional_operations": "global B sign via j->2j; global B,D sign via R->R^T",
        "three_class_statement_requires_Paley_core_for_B_sign_operation": True,
        "classes": classes,
    }


def build_order_716_family_report():
    classification = coefficient_classification_716()
    ds = [np.array(row["D"], dtype=np.int64) for row in quaternion_orientation_audit()]
    cases, survivors = [], []
    for b_index, coefficient in enumerate(classification["normalized_matrices"]):
        b = np.array(coefficient, dtype=np.int64)
        for d_index, d in enumerate(ds):
            remaining, local_records = [], []
            for sigma in required_sums_716(b, d):
                matrix, rhs = moment_arrays_716(b, d, sigma)
                if matrix is None:
                    local_records.append({"sigma": sigma.tolist(), "offset_nonintegral": True,
                                          "offset_residues": (b @ d @ sigma % PRIME).tolist()})
                    continue
                witness = first_modular_obstruction(matrix, rhs)
                local_records.append({"sigma": sigma.tolist(), "state_count": matrix.shape[1],
                                      "modular_obstruction": witness})
                if witness is None:
                    remaining.append(sigma)
            matrix = np.array([[int(s[i] * s[j]) for i, j in PAIRS] for s in remaining],
                              dtype=np.int64).T if remaining else np.zeros((10, 0), dtype=np.int64)
            target = np.array([4 * PRIME**2 if i == j else 0 for i, j in PAIRS], dtype=np.int64)
            witness = first_modular_obstruction(matrix, target)
            case = {
                "B_index": b_index, "D_index": d_index,
                "B": b.tolist(), "D": d.tolist(), "local_directions": local_records,
                "remaining_directions": [s.tolist() for s in remaining],
                "global_frame_moment_obstruction": witness,
                "full_frame_excluded_by_these_checks": witness is not None,
            }
            if witness is None:
                e = commuting_factor(b, d)
                case["integral_E_with_B_equals_D_E"] = e.tolist()
                case["E_first_row_norm_decomposition"] = [int(x * x) for x in e[0, 1:]]
                survivors.append(case)
            elif np.array_equal(b @ d, d @ b):
                raise ValueError("a commuting case was unexpectedly excluded")
            cases.append(case)
    if len(cases) != 352 or len(survivors) != 28:
        raise ValueError("the complete order-716 audit changed")
    return {
        "scope": "Complete order-716 constant tensor coefficient audit; survivors are not constructed frames",
        "prime": PRIME, "order": 4 * PRIME,
        "coefficient_classification": classification,
        "D_count": len(ds), "case_count": len(cases),
        "excluded_case_count": len(cases) - len(survivors), "surviving_case_count": len(survivors),
        "certificate_moduli": MODULI,
        "survivors_are_exactly_the_commuting_B_D_cases": True,
        "all_survivors_have_integral_skew_E_squared_minus_179I": True,
        "symmetric_regular_involution": "O=(E tensor Q_179+I_4 tensor J_179)/179",
        "factorization": "B=D E, C=179 (D tensor I_179) O, R=(I+(D tensor I_179)O)/2",
        "involution_positive_eigenspace_dimension": 360,
        "involution_negative_eigenspace_dimension": 356,
        "first_survivor_integer_local_moment_controls": integer_moment_controls(survivors[0]),
        "survivor_normal_forms": survivor_normal_forms(survivors),
        "cases": cases,
        "flat_frame_constructed": False,
        "new_real_hadamard_order_claimed": False,
    }

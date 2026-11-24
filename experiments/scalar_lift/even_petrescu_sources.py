"""Known Paley inputs and exact prerequisites for even Petrescu sources.

The normal corner D and the admissible X/Y split remain missing. These
checks do not construct either large Butson source or a new real order.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from functools import cache
from math import isqrt
from pathlib import Path

import numpy as np
from sympy import isprime

from experiments.scalar_lift.bordered_seed_assembly import check_sign_hadamard
from experiments.scalar_lift.commuting_parity import paley_core
from experiments.scalar_lift.dyadic_partial_lift import sylvester
from experiments.scalar_lift.fifth_source_assembly import exact_column_gram, root_vectors
from hadamard.cgw import is_cgw_exponent_matrix
from hadamard.coefficients import CyclotomicField


@cache
def paley_500():
    if not isprime(499) or 499 % 4 != 3:
        raise ValueError("the Paley parameter is not a prime congruent to3 modulo4")
    q = paley_core(499)
    h = np.ones((500, 500), dtype=np.int64)
    h[1:, 1:] = q - np.eye(499, dtype=np.int64)
    if np.any(h[0] != 1) or np.any(h[:, 0] != 1):
        raise ValueError("the Paley seed is not normalized")
    return h


def input_hadamard(order):
    if order not in (1000, 2000):
        raise ValueError("the known input orders are1000 and2000")
    return np.kron(sylvester(order // 500), paley_500())


def boundary_input(order):
    """Return all noninitial rows of the explicitly normalized input."""
    return input_hadamard(order)[1:]


def input_report(order):
    h = input_hadamard(order)
    h_check = check_sign_hadamard(h)
    if np.any(h[0] != 1) or np.any(h[:, 0] != 1):
        raise ValueError("the tensor input is not normalized")
    t = h[1:]
    # The row Gram is the principal submatrix of the fully verified H H^T.
    # Check every off-diagonal entry of T^T T independently by Hamming distance.
    words = [sum(1 << j for j, value in enumerate(column) if value == 1)
             for column in t.T]
    for i, word in enumerate(words):
        if any(2 * (word ^ words[j]).bit_count() != order for j in range(i)):
            raise ValueError("the exact column Gram T^T T=tI-J failed")
    if t.shape != (order - 1, order) or not np.all(np.abs(t) == 1):
        raise ValueError("the rectangular sign boundary has wrong shape or entries")
    return {
        "order": order, "construction": f"H{order // 500} tensor Paley H500",
        "complete_sign_Gram_check": h_check,
        "first_row_and_column_all_one": True,
        "T_shape": list(t.shape), "T_row_Gram": f"{order} I_{order - 1}",
        "T_column_Gram": f"{order} I_{order} - J_{order}",
        "T_row_Gram_inherited_from_full_verified_input": True,
        "T_column_Gram_all_off_diagonal_pairs_checked": order * (order - 1) // 2,
        "T_column_Gram_diagonal": order - 1,
        "T_binary_sha256": hashlib.sha256(((t + 1) // 2).astype(np.uint8).tobytes()).hexdigest(),
    }


def polynomial_element(field, coefficients, subfield_order):
    if field.conductor % subfield_order:
        raise ValueError("the subfield conductor does not divide the chosen phase")
    z = field.zeta() ** (field.conductor // subfield_order)
    return sum((coefficient * z ** j for j, coefficient in enumerate(coefficients)), field.zero())


def row_sum_witness(order, phase):
    if (order, phase) not in ((1000, 24), (1000, 48), (2000, 48)):
        raise ValueError("this input/phase route is not certified")
    if order == 1000:
        subfield, coefficients = 8, [21, 0, 21, 46]
    else:
        subfield, coefficients = 16, [7, -40, 50, 40, 7, 10, 0, 10]
    field = CyclotomicField(phase)
    lam = polynomial_element(field, coefficients, subfield)
    source_order = 3 * order - 2
    if lam * lam.conjugate() != field.element([source_order]):
        raise ValueError("the exact corner row-sum norm failed")
    counts = [0] * subfield
    for j, coefficient in enumerate(coefficients):
        counts[(j + (subfield // 2 if coefficient < 0 else 0)) % subfield] += abs(coefficient)
    original_letters = sum(counts)
    padding, remainder = divmod(order - original_letters, 2)
    if padding < 0 or remainder:
        raise ValueError("the corner row sum cannot be padded by cancelling pairs")
    counts[0] += padding
    counts[subfield // 2] += padding
    letter_sum = polynomial_element(field, counts, subfield)
    if sum(counts) != order or letter_sum != lam:
        raise ValueError("the exact root multiplicities do not reproduce the row sum")
    return {"subfield_conductor": subfield, "coefficients_ascending": coefficients,
            "exact_relative_norm": source_order, "unpadded_root_count": original_letters,
            "added_cancelling_pairs": padding, "padded_root_count": sum(counts),
            "subfield_root_exponent_counts": counts,
            "embedded_phase_root_exponents": [phase // subfield * j for j in range(subfield)],
            "all_letters_in_dyadic_part_of_restricted_alphabet": True,
            "full_row_sum_and_norm_verified": True}, lam


def determinant_witness(order, phase, lam):
    field = CyclotomicField(phase)
    z = field.zeta() ** (phase // 8)
    i, sqrt_minus_two = z ** 2, z + z ** 3
    if order == 1000:
        gamma = (1 + i) * (7 + 15 * sqrt_minus_two)
        formula = "(1+i)*(7+15*sqrt(-2))"
    else:
        gamma = (1 + i) * (1 + sqrt_minus_two) ** 3 * (6 + i)
        formula = "(1+i)*(1+sqrt(-2))^3*(6+i)"
    if (gamma * gamma.conjugate() != field.element([order - 2])
            or lam * lam.conjugate() != field.element([3 * order - 2])):
        raise ValueError("a factor of the corner determinant norm witness failed")
    if any(value.denominator != 1 for value in gamma.coefficients):
        raise ValueError("the corner determinant witness factor is not integral")
    return {"required_norm": f"{order - 2}^{order - 1}*{3 * order - 2}",
            "gamma_formula": formula,
            "gamma_coefficients_in_chosen_phase": [int(value) for value in gamma.coefficients],
            "gamma_norm": order - 2, "lambda_norm": 3 * order - 2,
            "factored_determinant_witness": f"gamma^{order - 1}*lambda",
            "factor_norms_verified_exactly": True,
            "scalar_integral_determinant_norm_obstruction_absent": True,
            "corner_matrix_constructed": False}


def two_phase_corner_audit(order):
    candidates = []
    for k in range(order + 1):
        if k * (k - 1) % (order - 1):
            continue
        lam = k * (k - 1) // (order - 1)
        weight = k - lam
        required = None if not weight else Fraction(order - 2, weight)
        possible = weight > 0 and 4 * weight >= order - 2
        candidates.append({"row_weight_k": k, "intersection_lambda": lam,
                           "k_minus_lambda": weight,
                           "required_squared_chord": None if required is None else
                           {"numerator": required.numerator, "denominator": required.denominator},
                           "compatible_with_squared_chord_at_most_four": possible})
    if any(record["compatible_with_squared_chord_at_most_four"] for record in candidates):
        raise ValueError("the claimed complete two-phase exclusion failed")
    return {"corner_order": order, "all_integer_k_values_checked": order + 1,
            "integer_symmetric_design_parameter_candidates": candidates,
            "maximum_k_minus_lambda": max(record["k_minus_lambda"] for record in candidates),
            "all_two_global_unit_phase_corners_excluded": True,
            "phases_need_not_be_roots_of_unity_for_this_exclusion": True,
            "corners_with_three_or_more_phases_excluded": False}


def exact_integer_ratio(numerator, denominator):
    """Test an integer ratio by every cyclotomic power-basis coefficient."""
    pivot = next((j for j, value in enumerate(denominator) if value), None)
    if pivot is None:
        raise ValueError("a phase chord denominator is zero")
    if numerator[pivot] % denominator[pivot]:
        return None
    quotient = int(numerator[pivot] // denominator[pivot])
    return quotient if np.array_equal(numerator, quotient * denominator) else None


def symmetric_three_phase_corner_audit(order, phase):
    if (order, phase) not in ((1000, 24), (1000, 48), (2000, 48)):
        raise ValueError("the exact symmetric-corner audits cover three specified routes")
    roots = root_vectors(phase)
    integral_cases = []
    valid_parameters = []
    for u in range(1, phase):
        for v in range(1, phase):
            if u == v:
                continue
            chord = 2 * roots[0] - roots[(u - v) % phase] - roots[(v - u) % phase]
            cross = (2 * roots[0] + roots[u] + roots[-u % phase]
                     - roots[v] - roots[-v % phase]
                     - roots[(u - v) % phase] - roots[(v - u) % phase])
            alpha = exact_integer_ratio((order - 4) * roots[0] + roots[v] + roots[-v % phase], chord)
            beta = exact_integer_ratio(-cross, chord)
            if alpha is None or beta is None:
                continue
            discriminant = (order + beta) ** 2 - 4 * alpha * (order - 1)
            square = discriminant >= 0 and isqrt(discriminant) ** 2 == discriminant
            candidates = []
            if square:
                for numerator in sorted({order + beta - isqrt(discriminant),
                                         order + beta + isqrt(discriminant)}):
                    if numerator % 2:
                        continue
                    k = numerator // 2
                    mu, lam = k - alpha, k - alpha + beta
                    valid = 0 < k < order - 1 and 0 <= lam < k and 0 <= mu <= k
                    candidate = {"k": k, "lambda": lam, "mu": mu, "valid": valid}
                    candidates.append(candidate)
                    if valid:
                        valid_parameters.append({"u_exponent": u, "v_exponent": v, **candidate})
            integral_cases.append({"u_exponent": u, "v_exponent": v,
                                   "alpha": alpha, "beta": beta,
                                   "k_quadratic_discriminant": discriminant,
                                   "discriminant_is_a_nonnegative_integer_square": square,
                                   "integer_k_candidates": candidates})
    expected_count = 16 if order == 2000 else 14
    if len(integral_cases) != expected_count or valid_parameters:
        raise ValueError("the complete symmetric three-phase corner census changed")
    return {
        "corner_order": order, "phase": phase,
        "normalized_diagonal_phase": 1,
        "all_ordered_distinct_nonidentity_off_diagonal_phase_pairs_checked": (phase - 1) * (phase - 2),
        "integral_alpha_beta_pair_count": len(integral_cases),
        "complete_integral_cases": integral_cases,
        "valid_strongly_regular_parameter_sets": valid_parameters,
        "symmetric_constant_diagonal_two_off_diagonal_phase_corners_excluded": True,
        "all_three_phase_corners_excluded": False,
        "phases_outside_the_certified_24_48_alphabets_excluded": False,
        "all_full_phase_roots_checked_before_any_lift_alphabet_restriction": True,
    }


def corner_formula_controls():
    # A regular graph with both edges and nonedges checks the complete
    # expansion even when its resulting corner does not satisfy Petrescu's Gram.
    order, phase, u, v, k = 5, 12, 2, 5, 2
    adjacency = np.array([[int((i - j) % 5 in (1, 4)) for j in range(5)]
                          for i in range(5)], dtype=np.int64)
    roots = root_vectors(phase)
    exponents = np.where(adjacency, u, v)
    np.fill_diagonal(exponents, 0)
    actual = np.array([[roots[(exponents[i] - exponents[j]) % phase].sum(axis=0)
                        for j in range(order)] for i in range(order)])
    chord = 2 * roots[0] - roots[(u - v) % phase] - roots[(v - u) % phase]
    diagonal = 2 * roots[0] - roots[v] - roots[-v % phase]
    cross = (2 * roots[0] + roots[u] + roots[-u % phase] - roots[v] - roots[-v % phase]
             - roots[(u - v) % phase] - roots[(v - u) % phase])
    constant = (order - 2) * roots[0] + roots[v] + roots[-v % phase] - k * chord
    predicted = (np.eye(order, dtype=np.int64)[:, :, None] * diagonal
                 + (adjacency @ adjacency)[:, :, None] * chord
                 + adjacency[:, :, None] * cross
                 + np.ones((order, order, 1), dtype=np.int64) * constant)
    target = ((order - 2) * np.eye(order, dtype=np.int64)
              + 2 * np.ones((order, order), dtype=np.int64))[:, :, None] * roots[0]
    if not np.array_equal(actual, predicted) or np.array_equal(actual, target):
        raise ValueError("the independent graph expansion control failed")

    # The published order19 example supplies a positive nontrivial corner.
    seed = json.loads(Path(__file__).with_name("bh19_6_seed.json").read_text())
    full = np.array([[int(value) for value in row] for row in seed["exponent_rows"]], dtype=np.int64)
    if not is_cgw_exponent_matrix(full.tolist(), 6, 19):
        raise ValueError("the published full BH19 control failed")
    corner = full[12:, 12:]
    expected = np.zeros((7, 7, len(root_vectors(6)[0])), dtype=np.int64)
    expected[:, :, 0] = 5 * np.eye(7, dtype=np.int64) + 2 * np.ones((7, 7), dtype=np.int64)
    if (not np.array_equal(exact_column_gram(corner, 6), expected)
            or not np.array_equal(exact_column_gram(corner.T, 6), expected)):
        raise ValueError("the published Petrescu corner lost its two exact Grams")
    return {
        "regular_graph_expansion": {
            "graph": "cycle on5 vertices", "phase": phase,
            "diagonal_exponent": 0, "edge_exponent": u, "nonedge_exponent": v,
            "complete_power_basis_Gram_expansion_verified": True,
            "nonzero_Petrescu_Gram_defect_verified": True},
        "published_Petrescu_corner": {
            "source_url": seed["source_url"], "source_location": seed["location"],
            "full_BH19_Gram_verified": True, "corner_order": 7, "phase": 6,
            "corner_exponents": corner.tolist(), "both_Grams_equal_5I_plus_2J": True,
            "constant_diagonal": bool(np.all(np.diag(corner) == corner[0, 0])),
            "symmetric": bool(np.array_equal(corner, corner.T)),
            "distinct_phase_count": len(set(int(value) for value in corner.flat)),
            "restricted_sign_lift_compatibility_asserted": False},
    }


def small_array_control():
    x, y = np.array([[1]], dtype=np.int64), np.array([[-1]], dtype=np.int64)
    t, d = np.array([[1, -1]], dtype=np.int64), np.ones((2, 2), dtype=np.int64)
    h = np.block([[x, y, t], [y, x, t], [t.T, t.T, d]])
    return {"s": 1, "t": 2, "source_order": 4,
            "actual_complete_Petrescu_array": h.tolist(),
            "full_sign_Gram_check": check_sign_hadamard(h), "known_order_control_only": True}


def build_report():
    base_check = check_sign_hadamard(paley_500())
    inputs = {order: input_report(order) for order in (1000, 2000)}
    routes = []
    for order, phase in ((1000, 24), (1000, 48), (2000, 48)):
        row_sum, lam = row_sum_witness(order, phase)
        n, s = 3 * order - 2, order - 1
        routes.append({
            "input_order_t": order, "core_order_s": s, "source_order": n, "phase": phase,
            "conditional_real_target_order": n * (phase // 3),
            "real_lift_width": phase // 3,
            "raw_alphabet": "mu_(phase/3) union zeta3*mu_(phase/3)",
            "raw_exponent_classes_mod3": [0, (phase // 3) % 3],
            "alignment_with_certified_sign_lift": "conjugate entire source" if phase == 24 else "direct",
            "corner_Gram_equation": f"D D*=D* D={order - 2} I_{order}+2 J_{order}",
            "corner_commutes_with_J": True,
            "E_formula": "E=X+Y=-(T D* T*)/t",
            "E_Gram_implied_by_corner_and_input": f"E E*={order - 2} I_{s}",
            "F_condition": f"F F*={n} I_{s}",
            "remaining_entry_condition": "X=(E+F)/2 and Y=(E-F)/2 lie in the admitted alphabet",
            "joint_incidence_required": "D and F must jointly satisfy the corner Gram and X/Y root conditions",
            "necessary_integrality_constraint_from_the_start": "T D* T* lies in t M_s(Z[zeta_q])",
            "input_order_divisible_by_5_cubed": order % 125 == 0,
            "corner_row_sum_witness": row_sum,
            "corner_determinant_norm_witness": determinant_witness(order, phase, lam),
            "corner_order_is_even_and_corner_is_not_Butson": True,
            "previous_odd_border_obstruction_applies": False,
            "reason_previous_odd_border_obstruction_does_not_apply": "the corner order is even",
            "corner_D_constructed": False, "X_Y_split_constructed": False,
            "large_source_constructed": False,
        })
    return {
        "schema": "even-petrescu-source-prerequisites-v2",
        "scope": "Known inputs and norms; exclusions of two-phase and symmetric constant-diagonal/two-off-phase corners",
        "primary_Petrescu_reference": "https://ajc.maths.uq.edu.au/pdf/55/ajc_v55_p031.pdf",
        "Paley_prime": 499, "Paley_base_complete_check": base_check,
        "known_input_certificates": [inputs[order] for order in (1000, 2000)],
        "routes": routes,
        "complete_two_phase_corner_audits": [two_phase_corner_audit(order) for order in (1000, 2000)],
        "symmetric_three_phase_corner_audits": [symmetric_three_phase_corner_audit(order, phase)
                                                for order, phase in ((1000, 24), (1000, 48), (2000, 48))],
        "independent_corner_formula_controls": corner_formula_controls(),
        "small_complete_array_control": small_array_control(),
        "native_solver_run": False, "new_real_Hadamard_order_constructed": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--output", type=Path)
    group.add_argument("--check", type=Path)
    args = parser.parse_args()
    report = build_report()
    if args.check:
        if json.loads(args.check.read_text()) != report:
            raise ValueError("the even Petrescu certificate differs from exact replay")
    else:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"known_input_orders": [1000, 2000], "conditional_source_orders": [2998, 5998],
                      "two_phase_corners_excluded": True, "large_source_constructed": False}))


if __name__ == "__main__":
    main()

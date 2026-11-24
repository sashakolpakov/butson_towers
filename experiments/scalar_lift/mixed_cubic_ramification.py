"""First ramification obstruction for the displayed dense mixed boundaries.

The source-order1499 and2999 F4 matrices are unitary, but their fixed
primitive boundary phases admit no root-valued lift modulo pi squared.
This excludes those residue/boundary pairings, not all Butson sources.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from experiments.scalar_lift.mixed_cubic_boundaries import mixed_frames, mixed_frames_2999

MULTIPLY = np.array([[0, 0, 0, 0], [0, 1, 2, 3], [0, 2, 3, 1], [0, 3, 1, 2]], dtype=np.int64)
CONJUGATE = np.array([0, 1, 3, 2], dtype=np.int64)


def product(left, right):
    left, right = np.asarray(left, dtype=np.int64), np.asarray(right, dtype=np.int64)
    output = np.zeros((left.shape[0], right.shape[1]), dtype=np.int64)
    for k in range(left.shape[1]):
        output ^= MULTIPLY[left[:, k, None], right[None, k, :]]
    return output


def derivative(h0, h1):
    return product(h1, CONJUGATE[h0].T) ^ product(h0, CONJUGATE[h1].T)


def core_derivative_formula(m, v):
    """Four binary coefficients of the two core Gram derivative blocks."""
    m, v = np.asarray(m, dtype=np.int64), np.asarray(v, dtype=np.int64)
    d = np.diag(v)
    x, y = m.sum(axis=1) % 2, v.sum(axis=1) % 2
    r = x ^ y
    diagonal_0 = r[:, None] ^ r[None, :] ^ v ^ d[None, :]
    diagonal_omega = v ^ v.T ^ d[:, None] ^ d[None, :]
    cross_0 = (x ^ d)[:, None] ^ y[None, :] ^ m.T
    cross_omega = m ^ m.T ^ (r ^ d)[:, None] ^ (r ^ d)[None, :]
    return diagonal_0, diagonal_omega, cross_0, cross_omega


def binary_rank(words):
    pivots = {}
    for original in words:
        word = original
        while word:
            pivot = word.bit_length() - 1
            if pivot in pivots:
                word ^= pivots[pivot]
            else:
                pivots[pivot] = word
                break
    return len(pivots)


def coefficient_oracle():
    """Independent F4 multiplication checks every basis of small linear maps.

    The arbitrary-size elimination is recorded in the proof. These
    checks catch orientation, conjugation and coefficient-basis mistakes.
    No small-size result is extrapolated as an existence assertion.
    """
    reports = []
    for s in (2, 4, 6):
        ones, eye = np.ones((s, s), dtype=np.int64), np.eye(s, dtype=np.int64)
        y0 = 2 * ones ^ 3 * eye
        constraint_rows = [0] * (2 * s * s)
        for index in range(2 * s * s):
            m, v = np.zeros((s, s), dtype=np.int64), np.zeros((s, s), dtype=np.int64)
            target = m if index < s * s else v
            target.flat[index % (s * s)] = 1
            y1 = 2 * v ^ 3 * np.diag(np.diag(v))
            direct_diagonal = (product(m, ones) ^ product(y1, CONJUGATE[y0].T)
                               ^ product(ones, m.T) ^ product(y0, CONJUGATE[y1].T))
            direct_cross = (product(m, CONJUGATE[y0].T) ^ product(y1, ones)
                            ^ product(ones, CONJUGATE[y1].T) ^ product(y0, m.T))
            a, b, c, d = core_derivative_formula(m, v)
            if not np.array_equal(direct_diagonal, a ^ (2 * b)) or not np.array_equal(direct_cross, c ^ (2 * d)):
                raise ValueError("the core coefficient formula disagrees with direct F4 multiplication")
            # Even the constant coefficients alone force the common
            # column-phase form. Verify their small-map ranks exactly.
            for j, value in enumerate(np.concatenate((a.ravel(), c.ravel()))):
                if value:
                    constraint_rows[j] |= 1 << index
        rank = binary_rank(constraint_rows)
        if rank != 2 * s * s - s:
            raise ValueError("the small core derivative kernel has the wrong dimension")
        for j in range(s):
            bits = np.zeros((s, s), dtype=np.int64)
            bits[:, j] = 1
            if any(np.any(block) for block in core_derivative_formula(bits, bits)):
                raise ValueError("a common column phase is absent from the core kernel")
        reports.append({"s": s, "all_linear_basis_vectors_checked": 2 * s * s,
                        "constant_coefficient_system_rank": rank,
                        "kernel_dimension": s, "kernel_is_common_column_phase_bits": True})

    # Audit the cross-block reduction in full small matrices, for every
    # column-phase vector, independently of the displayed formulas.
    cross_checks = 0
    for s in (2, 4):
        ones, eye = np.ones((s, s), dtype=np.int64), np.eye(s, dtype=np.int64)
        k = ones ^ eye
        active = np.array([1, 1] + [0] * (s - 2), dtype=np.int64)
        a = active[:, None] * np.array([[0, 1, 1]])
        u0, u1 = np.ones((s, 3), dtype=np.int64) ^ (3 * a), 2 * a
        d0 = np.ones((3, 3), dtype=np.int64) ^ (3 * np.eye(3, dtype=np.int64))
        y0 = 2 * ones ^ (3 * eye)
        h0 = np.block([[ones, y0, u0], [y0, ones, u0], [u0.T, u0.T, d0]])
        for mask in range(1 << s):
            d = np.array([(mask >> j) & 1 for j in range(s)], dtype=np.int64)
            m = np.ones((s, 1), dtype=np.int64) * d[None, :]
            y1 = MULTIPLY[y0, m]
            h1 = np.block([[m, y1, u1], [y1, m, u1], [u1.T, u1.T, np.zeros((3, 3), dtype=np.int64)]])
            direct = derivative(h0, h1)
            expected = (k @ np.diag(d) @ a % 2) ^ a ^ (3 * ((k @ d % 2)[:, None] * np.ones((1, 3), dtype=int)))
            if not np.array_equal(direct[:s, 2 * s:], expected):
                raise ValueError("the full boundary derivative disagrees with the reduced equation")
            if np.any(direct[:2 * s, :2 * s]):
                raise ValueError("the common-column core derivative unexpectedly has a residual")
            cross_checks += 1
    return {"exact_core_formula_controls": reports, "full_boundary_derivatives_checked": cross_checks}


def phase_bits(exponents, phase, cube_exponent):
    exponents = np.asarray(exponents, dtype=np.int64)
    colors = exponents % 3
    if np.any(colors == 2):
        raise ValueError("the entry is outside the permitted two-coset alphabet")
    difference = exponents - cube_exponent * colors
    if np.any(difference % 3):
        raise ValueError("the cubic/dyadic exponent did not split integrally")
    dyadic = (difference // 3) % (phase // 3)
    return colors, dyadic % 2


def target_obstruction(route):
    phase, s, t = route["phase"], route["s"], route["border_order"]
    a = np.asarray(route["cubic_support_A"], dtype=np.int64)
    cube = route["scalar_pair_action"]["cube_exponent"]
    u_colors, u_bits = phase_bits(route["U_exponents"], phase, cube)
    bottom_colors, bottom_bits = phase_bits(-np.asarray(route["V_exponents"]), phase, cube)
    c, d_bits = phase_bits(route["D_exponents"], phase, cube)
    if (not np.array_equal(u_colors, a) or not np.array_equal(u_bits, a)
            or not np.array_equal(bottom_colors, a) or not np.array_equal(bottom_bits, a)
            or np.any(d_bits)):
        raise ValueError("the actual boundary phase bits do not match the first-layer theorem")
    nonzero = np.flatnonzero(np.any(a, axis=1))
    if not len(nonzero):
        raise ValueError("this obstruction requires a nonzero primitive boundary bit")
    row = a[nonzero[0]]
    active = np.any(a, axis=1).astype(np.int64)
    if (s % 2 or t % 2 != 1 or not np.array_equal(a, active[:, None] * row[None, :])
            or row.sum() % 2 or active.sum() % 2
            or not np.array_equal(c @ c.T % 2, np.eye(t, dtype=np.int64))
            or not np.array_equal(a @ c.T % 2, a)):
        raise ValueError("the actual boundary support fails a required Gram or incidence identity")
    return {"source_order": 2 * s + t, "phase": phase, "dyadic_phase_order": phase // 3,
            "core_dimension": s, "border_order": t,
            "boundary_support_sha256": hashlib.sha256(a.astype(np.uint8).tobytes()).hexdigest(),
            "actual_upper_and_lower_boundary_first_bits_equal_A": True,
            "actual_border_first_bits_are_zero": True,
            "nonzero_A_entries": int(a.sum()),
            "A_rank_one_factor_weights": [int(active.sum()), int(row.sum())],
            "A_row_and_column_parities_even": True,
            "A_Atranspose_zero_and_A_Ctranspose_equals_A": True,
            "K_equals_J_plus_I_is_invertible_over_F2": True,
            "core_derivative_forces_bX_equals_bY_equals_ones_times_dtranspose": True,
            "remaining_cross_equation": "K diag(d) A + A + omega^2 (Kd) ones^T = 0",
            "independent_F2_coefficients_force_Kd_zero_then_d_zero_then_A_zero": True,
            "contradiction_coordinate": [int(nonzero[0]), int(np.flatnonzero(row)[0])],
            "no_root_phase_lift_modulo_pi_squared_for_this_fixed_boundary_and_core_support": True,
            "conditional_real_lift_order": route["conditional_real_lift_order"],
            "general_Butson_source_order_excluded": False}


def certificate():
    return {"schema": "mixed-cubic-first-ramification-v1",
            "coefficient_ring": "F4[pi]/(pi^2), pi=1-zeta_d, d in{16,32}",
            "entry_expansion": "omega^c zeta_d^b = omega^c(1+(b mod2)pi)",
            "conjugation": "omega->omega^2, pi->pi modulo pi^2",
            "known_dense_F4_residues_remain_unitary": True,
            "oracle": coefficient_oracle(),
            "targets": [target_obstruction(mixed_frames(96)),
                        target_obstruction(mixed_frames_2999(48)),
                        target_obstruction(mixed_frames_2999(96))],
            "row_Gram_alone_is_inconsistent": True,
            "column_Gram_follows_from_row_Gram_for_square_matrices_over_this_commutative_ring": True,
            "second_layer_search_needed": False, "native_probe_run": False,
            "larger_source_constructed": False,
            "scope": "The specified repeated X/Y dense cubic support with the saved primitive mixed boundaries only."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--output", type=Path)
    action.add_argument("--check", type=Path)
    args = parser.parse_args()
    report = certificate()
    if args.check:
        if json.loads(args.check.read_text()) != report:
            raise ValueError("the first-ramification certificate failed replay")
    else:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"excluded_fixed_residue_boundary_pairs": [[target["source_order"], target["phase"]]
                                                              for target in report["targets"]],
                      "obstruction_precision": "pi^2", "general_source_orders_excluded": False,
                      "native_probe_run": False}))


if __name__ == "__main__":
    main()

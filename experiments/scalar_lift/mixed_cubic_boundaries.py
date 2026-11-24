"""Actual mixed cubic partial frames at source orders1499 and2999.

Complete dense-core F4 residues are verified. The old monomial-core
residue remains reproducible as historical v1 data, but is now excluded
from characteristic-zero completion. No characteristic-zero X/Y core
is constructed. The adjoint-boundary obstruction does not apply to the
displayed transpose coupling.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from functools import cache
from itertools import product
from math import gcd
from pathlib import Path

import numpy as np

from experiments.scalar_lift.bordered_seed_assembly import (
    border_source,
    check_sign_hadamard,
    paley_748,
)
from experiments.scalar_lift.fifth_source_assembly import exact_column_gram, root_vectors
from hadamard.coefficients import CyclotomicField


def pair_action_audit(phase):
    cube = next(e for e in (phase // 3, 2 * phase // 3) if e % 3 == 1)
    eta = (cube + 3) % phase
    if gcd(eta, phase) != 1:
        raise ValueError("the pair parameter is not a primitive full-phase root")
    roots = root_vectors(phase)
    entries = [[[(0, 1), (eta, 1)], [(0, 1), (eta, -1)]],
               [[(0, 1), (eta, -1)], [(0, 1), (eta, 1)]]]
    gram = np.zeros((2, 2, roots.shape[1]), dtype=np.int64)
    for i, j, k in product(range(2), repeat=3):
        for e, sign in entries[i][k]:
            for f, other_sign in entries[j][k]:
                gram[i, j] += sign * other_sign * roots[(e - f) % phase]
    expected = np.zeros_like(gram)
    expected[0, 0, 0] = expected[1, 1, 0] = 4
    if not np.array_equal(gram, expected):
        raise ValueError("the exact2M Gram is not4I")
    for first, second in product((-1, 1), repeat=2):
        output = [sum(value * sum(sign * roots[e] for e, sign in entries[i][j])
                      for j, value in enumerate((first, second))) for i in range(2)]
        expected_root = roots[0 if first == second else eta]
        if any(not np.array_equal(output[i], 2 * value * expected_root)
               for i, value in enumerate((first, second))):
            raise ValueError("the pair action does not map sign entries to the required roots")
    return {"phase": phase, "cube_exponent": cube, "eta_exponent": eta,
            "twice_pair_action_Gram_equals_4I": True,
            "all_four_sign_pair_images_verified": True}


@cache
def real_frame():
    h = paley_748()
    h = h * h[:, :1]
    h = h * h[:1]
    check = check_sign_hadamard(h)
    t = h[:, :3]
    buckets = {pattern: np.flatnonzero(np.all(t == pattern, axis=1)).tolist()
               for pattern in ((1, 1, 1), (1, -1, -1), (1, 1, -1), (1, -1, 1))}
    if any(len(rows) != 187 for rows in buckets.values()):
        raise ValueError("the normalized three-column pattern counts changed")
    pairs = []
    for first_pattern, second_pattern in (((1, 1, 1), (1, -1, -1)), ((1, 1, -1), (1, -1, 1))):
        first, second = buckets[first_pattern], buckets[second_pattern]
        # Mix93 unequal-pattern pairs, then pair the94 leftovers within
        # each pattern.  This keeps both active and inactive rows in every
        # pattern, so the partial-frame deformation is not a phase gauge.
        pairs.extend(zip(first[:93], second[:93], strict=True))
        for remaining in (first[93:], second[93:]):
            pairs.extend(zip(remaining[::2], remaining[1::2], strict=True))
    if sorted(index for pair in pairs for index in pair) != list(range(748)):
        raise ValueError("the row pairs do not partition the real input")
    support = np.zeros_like(t)
    for first, second in pairs:
        support[first] = support[second] = (t[first] != t[second]).astype(int)
    if (not np.all(support[:, 0] == 0) or not np.array_equal(support[:, 1], support[:, 2])
            or np.count_nonzero(support[:, 1]) != 372):
        raise ValueError("the mixed boundary support is not the intended even-row pattern")
    return h, t, support, [list(pair) for pair in pairs], check


def mixed_frames(phase=96):
    if phase not in (24, 48, 96):
        raise ValueError("the audited dyadic cubic phases are24,48,96")
    h, t, support, pairs, real_check = real_frame()
    scalar = pair_action_audit(phase)
    eta, cube = scalar["eta_exponent"], scalar["cube_exponent"]
    u = ((t == -1) * (phase // 2) + eta * support) % phase
    v = (-u) % phase
    expected = np.zeros((3, 3, root_vectors(phase).shape[1]), dtype=np.int64)
    expected[np.arange(3), np.arange(3), 0] = 748
    for frame in (u, v):
        if not np.array_equal(exact_column_gram(frame, phase), expected):
            raise ValueError("a mixed partial frame failed its exact column Gram")
    d = np.array([[cube if i == j else 0 for j in range(3)] for i in range(3)])
    expected_d = np.zeros_like(expected)
    expected_d[np.arange(3), np.arange(3), 0] = 3
    if not np.array_equal(exact_column_gram(d, phase), expected_d):
        raise ValueError("the order-three border failed its exact Gram")
    if np.any(u % 3 == 2) or np.any((-v) % 3 == 2) or np.array_equal(u, v):
        raise ValueError("the transpose boundary or permitted alphabets failed")
    active = next(i for i in range(748) if support[i, 1] and np.array_equal(t[i], (1, 1, 1)))
    inactive = next(i for i in range(748) if not support[i, 1] and np.array_equal(t[i], (1, 1, 1)))
    cross = int((u[active, 1] + u[inactive, 0] - u[active, 0] - u[inactive, 1]) % phase)
    if cross != eta:
        raise ValueError("the genuine cyclotomic frame cross ratio disappeared")
    full_exponents = (h == -1).astype(np.int64) * (phase // 2)
    for first, second in pairs:
        changed = h[first] != h[second]
        full_exponents[first, changed] += eta
        full_exponents[second, changed] += eta
    full_exponents %= phase
    if not np.array_equal(full_exponents[:, :3], u):
        raise ValueError("the mixed frame is not the first three columns of the full constructed Butson")
    return {"phase": phase, "source_order_sought": 1499, "s": 748, "border_order": 3,
            "conditional_real_lift_order": 1499 * (phase // 3), "new_real_order_claimed": False,
            "real_input": real_check, "scalar_pair_action": scalar,
            "row_pair_indices": pairs, "T_signs": t.tolist(),
            "U_exponents": u.tolist(), "V_exponents": v.tolist(),
            "cubic_support_A": support.tolist(), "active_rows": int(support[:, 1].sum()),
            "D_exponents": d.tolist(), "Dstar_D_and_D_Dstar_equal_3I": True,
            "Ustar_U_and_Vstar_V_equal_748I": True,
            "stacked_boundaries_have_column_Gram_1496I": True,
            "both_source_boundary_blocks_in_permitted_two_cosets": True,
            "V_equals_conjugate_U_and_bottom_block_is_U_transpose": True,
            "adjoint_boundary_U_equals_V_is_false": True,
            "cross_ratio": {"rows": [active, inactive], "columns": [1, 0], "exponent": cross,
                            "primitive_full_phase_order": phase,
                            "excludes_gauge_to_any_lower_dyadic_level": True},
            "full_BH748_constructed_by_exact_unitary_pair_action": True,
            "full_BH748_exponents_sha256": hashlib.sha256(full_exponents.astype(np.uint8).tobytes()).hexdigest(),
            "X_Y_core_constructed": False}


def full_residue_audit(support, border_support=None, *, dense_core=False):
    a = np.asarray(support, dtype=np.int64)
    s, border = a.shape
    c = np.eye(border, dtype=np.int64) if border_support is None else np.asarray(border_support, dtype=np.int64)
    if c.shape != (border, border) or not np.all((c == 0) | (c == 1)):
        raise ValueError("the border cubic support has the wrong shape or alphabet")
    n = 2 * s + border
    rows = []
    boundary_words = [sum(int(value) << j for j, value in enumerate(row)) for row in a]
    core_words = [(((1 << s) - 1) ^ (1 << i)) if dense_core else (1 << i) for i in range(s)]
    rows.extend((core_words[i] << s) | (boundary_words[i] << (2 * s)) for i in range(s))
    rows.extend(core_words[i] | (boundary_words[i] << (2 * s)) for i in range(s))
    for j in range(border):
        row = sum(int(value) << (2 * s + k) for k, value in enumerate(c[j]))
        for i in range(s):
            if a[i, j]:
                row |= (1 << i) | (1 << (s + i))
        rows.append(row)
    if any(row.bit_count() % 2 != 1 for row in rows):
        raise ValueError("the full cubic support does not have row parity1")
    checks = 0
    for i, row in enumerate(rows):
        for j in range(i + 1):
            other = rows[j]
            n11 = (row & other).bit_count() % 2
            n10 = (row & ~other).bit_count() % 2
            n01 = (other & ~row).bit_count() % 2
            n00 = (n - (row | other).bit_count()) % 2
            # Codes1,2,3 denote1,omega,omega^2 in F4.
            value = (n00 ^ n11) ^ (2 * n10) ^ (3 * n01)
            if n11 != int(i == j) or value != int(i == j):
                raise ValueError("the full source binary or F4 row Gram failed")
            checks += 1
    width = (n + 7) // 8
    digest = hashlib.sha256(b"".join(row.to_bytes(width, "little") for row in rows)).hexdigest()
    report = {"source_residue_order": n, "support_formula": "P=[[0,I,A],[I,0,A],[A^T,A^T,I3]]",
            "core_residue_choices": "X=J, Y=J+(omega-1)I",
            "binary_support_row_parity_one": True, "binary_support_Gram_equals_I": True,
            "complete_F4_source_row_Gram_equals_I": True,
            "all_row_pairs_including_diagonal_checked": checks,
            "packed_support_sha256": digest,
            "literal_1_omega_core_does_not_lift": "Its F=X-Y is (1-omega)I, with Gram3I rather than1499I.",
            "other_dyadic_phase_choices_in_the_same_support_are_not_excluded": True}
    if dense_core:
        report.update({"support_formula": "P=[[0,K,A],[K,0,A],[A^T,A^T,C]], K=J_s+I_s overF2",
                       "core_residue_choices": "X=J, Y=J+(omega-1)(J-I)",
                       "core_cubic_support_is_connected": s > 2,
                       "border_cubic_support": c.tolist(),
                       "literal_1_omega_core_does_not_lift":
                       f"F=(1-omega)(J-I) has Gram3((s-2)J+I), not{n}I.",
                       "only_a_finite_field_completion_is_claimed": True})
        del report["other_dyadic_phase_choices_in_the_same_support_are_not_excluded"]
    return report


def seven_border(phase):
    """Permute a known border so its cubic support fixes a weight-four row."""
    original = border_source(7, phase)
    original_support = (original % 3 == 1).astype(np.int64)
    for mask in range(1, 128):
        if mask.bit_count() != 4:
            continue
        a = np.array([(mask >> j) & 1 for j in range(7)], dtype=np.int64)
        image = (a @ original_support) % 2
        if image.sum() != 4:
            continue
        choices = {value: iter(np.flatnonzero(image == value).tolist()) for value in (0, 1)}
        permutation = [next(choices[int(value)]) for value in a]
        d = original[:, permutation]
        c = (d % 3 == 1).astype(np.int64)
        if not np.array_equal(a @ c % 2, a) or not np.array_equal(c @ c.T % 2, np.eye(7, dtype=int)):
            raise ValueError("the permuted order-seven support does not have the required fixed vector")
        return d, c, a, permutation
    raise ValueError("no weight-four vector retained its weight under the border support")


@cache
def doubled_real_frame():
    base = real_frame()[0]
    h = np.block([[base, base], [base, -base]])
    _, c, a, permutation = seven_border(48)
    columns = [j + 748 * int(a[j]) for j in range(7)]
    t = h[:, columns]
    pairs = [[i, i + 748] for i in range(748)]
    # A proper subset of the opposite-Walsh-row pairs receives the
    # primitive parameter. Both rows in each active pair have support a.
    active_pairs = list(range(374))
    support = np.zeros_like(t)
    for i in active_pairs:
        support[i] = support[i + 748] = a
    if (not np.array_equal(t.T @ t, 1496 * np.eye(7, dtype=int))
            or np.any(support.sum(axis=0) % 2)
            or not np.array_equal(support @ c.T % 2, support)):
        raise ValueError("the doubled real partial frame or its cubic incidence failed")
    return h, t, support, pairs, active_pairs, columns, a, permutation, check_sign_hadamard(h)


def mixed_frames_2999(phase=48):
    if phase not in (48, 96):
        raise ValueError("the audited order2999 phases are48 and96")
    h, t, support, pairs, active_pairs, columns, a, permutation, real_check = doubled_real_frame()
    scalar = pair_action_audit(phase)
    eta = scalar["eta_exponent"]
    u = ((t == -1) * (phase // 2) + eta * support) % phase
    v = (-u) % phase
    d, c, fixed, border_permutation = seven_border(phase)
    if not np.array_equal(fixed, a) or border_permutation != permutation:
        raise ValueError("the two phases use incompatible border column permutations")
    expected = np.zeros((7, 7, root_vectors(phase).shape[1]), dtype=np.int64)
    expected[np.arange(7), np.arange(7), 0] = 1496
    if any(not np.array_equal(exact_column_gram(frame, phase), expected) for frame in (u, v)):
        raise ValueError("a1496-by7 frame failed its exact column Gram")
    expected[np.arange(7), np.arange(7), 0] = 7
    if not np.array_equal(exact_column_gram(d, phase), expected):
        raise ValueError("the permuted order-seven border failed its exact Gram")
    if np.any(u % 3 == 2) or np.any((-v) % 3 == 2) or np.array_equal(u, v):
        raise ValueError("the new mixed boundary failed its permitted alphabet or transpose condition")
    full_exponents = (h == -1).astype(np.int64) * (phase // 2)
    for index in active_pairs:
        first, second = pairs[index]
        changed = h[first] != h[second]
        full_exponents[first, changed] += eta
        full_exponents[second, changed] += eta
    full_exponents %= phase
    if not np.array_equal(full_exponents[:, columns], u):
        raise ValueError("the selected columns are not those of the full constructed Butson")
    active, inactive = active_pairs[0], next(i for i in range(748) if i not in active_pairs)
    one, zero = int(np.flatnonzero(a)[0]), int(np.flatnonzero(a == 0)[0])
    cross = int((u[active, one] + u[inactive, zero] - u[active, zero] - u[inactive, one]) % phase)
    if gcd(cross, phase) != 1:
        raise ValueError("the mixed partial frame has lost the primitive cyclotomic cross ratio")
    return {"phase": phase, "source_order_sought": 2999, "s": 1496, "border_order": 7,
            "conditional_real_lift_order": 2999 * (phase // 3), "new_real_order_claimed": False,
            "real_input": real_check, "scalar_pair_action": scalar,
            "row_pair_indices": pairs, "active_row_pair_indices": active_pairs,
            "selected_real_column_indices": columns, "T_signs": t.tolist(),
            "U_exponents": u.tolist(), "V_exponents": v.tolist(),
            "cubic_support_A": support.tolist(), "active_rows": 2 * len(active_pairs),
            "D_exponents": d.tolist(), "border_cubic_support_C": c.tolist(),
            "border_column_permutation": permutation, "even_fixed_vector": a.tolist(),
            "Dstar_D_and_D_Dstar_equal_7I": True,
            "Ustar_U_and_Vstar_V_equal_1496I": True,
            "stacked_boundaries_have_column_Gram_2992I": True,
            "both_source_boundary_blocks_in_permitted_two_cosets": True,
            "V_equals_conjugate_U_and_bottom_block_is_U_transpose": True,
            "adjoint_boundary_U_equals_V_is_false": True,
            "cross_ratio": {"rows": [active, inactive], "columns": [one, zero], "exponent": cross,
                            "primitive_full_phase_order": phase,
                            "excludes_gauge_to_any_lower_dyadic_level": True},
            "full_BH1496_constructed_by_exact_unitary_pair_action": True,
            "full_BH1496_exponents_sha256": hashlib.sha256(full_exponents.astype(np.uint8).tobytes()).hexdigest(),
            "X_Y_core_constructed": False}


def norm_2999_control(phase):
    """Verify an explicit integral scalar norm, without PARI assumptions."""
    coefficients = [8, -32, 10, -7, 20, 14, -12, -10, 12, 14, -20, -7, -10, -32, -8, 0]
    field = CyclotomicField(phase)
    root = field.zeta()
    z = sum((coefficient * root ** (j * (phase // 48)) for j, coefficient in enumerate(coefficients)), field.zero())
    if z * z.conjugate() != field.coerce(2999):
        raise ValueError("the exact scalar norm2999 witness failed")
    return {"phase": phase, "base_phase": 48, "coefficients_ascending_in_zeta48": coefficients,
            "z_conjugate_z_equals_2999": True,
            "L_equals_z_times_I1489_solves_metric_norm_only": True,
            "does_not_supply_a_root_valued_X_Y_split": True}


def dense_certificate():
    old_route = mixed_frames(96)
    new_routes = [mixed_frames_2999(phase) for phase in (48, 96)]
    new_route = new_routes[0]
    return {"schema": "mixed-cubic-boundaries-dense-v2",
            "followup_first_ramification_audit": {
                "fixed_dense_residue_primitive_boundary_pairings_have_no_pi_squared_root_lift": True,
                "finite_field_Gram_and_known_mixed_frames_remain_valid": True,
                "proof": "docs/MIXED_CUBIC_RAMIFICATION.md",
                "certificate": "certificates/scalar-lift/mixed-cubic-ramification.json",
                "general_Butson_source_orders_are_not_excluded": True},
            "historical_v1": {"certificate": "mixed-cubic-boundaries.json",
                              "finite_field_computation_remains_correct": True,
                              "old_PX_zero_PY_identity_support_cannot_lift_at_any_dyadic_level": True,
                              "reason": "Its cubic coefficient is monomial; the connected-support obstruction excludes it.",
                              "proof": "docs/CUBIC_CORE_CONNECTIVITY.md",
                              "old_not_excluded_field_was_limited_to_the_residue_check_and_is_superseded": True},
            "routes": [old_route, *new_routes],
            "dense_full_source_residues": [
                full_residue_audit(old_route["cubic_support_A"], dense_core=True),
                full_residue_audit(new_route["cubic_support_A"], new_route["border_cubic_support_C"], dense_core=True)],
            "exact_norm_2999_controls": [norm_2999_control(phase) for phase in (48, 96)],
            "exact_remaining_completion_2999": {
                "E": "X+Y", "F": "X-Y",
                "equations": ["FF*=2999I1496", "EE*=2999I1496-2UU*", "E conjugate(U)=-UD*"],
                "root_condition": "(E+F)/2 and (E-F)/2 belong to the admitted two cubic cosets",
                "known_full_frame": "G=[U Z] is the constructed BH1496,q after a column permutation; Z has1489 columns",
                "metric_factorization": "E=(-U D* U^T+Z L Z^T)/1496, LL*=2999I1489",
                "metric_factorization_solves_boundary_and_E_Gram": True,
                "metric_factorization_does_not_guarantee_integrality_or_the_root_split": True},
            "independent_exact_metric_control": metric_control(),
            "characteristic_zero_large_source_constructed": False,
            "new_real_Hadamard_constructed": False, "native_probe_run": False}


def metric_control():
    """Independent exact rational-cyclotomic test of the metric placement."""
    field = CyclotomicField(24)
    root = field.zeta()
    h = np.array([[1, 1, 1, 1], [1, -1, 1, -1],
                  [1, 1, -1, -1], [1, -1, -1, 1]])
    scalar = pair_action_audit(24)
    exponents = (h == -1).astype(int) * 12
    for first, second in ((0, 1), (2, 3)):
        changed = h[first] != h[second]
        exponents[first, changed] += scalar["eta_exponent"]
        exponents[second, changed] += scalar["eta_exponent"]
    exponents %= 24
    g = np.array([[root ** int(e) for e in row] for row in exponents], dtype=object)
    u, z = g[:, :3], g[:, 3:]
    d = np.array([[root ** scalar["cube_exponent"] if i == j else field.one()
                   for j in range(3)] for i in range(3)], dtype=object)

    def adjoint(matrix):
        return np.array([[value.conjugate() for value in row] for row in matrix], dtype=object).T

    scalar_norm = 3 + root**3 + root**9  # 3+sqrt(-2), exactly inQ(zeta24).
    if scalar_norm * scalar_norm.conjugate() != field.coerce(11):
        raise ValueError("the independent metric-control scalar norm failed")
    e = (-u @ adjoint(d) @ u.T + z @ np.array([[scalar_norm]], dtype=object) @ z.T) / 4
    identity = np.array([[field.coerce(11 if i == j else 0) for j in range(4)] for i in range(4)])
    if (not np.array_equal(e @ np.array([[value.conjugate() for value in row] for row in u]), -u @ adjoint(d))
            or not np.array_equal(e @ adjoint(e), identity - 2 * u @ adjoint(u))):
        raise ValueError("the exact metric factorization has a misplaced transpose or conjugate")
    return {"s": 4, "border_order": 3, "formal_source_order": 11, "phase": 24,
            "full_frame_exponents": exponents.tolist(),
            "L": "3+zeta24^3+zeta24^9=3+sqrt(-2)", "L_norm": 11,
            "E_coefficients_ascending": [[[str(c) for c in value.coefficients] for value in row] for row in e],
            "E_conjugate_U_equals_minus_U_Dstar_verified": True,
            "E_Estar_equals_11I_minus_2U_Ustar_verified": True,
            "F_or_BH11_constructed": False}


def certificate():
    route = mixed_frames(96)
    return {"schema": "mixed-cubic-boundaries-v1", "route": route,
            "full_source_residue": full_residue_audit(route["cubic_support_A"]),
            "independent_exact_metric_control": metric_control(),
            "exact_remaining_completion": {
                "E": "X+Y", "F": "X-Y",
                "equations": ["FF*=1499I748", "EE*=1499I748-2UU*", "EV=-UD*"],
                "entry_condition": "(E+F)/2 and (E-F)/2 lie in mu32 union omega*mu32",
                "known_full_frame": "G=[U Z] is the constructed BH748,96, with G*G=748I",
                "metric_factorization": "E=(-U D* U^T+Z L Z^T)/748, LL*=1499I745",
                "metric_factorization_solves_boundary_and_E_Gram": True,
                "metric_factorization_does_not_guarantee_integrality_or_the_root_split": True},
            "characteristic_zero_large_source_constructed": False,
            "new_real_Hadamard_constructed": False, "native_probe_run": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dense", action="store_true", help="replay the current dense-core1499/2999 certificate")
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--output", type=Path)
    action.add_argument("--check", type=Path)
    args = parser.parse_args()
    report = dense_certificate() if args.dense else certificate()
    if args.check:
        if json.loads(args.check.read_text()) != report:
            raise ValueError("the mixed-boundary certificate failed replay")
    else:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"known_full_input_orders": [748, 1496] if args.dense else [748],
                      "full_source_residue_orders": [1499, 2999] if args.dense else [1499],
                      "core_support": "dense" if args.dense else "historical monomial; lift excluded",
                      "characteristic_zero_core_constructed": False, "native_probe_run": False}))


if __name__ == "__main__":
    main()

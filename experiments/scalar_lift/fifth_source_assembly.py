"""Exact smaller-Butson inputs for fifth-coset bordered source completions.

This certifies D and T, not the missing square blocks X and Y.  The partial
column stack is not itself a new Butson matrix or a real Hadamard matrix.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from functools import cache
from pathlib import Path

import numpy as np

from experiments.scalar_lift.bordered_seed_assembly import check_sign_hadamard
from experiments.scalar_lift.fifth_dyadic_tower import phase_blocks
from experiments.scalar_lift.fifth_gauged_seed_escape import gauged_five
from hadamard.coefficients import CyclotomicField


@cache
def root_vectors(phase):
    field = CyclotomicField(phase)
    return np.array([[int(c) for c in (field.zeta() ** e).coefficients]
                     for e in range(phase)], dtype=np.int64)


def exact_column_gram(exponents, phase):
    """Return the complete column Gram in the cyclotomic power basis."""
    source = np.asarray(exponents, dtype=np.int64) % phase
    roots = root_vectors(phase)
    columns = source.shape[1]
    return np.array([[roots[(source[:, j] - source[:, i]) % phase].sum(axis=0)
                      for j in range(columns)] for i in range(columns)])


def check_columns(exponents, phase):
    source = np.asarray(exponents, dtype=np.int64)
    rows, columns = source.shape
    if np.any(source % 5 == 0):
        raise ValueError("a source entry lies in the forbidden dyadic subgroup")
    gram = exact_column_gram(source, phase)
    expected = np.zeros_like(gram)
    expected[np.arange(columns), np.arange(columns), 0] = rows
    if not np.array_equal(gram, expected):
        raise ValueError("the exact cyclotomic column Gram failed")
    return {"shape": [rows, columns], "column_Gram_equals_rows_I": True,
            "all_entries_avoid_dyadic_subgroup": True,
            "exponent_sha256": hashlib.sha256(
                json.dumps(source.tolist(), separators=(",", ":")).encode()
            ).hexdigest()}


def fifth_seed(phase):
    if phase not in (40, 80):
        raise ValueError("the certified source phases are 40 and 80")
    return np.array([[(phase // 5) * ((3 * (i + j) ** 2 + 1) % 5)
                      for j in range(5)] for i in range(5)], dtype=np.int64)


def extended_eight_seed(phase):
    """zeta5 times a Dita BH(8,d), d=phase/5, with primitive d cross ratio."""
    fifth_seed(phase)
    return np.array([[(phase // 5 + (phase // 2) * a * b
                       + (phase // 4) * r * s + 5 * b * r) % phase
                      for b in range(2) for s in range(4)]
                     for a in range(2) for r in range(4)], dtype=np.int64)


def seed_assembly(source_order, phase, *, rotate_one_block=True,
                  rotate_one_eight=False, use_gauged_five=True):
    if (source_order, phase) not in ((1499, 40), (1499, 80), (2999, 80)):
        raise ValueError("this pair is not one of the audited target routes")
    s = (source_order - 5) // 2
    copies_five, remainder = divmod(s - 4 * 8, 5)
    if remainder or copies_five < 0:
        raise ValueError("the smaller-source decomposition failed")
    d = fifth_seed(phase)
    b = extended_eight_seed(phase)
    blocks = [d] * copies_five
    if rotate_one_block:
        # Either rotation alone is excluded by a joint trace/unitarity
        # condition. Both rotations pass the necessary compression test;
        # they do not solve the simultaneous square completion.
        blocks[-1] = (d + 3 * (phase // 5)) % phase
    if use_gauged_five:
        blocks[-2] = (phase // 5) * np.array(gauged_five()[0], dtype=np.int64)
    eight_blocks = [b[:, :5]] * 4
    if rotate_one_eight:
        eight_blocks[-1] = (b[:, :5] + phase // 5) % phase
    t = np.vstack(blocks + eight_blocks)
    return t, d, b


def complete_tensor_source(phase):
    """A verified known-order control using the smaller BH5 and BH8 seeds."""
    d, b = fifth_seed(phase), extended_eight_seed(phase)
    # Remove B8's global fifth root before tensoring, to avoid cancelling
    # a fifth coset of D. The remaining B8 entries are dyadic roots.
    dyadic = (b - phase // 5) % phase
    source = np.array([[(int(e) + int(f)) % phase for e in row_d for f in row_b]
                       for row_d in d for row_b in dyadic], dtype=np.int64)
    check_columns(source, phase)
    return source


def completion_exponents(x, y, t, d, phase):
    """Assemble a candidate only; callers must separately verify its Gram."""
    x, y, t, d = (np.asarray(v, dtype=np.int64) for v in (x, y, t, d))
    s = len(x)
    if (x.shape != (s, s) or y.shape != (s, s)
            or t.shape != (s, 5) or d.shape != (5, 5)):
        raise ValueError("incompatible source block dimensions")
    result = np.block([[x, y, t], [y, x, t], [-t.T, -t.T, d]]) % phase
    if np.any(result % 5 == 0):
        raise ValueError("the candidate leaves the admitted fifth-root alphabet")
    return result


def build_report():
    routes = []
    for order, phase in ((1499, 40), (1499, 80), (2999, 80)):
        t, d, b = seed_assembly(order, phase)
        s = len(t)
        checks = {"D": check_columns(d, phase), "B8": check_columns(b, phase),
                  "M": check_columns((phase // 5) * np.array(gauged_five()[0]), phase),
                  "T": check_columns(t, phase),
                  "conjugate_T": check_columns(-t, phase)}
        cross_ratio = int((b[0, 0] + b[1, 4] - b[0, 4] - b[1, 0]) % phase)
        if cross_ratio != 5:
            raise ValueError("the primitive dyadic cross ratio disappeared")
        routes.append({
            "source_order": order, "phase": phase, "real_lift_width": 2 * phase // 5,
            "real_target_order": order * (2 * phase // 5), "s": s,
            "small_seed_stack": {"BH5_copies": (s - 32) // 5,
                                 "ordinary_D_copies": (s - 32) // 5 - 2,
                                 "BH5_copies_rotated_by_zeta5_cubed": 1,
                                 "gauged_M_copies": 1,
                                 "partial_BH8_copies": 4,
                                 "partial_BH8_fifth_cosets": [1, 1, 1, 1], "columns": 5},
            "D_exponents": d.tolist(), "B8_exponents": b.tolist(),
            "M_exponents": ((phase // 5) * np.array(gauged_five()[0])).tolist(),
            "exact_checks": checks,
            "T_cross_ratio": {"B8_rows": [0, 1], "B8_columns": [0, 4],
                              "exponent": cross_ratio, "root_order": phase // 5,
                              "excludes_gauge_to_phase_dividing_half_phase": True},
            "bottom_Gram_identity_2T_star_T_plus_D_D_star_equals_NI": True,
            "initial_unrotated_stack_excluded_by_mod_two_trace": True,
            "single_rotation_variants_excluded_by_joint_trace_and_unitarity": True,
            "single_rotation_exclusion_certificate": "certificates/scalar-lift/fifth-global-residue.json",
            "rotated_stack_rowwise_mod_two_preflight": "surjective; not full Gram feasibility",
            "rowwise_mod_two_boundary_and_unit_norm_compatible": True,
            "mod_two_preflight_certificate": "certificates/scalar-lift/fifth-source-boundary-audit.json",
            "historical_simultaneous_rotation_compression_passed_but_sparse_witness_unrealizable": True,
            "current_gauged_seed_compression_rank": 100,
            "current_gauged_seed_compression_surjective": True,
            "compression_certificate": "certificates/scalar-lift/fifth-gauged-seed-escape.json",
            "necessary_repeated_group_E_core_order": 27 if order == 1499 else 29,
            "E_core_factorization_certificate": "certificates/scalar-lift/fifth-grouped-completion.json",
            "X_Y_split_still_required_after_E_core": True,
            "characteristic_zero_lift_must_vary_phases_across_repetition_labels": True,
            "source_completion_status": "UNSOLVED; necessary compression is surjective",
            "completion_conditions": ["F F*=N I_s", "E E*=N I_s-2 T T*",
                                      "E T=-T D*", "X=(E+F)/2 and Y=(E-F)/2 in L"],
            "X_and_Y_constructed": False,
            "source_completion_constructed": False,
        })
    controls = []
    for phase in (40, 80):
        source = complete_tensor_source(phase)
        blocks = phase_blocks(phase)
        real = np.block([[blocks[int(e)] for e in row] for row in source])
        controls.append({"phase": phase, "source_order": 40,
                         "source_exponents": source.tolist(),
                         "exact_source_check": check_columns(source, phase),
                         "implied_real_order": 40 * (2 * phase // 5),
                         "full_real_check": check_sign_hadamard(real),
                         "known_real_order_only": True})
    return {"scope": "Exact small-source assembly; high-order square completion remains unsolved",
            "alphabet": "L=mu_(5d) minus mu_d; equivalently exponent not divisible by5",
            "alphabet_is_closed_under_conjugation": True,
            "complete_source_controls": controls,
            "routes": routes, "new_real_Hadamard_order_constructed": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", type=Path)
    args = parser.parse_args()
    report = build_report()
    if args.check and json.loads(args.check.read_text()) != report:
        raise SystemExit("fifth source assembly certificate differs from replay")
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"verified_seed_routes": [r["real_target_order"] for r in report["routes"]],
                      "source_completions_constructed": False}))


if __name__ == "__main__":
    main()

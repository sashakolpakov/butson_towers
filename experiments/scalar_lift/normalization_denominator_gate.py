"""Necessary odd denominator primes for normalized scalar sign blocks.

This gate concerns an orthogonal rational field action and a Hadamard
normalization frame of width4r, r odd>1. Arbitrary weighted restrictions
or nonuniform source-coordinate changes are outside its scope.
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np

from experiments.scalar_lift.inverse_blocks import known_order_twelve_pair


def denominator_gate(width, localization_denominator, primitive_ratio_order):
    """Assume all scalar-action entries lie in Z[1/localization_denominator].

    The ratio must be an actually used root ratio, and its action must
    satisfy its cyclotomic polynomial, not just have that finite order.
    """
    for value in (width, localization_denominator, primitive_ratio_order):
        if type(value) is not int or value < 1:
            raise ValueError("positive integer parameters are required")
    odd_part_root = primitive_ratio_order
    while odd_part_root % 2 == 0:
        odd_part_root //= 2
    in_width_scope = width > 4 and width % 8 == 4
    odd_part_width = width//4 if in_width_scope else None
    excluded = (in_width_scope and odd_part_root > 1
                and math.gcd(localization_denominator, odd_part_width) == 1)
    return {
        "width": width, "primitive_used_ratio_order": primitive_ratio_order,
        "localization_denominator": localization_denominator,
        "width_is_four_times_an_odd_integer_greater_than_one": in_width_scope,
        "odd_width_factor": odd_part_width, "ratio_has_nontrivial_odd_part": odd_part_root > 1,
        "excluded": excluded,
        "required_odd_denominator_prime_must_divide": odd_part_width,
        "scope": "Uniform rational orthogonal FIELD action; every used image is a sign block with Gram width*I",
        "faithful_cyclic_group_action_alone_is_sufficient": False,
        "general_weighted_or_nonuniform_maps_excluded": False,
    }


def _pair_control(first, second, expected_denominator):
    d = len(first)
    unit = np.eye(d, dtype=np.int64)
    if not all(np.all(np.abs(h) == 1) and np.array_equal(h @ h.T, d*unit)
               for h in (first, second)):
        raise ArithmeticError("a control is not a pair of sign Hadamards")
    cross = second @ first.T
    if not np.array_equal(cross @ first, d*second):
        raise ArithmeticError("the control scalar action does not map its frame")
    if not np.array_equal(cross @ cross.T, d*d*unit):
        raise ArithmeticError("the control rational scalar action is not orthogonal")
    if np.any(cross @ cross-d*cross+d*d*unit):
        raise ArithmeticError("the primitive sixth-root polynomial failed")
    denominators = sorted({Fraction(int(x), d).denominator for x in cross.flat})
    common = math.lcm(*denominators)
    if common != expected_denominator:
        raise ArithmeticError("the exact control denominator changed")
    nonmonomial = bool(np.any(np.count_nonzero(cross, axis=1) != 1))
    if not nonmonomial:
        raise ArithmeticError("the intended nonmonomial control became monomial")
    return {"width": d, "both_complete_sign_Grams_verified": True,
            "primitive_sixth_root_action_verified": True,
            "exact_action_entry_denominators": denominators,
            "exact_common_action_denominator": common,
            "nonmonomial": nonmonomial,
            "gate": denominator_gate(d, common, 6)}


def certificate():
    h2 = np.array([[1, 1], [1, -1]], dtype=np.int64)
    h4 = np.kron(h2, h2)
    twice_r4 = np.array([[1, -1, -1, -1], [1, 1, -1, 1],
                         [1, 1, 1, -1], [1, -1, 1, 1]], dtype=np.int64)
    g4 = (twice_r4 @ h4)//2
    h8, g8 = np.kron(h4, h2), np.kron(g4, h2)
    h12, g12 = (np.array(h, dtype=np.int64) for h in known_order_twelve_pair())
    controls = [_pair_control(h4, g4, 2), _pair_control(h8, g8, 2),
                _pair_control(h12, g12, 6)]
    if any(control["gate"]["excluded"] for control in controls):
        raise ArithmeticError("a true positive control was incorrectly excluded")
    rejected = [denominator_gate(d, m, root) for d in (668, 716, 764, 892, 1132, 1388, 1436)
                for m in (2, 4, 8, 16, 32, 3, 5, 6, 10, 30)
                for root in (3, 6, 12)]
    if not all(case["excluded"] for case in rejected):
        raise ArithmeticError("a denominator lacking every width prime escaped the gate")
    admitted = [denominator_gate(4*r, 2*r, root) for r in (3, 167, 179, 191, 223)
                for root in (3, 6, 12)]
    if any(case["excluded"] for case in admitted):
        raise ArithmeticError("a width-prime denominator was incorrectly excluded")
    return {"schema": "normalization-denominator-gate-v1",
            "theorem": "At width4r, oddr>1, a nonmonomial rational orthogonal map between sign Hadamards needs an odd denominator prime dividingr",
            "signed_permutation_exclusion": "Every signed cycle has an eigenvalue of power-of-two order, incompatible with a cyclotomic polynomial of nontrivial odd part",
            "cubic_refinement": "A primitive sixth-root map has every entry odd/(2r); ifr prime its exact common denominator is2r",
            "nontrivial_positive_controls": controls,
            "rejected_localization_cases": rejected,
            "width_prime_localizations_not_excluded": admitted,
            "broad_dyadic_chart_excluded": True,
            "all_rational_basis_changes_excluded": False,
            "new_Hadamard_order_claimed": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--report", type=Path)
    group.add_argument("--check", type=Path)
    args = parser.parse_args()
    report = certificate()
    if args.report is not None:
        args.report.write_text(json.dumps(report, indent=2)+"\n")
    elif json.loads(args.check.read_text()) != report:
        raise ValueError("the denominator-gate certificate differs from exact replay")
    print(json.dumps({"rejected_localizations": len(report["rejected_localization_cases"]),
                      "positive_controls": 3, "all_rational_bases_excluded": False}))


if __name__ == "__main__":
    main()

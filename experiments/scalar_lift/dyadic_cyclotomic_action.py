"""Exact minimal-degree rational orthogonal actions for phases 24, 48, 96.

These certify the scalar/metric component only. Neither a Butson source
nor a compatible sign frame is supplied at the large target orders.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from sympy import Poly, Symbol, cyclotomic_poly


def action_numerator(negacyclic_size):
    """Return M with S=M/2, using Norm(1+sqrt(-2))=3 inside Q(zeta_8)."""
    n = negacyclic_size
    if n not in (4, 8, 16):
        raise ValueError("the certified sizes are 4, 8, and 16")
    r = np.zeros((n, n), dtype=np.int64)
    r[np.arange(n - 1), np.arange(1, n)] = 1
    r[-1, 0] = -1
    eye = np.eye(n, dtype=np.int64)
    j = np.linalg.matrix_power(r, n // 4) + np.linalg.matrix_power(r, 3 * n // 4)
    a = eye + j
    v = np.block([[eye, a], [-a.T, eye]])
    zero = np.zeros_like(r)
    diagonal_r = np.block([[r, zero], [zero, r]])
    if not np.array_equal(j.T, -j) or not np.array_equal(j @ j, -2 * eye):
        raise ValueError("the quadratic subfield identity failed")
    if not np.array_equal(a @ a.T, 3 * eye):
        raise ValueError("the relative norm identity failed")
    if not np.array_equal(diagonal_r @ v, v @ diagonal_r):
        raise ValueError("the scalar factors do not commute")
    return diagonal_r @ v


def verify_action(matrix_numerator, phase):
    """Use integer Gram and homogenized cyclotomic identities, without floats."""
    if phase not in (24, 48, 96):
        raise ValueError("the certified phases are 24, 48, and 96")
    entries = np.asarray(matrix_numerator)
    d = phase // 3
    if entries.shape != (d, d) or entries.dtype.kind not in "iu":
        raise ValueError("an integral square matrix of the prescribed size is required")
    if np.any(np.abs(entries.astype(object)) > 2):
        raise ValueError("the matrix entries exceed the certified numerator bound")
    # Object integers keep the independent verifier safe under arbitrary input.
    m = entries.astype(object)
    eye = np.eye(d, dtype=object)
    if not np.array_equal(m @ m.T, 4 * eye):
        raise ValueError("the rational action does not preserve the identity metric")
    n = d // 2
    residual = (np.linalg.matrix_power(m, 2 * n)
                - 2**n * np.linalg.matrix_power(m, n) + 2**(2 * n) * eye)
    if np.any(residual):
        raise ValueError("the cyclotomic identity failed")
    x = Symbol("x")
    polynomial = Poly(cyclotomic_poly(phase, x), x)
    if polynomial != Poly(x**d - x**n + 1, x):
        raise ValueError("the homogenized polynomial is not the cyclotomic polynomial")
    # Phi_q is irreducible over Q; an annihilating Phi_q of degree d in
    # dimension d is necessarily the minimal and characteristic polynomial.
    return {
        "phase": phase, "dimension": d, "denominator": 2,
        "minimal_and_characteristic_polynomial": str(polynomial.as_expr()),
        "faithful_field_degree": d,
        "rational_identity_metric_verified": True,
        "numerator": entries.tolist(),
    }


def build_report():
    actions = [verify_action(action_numerator(n), 6 * n) for n in (4, 8, 16)]
    m = np.asarray(actions[-1]["numerator"], dtype=object)
    eye = np.eye(32, dtype=object)
    restriction_tower = []
    for phase in (3, 12, 24, 48, 96):
        step = 96 // phase
        x = Symbol("x")
        polynomial = Poly(cyclotomic_poly(phase, x), x)
        degree = polynomial.degree()
        total_degree = step * degree
        residual = sum((
            int(polynomial.nth(k)) * 2**(total_degree - k * step)
            * np.linalg.matrix_power(m, k * step)
            for k in range(degree + 1) if polynomial.nth(k)
        ), np.zeros_like(eye))
        if np.any(residual):
            raise ValueError("a nested field restriction failed")
        restriction_tower.append({
            "phase": phase, "generator_power_of_S96": step,
            "field_degree": degree, "module_rank": 32 // degree,
            "rational_module_dimension": 32,
            "same_identity_metric_retained": True,
            "cyclotomic_annihilator_verified": True,
        })
    return {
        "scope": "Constructive scalar/metric inputs; no source or sign-frame witness",
        "quadratic_norm_identity": "(1+sqrt(-2))(1-sqrt(-2))=3",
        "tower_degrees": [2, 4, 8, 16, 32],
        "tower_fields": ["Q(zeta3)", "Q(zeta12)", "Q(zeta24)",
                         "Q(zeta48)", "Q(zeta96)"],
        "actions": actions,
        "nested_restrictions_of_one_action": restriction_tower,
        "conditional_target_examples": [
            {"order": 47984, "source_order": 2999, "phase": 48, "block_order": 16},
            {"order": 95968, "source_order": 2999, "phase": 96, "block_order": 32},
            {"order": 47968, "source_order": 1499, "phase": 96, "block_order": 32},
        ],
        "source_matrices_constructed": False,
        "compatible_sign_bases_constructed": False,
        "new_real_Hadamard_order_constructed": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", type=Path)
    args = parser.parse_args()
    report = build_report()
    if args.check and json.loads(args.check.read_text()) != report:
        raise SystemExit("the scalar-action certificate differs from its replay")
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"verified_phase_dimension_pairs": [
        [a["phase"], a["dimension"]] for a in report["actions"]
    ], "new_real_Hadamard_order_constructed": False}))


if __name__ == "__main__":
    main()

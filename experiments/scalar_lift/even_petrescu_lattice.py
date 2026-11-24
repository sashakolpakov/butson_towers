"""Exact 5-adic lattice factorization for the known even Petrescu inputs.

No corner or source is constructed. The congruence is a necessary part
of the joint source incidence, independent of any permutation symmetry.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from functools import cache
from pathlib import Path

import numpy as np

from experiments.scalar_lift.commuting_parity import paley_core
from experiments.scalar_lift.dyadic_partial_lift import sylvester
from experiments.scalar_lift.even_petrescu_sources import paley_500
from experiments.scalar_lift.fifth_source_assembly import root_vectors
from experiments.scalar_lift.paley_sextic_operator import rank_mod_prime

MODULUS = 125


@cache
def lattice_data():
    k = paley_500()
    r = np.ones(500, dtype=np.int64)
    r[1:] = -1
    s = r[:, None] * k
    identity = np.eye(500, dtype=np.int64)
    if (not np.array_equal(s + s.T, 2 * identity)
            or not np.array_equal(s @ s.T, 500 * identity)):
        raise ValueError("the skew Paley identities failed")
    p = s * pow(2, -1, MODULUS) % MODULUS
    complement = (identity - p) % MODULUS
    if (np.any((p @ p - p) % MODULUS)
            or np.any((p.T - complement) % MODULUS)):
        raise ValueError("the complementary idempotent identities failed")
    return k, r, s, p, complement


def block_projection(matrix):
    """P M (I-P) modulo 125, for a coefficient block M of D*."""
    _, _, _, p, complement = lattice_data()
    matrix = np.asarray(matrix, dtype=np.int64)
    if matrix.shape != (500, 500):
        raise ValueError("a Paley coefficient block must have order 500")
    return ((p @ (matrix % MODULUS)) % MODULUS) @ complement % MODULUS


def dense_controls():
    k, _, _, p, complement = lattice_data()
    for matrix in (k, k.T):
        if np.any(block_projection(matrix)):
            raise ValueError("a dense Hadamard control violates the lattice")
    # One full product verifies K K^T=500I. Associativity then proves both
    # K K K^T/500=K and K K^T K^T/500=K^T without extra large products.
    if not np.array_equal(k @ k.T, 500 * np.eye(500, dtype=np.int64)):
        raise ValueError("the dense-control exact Gram failed")
    changed = k.copy()
    changed[7, 13] *= -1
    delta = int(changed[7, 13] - k[7, 13])
    predicted = delta * np.outer(p[:, 7], complement[13, :]) % MODULUS
    actual = block_projection(changed)
    if not np.array_equal(actual, predicted) or not np.all(actual):
        raise ValueError("the full one-sign-flip negative control failed")
    return {
        "dense_coefficient_blocks_M": ["K", "K^T"],
        "both_lattice_projections_zero": True,
        "full_Gram_identity_K_Kt_equals_500I": True,
        "integral_full_transform_identities": [
            "K K K^T/500=K", "K K^T K^T/500=K^T"],
        "adjoint_convention": "M=D*, so D=K gives the transformed matrix K^T",
        "one_sign_flip_position": [7, 13],
        "negative_projection_equals_exact_rank_one_update": True,
        "nonzero_negative_projection_entries": int(np.count_nonzero(actual)),
        "Petrescu_corner_Gram_claimed": False,
    }


def root_support_audit(phase):
    roots = root_vectors(phase)
    differences = roots[:, None, :] - roots[None, :, :]
    root_bound = int(np.max(np.abs(roots)))
    difference_bound = int(np.max(np.abs(differences)))
    if root_bound != 1 or difference_bound != 2:
        raise ValueError("the cyclotomic coefficient bounds changed")
    if any(not np.any((roots[e] - roots[0]) % MODULUS) for e in range(1, phase)):
        raise ValueError("a distinct root equals one modulo 125")
    return {
        "phase": phase, "power_basis_degree": roots.shape[1],
        "all_root_coefficients_in_minus_one_zero_one": True,
        "ordered_root_difference_pairs_checked": phase * phase,
        "maximum_root_difference_coefficient": difference_bound,
        "all_nontrivial_single_entry_phase_changes_rejected": phase - 1,
        "minimum_support_of_nonzero_single_row_change": 63,
        "minimum_support_of_nonzero_single_row_sign_flips": 125,
        "bound_scope": "A change confined to one row (or one column) of one 500 by 500 block",
        "multirow_trade_lower_bound_claimed": False,
    }


def small_walsh_factor_controls():
    """Direct order-40/80 checks of the same factorization over Z/5Z."""
    modulus, m = 5, 20
    k = np.ones((m, m), dtype=np.int64)
    k[1:, 1:] = paley_core(m - 1) - np.eye(m - 1, dtype=np.int64)
    gauge = np.ones(m, dtype=np.int64)
    gauge[1:] = -1
    s = gauge[:, None] * k
    p = s * pow(2, -1, modulus) % modulus
    complement = (np.eye(m, dtype=np.int64) - p) % modulus
    if np.any((p @ p - p) % modulus) or np.any((p.T - complement) % modulus):
        raise ValueError("the small Paley projectors failed")
    records = []
    for r in (2, 4):
        walsh = sylvester(r)
        full_k = np.kron(walsh, k)
        w = np.kron(walsh, np.eye(m, dtype=np.int64))
        w_inverse = w * pow(r, -1, modulus) % modulus
        full_gauge = np.tile(gauge, r)
        for flipped in (False, True):
            matrix = np.kron(np.ones((r, r), dtype=np.int64), k)
            if flipped:
                matrix[3, m + 7] *= -1
            projected = np.zeros_like(matrix)
            for i in range(r):
                for j in range(r):
                    block = matrix[i * m:(i + 1) * m, j * m:(j + 1) * m]
                    projected[i * m:(i + 1) * m, j * m:(j + 1) * m] = (
                        p @ block @ complement) % modulus
            actual = full_k @ matrix @ full_k.T % modulus
            mixed = 4 * (w @ projected @ w.T)
            predicted = full_gauge[:, None] * mixed * full_gauge[None, :] % modulus
            recovered = (pow(4, -1, modulus) * w_inverse
                         @ (full_gauge[:, None] * actual * full_gauge[None, :])
                         @ w_inverse.T) % modulus
            if not np.array_equal(actual, predicted) or not np.array_equal(recovered, projected):
                raise ValueError("the direct Walsh factorization or its inverse failed")
            if bool(np.any(actual)) != flipped or bool(np.any(projected)) != flipped:
                raise ValueError("the direct zero-congruence equivalence failed")
            records.append({"Walsh_order": r, "Paley_order": m,
                            "full_order": r * m, "modulus": modulus,
                            "one_entry_flipped": flipped,
                            "forward_factorization_verified": True,
                            "inverse_recovers_all_block_projections": True,
                            "both_zero_congruences_hold": not flipped})
    return records


@cache
def build_report():
    k, _, s, p, complement = lattice_data()
    rank = rank_mod_prime(p, 5)
    if rank != 250 or not np.all(p % 5) or not np.all(complement % 5):
        raise ValueError("the projector rank or unit-entry property failed")
    return {
        "schema": "even-petrescu-lattice-v1",
        "scope": "Exact 5^3-part of the joint integrality condition; no symmetry ansatz",
        "Paley_order": 500, "modulus": MODULUS,
        "K_binary_sha256": hashlib.sha256(((k + 1) // 2).astype(np.uint8).tobytes()).hexdigest(),
        "S_sign_sha256": hashlib.sha256(((s + 1) // 2).astype(np.uint8).tobytes()).hexdigest(),
        "S_definition": "S=R K, R=diag(1,-1,...,-1)",
        "S_plus_St_equals_2I": True, "S_St_equals_500I": True,
        "P_definition": "P=S/2 modulo 125",
        "P_squared_equals_P": True, "Pt_equals_I_minus_P": True,
        "rank_P_modulo_5_by_full_elimination": rank,
        "all_entries_of_P_and_I_minus_P_are_units_modulo_125": True,
        "input_orders": [1000, 2000],
        "block_condition": "P (D*)_ij (I-P)=0 modulo 125 in every power-basis coefficient",
        "equivalent_to_5_cubed_part_of_E_integrality_for_regular_D": True,
        "powers_of_two_integrality_and_root_sum_constraints_remain": True,
        "dense_controls": dense_controls(),
        "root_support_audits": [root_support_audit(phase) for phase in (24, 48)],
        "direct_small_Walsh_factor_controls": small_walsh_factor_controls(),
        "monomial_or_permutation_symmetry_required": False,
        "native_solver_run": False, "large_source_constructed": False,
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
            raise ValueError("the lattice certificate differs from exact replay")
    else:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"projector_rank_mod_5": 250,
                      "minimum_single_row_change_support": 63,
                      "minimum_single_row_sign_flip_support": 125,
                      "large_source_constructed": False}))


if __name__ == "__main__":
    main()

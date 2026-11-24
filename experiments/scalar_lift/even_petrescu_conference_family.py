"""Exact cyclotomic boundary sources from the Paley conference of order 500.

This constructs complete known-order BH(1000,8/16) and BH(2000,8/16)
inputs. A related symmetric array is only a partial Petrescu corner.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from functools import cache
from math import gcd
from pathlib import Path

import numpy as np

from experiments.scalar_lift.commuting_parity import paley_core
from experiments.scalar_lift.even_petrescu_lattice import block_projection, lattice_data
from experiments.scalar_lift.paley_sextic_operator import rank_mod_prime
from hadamard.cgw import is_cgw_exponent_matrix
from hadamard.coefficients import CyclotomicField


def multiply(first, second, weight):
    """Multiply a+bC using C^2=-weight I."""
    a, b = first
    c, d = second
    return a * c - weight * b * d, a * d + b * c


def adjoint(value):
    a, b = value
    return a.conjugate(), -b.conjugate()


def add(first, second):
    return tuple(a + b for a, b in zip(first, second, strict=True))


def gram(blocks, weight, field):
    zero = (field.zero(), field.zero())
    return [[sum_pairs((multiply(blocks[i][k], adjoint(blocks[j][k]), weight)
                        for k in range(len(blocks))), zero)
             for j in range(len(blocks))] for i in range(len(blocks))]


def sum_pairs(values, initial):
    result = initial
    for value in values:
        result = add(result, value)
    return result


def coefficients(value):
    return [[int(coefficient) for coefficient in entry.coefficients] for entry in value]


def block_family(phase):
    if phase not in (8, 16):
        raise ValueError("the saved boundary controls use phases 8 and 16")
    field = CyclotomicField(phase)
    z, one = field.zeta(), field.one()
    imaginary = z ** (phase // 4)
    a, b = (z, one), (imaginary * z.conjugate(), imaginary)
    minus_b_star = tuple(-entry for entry in adjoint(b))
    return field, [[a, b], [minus_b_star, adjoint(a)]]


def formal_family_certificate(phase, weight=499):
    field, blocks = block_family(phase)
    result = gram(blocks, weight, field)
    for i in range(2):
        for j in range(2):
            expected = (field.element([2 * (weight + 1) if i == j else 0]), field.zero())
            if result[i][j] != expected:
                raise ValueError("the complete formal conference-family Gram failed")
    a, b = blocks[0]
    complementary = add(multiply(a, adjoint(a), weight), multiply(b, adjoint(b), weight))
    if complementary != (field.element([2 * (weight + 1)]), field.zero()):
        raise ValueError("the complementary block Gram failed")
    return {
        "phase": phase, "conference_weight": weight,
        "block_coefficients_in_I_C": [[coefficients(value) for value in row] for row in blocks],
        "all_four_full_Gram_blocks": [[coefficients(value) for value in row] for row in result],
        "complete_Gram_equals_twice_weight_plus_one_I": True,
        "proof_uses_exact_C_star_minus_C_and_C_squared_minus_weight_I": True,
        "complementary_pair_Gram_verified": True,
        "root_entry_proof": "C has zero diagonal and off-diagonal signs; every diagonal and off-diagonal coefficient is a root",
    }


def zero_modulo(value, modulus=125):
    return all(coefficient.denominator == 1 and int(coefficient) % modulus == 0
               for entry in value for coefficient in entry.coefficients)


def kernel_certificate(phase):
    """Exact graph kernel and block factorization over the cyclotomic ring mod125."""
    field, blocks = block_family(phase)
    z, one, zero = field.zeta(), field.one(), field.zero()
    imaginary = z ** (phase // 4)
    graph = (zero, -imaginary * z.conjugate())
    identity = (one, zero)
    graph_image = [add(multiply(row[0], graph, 499), row[1]) for row in blocks]
    expected = [(500 * imaginary * z.conjugate(), zero),
                (500 * z.conjugate(), zero)]
    if graph_image != expected or not all(zero_modulo(value) for value in graph_image):
        raise ValueError("the full graph-kernel identity failed")
    if not zero_modulo(add(multiply(graph, adjoint(graph), 499), identity)):
        raise ValueError("the graph is not totally isotropic modulo125")

    channels = []
    for epsilon in (1, -1):
        channel = [[a + epsilon * b for a, b in row] for row in blocks]
        determinant = channel[0][0] * channel[1][1] - channel[0][1] * channel[1][0]
        kernel = [-imaginary * epsilon * z.conjugate(), one]
        if determinant != zero or any(sum((row[j] * kernel[j] for j in range(2)), zero)
                                      != zero for row in channel):
            raise ValueError("a two-dimensional channel or its exact kernel failed")
        inverse_pivot = (63 * sum(((-z) ** j for j in range(phase // 2)), zero)
                         if epsilon == 1 else
                         -63 * sum((z ** j for j in range(phase // 2)), zero))
        if not zero_modulo((channel[0][0] * inverse_pivot - one, zero)):
            raise ValueError("a channel pivot is not a unit modulo125")
        channels.append({"C_eigenvalue_mod_125": epsilon,
                         "channel_coefficients": [[coefficients((entry, zero))[0]
                                                   for entry in row] for row in channel],
                         "determinant_zero_exactly": True,
                         "unit_pivot_inverse_coefficients": coefficients((inverse_pivot, zero))[0],
                         "kernel_vector": [coefficients((entry, zero))[0] for entry in kernel],
                         "rank_over_cyclotomic_ring_mod_125": 1,
                         "Paley_summand_rank": 250})

    # H = U [I,-G] mod125, with U its first block column. Its top block
    # C+zI is invertible, proving the lossless coupled congruence.
    minus_graph = tuple(-entry for entry in graph)
    factor_defects = [tuple(row[1][j] - multiply(row[0], minus_graph, 499)[j]
                           for j in range(2)) for row in blocks]
    if not all(zero_modulo(value) for value in factor_defects):
        raise ValueError("the graph block factorization failed modulo125")
    inverse_denominator = -63 * sum((z ** (2 * j) for j in range(phase // 4)), zero)
    inverse_a = (z * inverse_denominator, -inverse_denominator)
    product = multiply(blocks[0][0], inverse_a, 499)
    if not zero_modulo((product[0] - one, product[1])):
        raise ValueError("the top block of the graph factorization is not invertible")
    _, _, _, p, _ = lattice_data()
    if rank_mod_prime(p, 5) != 250:
        raise ValueError("the two Paley summands do not have the claimed rank")
    return {
        "phase": phase, "modulus": 125,
        "channels": channels,
        "graph_operator": "G_z=-i*conjugate(z)*C",
        "graph_operator_coefficients_in_I_C": coefficients(graph),
        "full_kernel": "{(G_z*w,w): w in R^500}, R=Z[zeta_d]/125",
        "exact_graph_image_before_reduction": [coefficients(value) for value in graph_image],
        "exact_graph_image_formula": "H(z)[G_z;I]=500*conjugate(z)*[iI;I]",
        "kernel_rank": 500,
        "each_Paley_summand_contributes_rank": 250,
        "different_from_two_copies_of_old_Paley_kernel": True,
        "graph_is_totally_isotropic_modulo_125": True,
        "lossless_factorization": "H(z)=U[I,-G_z] modulo125; U is the first block column and has an invertible top block",
        "inverse_top_block_coefficients_in_I_C": coefficients(inverse_a),
        "new_coupled_congruence_for_M_equals_Dstar": "M11*G_z+M12=G_z*(M21*G_z+M22) modulo125",
        "equivalent_reduced_Gram_congruence": "[I,-G_z] M [I,-G_z]*=0 modulo125",
        "normalized_kernel_transport": "C0* ker(H(z)), C0=diag(conjugate(raw first row)*raw[0,0])",
        "normalization_star_means_adjoint": True,
        "order_2000_kernel_before_normalization": "two copies of this graph, since the outer H2 is invertible modulo125",
    }


def conference_exponents(conference, phase):
    """Build the full two-by-two family from a real skew conference."""
    c = np.asarray(conference, dtype=np.int64)
    m = len(c)
    if c.shape != (m, m) or np.any(np.diag(c)):
        raise ValueError("expected a square conference matrix with zero diagonal")
    mask = ~np.eye(m, dtype=bool)
    if not np.all(np.abs(c[mask]) == 1):
        raise ValueError("the conference off-diagonal entries must be signs")
    base = np.where(c == -1, phase // 2, 0).astype(np.int16)
    upper_left = base.copy()
    upper_right = (base + phase // 4) % phase
    lower_left = (base - phase // 4) % phase
    lower_right = (base + phase // 2) % phase
    np.fill_diagonal(upper_left, 1)
    np.fill_diagonal(upper_right, phase // 4 - 1)
    np.fill_diagonal(lower_left, phase // 4 + 1)
    np.fill_diagonal(lower_right, phase - 1)
    return np.block([[upper_left, upper_right], [lower_left, lower_right]]) % phase


def source_exponents(order, phase, normalized=True):
    if order not in (1000, 2000):
        raise ValueError("the saved known-order sources have orders 1000 and 2000")
    _, _, s, _, _ = lattice_data()
    c = s - np.eye(500, dtype=np.int64)
    source = conference_exponents(c, phase)
    if order == 2000:
        source = np.block([[source, source], [source, (source + phase // 2) % phase]])
    if normalized:
        source = (source - source[:, :1] - source[:1, :] + source[0, 0]) % phase
    return source.astype(np.int16)


def source_report(order, phase):
    source = source_exponents(order, phase)
    if np.any(source[0]) or np.any(source[:, 0]) or source[1, 2] != 1:
        raise ValueError("normalization or the primitive cross ratio failed")
    counts = np.bincount(source.ravel(), minlength=phase)
    return {
        "order": order, "phase": phase,
        "construction": "conference family" if order == 1000 else "H2 tensor conference family",
        "normalized_first_row_and_column_all_one": True,
        "all_entries_exact_dyadic_roots": True,
        "root_exponent_counts": [int(value) for value in counts],
        "normalized_exponent_sha256": hashlib.sha256(source.astype(np.uint8).tobytes()).hexdigest(),
        "primitive_cross_ratio_rows": [0, 1], "primitive_cross_ratio_columns": [0, 2],
        "primitive_cross_ratio_exponent": 1,
        "cross_ratio_order": phase // gcd(phase, 1),
        "full_Gram_certification": "exact order-500 conference identity and every formal two-by-two Gram block; tensor identity for order 2000",
        "T_shape_after_deleting_initial_row": [order - 1, order],
        "T_row_Gram": f"{order} I", "T_column_Gram": f"{order} I-J",
        "compatible_restricted_cubic_alphabet": f"mu_{phase} subset mu_{phase} union omega*mu_{phase}",
        "known_order_control_only": True,
        "higher_Petrescu_source_completed": False,
    }


def partial_corner_certificate():
    phase, weight = 8, 499
    field, blocks = block_family(phase)
    a, b = blocks[0]
    corner = [[a, b], [b, a]]
    result = gram(corner, weight, field)
    for i in range(2):
        for j in range(2):
            expected = (field.element([1000 if i == j else 2]), field.zero())
            if result[i][j] != expected:
                raise ValueError("the exact partial-corner Gram failed")
    difference = tuple(first - second for first, second in zip(a, b, strict=True))
    total = add(a, b)
    if multiply(difference, adjoint(difference), weight) != (field.element([998]), field.zero()):
        raise ValueError("the required difference-sector Gram failed")
    if multiply(total, adjoint(total), weight) != (field.element([1002]), field.zero()):
        raise ValueError("the current sum-sector Gram failed")
    _, _, s, _, _ = lattice_data()
    if np.any(block_projection(s - np.eye(500, dtype=np.int64))):
        raise ValueError("the conference generator violates the Paley lattice")
    return {
        "order": 1000, "phase": 8,
        "formula": "D0=[[A,B],[B,A]], A=C+zeta8 I, B=i(C+conjugate(zeta8) I)",
        "all_entries_are_roots": True,
        "full_Gram": "1000 I_1000 + 2 ([[0,1],[1,0]] tensor I_500)",
        "difference_sector_Gram": "998 I_500",
        "sum_sector_Gram": "1002 I_500",
        "required_sum_sector_Gram_not_yet_met": "998 I_500+4 J_500",
        "all_blocks_preserve_old_Paley_lattice_modulo_125": True,
        "lattice_proof": "every block and its adjoint is a linear combination of I,C; C=S-I commutes with P",
        "desired_Petrescu_corner_constructed": False,
    }


def small_direct_controls():
    q = paley_core(3)
    c = np.zeros((4, 4), dtype=np.int64)
    c[0, 1:] = 1
    c[1:, 0] = -1
    c[1:, 1:] = -q
    if not np.array_equal(c @ c.T, 3 * np.eye(4, dtype=np.int64)):
        raise ValueError("the small skew conference control failed")
    records = []
    for phase in (8, 16):
        source = conference_exponents(c, phase)
        for tensor in (False, True):
            actual = np.block([[source, source], [source, (source + phase // 2) % phase]]) if tensor else source
            n = len(actual)
            if not is_cgw_exponent_matrix(actual.tolist(), phase, n):
                raise ValueError("a complete direct small cyclotomic Gram failed")
            records.append({"order": n, "phase": phase, "H2_tensor": tensor,
                            "complete_entrywise_cyclotomic_Gram_checked": True})
    return records


@cache
def build_report():
    return {
        "schema": "even-petrescu-conference-boundary-v1",
        "scope": "Complete known-order cyclotomic boundary sources and one partial corner; no new real order",
        "formal_full_Gram_certificates": [formal_family_certificate(phase) for phase in (8, 16)],
        "new_boundary_kernel_certificates": [kernel_certificate(phase) for phase in (8, 16)],
        "complete_known_sources": [source_report(order, phase)
                                   for order in (1000, 2000) for phase in (8, 16)],
        "partial_corner": partial_corner_certificate(),
        "independent_small_direct_Gram_controls": small_direct_controls(),
        "new_normalized_boundary_requires_its_own_joint_lattice_equations": True,
        "parameter_need_not_be_a_root_for_unitarity": True,
        "novelty_of_the_conference_construction_claimed": False,
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
            raise ValueError("the conference family certificate differs from exact replay")
    else:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"complete_known_source_orders": [1000, 2000],
                      "genuine_phases": [8, 16], "partial_corner_order": 1000,
                      "higher_Petrescu_source_completed": False}))


if __name__ == "__main__":
    main()

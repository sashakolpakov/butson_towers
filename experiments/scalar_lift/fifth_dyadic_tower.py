"""Faithful phase20/40/80 partial scalar lifts with exact sign outputs.

These extend the verified fifth-root/Gaussian lift. Large compatible
Butson sources are still required; the included controls are known orders.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from functools import cache
from itertools import product
from pathlib import Path

import numpy as np
from sympy import Poly, Symbol, cyclotomic_poly

from experiments.scalar_lift.cyclotomic_towers import source_obstructions
from hadamard.cgw import is_cgw_exponent_matrix
from hadamard.dense_lift import FIFTH_POWER_QUATERNARY_BLOCKS


def sylvester(n):
    result = np.ones((1, 1), dtype=np.int64)
    while len(result) < n:
        result = np.block([[result, result], [result, -result]])
    if len(result) != n:
        raise ValueError("a Sylvester order must be a power of two")
    return result


def fifth_numerators():
    j = np.array([[0, -1], [1, 0]], dtype=np.int64)
    gaussian = (np.eye(2, dtype=np.int64), j, -np.eye(2, dtype=np.int64), -j)
    matrices = [np.block([[gaussian[e] for e in row] for row in exponents])
                for exponents in FIFTH_POWER_QUATERNARY_BLOCKS]
    identity = np.eye(8, dtype=np.int64)
    if not np.array_equal(matrices[0] @ matrices[0].T, 4 * identity):
        raise ValueError("the realified fifth-root action is not orthogonal")
    for r, s in product(range(1, 5), repeat=2):
        power = (r + s) % 5
        expected = 4 * identity if power == 0 else 2 * matrices[power - 1]
        if not np.array_equal(matrices[r - 1] @ matrices[s - 1], expected):
            raise ValueError("the stored fifth powers do not multiply exactly")
    if np.any(2 * identity + sum(matrices)):
        raise ValueError("Phi5 does not annihilate the realified action")
    return matrices


def matchings(indices):
    if not indices:
        yield []
        return
    first = indices[0]
    for second in indices[1:]:
        for rest in matchings([i for i in indices[1:] if i != second]):
            yield [(first, second), *rest]


def complex_structure_audit():
    action = fifth_numerators()[0]
    identity = np.eye(8, dtype=np.int64)
    tested, commuters = 0, []
    for matching in matchings(list(range(8))):
        for signs in product((-1, 1), repeat=4):
            j = np.zeros((8, 8), dtype=np.int64)
            for (a, b), sign in zip(matching, signs, strict=True):
                j[a, b], j[b, a] = sign, -sign
            tested += 1
            if not np.array_equal(j @ j, -identity):
                raise ValueError("a purported complex structure has the wrong square")
            if np.array_equal(j @ action, action @ j):
                encoding = []
                for row in j:
                    col = int(np.flatnonzero(row)[0])
                    encoding.append((col + 1) * int(row[col]))
                commuters.append(encoding)
    if tested != 1680 or len(commuters) != 6:
        raise ValueError("the complete signed-permutation complex-structure audit changed")
    return {"perfect_matchings": 105, "orientations_per_matching": 16,
            "total_complex_structures_tested": tested,
            "commuting_complex_structures": commuters,
            "commuting_count": len(commuters),
            "exhaustive_for_this_fixed_fifth_root_action": True}


@cache
def tower_level(phase):
    if type(phase) is not int or phase not in (20, 40, 80):
        raise ValueError("the verified tower levels have phase20,40,80")
    dyadic = phase // 5
    copies = dyadic // 4
    dimension = 2 * dyadic
    fifth = [np.kron(np.eye(copies, dtype=np.int64), a) for a in fifth_numerators()]
    r = np.kron(np.eye(4, dtype=np.int64), np.array([[0, -1], [1, 0]], dtype=np.int64))
    for _ in range(copies.bit_length() - 1):
        identity = np.eye(len(r), dtype=np.int64)
        zero = np.zeros_like(r)
        r = np.block([[zero, r], [identity, zero]])
    frame0 = np.kron(sylvester(copies), 2 * np.kron(np.eye(4, dtype=np.int64), sylvester(2)))
    numerator = fifth[0] @ frame0
    if np.any(numerator % 2):
        raise ValueError("the shifted flat frame is not integral")
    flat_frame = numerator // 2
    a, b = pow(5, -1, dyadic), pow(dyadic, -1, 5)
    generator = np.linalg.matrix_power(r, a) @ fifth[b - 1]
    identity = np.eye(dimension, dtype=np.int64)
    if (not np.array_equal(r @ r.T, identity)
            or not np.array_equal(np.linalg.matrix_power(r, dyadic // 2), -identity)
            or not np.array_equal(r @ fifth[0], fifth[0] @ r)
            or not np.array_equal(generator @ generator.T, 4 * identity)):
        raise ValueError("the exact dyadic/fifth commuting actions failed")
    if phase > 20:
        lower = tower_level(phase // 2)
        embedded = np.kron(np.eye(2, dtype=np.int64), lower["generator_numerator"])
        if not np.array_equal(generator @ generator, 2 * embedded):
            raise ValueError("the primitive generator does not square to the lower scalar action")
    if (not np.array_equal(frame0 @ frame0.T, dimension * identity)
            or not np.array_equal(flat_frame @ flat_frame.T, dimension * identity)
            or not np.all(np.abs(flat_frame) == 1)):
        raise ValueError("a frame has the wrong Gram or flatness")
    return {"phase": phase, "dyadic_order": dyadic, "dimension": dimension,
            "dyadic_action": r, "fifth_numerators": fifth,
            "generator_numerator": generator, "generator_denominator": 2,
            "dyadic_power": a, "fifth_power": b,
            "sparse_frame": frame0, "flat_frame": flat_frame}


def action_power(level, exponent):
    """Exact rho(zeta_phase^exponent) as an integer numerator and denominator."""
    r = np.linalg.matrix_power(level["dyadic_action"],
                               level["dyadic_power"] * exponent % level["dyadic_order"])
    fifth_power = level["fifth_power"] * exponent % 5
    if fifth_power == 0:
        return r, 1
    return r @ level["fifth_numerators"][fifth_power - 1], 2


@cache
def phase_blocks(phase, shifted_flat_frame=False):
    level = tower_level(phase)
    frame = level["flat_frame"] if shifted_flat_frame else level["sparse_frame"]
    blocks = []
    for e in range(phase):
        numerator, denominator = action_power(level, e)
        product_matrix = numerator @ frame
        if np.any(product_matrix % denominator):
            raise ValueError("a phase block is not integral")
        blocks.append(product_matrix // denominator)
    return blocks


def lift_source(source, phase, shifted_flat_frame=False):
    rows = [list(row) for row in source]
    if not rows or any(len(row) != len(rows) for row in rows):
        raise ValueError("a nonempty square source is required")
    if any(type(e) is not int for row in rows for e in row):
        raise ValueError("source exponents must be exact integers")
    blocks = phase_blocks(phase, shifted_flat_frame)
    rows = [[e % phase for e in row] for row in rows]
    if any(not np.all(np.abs(blocks[e]) == 1) for row in rows for e in row):
        raise ValueError("the source uses an excluded fifth-root coset")
    if not is_cgw_exponent_matrix(rows, phase, len(rows)):
        raise ValueError("the source does not satisfy its exact cyclotomic Gram identity")
    return np.block([[blocks[e] for e in row] for row in rows])


def level_certificate(phase):
    level = tower_level(phase)
    n, c = level["dimension"], level["dyadic_order"]
    identity = np.eye(n, dtype=np.int64)
    x = Symbol("x")
    polynomial = Poly(cyclotomic_poly(phase, x), x)
    total = np.zeros((n, n), dtype=np.int64)
    for (power,), coefficient in polynomial.terms():
        numerator, denominator = action_power(level, power)
        total += int(coefficient) * (2 // denominator) * numerator
    if np.any(total) or polynomial.degree() != n:
        raise ValueError("the full-degree cyclotomic polynomial identity failed")
    numerator, denominator = action_power(level, 5)
    if not np.array_equal(numerator, denominator * level["dyadic_action"]):
        raise ValueError("the generator's fifth power is not the dyadic scalar")
    numerator, denominator = action_power(level, c)
    if not np.array_equal(2 * numerator, denominator * level["fifth_numerators"][0]):
        raise ValueError("the generator does not recover the fifth-root scalar")
    blocks = phase_blocks(phase)
    shifted = phase_blocks(phase, True)
    allowed = [e for e, block in enumerate(blocks) if np.all(np.abs(block) == 1)]
    flat_allowed = [e for e, block in enumerate(shifted) if np.all(np.abs(block) == 1)]
    if allowed != [e for e in range(phase) if e % 5]:
        raise ValueError("the raw allowed alphabet is not exactly the four nontrivial fifth cosets")
    if flat_allowed != [e for e in range(phase) if (e + c) % 5]:
        raise ValueError("the shifted flat-frame alphabet has the wrong omitted coset")
    for e, block in enumerate(blocks):
        if not np.array_equal(block @ block.T, n * identity):
            raise ValueError("a scalar block has the wrong Gram matrix")
        if not np.array_equal(shifted[e], blocks[(e + c) % phase]):
            raise ValueError("the sparse and flat frames have inconsistent phase shifts")
        if not np.array_equal(level["generator_numerator"] @ block,
                              2 * blocks[(e + 1) % phase]):
            raise ValueError("successive powers do not follow the generator action")
        for f, other in enumerate(blocks):
            if not np.array_equal(block @ other.T, blocks[(e - f) % phase] @ blocks[0].T):
                raise ValueError("the scalar quotient-to-cross-Gram identity failed")
    # A known source, gauged by genuinely higher dyadic phases, exercises
    # both the fifth-root and every available dyadic scalar mechanism.
    source = [[(c * ((3 * (i + j)**2 + 1) % 5) + 5 * (i + j * j)) % phase
               for j in range(5)] for i in range(5)]
    sign_matrix = lift_source(source, phase)
    shifted_source = [[(e - c) % phase for e in row] for row in source]
    if not np.array_equal(lift_source(shifted_source, phase, True), sign_matrix):
        raise ValueError("the shifted flat frame changed the constructed matrix")
    order = 5 * n
    if (not np.all(np.abs(sign_matrix) == 1)
            or not np.array_equal(sign_matrix @ sign_matrix.T, order * np.eye(order, dtype=np.int64))
            or not np.array_equal(sign_matrix.T @ sign_matrix, order * np.eye(order, dtype=np.int64))):
        raise ValueError("the exact source-to-sign control failed")
    digest = hashlib.sha256(sign_matrix.astype(np.int8).tobytes()).hexdigest()
    return {"phase": phase, "dimension": n, "cyclotomic_degree": polynomial.degree(),
            "generator_numerator": level["generator_numerator"].tolist(),
            "generator_denominator": 2, "dyadic_action": level["dyadic_action"].tolist(),
            "sparse_frame": level["sparse_frame"].tolist(),
            "flat_frame": level["flat_frame"].tolist(),
            "Phi_phase_and_orthogonal_star_action_verified": True,
            "generator_fifth_power_is_dyadic_action": True,
            "generator_dyadic_power_is_fifth_action": True,
            "generator_square_is_two_copies_of_lower_generator": phase > 20,
            "allowed_exponents_sparse_frame": allowed,
            "allowed_exponents_flat_frame": flat_allowed,
            "flat_frame_omitted_residue_modulo5": -c % 5,
            "all_phase_block_Grams_verified": True,
            "all_cross_Gram_identities_verified": phase * phase,
            "control": {"source_order": 5, "source_exponents": source,
                        "real_order": order, "both_exact_integer_Grams_verified": True,
                        "same_output_with_shifted_flat_frame": True,
                        "sign_matrix_sha256": digest, "known_order_only": True}}


def build_report():
    return {"schema": "hadamard.fifth-dyadic-scalar-tower.v1",
            "complex_structure_audit": complex_structure_audit(),
            "levels": [level_certificate(q) for q in (20, 40, 80)],
            "high_source_route_audit": [{"source_order": p, "source_phase": q,
                                         "block_width": 2 * (q // 5),
                                         "real_order": 2 * (q // 5) * p,
                                         "source_witness_available": False,
                                         "source_obstructions": source_obstructions(p, q),
                                         "passes_checked_necessary_source_filters": not source_obstructions(p, q)}
                                        for p in (1499, 2999) for q in (20, 40, 80)],
            "large_Hadamard_constructed": False,
            "scope": "Explicit scalar actions and flat alphabets; missing large Butson sources remain necessary"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", type=Path)
    args = parser.parse_args()
    report = build_report()
    if args.check and report != json.loads(args.check.read_text()):
        raise ValueError("the exact fifth-dyadic tower certificate did not replay")
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"commuting_signed_complex_structures": report["complex_structure_audit"]["commuting_count"],
                      "levels": [{"phase": row["phase"], "width": row["dimension"],
                                  "flat_alphabet_size": len(row["allowed_exponents_flat_frame"]),
                                  "control_order": row["control"]["real_order"]}
                                 for row in report["levels"]]}))


if __name__ == "__main__":
    main()

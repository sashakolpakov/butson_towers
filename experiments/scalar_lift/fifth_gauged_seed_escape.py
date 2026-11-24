"""A known BH5 gauge giving a surjective binary boundary-compression map.

Only the small sources and necessary finite-field compression are solved.
The binary witness is not asserted to be E+conjugate(E) for a unitary E.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from itertools import permutations
from pathlib import Path

import numpy as np

from experiments.scalar_lift.fifth_source_boundary_audit import (
    FIFTH_ROOTS,
    MULTIPLY,
    binary_rank,
    chirped_five,
    dagger,
    field_power,
    matrix_multiply,
)
from hadamard.cgw import is_cgw_exponent_matrix


def gauged_five(first_row_phase=1):
    if first_row_phase not in (1, 2, 3):
        raise ValueError("this explicit family uses first-row phases1,2,3")
    row_phases = (first_row_phase, 1, 2, 3, 4)
    column_phases = (0, 0, 0, 0, 1)
    exponents = [[(row_phases[i] + i * j + column_phases[j]) % 5 for j in range(5)]
                 for i in range(5)]
    if any(e == 0 for row in exponents for e in row) or not is_cgw_exponent_matrix(exponents, 5, 5):
        raise ValueError("the exact admitted BH5 gauge failed")
    reduction = np.array([[FIFTH_ROOTS[e] for e in row] for row in exponents], dtype=np.uint8)
    return exponents, reduction


def pack(matrix):
    return sum(int(value) << (4 * j) for j, value in enumerate(np.asarray(matrix).flat))


def compression_basis(rows):
    rows = np.asarray(rows, dtype=np.uint8)
    pivots, independent = {}, []
    for i, a in enumerate(rows):
        for j, b in enumerate(rows):
            value = pack(MULTIPLY[a[:, None], b[None, :]])
            proof = 1 << (i * len(rows) + j)
            while value:
                leading = value.bit_length() - 1
                if leading in pivots:
                    old_value, old_proof = pivots[leading]
                    value ^= old_value
                    proof ^= old_proof
                else:
                    pivots[leading] = (value, proof)
                    independent.append([i, j])
                    break
    return pivots, independent


def compression_witness(rows, target):
    pivots, independent = compression_basis(rows)
    value, proof = pack(target), 0
    while value:
        leading = value.bit_length() - 1
        if leading not in pivots:
            return None, independent
        old_value, old_proof = pivots[leading]
        value ^= old_value
        proof ^= old_proof
    witness = np.array([[(proof >> (i * len(rows) + j)) & 1 for j in range(len(rows))]
                        for i in range(len(rows))], dtype=np.uint8)
    actual = matrix_multiply(matrix_multiply(rows.T, witness), rows)
    if not np.array_equal(actual, target):
        raise ValueError("the recovered binary compression witness failed")
    return witness, independent


def symmetric_completion_audit(rows, d):
    """A direct dual functional excludes BT=T(D+D*) for binary B."""
    images = matrix_multiply(rows, d ^ dagger(d))
    packed_rows = [pack(row) for row in rows]
    rank = binary_rank(packed_rows)
    augmented_rank = binary_rank(packed_rows + [pack(row) for row in images])
    bad_index = next(i for i, row in enumerate(images)
                     if binary_rank(packed_rows + [pack(row)]) > rank)
    equations = [(value, 0) for value in packed_rows] + [(pack(images[bad_index]), 1)]
    pivots = {}
    for mask, rhs in equations:
        while mask:
            leading = mask.bit_length() - 1
            if leading in pivots:
                old_mask, old_rhs = pivots[leading]
                mask ^= old_mask
                rhs ^= old_rhs
            else:
                pivots[leading] = (mask, rhs)
                break
        if not mask and rhs:
            raise ValueError("the requested separating functional was inconsistent")
    dual = 0
    for leading in sorted(pivots):
        mask, rhs = pivots[leading]
        if rhs ^ ((mask & dual).bit_count() % 2):
            dual |= 1 << leading
    if (any((row & dual).bit_count() % 2 for row in packed_rows)
            or (pack(images[bad_index]) & dual).bit_count() % 2 != 1):
        raise ValueError("the direct symmetry-exclusion dual certificate failed")
    return {"binary_row_space_dimension": rank,
            "dimension_after_adjoining_T_D_plus_D_star": augmented_rank,
            "witness_representative_row": bad_index,
            "witness_image_row": images[bad_index].tolist(),
            "dual_binary_mask": dual,
            "dual_power_basis_masks_per_field_coordinate": [(dual >> (4 * j)) & 15
                                                             for j in range(5)],
            "dual_vanishes_on_every_T_row": True,
            "dual_equals_one_on_witness_image": True,
            "symmetric_E_completion_excluded": True,
            "nonsymmetric_E_completion_excluded": False}


def one_gauge_report(first_row_phase):
    exponents, m = gauged_five(first_row_phase)
    d, gamma, rho = chirped_five(), FIFTH_ROOTS[3], FIFTH_ROOTS[1]
    d2 = matrix_multiply(d, d)
    # 141/291 ordinary copies: odd; one gamma D; one M; 32 B8 rows.
    gram = MULTIPLY[1 ^ field_power(gamma, 2), d2] ^ matrix_multiply(m.T, m)
    target = matrix_multiply(d, gram) ^ matrix_multiply(gram, dagger(d))
    rows = np.vstack([d, MULTIPLY[gamma, d], m, np.full((1, 5), rho, dtype=np.uint8)])
    witness, independent = compression_witness(rows, target)
    if witness is None or len(independent) != 100:
        raise ValueError("the explicit gauge did not give a surjective compression")
    return {"first_row_phase": first_row_phase, "M_exponents_mod5": exponents,
            "M_exact_BH5_Gram_verified": True, "M_avoids_fifth_root_one": True,
            "M_row_sums_mod2": np.bitwise_xor.reduce(m, axis=1).tolist(),
            "representative_rows": rows.tolist(), "ordinary_boundary_Gram": gram.tolist(),
            "compression_target": target.tolist(),
            "compression_F2_rank": len(independent), "compression_F2_codomain_dimension": 100,
            "independent_generator_row_pairs": independent,
            "binary_compression_witness": witness.tolist(),
            "full_compression_equation_verified": True,
            "symmetric_completion_preflight": symmetric_completion_audit(rows, d),
            "witness_realized_as_relative_trace_of_unitary_E": False}


def column_permutation_audit():
    d, gamma, rho = chirped_five(), FIFTH_ROOTS[3], FIFTH_ROOTS[1]
    d2 = matrix_multiply(d, d)
    counts = Counter()
    for permutation in permutations(range(5)):
        p = np.eye(5, dtype=np.uint8)[:, permutation]
        m = matrix_multiply(d, p)
        rows = np.vstack([d, MULTIPLY[gamma, d], m, np.full((1, 5), rho, dtype=np.uint8)])
        gram = MULTIPLY[1 ^ field_power(gamma, 2), d2] ^ matrix_multiply(m.T, m)
        target = matrix_multiply(d, gram) ^ matrix_multiply(gram, dagger(d))
        witness, independent = compression_witness(rows, target)
        if witness is not None:
            raise ValueError("a pure column permutation escaped the common row-sum obstruction")
        counts[len(independent)] += 1
    if dict(counts) != {75: 10, 91: 50, 99: 60}:
        raise ValueError("the complete column-permutation audit changed")
    return {"permutations_checked": 120, "all_fail_binary_compression": True,
            "compression_rank_counts": {str(k): v for k, v in sorted(counts.items())}}


def build_report():
    return {
        "scope": "Known small-source gauges and a necessary binary compression only",
        "source_boundary_recipe": {
            "source1499": {"ordinary_D_copies": 141, "gamma_D_copies": 1,
                           "M_copies": 1, "partial_BH8_copies": 4},
            "source2999": {"ordinary_D_copies": 291, "gamma_D_copies": 1,
                           "M_copies": 1, "partial_BH8_copies": 4},
            "partial_BH8_fifth_cosets": [1, 1, 1, 1],
            "all_column_Grams_preserved_by_replacing_D_with_M": True,
        },
        "explicit_passing_gauge_family": [one_gauge_report(r) for r in (1, 2, 3)],
        "pure_column_permutation_control": column_permutation_audit(),
        "simultaneous_unitary_residue_completion_constructed": False,
        "characteristic_zero_source_constructed": False,
        "new_real_Hadamard_matrix_constructed": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", type=Path)
    args = parser.parse_args()
    report = build_report()
    if args.check and json.loads(args.check.read_text()) != report:
        raise SystemExit("gauged-seed certificate differs from replay")
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"explicit_surjective_gauges": 3, "compression_rank": 100,
                      "unitary_completion_constructed": False}))


if __name__ == "__main__":
    main()

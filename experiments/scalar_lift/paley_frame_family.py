"""Exact modular moment certificates for the order-668 tensor-operator family.

Every symmetric integral B=I mod2 with B^2=501I is covered, together with
every skew sign-conference matrix D of order four. The proof requires only
the skew conference-core identities of Q, not its cyclic Paley realization.
"""

from __future__ import annotations

from collections import Counter
from itertools import product

import numpy as np

from experiments.scalar_lift.paley_coefficient_forms import normalized_coefficient_forms
from experiments.scalar_lift.paley_sextic_operator import quaternion_orientation_audit

PRIME = 167
CERTIFICATE_PRIME = 1009
SIGNS = np.array(list(product((-1, 1), repeat=4)), dtype=np.int64)
PAIRS = tuple((i, j) for i in range(4) for j in range(i, 4))


def required_antipodal_sums(b, d):
    """All possible frame sums, selecting the positive first coordinate."""
    sums = []
    for large in range(4):
        for signs in SIGNS:
            k = signs.copy()
            k[large] *= 3
            numerator = b @ k
            if np.any(numerator % 3):
                continue
            sigma = numerator // 3
            if sigma[0] < 0 or np.any((sigma + d @ sigma) % 4 != 2):
                continue
            if int(sigma @ sigma) != 4 * PRIME or np.any(sigma % 2 != 1):
                raise ValueError("the extremal block sum has an incorrect norm or parity")
            sums.append(sigma)
    return sums


def local_count_system(b, d, sigma):
    """Independent integer-array assembly of the 33 necessary moments.

    Rows: count; four v sums; four y sums; four Qv sums; ten symmetric
    v/Qv mixed moments; ten Qv Gram moments. No positivity assumption on
    the unknown counts is needed for the obstruction.
    """
    p = PRIME
    offset = b @ d @ sigma
    if np.any(offset % p):
        raise ValueError("a required sum has nonintegral BD sigma / p")
    offset //= p
    columns = []
    for v in SIGNS:
        for y in SIGNS:
            z = b @ (2 * y - v) - offset
            if np.any(z % 6):
                continue
            z //= 3
            columns.append(
                [1]
                + v.tolist()
                + y.tolist()
                + z.tolist()
                + [int(v[i] * z[j] + v[j] * z[i]) for i, j in PAIRS]
                + [int(z[i] * z[j] - p * v[i] * v[j]) for i, j in PAIRS]
            )
    rhs = (
        [p]
        + sigma.tolist()
        + ((sigma + d @ sigma) // 2).tolist()
        + [0] * 14
        + [-int(sigma[i] * sigma[j]) for i, j in PAIRS]
    )
    matrix = np.array(columns, dtype=np.int64).T
    if not columns:
        matrix = np.zeros((33, 0), dtype=np.int64)
    return matrix, np.array(rhs, dtype=np.int64)


def modular_inconsistency(matrix, rhs, modulus=CERTIFICATE_PRIME):
    """Return an exactly checked left-null witness, or None if consistent."""
    matrix = np.array(matrix, dtype=np.int64)
    rhs = np.array(rhs, dtype=np.int64)
    rows, columns = matrix.shape
    augmented = np.column_stack((matrix, rhs)) % modulus
    transform = np.eye(rows, dtype=np.int64)
    rank = 0
    for column in range(columns):
        nonzero = np.flatnonzero(augmented[rank:, column])
        if not len(nonzero):
            continue
        pivot = rank + int(nonzero[0])
        augmented[[rank, pivot]] = augmented[[pivot, rank]]
        transform[[rank, pivot]] = transform[[pivot, rank]]
        inverse = pow(int(augmented[rank, column]), -1, modulus)
        augmented[rank] = augmented[rank] * inverse % modulus
        transform[rank] = transform[rank] * inverse % modulus
        multipliers = augmented[rank + 1 :, column].copy()
        augmented[rank + 1 :] = (
            augmented[rank + 1 :] - multipliers[:, None] * augmented[rank]
        ) % modulus
        transform[rank + 1 :] = (
            transform[rank + 1 :] - multipliers[:, None] * transform[rank]
        ) % modulus
        rank += 1
        if rank == rows:
            break
    for row in range(rank, rows):
        if augmented[row, -1]:
            witness = transform[row]
            residue = int(witness @ rhs % modulus)
            if np.any(witness @ matrix % modulus) or residue == 0:
                raise ValueError("the modular elimination produced an invalid certificate")
            return rank, {"left_null_vector": witness.tolist(), "rhs_residue": residue}
    return rank, None


def audit_one_operator(b, d):
    b = np.array(b, dtype=np.int64)
    if not np.array_equal(b @ b, 3 * PRIME * np.eye(4, dtype=np.int64)):
        raise ValueError("the Clifford square failed")
    sigma_values = required_antipodal_sums(b, d)
    local_records, survivors = [], []
    for sigma in sigma_values:
        offset_residues = b @ d @ sigma % PRIME
        if np.any(offset_residues):
            local_records.append(
                {
                    "sigma": sigma.tolist(),
                    "excluded_by_nonintegral_BD_sigma_over_p": True,
                    "BD_sigma_residues_modulo_167": offset_residues.tolist(),
                }
            )
            continue
        matrix, rhs = local_count_system(b, d, sigma)
        rank, witness = modular_inconsistency(matrix, rhs)
        local_records.append(
            {
                "sigma": sigma.tolist(),
                "local_state_count": matrix.shape[1],
                "coefficient_rank_modulo_1009": rank,
                "inconsistency_certificate": witness,
            }
        )
        if witness is None:
            survivors.append(sigma)
    # Counts of antipodal directions must satisfy the full-frame Gram moment.
    moments = np.array(
        [[int(sigma[i] * sigma[j]) for i, j in PAIRS] for sigma in survivors],
        dtype=np.int64,
    ).T
    if not survivors:
        moments = np.zeros((10, 0), dtype=np.int64)
    target = np.array([4 * PRIME**2 if i == j else 0 for i, j in PAIRS], dtype=np.int64)
    _, global_witness = modular_inconsistency(moments, target)
    return {
        "B": b.tolist(),
        "D_upper_triangle": [int(d[i, j]) for i in range(4) for j in range(i + 1, 4)],
        "local_directions": local_records,
        "remaining_antipodal_directions": [sigma.tolist() for sigma in survivors],
        "global_frame_moment_inconsistency_certificate": global_witness,
        "full_frame_excluded": global_witness is not None,
    }


def build_paley_frame_family_report():
    ds = [np.array(row["D"], dtype=np.int64) for row in quaternion_orientation_audit()]
    classification = normalized_coefficient_forms()
    cases = [audit_one_operator(b, d) for b in classification["matrices"] for d in ds]
    if len(cases) != 128 or not all(case["full_frame_excluded"] for case in cases):
        raise ValueError("the stated exhaustive tensor-operator exclusion failed")
    histogram = Counter(
        (
            row["local_state_count"],
            row["coefficient_rank_modulo_1009"],
            row["inconsistency_certificate"] is not None,
        )
        for case in cases
        for row in case["local_directions"]
        if "local_state_count" in row
    )
    return {
        "scope": "C=B tensor Q+D tensor J at p=167, integral symmetric B=I mod2, B^2=501I, skew sign-conference D",
        "prime": PRIME,
        "certificate_modulus": CERTIFICATE_PRIME,
        "coefficient_classification": classification,
        "all_integral_symmetric_B_in_scope_included_up_to_signed_permutation": True,
        "rational_coefficients_of_a_flat_frame_are_forced_into_this_integral_family": True,
        "all_skew_sign_conference_D_included": True,
        "normalized_B_count": len(classification["matrices"]),
        "D_count": len(ds),
        "operator_count": len(cases),
        "local_antipodal_direction_count": sum(len(case["local_directions"]) for case in cases),
        "directions_excluded_by_nonintegral_BD_sigma_over_p": sum(
            "excluded_by_nonintegral_BD_sigma_over_p" in row
            for case in cases
            for row in case["local_directions"]
        ),
        "local_moment_histogram": [
            {"state_count": key[0], "rank": key[1], "inconsistent": key[2], "count": value}
            for key, value in sorted(histogram.items())
        ],
        "all_operators_excluded_by_exact_modular_certificates": True,
        "all_skew_sign_conference_cores_Q_with_Q_squared_J_minus_167I_included": True,
        "other_scalar_operator_templates_excluded": False,
        "cases": cases,
        "new_real_hadamard_order_claimed": False,
    }

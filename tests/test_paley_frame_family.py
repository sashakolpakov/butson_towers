from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from experiments.scalar_lift.paley_frame_family import (
    CERTIFICATE_PRIME,
    PAIRS,
    PRIME,
    build_paley_frame_family_report,
    local_count_system,
)
from experiments.scalar_lift.paley_frame_moments import moment_system
from experiments.scalar_lift.paley_sextic_operator import D_COMPATIBLE, B


def test_independent_numpy_and_symbolic_moment_assemblies_agree():
    b, d = np.array(B), np.array(D_COMPATIBLE)
    sigma = np.array([7, -3, 23, -9])
    matrix, rhs = local_count_system(b, d, sigma)
    reference = moment_system(PRIME, B, D_COMPATIBLE, sigma)
    permutation = list(range(13)) + list(range(14, 33, 2)) + list(range(13, 33, 2))
    assert np.array_equal(matrix, np.array(reference["coefficient_matrix"])[permutation])
    assert np.array_equal(rhs, np.array(reference["right_hand_side"])[permutation])


def test_saved_certificates_directly_exclude_every_normalized_operator():
    path = Path(__file__).resolve().parents[1] / "certificates/scalar-lift/paley-frame-family-668.json"
    report = json.loads(path.read_text())
    replay = json.loads(json.dumps(build_paley_frame_family_report()))
    assert replay == report
    local_obstructions = 0
    for case in report["cases"]:
        b = np.array(case["B"], dtype=np.int64)
        d = np.zeros((4, 4), dtype=np.int64)
        pairs = [(i, j) for i in range(4) for j in range(i + 1, 4)]
        for (i, j), sign in zip(pairs, case["D_upper_triangle"], strict=True):
            d[i, j], d[j, i] = sign, -sign
        survivors = []
        for row in case["local_directions"]:
            sigma = np.array(row["sigma"], dtype=np.int64)
            if row.get("excluded_by_nonintegral_BD_sigma_over_p"):
                assert np.any(b @ d @ sigma % PRIME)
                continue
            matrix, rhs = local_count_system(b, d, sigma)
            certificate = row["inconsistency_certificate"]
            if certificate is None:
                survivors.append(sigma)
            else:
                witness = np.array(certificate["left_null_vector"], dtype=np.int64)
                assert np.all(witness @ matrix % CERTIFICATE_PRIME == 0)
                assert int(witness @ rhs % CERTIFICATE_PRIME) == certificate["rhs_residue"] != 0
                local_obstructions += 1
        witness = np.array(
            case["global_frame_moment_inconsistency_certificate"]["left_null_vector"],
            dtype=np.int64,
        )
        for sigma in survivors:
            moment = np.array([int(sigma[i] * sigma[j]) for i, j in PAIRS])
            assert witness @ moment % CERTIFICATE_PRIME == 0
        target = np.array([4 * PRIME**2 if i == j else 0 for i, j in PAIRS])
        assert witness @ target % CERTIFICATE_PRIME != 0
    assert len(report["cases"]) == 128 and local_obstructions == 64

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from experiments.scalar_lift.order_716_family import (
    PAIRS,
    build_order_716_family_report,
    coefficient_classification_716,
    moment_arrays_716,
)
from experiments.scalar_lift.paley_frame_moments import moment_system
from experiments.scalar_lift.paley_sextic_operator import quaternion_orientation_audit


def test_generalized_integer_moments_match_independent_symbolic_builder():
    b = np.array(coefficient_classification_716()["normalized_matrices"][0])
    d = np.array(quaternion_orientation_audit()[3]["D"])
    sigma = np.array([1, -7, 21, 15])
    matrix, rhs = moment_arrays_716(b, d, sigma)
    independent = moment_system(179, b.tolist(), d.tolist(), sigma.tolist())
    labels = independent["row_labels"]
    order = list(range(13)) + [labels.index(f"skew{i}{j}") for i, j in PAIRS]
    order += [labels.index(f"gram{i}{j}") for i, j in PAIRS]
    assert np.array_equal(matrix, np.array(independent["coefficient_matrix"])[order])
    assert np.array_equal(rhs, np.array(independent["right_hand_side"])[order])


def test_full_coefficient_audit_replays_certificates_and_preserves_survivors():
    report = build_order_716_family_report()
    root = Path(__file__).resolve().parents[1]
    saved = json.loads((root / "certificates/scalar-lift/order-716-family.json").read_text())
    assert saved == json.loads(json.dumps(report))
    assert report["case_count"] == 352
    assert report["excluded_case_count"] == 324
    assert report["surviving_case_count"] == 28
    survivors = [case for case in report["cases"] if not case["full_frame_excluded_by_these_checks"]]
    assert len(survivors) == 28
    for case in survivors:
        b, d, e = map(np.array, (case["B"], case["D"], case["integral_E_with_B_equals_D_E"]))
        assert np.array_equal(b, d @ e)
        assert np.array_equal(d @ e, e @ d)
        assert np.array_equal(e @ e, -179 * np.eye(4, dtype=np.int64))
    forms = report["survivor_normal_forms"]
    assert forms["simultaneous_signed_permutation_orbit_count"] == 10
    assert forms["Paley_frame_existence_class_count"] == 3
    for group in forms["classes"]:
        target = survivors[group["representative_survivor_index"]]
        for witness in group["members_and_equivalence_witnesses"]:
            source = survivors[witness["survivor_index"]]
            permutation, signs = witness["permutation"], witness["coordinate_signs"]
            factors = np.outer(signs, signs)
            for name in ("B", "D"):
                transformed = np.array(source[name])[np.ix_(permutation, permutation)] * factors
                assert np.array_equal(witness[f"global_{name}_sign"] * transformed, target[name])
    assert not report["flat_frame_constructed"]

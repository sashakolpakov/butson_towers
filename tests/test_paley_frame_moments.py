from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import numpy as np
import pytest

from experiments.scalar_lift.paley_frame_moments import (
    LEFT_NULL_CERTIFICATE,
    build_paley_frame_moment_report,
    moment_system,
    verify_left_null_certificate,
)
from experiments.scalar_lift.paley_sextic_operator import D_COMPATIBLE, B


def test_integer_certificate_rejects_the_required_direction_and_detects_corruption():
    system = moment_system(167, B, D_COMPATIBLE, [7, -3, 23, -9])
    assert verify_left_null_certificate(system, LEFT_NULL_CERTIFICATE) == 16032
    corrupted = list(LEFT_NULL_CERTIFICATE)
    corrupted[0] += 1
    with pytest.raises(ValueError, match="annihilate"):
        verify_left_null_certificate(system, corrupted)


def test_actual_nonconstant_flat_column_satisfies_all_small_prime_moments():
    # A direct 12-dimensional positive control, with no moment solver involved.
    p = 3
    b = np.array([[1, 2, 2, 0], [2, -1, 0, 2], [2, 0, -1, -2], [0, 2, -2, 1]])
    d = np.array(D_COMPATIBLE)
    q = np.array([[0, 1, -1], [-1, 0, 1], [1, -1, 0]])
    c = np.kron(b, q) + np.kron(d, np.ones((p, p), dtype=int))
    v = np.array([-1, -1, -1, -1, -1, 1, -1, 1, 1, 1, 1, -1])
    numerator = p * v + c @ v
    assert np.all(np.abs(numerator) == 2 * p)
    x, y = v.reshape(4, p), (numerator // (2 * p)).reshape(4, p)
    z = (q @ x.T).T
    sigma = x.sum(axis=1)
    assert np.any(np.abs(sigma) != p)
    system = moment_system(p, b.tolist(), d.tolist(), sigma.tolist())
    observed = Counter(
        tuple(tuple(int(a) for a in block[:, j]) for block in (x, y, z)) for j in range(p)
    )
    assert set(observed).issubset(set(system["local_states"]))
    counts = np.array([observed[state] for state in system["local_states"]])
    assert np.array_equal(
        np.array(system["coefficient_matrix"]) @ counts, system["right_hand_side"]
    )


def test_saved_moment_certificate_replays_exactly():
    root = Path(__file__).resolve().parents[1]
    path = root / "certificates/scalar-lift/paley-frame-moments-668.json"
    assert json.loads(json.dumps(build_paley_frame_moment_report())) == json.loads(path.read_text())

"""Exact sign criterion for a commuting quaternion/Paley scalar tower.

This characterizes individual flat columns. It does not construct an
orthogonal frame, or assert that the required columns exist.
"""

from __future__ import annotations

from itertools import product

import numpy as np

FIRST_B = ((1, 4, 14, 18), (4, 21, 4, -8), (14, 4, -15, 10), (18, -8, 10, -7))
FIRST_D = ((0, -1, -1, 1), (1, 0, 1, 1), (1, -1, 0, -1), (-1, -1, 1, 0))


def paley_core(prime):
    """Build Q_ij=chi(j-i); verify the identities used in the argument."""
    character = np.array(
        [0] + [1 if pow(a, (prime - 1) // 2, prime) == 1 else -1 for a in range(1, prime)],
        dtype=np.int64,
    )
    indices = np.arange(prime)
    q = character[(indices[None, :] - indices[:, None]) % prime]
    if not np.array_equal(q.T, -q) or not np.array_equal(
        q @ q, np.ones((prime, prime), dtype=np.int64) - prime * np.eye(prime, dtype=np.int64)
    ) or np.any(q.sum(axis=0)):
        raise ValueError("the proposed core fails the skew conference identities")
    return q


def cube_switch(signs, d):
    """T(v)_ij=-kappa*product(v_kj:k!=i), for a skew sign D of order four."""
    signs, d = np.asarray(signs, dtype=np.int64), np.asarray(d, dtype=np.int64)
    if signs.ndim != 2 or signs.shape[0] != 4 or not np.all(np.abs(signs) == 1):
        raise ValueError("expected four rows of signs")
    if d.shape != (4, 4) or not np.array_equal(d.T, -d) or not np.array_equal(
        d @ d, -3 * np.eye(4, dtype=np.int64)
    ) or not np.all(np.abs(d[~np.eye(4, dtype=bool)]) == 1):
        raise ValueError("D must be a skew sign conference matrix of order four")
    kappas = [int(np.prod(np.delete(d[i], i))) for i in range(4)]
    if len(set(kappas)) != 1:
        raise ValueError("the row-product signs disagree")
    return -kappas[0] * signs.prod(axis=0)[None, :] * signs


def commuting_factors():
    b, d = np.array(FIRST_B, dtype=np.int64), np.array(FIRST_D, dtype=np.int64)
    numerator = -d @ b
    if np.any(numerator % 3):
        raise ValueError("the quaternion factor is not integral")
    e = numerator // 3
    unit = np.eye(4, dtype=np.int64)
    if not np.array_equal(e.T, -e) or not np.array_equal(e @ e, -179 * unit):
        raise ValueError("the quaternion norm identity failed")
    if not np.array_equal(d @ e, b) or not np.array_equal(e @ d, b):
        raise ValueError("the factors do not commute")
    return b, d, e


def involution_numerator(signs, e, q):
    """Compute p O v, with four blocks stored as rows of a 4 by p array."""
    signs = np.asarray(signs, dtype=np.int64)
    return e @ (signs @ q.T) + signs.sum(axis=1)[:, None]


def reconstruct_flat_image(signs, e, d, q):
    """Return the exact sign image when the cubic criterion holds, else None.

    The caller must verify the commuting factor identities E^2=-p I,
    D^2=-3 I, [D,E]=0 and the skew conference identities for Q. Necessity
    uses p an odd prime other than three. The sufficiency check is direct.
    """
    signs = np.asarray(signs, dtype=np.int64)
    switched = cube_switch(signs, d)
    prime = len(q)
    if signs.shape != (4, prime):
        raise ValueError("the sign blocks have the wrong length")
    if not np.array_equal(involution_numerator(signs, e, q), prime * switched):
        return None
    numerator = signs + d @ switched
    if np.any(numerator % 2):
        raise ValueError("the image is not integral")
    image = numerator // 2
    if not np.all(np.abs(image) == 1):
        raise ValueError("the local cube identity failed")
    return image


def build_commuting_parity_report():
    """Replay all local cube states and the first order-716 factorization."""
    from experiments.scalar_lift.order_716_family import required_sums_716
    from experiments.scalar_lift.paley_sextic_operator import quaternion_orientation_audit

    signs = np.array(list(product((-1, 1), repeat=4)), dtype=np.int64).T
    orientation_checks = []
    for record in quaternion_orientation_audit():
        d = np.array(record["D"], dtype=np.int64)
        switched = cube_switch(signs, d)
        numerator = signs + d @ switched
        if not np.all(np.abs(numerator) == 2):
            raise ValueError("a proposed local sign image failed")
        # Uniqueness: for a fixed w in the cube, exactly one v in the cube
        # makes (v+D w)/2 a sign vector.
        for w in signs.T:
            valid = [v.tolist() for v in signs.T if np.all(np.abs(v + d @ w) == 2)]
            if valid != [cube_switch(w[:, None], d)[:, 0].tolist()]:
                raise ValueError("the local parity switch is not the unique solution")
        orientation_checks.append({"D": d.tolist(), "verified_local_sign_states": 16})
    b, d, e = commuting_factors()
    q = paley_core(179)
    controls = []
    for local in signs.T:
        v = np.repeat(local[:, None], 179, axis=1)
        image = reconstruct_flat_image(v, e, d, q)
        if image is not None:
            controls.append({"constant_input": local.tolist(), "constant_image": image[:, 0].tolist()})
    if len(controls) != 8:
        raise ValueError("the constant-column control changed")
    directions = []
    for sigma in required_sums_716(b, d):
        numerator = e @ sigma
        if np.any(numerator % 179):
            raise ValueError("the required direction has no integral border")
        h = numerator // 179
        if not np.all(np.abs(h) == 1) or not np.array_equal(sigma, -e @ h):
            raise ValueError("the required border is not a sign vector")
        directions.append({"sigma": sigma.tolist(), "border_signs_h": h.tolist()})
    return {
        "scope": "Exact individual-column characterization for a commuting scalar tower",
        "prime": 179,
        "B": b.tolist(), "D": d.tolist(), "E": e.tolist(),
        "local_orientation_checks": orientation_checks,
        "criterion": "O v=T(v); T_i(v_j)=-product(v_kj:k!=i) for this D",
        "equivalent_required_direction_equation": "Q v=h-E T(v), h=E sigma/179",
        "required_directions": directions,
        "constant_column_controls": controls,
        "constant_controls_have_required_extremal_block_sums": False,
        "required_flat_column_constructed": False,
        "full_frame_constructed": False,
    }

"""Exact moment obstruction to the previously surviving Paley sextic frame.

The certificate is an integer linear dependence, not a numerical solver
status.  It excludes even one flat column with a block sum that any full
frame for the fixed operator is required to use.
"""

from __future__ import annotations

from itertools import combinations_with_replacement, product

from sympy import Matrix

from experiments.scalar_lift.paley_sextic_operator import (
    D_COMPATIBLE,
    PRIME,
    B,
    necessary_regular_block_sum_audit,
)

# Row labels are generated below; zero coefficients retain the exact ordering.
LEFT_NULL_CERTIFICATE = (
    -777,
    288,
    80,
    276,
    12,
    2976,
    -9952,
    -512,
    5140,
    0,
    0,
    0,
    0,
    178,
    -24,
    -10,
    630,
    -6,
    2590,
    8,
    -160,
    186,
    -102,
    6,
    54,
    -16,
    2700,
    177,
    -41,
    -6,
    -514,
    184,
    0,
)


def moment_system(prime, b, d, sigma):
    """Necessary exact count equations for signs v and Rv with block sums sigma.

    A local state is (a,y,z), where a and y are the four signs at one
    cyclic coordinate, and z contains the corresponding four Qv entries.
    Each Qv entry is even because Q has prime-1 nonzero sign entries.
    """
    b, d, sigma = Matrix(b), Matrix(d), Matrix(sigma)
    unit = Matrix.eye(4)
    if b != b.T or b * b != 3 * prime * unit:
        raise ValueError("B must be symmetric with square 3p I")
    if d != -d.T or d * d != -3 * unit:
        raise ValueError("D must be skew with square -3 I")
    if sigma.shape != (4, 1) or any(int(x) % 2 != 1 for x in sigma):
        raise ValueError("four odd block sums are required")
    t_numerator = b * d * sigma
    tau_numerator = (unit + d) * sigma
    if any(x % prime for x in t_numerator) or any(x % 2 for x in tau_numerator):
        raise ValueError("the proposed block sums do not give integral t and tau")
    t, tau = t_numerator / prime, tau_numerator / 2
    sign_states = [Matrix(signs) for signs in product((-1, 1), repeat=4)]
    states = []
    for a in sign_states:
        for y in sign_states:
            numerator = b * (2 * y - a) - t
            if any(x % 6 for x in numerator):
                continue
            z = numerator / 3
            states.append(tuple(tuple(int(x) for x in block) for block in (a, y, z)))

    rows = [[1] * len(states)]
    rhs, labels = [prime], ["count"]
    for block, totals, name in ((0, sigma, "x"), (1, tau, "y"), (2, [0] * 4, "z")):
        for i in range(4):
            rows.append([state[block][i] for state in states])
            rhs.append(int(totals[i]))
            labels.append(f"{name}{i}")
    for i, j in combinations_with_replacement(range(4), 2):
        rows.append([z[i] * z[j] - prime * a[i] * a[j] for a, y, z in states])
        rhs.append(-int(sigma[i] * sigma[j]))
        labels.append(f"gram{i}{j}")
        rows.append([a[i] * z[j] + a[j] * z[i] for a, y, z in states])
        rhs.append(0)
        labels.append(f"skew{i}{j}")
    return {
        "prime": prime,
        "sigma": list(map(int, sigma)),
        "tau": list(map(int, tau)),
        "t": list(map(int, t)),
        "local_states": states,
        "row_labels": labels,
        "coefficient_matrix": rows,
        "right_hand_side": rhs,
    }


def verify_left_null_certificate(system, certificate):
    """Replay only integer arithmetic: f^T A=0 and f^T rhs is nonzero."""
    rows, rhs = system["coefficient_matrix"], system["right_hand_side"]
    if len(certificate) != len(rows):
        raise ValueError("the certificate has the wrong number of coefficients")
    products = [
        sum(f * row[j] for f, row in zip(certificate, rows, strict=True))
        for j in range(len(rows[0]))
    ]
    if any(products):
        raise ValueError("the alleged certificate does not annihilate every local state")
    contradiction = sum(f * value for f, value in zip(certificate, rhs, strict=True))
    if contradiction == 0:
        raise ValueError("the certificate has no nonzero contradiction")
    return contradiction


def build_paley_frame_moment_report():
    sums = necessary_regular_block_sum_audit()["required_block_sum_matrix_M"]
    sigma = [row[0] for row in sums]
    system = moment_system(PRIME, B, D_COMPATIBLE, sigma)
    contradiction = verify_left_null_certificate(system, LEFT_NULL_CERTIFICATE)
    if len(system["local_states"]) != 28 or contradiction != 16032:
        raise ValueError("the exact obstruction has changed")
    matrix = Matrix(system["coefficient_matrix"])
    rhs = Matrix(system["right_hand_side"])
    if matrix.rank() != 27 or matrix.row_join(rhs).rank() != 28:
        raise ValueError("the independent rational rank check disagrees")
    return {
        "scope": "No full sign frame for the fixed B,D_COMPATIBLE Paley sextic operator",
        "order": 4 * PRIME,
        "B": B,
        "D": D_COMPATIBLE,
        "required_regular_block_sum_matrix_M": sums,
        "required_columns_with_this_antipodal_block_sum_direction": PRIME,
        "moment_system": system,
        "left_null_certificate": LEFT_NULL_CERTIFICATE,
        "left_null_product_zero_for_every_local_state": True,
        "left_null_product_with_rhs": contradiction,
        "coefficient_matrix_rational_rank": 27,
        "augmented_matrix_rational_rank": 28,
        "a_flat_column_with_this_sigma_is_impossible": True,
        "a_flat_column_with_negative_sigma_is_impossible_by_sign_reversal": True,
        "full_sign_frame_excluded_for_this_fixed_operator": True,
        "other_operators_or_B_matrices_are_not_excluded": True,
        "new_real_hadamard_order_claimed": False,
    }

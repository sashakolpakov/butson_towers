"""Exact dense sextic operators and obstructions to their sign frames.

The operator meets the denominator and dyadic rank conditions for a real
order-668 scalar lift.  Its flat vectors span only four dimensions, so this
particular Paley/Clifford construction cannot provide the lift. The opposite
orientation is excluded by a separate, replayed integer moment certificate.
"""

from __future__ import annotations

from itertools import product

import numpy as np

PRIME = 167
B = ((1, 4, 22, 0), (4, -1, 0, 22), (22, 0, -1, -4), (0, 22, -4, 1))
D = ((0, 1, 1, 1), (-1, 0, 1, -1), (-1, -1, 0, 1), (-1, 1, -1, 0))
D_COMPATIBLE = ((0, 1, 1, 1), (-1, 0, -1, 1), (-1, 1, 0, -1), (-1, -1, 1, 0))


def rank_mod_prime(matrix, prime):
    """Exact row elimination; products fit int64 for the fixed prime 167."""
    a = np.array(matrix, dtype=np.int64) % prime
    rank = 0
    for column in range(a.shape[1]):
        nonzero = np.flatnonzero(a[rank:, column])
        if not len(nonzero):
            continue
        pivot = rank + int(nonzero[0])
        a[[rank, pivot]] = a[[pivot, rank]]
        a[rank] = a[rank] * pow(int(a[rank, column]), -1, prime) % prime
        a[rank + 1 :] = (
            a[rank + 1 :] - a[rank + 1 :, column, None] * a[rank]
        ) % prime
        rank += 1
        if rank == len(a):
            break
    return rank


def jacobsthal_matrix():
    p = PRIME
    character = np.array(
        [0] + [1 if pow(a, (p - 1) // 2, p) == 1 else -1 for a in range(1, p)],
        dtype=np.int64,
    )
    indices = np.arange(p)
    return character[(indices[None, :] - indices[:, None]) % p]


def dense_sextic_operator_numerator(d=D):
    """Return C with C^2=-3*167^2 I for a skew sign-conference matrix d."""
    return np.kron(np.array(B, dtype=np.int64), jacobsthal_matrix()) + np.kron(
        np.array(d, dtype=np.int64), np.ones((PRIME, PRIME), dtype=np.int64)
    )


def flat_sign_images(c, signs, prime=PRIME):
    """Test Rv in {+/-1} exactly, via odd integrality and its unit norm.

    Requires the caller's verified C^T=-C, C^2=-3p^2 I and sign vector v.
    Congruence modulo 4p makes (pI+C)v/(2p) odd integral.  Orthogonality
    then forces every entry to have magnitude one.
    """
    signs = np.array(signs, dtype=np.int64)
    if signs.shape != (len(c),) or not np.all(np.abs(signs) == 1):
        raise ValueError("expected a sign vector of the operator dimension")
    numerator = prime * signs + c @ signs
    return bool(np.all(numerator % (4 * prime) == 2 * prime))


def quaternion_orientation_audit():
    """Enumerate all sixteen skew sign-conference matrices of order four."""
    b = np.array(B, dtype=np.int64)
    pairs = [(i, j) for i in range(4) for j in range(i + 1, 4)]
    result = []
    for signs in product((-1, 1), repeat=6):
        d = np.zeros((4, 4), dtype=np.int64)
        for (i, j), sign in zip(pairs, signs, strict=True):
            d[i, j], d[j, i] = sign, -sign
        if np.array_equal(d @ d, -3 * np.eye(4, dtype=np.int64)):
            result.append(
                {
                    "D_upper_triangle": list(signs),
                    "D": d.tolist(),
                    "rank_stacked_B_BD_modulo_167": rank_mod_prime(np.vstack((b, b @ d)), PRIME),
                }
            )
    return result


def known_hadamard_frame_audit():
    """Check the orientations passing the first test on sixteen Hall masks.

    Two prescribed row orders are tested; independent column equivalence
    operations cannot change whether every column is flat. Arbitrary row
    operations and other Hadamard matrices are not excluded by this audit.
    """
    from experiments.scalar_lift.inverse_blocks import hall_variant, known_h668

    p, order = PRIME, 4 * PRIME
    b, q = np.array(B, dtype=np.int64), jacobsthal_matrix()
    live = [
        row for row in quaternion_orientation_audit() if row["rank_stacked_B_BD_modulo_167"] == 2
    ]
    result = [
        {
            "D_upper_triangle": row["D_upper_triangle"],
            "native_flat_columns_by_Hall_mask": [],
            "border_grouped_flat_columns_by_Hall_mask": [],
        }
        for row in live
    ]
    grouped = np.array(
        [j for block in range(4) for j in [block] + list(range(4 + 166 * block, 4 + 166 * (block + 1)))]
    )
    base = known_h668()
    for mask in range(16):
        matrix = np.array(hall_variant(base, mask), dtype=np.int64)
        for mode, h in (("native", matrix), ("border_grouped", matrix[grouped])):
            blocks = h.reshape(4, p, order)
            q_products = (q @ blocks.transpose(1, 0, 2).reshape(p, 4 * order)).reshape(
                p, 4, order
            ).transpose(1, 0, 2)
            bq_products = np.einsum("ab,bij->aij", b, q_products)
            sums = blocks.sum(axis=1)
            for row, outcome in zip(live, result, strict=True):
                d = np.array(row["D"], dtype=np.int64)
                numerator = p * blocks + bq_products + (d @ sums)[:, None, :]
                flat_columns = np.all(np.abs(numerator) == 2 * p, axis=(0, 1))
                outcome[f"{mode}_flat_columns_by_Hall_mask"].append(int(flat_columns.sum()))
    return result


def necessary_regular_block_sum_audit():
    """Enumerate the extremal block sums forced by a full compatible frame.

    The proof in docs/TOWER_LATTICE_INDEX.md shows every frame column has
    ||sigma||^2=4p and B*sigma=p*k, with k odd and ||k||^2=12.
    This is a necessary condition on an entire Hadamard frame, not on an
    arbitrary individual sign vector.
    """
    b = np.array(B, dtype=np.int64)
    d = np.array(D_COMPATIBLE, dtype=np.int64)
    integral_sums, odd_image_sums = [], []
    for large in range(4):
        for signs in product((-1, 1), repeat=4):
            k = np.array(signs, dtype=np.int64)
            k[large] *= 3
            numerator = b @ k
            if np.any(numerator % 3):
                continue
            sigma = numerator // 3
            if int(sigma @ sigma) != 4 * PRIME or np.any(sigma % 2 != 1):
                raise ValueError("an enumerated extremal sum failed its norm or parity")
            integral_sums.append(sigma.tolist())
            image_numerator = sigma + d @ sigma
            if np.all(image_numerator % 4 == 2):
                odd_image_sums.append(sigma.tolist())
    m = np.array(
        [[7, 3, 23, 9], [-3, 7, -9, 23], [23, -9, -7, 3], [-9, -23, 3, 7]],
        dtype=np.int64,
    )
    if not np.array_equal(m.T @ m, 4 * PRIME * np.eye(4, dtype=np.int64)):
        raise ValueError("the required block-sum frame is not exactly orthogonal")
    expected = {tuple(sign * m[:, j]) for j in range(4) for sign in (-1, 1)}
    if len(integral_sums) != 16 or set(map(tuple, odd_image_sums)) != expected:
        raise ValueError("the extremal enumeration differs from the stated four directions")
    return {
        "scope": "Necessary block sums if a full sign Hadamard frame exists for D_COMPATIBLE",
        "odd_k_norm_squared": 12,
        "integral_sigma_candidates": integral_sums,
        "sigma_candidates_with_odd_R_image_block_sums": odd_image_sums,
        "required_block_sum_matrix_M": m.tolist(),
        "M_transpose_M_exactly_668_I": True,
        "required_column_multiplicity_per_antipodal_direction": PRIME,
        "all_16_blocks_must_have_constant_row_and_column_sums_M": True,
        "circulant_blocks_assumed": False,
        "frame_constructed": False,
    }


def build_paley_sextic_operator_report():
    p, order = PRIME, 4 * PRIME
    b, d = np.array(B, dtype=np.int64), np.array(D, dtype=np.int64)
    c = dense_sextic_operator_numerator()
    identity = np.eye(order, dtype=np.int64)
    if not np.array_equal(b @ b, 3 * p * np.eye(4, dtype=np.int64)):
        raise ValueError("the Clifford square failed")
    if not np.array_equal(d @ d, -3 * np.eye(4, dtype=np.int64)):
        raise ValueError("the quaternion square failed")
    if not np.array_equal(c.T, -c) or not np.array_equal(c @ c, -3 * p * p * identity):
        raise ValueError("the full integer operator identity failed")
    numerator = p * identity + c
    if not np.all(numerator % 2 == 1):
        raise ValueError("the dyadic rank-one residue failed")
    off_diagonal = c[~np.eye(order, dtype=bool)]
    if not np.all(off_diagonal % 2 == 1) or np.max(np.abs(off_diagonal)) >= p:
        raise ValueError("the claimed denominator pattern failed")
    if not np.all((numerator @ np.ones(order, dtype=np.int64)) % (4 * p) == 0):
        raise ValueError("the binary-kernel simplification is unavailable")
    stacked = np.vstack((b, b @ d))
    if rank_mod_prime(stacked, p) != 4:
        raise ValueError("block sums are not forced to zero modulo p")
    # The block-sum proof reduces ALL sign vectors to these sixteen choices.
    constant_flat = []
    for block_signs in product((-1, 1), repeat=4):
        signs = np.repeat(block_signs, p)
        if flat_sign_images(c, signs):
            image = (p * signs + c @ signs) // (2 * p)
            if not np.all(np.abs(image) == 1):
                raise ValueError("the congruence did not yield actual signs")
            constant_flat.append(list(block_signs))
    flat_span = rank_mod_prime(constant_flat, p)
    if len(constant_flat) != 8 or flat_span != 4:
        raise ValueError("the exhaustive constant-block calculation changed")
    orientations = quaternion_orientation_audit()
    compatible_c = dense_sextic_operator_numerator(D_COMPATIBLE)
    if not np.array_equal(compatible_c @ compatible_c, -3 * p * p * identity):
        raise ValueError("the opposite-chirality operator square failed")
    # The moment builder uses definitions here, but never calls this report.
    from experiments.scalar_lift.paley_frame_moments import build_paley_frame_moment_report

    moment_obstruction = build_paley_frame_moment_report()
    return {
        "scope": "Explicit dense rational sextic operator; complete exclusion of its sign frame",
        "order": order,
        "prime": p,
        "B": [list(row) for row in B],
        "D": [list(row) for row in D],
        "construction": "C=B tensor Q_167 + D tensor J_167; R=(I+C/167)/2",
        "B_square": 501,
        "D_square": -3,
        "C_skew_and_square_verified_over_integers": True,
        "C_square": -3 * p * p,
        "maximum_absolute_off_diagonal_C_entry": int(np.max(np.abs(off_diagonal))),
        "R_reduced_denominators_per_row": {"2": 1, "334": order - 1},
        "twice_R_modulo_two_is_all_ones": True,
        "rank_C_modulo_167": rank_mod_prime(c, p),
        "rank_stacked_B_BD_modulo_167": rank_mod_prime(stacked, p),
        "BD_plus_DB": (b @ d + d @ b).tolist(),
        "flat_vector_criterion": "v=1-2x: Cx=0 mod167 and weight(x) odd",
        "block_sum_proof": "Bt=BDt=0 mod167 forces each binary block count t_i in {0,167}",
        "all_flat_vectors_constant_on_each_of_four_blocks": True,
        "flat_block_sign_vectors": constant_flat,
        "flat_vector_count": len(constant_flat),
        "dimension_of_flat_vector_span": flat_span,
        "full_sign_frame_excluded_for_this_operator": True,
        "quaternion_orientation_audit": orientations,
        "compatible_opposite_chirality_candidate": {
            "D": [list(row) for row in D_COMPATIBLE],
            "C_skew_and_square_verified_over_integers": bool(
                np.array_equal(compatible_c.T, -compatible_c)
            ),
            "rank_C_modulo_167": rank_mod_prime(compatible_c, p),
            "C_row_sums_by_block": compatible_c.sum(axis=1).reshape(4, p)[:, 0].tolist(),
            "rank_stacked_B_BD_modulo_167": rank_mod_prime(
                np.vstack((b, b @ np.array(D_COMPATIBLE, dtype=np.int64))), p
            ),
            "full_sign_frame_constructed": False,
            "full_sign_frame_excluded": True,
            "exclusion": "Exact local moment certificate for a required block-sum direction",
            "integer_left_null_contradiction": moment_obstruction["left_null_product_with_rhs"],
            "certificate": "certificates/scalar-lift/paley-frame-moments-668.json",
        },
        "known_hadamard_frame_audit": known_hadamard_frame_audit(),
        "compatible_candidate_regular_block_sum_audit": necessary_regular_block_sum_audit(),
        "new_real_hadamard_order_claimed": False,
    }

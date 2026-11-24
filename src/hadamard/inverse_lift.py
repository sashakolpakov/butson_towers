"""Exact inverse conditions for dense sixth-root scalar lifts.

Two real Hadamard blocks P,Q of order d realize the nonreal sixth roots
if and only if PQ^T+QP^T=dI. The source Butson matrix and the block pair
are independently verified before any output matrix is assembled.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Sequence
from operator import index

from .cgw import is_cgw_exponent_matrix
from .dense_lift import SignMatrix
from .scalar_lift import sixth_root_exponents


def _sign_rows(matrix: Sequence[Sequence[int]]) -> SignMatrix:
    rows = tuple(tuple(row) for row in matrix)
    if not rows or any(len(row) != len(rows) for row in rows):
        raise ValueError("expected a nonempty square sign matrix")
    if any(isinstance(x, bool) for row in rows for x in row):
        raise ValueError("every matrix entry must be an integer sign")
    try:
        rows = tuple(tuple(index(x) for x in row) for row in rows)
    except TypeError as exc:
        raise ValueError("every matrix entry must be an integer sign") from exc
    if any(x not in (-1, 1) for row in rows for x in row):
        raise ValueError("every matrix entry must be an integer sign")
    return rows


def _row_bits(rows: SignMatrix) -> tuple[int, ...]:
    return tuple(sum((x < 0) << j for j, x in enumerate(row)) for row in rows)


def _hadamard_bits(bits: tuple[int, ...]) -> bool:
    n = len(bits)
    return all(2 * (bits[i] ^ bits[j]).bit_count() == n for i in range(n) for j in range(i))


def is_sign_hadamard(matrix: Sequence[Sequence[int]]) -> bool:
    """Check the full integer Gram identity by exact binary row distances."""
    return _hadamard_bits(_row_bits(_sign_rows(matrix)))


def sixth_root_block_pair_conditions(
    first: Sequence[Sequence[int]], second: Sequence[Sequence[int]]
) -> dict:
    """Verify the two Hadamard identities and PQ^T+QP^T=dI exactly."""
    p, q = _sign_rows(first), _sign_rows(second)
    if len(p) != len(q):
        raise ValueError("the two block orders must agree")
    d = len(p)
    pb, qb = _row_bits(p), _row_bits(q)
    distances = tuple((a ^ b).bit_count() for a, b in zip(pb, qb, strict=True))
    diagonal = all(4 * distance == d for distance in distances)
    mixed = diagonal and all(
        (pb[i] ^ qb[j]).bit_count() + (qb[i] ^ pb[j]).bit_count() == d
        for i in range(d)
        for j in range(i)
    )
    first_valid, second_valid = _hadamard_bits(pb), _hadamard_bits(qb)
    return {
        "block_order": d,
        "first_is_hadamard": first_valid,
        "second_is_hadamard": second_valid,
        "diagonal_angle_condition": diagonal,
        "cross_symmetric_identity": mixed,
        "required_diagonal_row_distance": f"{d}/4",
        "observed_diagonal_row_distance_counts": dict(sorted(Counter(distances).items())),
        "valid": first_valid and second_valid and mixed,
    }


def lift_unreal_sixth_from_pair(
    matrix: Sequence[Sequence[int]],
    first: Sequence[Sequence[int]],
    second: Sequence[Sequence[int]],
    phase: int = 6,
) -> SignMatrix:
    """Build a dense real H(dn) from a verified unreal source and block pair.

    The substitution is zeta_6 -> P, zeta_6^2 -> Q, with their negatives.
    Its star representation is rho(zeta_6)=QP^T/d and frame Z=P-Q.
    A pair that fails the exact inverse identity is never used.
    """
    source = sixth_root_exponents(matrix, phase)
    if any(e is None or e % 3 == 0 for row in source for e in row):
        raise ValueError("the source must be dense and have no real sixth-root entries")
    if not is_cgw_exponent_matrix(source, 6, len(source)):
        raise ValueError("the source is not an exact Hadamard matrix")
    p, q = _sign_rows(first), _sign_rows(second)
    if not sixth_root_block_pair_conditions(p, q)["valid"]:
        raise ValueError("the block pair fails its exact Hadamard or angle identity")
    negative_p = tuple(tuple(-x for x in row) for row in p)
    negative_q = tuple(tuple(-x for x in row) for row in q)
    blocks = {1: p, 2: q, 4: negative_p, 5: negative_q}
    d = len(p)
    return tuple(
        tuple(blocks[e][r][c] for e in row for c in range(d)) for row in source for r in range(d)
    )

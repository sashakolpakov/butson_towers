"""Exact blocks and operator arithmetic for the inverse dense-lift audit."""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from math import gcd
from pathlib import Path

from sympy import isprime

from hadamard.cgw import is_cgw_exponent_matrix
from hadamard.inverse_lift import is_sign_hadamard, sixth_root_block_pair_conditions
from hadamard.scalar_lift import unreal_components

SEED = Path(__file__).with_name("h668_seed.json")
BH13_SEED = Path(__file__).with_name("bh13_6_seed.json")


def known_order_thirteen_source():
    """Verify the catalogue point and its explicit gauge into the unreal domain."""
    data = json.loads(BH13_SEED.read_text())
    source = data["exponents"]
    if not is_cgw_exponent_matrix(source, 6, 13):
        raise ValueError("the catalogue order-thirteen source fails its exact Gram identity")
    rows = data["unreal_row_phase_shifts"]
    columns = data["unreal_column_phase_shifts"]
    gauged = tuple(
        tuple((source[i][j] + rows[i] + columns[j]) % 6 for j in range(13)) for i in range(13)
    )
    if any(e % 3 == 0 for row in gauged for e in row):
        raise ValueError("the order-thirteen gauge is not unreal")
    if not is_cgw_exponent_matrix(gauged, 6, 13):
        raise ValueError("the order-thirteen gauge changes its Gram identity")
    return gauged


def triple_unreal_source(source):
    """The known tripling closure, implemented with a dense signed component proof.

    For H=(A+sqrt(-3)B)/2 let A'=B tensor J_3 and
    B'=A tensor I_3+B tensor C_3, where C_3 is the oriented 3-cycle.
    Both components remain dense signs, and their weighted Gram and
    amicability identities give an unreal BH(3n,6).
    """
    n = len(source)
    a, b = unreal_components(source)
    if any(x not in (-1, 1) for matrix in (a, b) for row in matrix for x in row):
        raise ValueError("tripling requires a dense unreal source")
    if not is_cgw_exponent_matrix(source, 6, n):
        raise ValueError("tripling requires an exact Hadamard source")
    c = ((0, 1, -1), (-1, 0, 1), (1, -1, 0))
    exponents = {(1, 1): 1, (-1, 1): 2, (-1, -1): 4, (1, -1): 5}
    result = tuple(
        tuple(
            exponents[b[i][j], a[i][j] if r == s else b[i][j] * c[r][s]]
            for j in range(n)
            for s in range(3)
        )
        for i in range(n)
        for r in range(3)
    )
    if not is_cgw_exponent_matrix(result, 6, 3 * n):
        raise ValueError("the tripled source failed its exact Gram identity")
    return result


def prime_block_operator_requirements(block_order):
    """Necessary operator arithmetic at d=4p, p an odd prime.

    This is a proved necessary condition, not a construction of the pair.
    Every cross inner product is twice an odd integer. Hence R=QP^T/d
    is dense with reduced denominators 2 or 2p. A row can contain at most
    three entries with denominator 2, since each has absolute value 1/2
    and all other entries are nonzero. The diagonal already uses one.
    """
    p, remainder = divmod(block_order, 4)
    if remainder or p <= 2 or not isprime(p):
        raise ValueError("this arithmetic bound requires d=4p with p an odd prime")
    return {
        "block_order": block_order,
        "odd_prime": p,
        "operator_must_be_dense": True,
        "allowed_reduced_denominators": [2, 2 * p],
        "maximum_denominator_two_entries_per_row": 3,
        "minimum_denominator_twice_prime_entries_per_row": block_order - 3,
        "half_integral_operator_excluded": True,
        "necessary_only": True,
    }


def pair_operator_arithmetic(first, second):
    """Inspect a verified pair using exact integer cross inner products."""
    if not sixth_root_block_pair_conditions(first, second)["valid"]:
        raise ValueError("operator arithmetic requires a valid block pair")
    d = len(first)
    products = tuple(
        tuple(sum(a * b for a, b in zip(q, p, strict=True)) for p in first) for q in second
    )
    denominators = tuple(Counter(d // gcd(abs(x), d) for x in row) for row in products)
    return {
        "block_order": d,
        "operator": "Q P^T / d",
        "zero_entries": sum(x == 0 for row in products for x in row),
        "reduced_denominator_counts_per_row": [
            {str(k): int(v) for k, v in sorted(row.items())} for row in denominators
        ],
    }


def quaternion_order_twelve_array(x, y, z, w):
    """The classical OD(12;3,3,3,3), with a directly checkable Gram."""
    rows = ((x, y, z), (-y, x, w), (-z, -w, x), (-w, z, -y))
    layout = (
        ((0, 0, 1), (1, -1, 1), (2, -1, 1), (3, -1, 1)),
        ((1, -1, -1), (0, 0, 1), (3, 1, 1), (2, 1, -1)),
        ((2, -1, -1), (3, 1, -1), (0, 0, 1), (1, 1, 1)),
        ((3, -1, -1), (2, 1, 1), (1, 1, -1), (0, 0, 1)),
    )

    def entry(i, j):
        r, a = divmod(i, 3)
        c, b = divmod(j, 3)
        sequence, mode, sign = layout[r][c]
        index = b - a if mode == 0 else mode * (a + b + 1)
        return sign * rows[sequence][index % 3]

    return tuple(tuple(entry(i, j) for j in range(12)) for i in range(12))


def known_order_twelve_pair():
    """A dense sextic block pair whose operator genuinely has denominator 6."""
    return (
        quaternion_order_twelve_array(1, 1, 1, 1),
        quaternion_order_twelve_array(-1, 1, 1, 1),
    )


def binary_digest(matrix):
    return hashlib.sha256(bytes((x + 1) // 2 for row in matrix for x in row)).hexdigest()


def known_h668():
    data = json.loads(SEED.read_text())
    v = data["cyclic_order"]
    sets = tuple(set(s) for s in data["negative_supports"])
    if tuple(map(len, sets)) != (82, 83, 83, 83):
        raise ValueError("incorrect source support sizes")
    if any(sum(sum((a + t) % v in s for a in s) for s in sets) != 164 for t in range(1, v)):
        raise ValueError("the source sets are not the stated cyclic difference family")
    rows = tuple(tuple(-1 if j in s else 1 for j in range(v)) for s in sets)
    # (sequence, index of its cyclic entry, sign) in the four-by-four core.
    layout = (
        ((0, 0, 1), (1, -1, 1), (2, -1, 1), (3, -1, 1)),
        ((1, -1, -1), (0, 0, 1), (3, 1, 1), (2, 1, -1)),
        ((2, -1, -1), (3, 1, -1), (0, 0, 1), (1, 1, 1)),
        ((3, -1, -1), (2, 1, 1), (1, 1, -1), (0, 0, 1)),
    )

    def entry(i, j):
        if i < 4:
            return data["corner"][i][j] if j < 4 else data["top_border"][i][(j - 4) // v]
        if j < 4:
            return data["left_border"][(i - 4) // v][j]
        r, a = divmod(i - 4, v)
        c, b = divmod(j - 4, v)
        sequence, mode, sign = layout[r][c]
        index = b - a if mode == 0 else mode * (a + b + 1)
        return sign * rows[sequence][index % v]

    matrix = tuple(tuple(entry(i, j) for j in range(668)) for i in range(668))
    if binary_digest(matrix) != data["base_matrix_binary_sha256"]:
        raise ValueError("the reconstructed order-668 matrix differs from the published hash")
    if not is_sign_hadamard(matrix):
        raise ValueError("the reconstructed order-668 matrix fails exact orthogonality")
    return matrix


def hall_variant(base, mask):
    if type(mask) is not int or not 0 <= mask < 16:
        raise ValueError("a Hall switch mask must be an integer from 0 through 15")
    return tuple(
        tuple(
            -x
            if (
                (i < 4 <= j and mask >> ((j - 4) // 166) & 1)
                or (j < 4 <= i and mask >> ((i - 4) // 166) & 1)
            )
            else x
            for j, x in enumerate(row)
        )
        for i, row in enumerate(base)
    )

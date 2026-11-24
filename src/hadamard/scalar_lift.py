"""Exact restriction of Eisenstein scalars and partial real weighing lifts.

Restriction in the integral basis ``(1, omega)`` has entries in ``{0,+/-1}``
but preserves an A2 metric, not the Euclidean metric.  A four by four
representation gives a Euclidean real weighing matrix for every sixth-root
source.  A denser lift requires all nonzero sixth roots to be nonreal; this
is the classical unreal-sixth-root construction, not a complete morphism
from all sixth roots to dense real Hadamard matrices.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from .butson import prime_divisors
from .cgw import ExponentMatrix, is_cgw_exponent_matrix

IntegerMatrix = tuple[tuple[int, ...], ...]

# Multiplication by zeta_6 = 1 + omega in the integral basis (1, omega).
_POWER_BASIS_BLOCKS: tuple[IntegerMatrix, ...] = (
    ((1, 0), (0, 1)),
    ((1, -1), (1, 0)),
    ((0, -1), (1, -1)),
    ((-1, 0), (0, -1)),
    ((-1, 1), (-1, 0)),
    ((0, 1), (-1, 1)),
)

# C^T = -C and C C^T = 3 I; its support complements that of I_4.
SKEW_CONFERENCE_4: IntegerMatrix = (
    (0, 1, 1, 1),
    (-1, 0, 1, -1),
    (-1, -1, 0, 1),
    (-1, 1, -1, 0),
)

_UNREAL_COMPONENTS = {1: (1, 1), 2: (-1, 1), 4: (-1, -1), 5: (1, -1)}

# U = (I_4 + C)/2 is orthogonal, Phi_6(U)=0.  These are U^e Z,
# where Z = diag(H_2,H_2) and Z Z^T = 2 I.  Every block has weight two.
EISENSTEIN_WEIGHING_BLOCKS: tuple[IntegerMatrix, ...] = (
    ((1, 1, 0, 0), (1, -1, 0, 0), (0, 0, 1, 1), (0, 0, 1, -1)),
    ((1, 0, 1, 0), (0, -1, 0, 1), (-1, 0, 1, 0), (0, -1, 0, -1)),
    ((0, -1, 1, 0), (-1, 0, 0, 1), (-1, 0, 0, -1), (0, -1, -1, 0)),
    ((-1, -1, 0, 0), (-1, 1, 0, 0), (0, 0, -1, -1), (0, 0, -1, 1)),
    ((-1, 0, -1, 0), (0, 1, 0, -1), (1, 0, -1, 0), (0, 1, 0, 1)),
    ((0, 1, -1, 0), (1, 0, 0, -1), (1, 0, 0, 1), (0, 1, 1, 0)),
)


def _square_integer_matrix(matrix: Sequence[Sequence[int]]) -> IntegerMatrix:
    rows = tuple(tuple(row) for row in matrix)
    if not rows or any(len(row) != len(rows) for row in rows):
        raise ValueError("expected a nonempty square integer matrix")
    if any(type(entry) is not int for row in rows for entry in row):
        raise ValueError("matrix entries must be integers")
    return rows


def _parameters(order: int, weight: int) -> None:
    if type(order) is not int or type(weight) is not int or not 1 <= weight <= order:
        raise ValueError("expected integers with 1 <= weight <= order")


def sixth_root_exponents(matrix: ExponentMatrix, phase: int = 6) -> ExponentMatrix:
    """Embed third or sixth roots in mu_6, retaining ``None`` for zeros."""
    if type(phase) is not int or phase not in (3, 6):
        raise ValueError("the source alphabet must be third or sixth roots")
    rows = tuple(tuple(row) for row in matrix)
    if not rows or any(len(row) != len(rows) for row in rows):
        raise ValueError("expected a nonempty square exponent matrix")
    if any(entry is not None and type(entry) is not int for row in rows for entry in row):
        raise ValueError("an exponent must be an integer or None")
    return tuple(
        tuple(None if entry is None else (6 // phase * entry) % 6 for entry in row) for row in rows
    )


def restrict_eisenstein_scalars(matrix: ExponentMatrix, phase: int = 6) -> IntegerMatrix:
    """Return the integral 2n by 2n representation, without claiming Euclidean orthogonality.

    If the source satisfies ``H H* = w I``, the result L satisfies
    ``L (I_n tensor Q) L^T = w (I_n tensor Q)`` for Q = [[2,1],[1,2]].
    The dual trace form Q is essential; in general ``L L^T != w I``.
    """
    source = sixth_root_exponents(matrix, phase)
    zero: IntegerMatrix = ((0, 0), (0, 0))
    return tuple(
        tuple(
            (zero if entry is None else _POWER_BASIS_BLOCKS[entry])[inside_row][inside_col]
            for entry in row
            for inside_col in range(2)
        )
        for row in source
        for inside_row in range(2)
    )


def lift_eisenstein_weighing(matrix: ExponentMatrix, weight: int, phase: int = 6) -> IntegerMatrix:
    """Lift any CGW(n,w;6) to a real W(4n,2w), with no unreal restriction.

    The output has zeros, even when the source is dense.  No Hadamard
    completion is implied.  Source orthogonality is checked exactly.
    """
    source = sixth_root_exponents(matrix, phase)
    _parameters(len(source), weight)
    if not is_cgw_exponent_matrix(source, 6, weight):
        raise ValueError("the source is not an exact sixth-root weighing matrix")
    zero: IntegerMatrix = ((0,) * 4,) * 4
    return tuple(
        tuple(
            (zero if entry is None else EISENSTEIN_WEIGHING_BLOCKS[entry])[inside_row][inside_col]
            for entry in row
            for inside_col in range(4)
        )
        for row in source
        for inside_row in range(4)
    )


@dataclass(frozen=True)
class UnrealGaugeResult:
    """A bounded F_3 constraint search; only ``exhausted`` proves no gauge exists."""

    status: str
    nodes_visited: int
    row_shifts: tuple[int, ...] | None = None
    column_shifts: tuple[int, ...] | None = None
    matrix: ExponentMatrix | None = None


def find_unreal_gauge(
    matrix: ExponentMatrix, phase: int = 6, *, max_nodes: int = 100_000
) -> UnrealGaugeResult:
    """Find mu_6 row/column phases avoiding real nonzero entries.

    Realness only depends on the exponent modulo three.  For each supported
    entry solve ``e_ij + r_i + c_j != 0 (mod 3)``.  Fixing one row shift to
    zero removes a global gauge, while bit masks propagate each column's
    remaining shifts.  No floating point or modular-to-integer inference is
    used.  A node limit returns ``node_limit``, never ``exhausted``.
    """
    source = sixth_root_exponents(matrix, phase)
    if type(max_nodes) is not int or max_nodes < 1:
        raise ValueError("max_nodes must be a positive integer")
    order = len(source)
    ordering = sorted(range(order), key=lambda i: -sum(e is not None for e in source[i]))
    stack = [(0, (7,) * order, ())]
    visited = 0
    while stack:
        if visited == max_nodes:
            return UnrealGaugeResult("node_limit", visited)
        depth, masks, shifts = stack.pop()
        visited += 1
        if depth == order:
            row_shifts = [0] * order
            for row, shift in zip(ordering, shifts, strict=True):
                row_shifts[row] = shift
            column_shifts = tuple((mask & -mask).bit_length() - 1 for mask in masks)
            adjusted = tuple(
                tuple(
                    None if entry is None else (entry + row_shifts[i] + column_shifts[j]) % 6
                    for j, entry in enumerate(row)
                )
                for i, row in enumerate(source)
            )
            return UnrealGaugeResult("found", visited, tuple(row_shifts), column_shifts, adjusted)
        row = source[ordering[depth]]
        for shift in (0,) if depth == 0 else (2, 1, 0):
            updated = tuple(
                mask if entry is None else mask & ~(1 << ((-entry - shift) % 3))
                for mask, entry in zip(masks, row, strict=True)
            )
            if all(updated):
                stack.append((depth + 1, updated, shifts + (shift,)))
    return UnrealGaugeResult("exhausted", visited)


def unreal_components(
    matrix: ExponentMatrix, phase: int = 6
) -> tuple[IntegerMatrix, IntegerMatrix]:
    """Write an unreal source as ``(A + sqrt(-3) B)/2`` with shared support."""
    source = sixth_root_exponents(matrix, phase)
    if any(entry in (0, 3) for row in source for entry in row):
        raise ValueError("the source has real nonzero entries; find an unreal gauge first")
    pairs = tuple(
        tuple((0, 0) if entry is None else _UNREAL_COMPONENTS[entry] for entry in row)
        for row in source
    )
    return tuple(tuple(pair[0] for pair in row) for row in pairs), tuple(
        tuple(pair[1] for pair in row) for row in pairs
    )


def unreal_pair_conditions(
    first: Sequence[Sequence[int]], second: Sequence[Sequence[int]], weight: int
) -> dict[str, bool]:
    """Check support, ``AA^T+3BB^T=4wI``, and ``AB^T=BA^T`` over the integers."""
    a = _square_integer_matrix(first)
    b = _square_integer_matrix(second)
    if len(a) != len(b):
        raise ValueError("the two component matrices must have the same order")
    _parameters(len(a), weight)
    support = all(
        (x == y == 0) or (abs(x) == abs(y) == 1)
        for ar, br in zip(a, b, strict=True)
        for x, y in zip(ar, br, strict=True)
    )
    gram = all(
        sum(x * y + 3 * u * v for x, y, u, v in zip(ar, ac, br, bc, strict=True))
        == (4 * weight if i == j else 0)
        for i, (ar, br) in enumerate(zip(a, b, strict=True))
        for j, (ac, bc) in enumerate(zip(a, b, strict=True))
    )
    amicable = all(
        sum(x * v - u * y for x, y, u, v in zip(ar, ac, br, bc, strict=True)) == 0
        for ar, br in zip(a, b, strict=True)
        for ac, bc in zip(a, b, strict=True)
    )
    return {
        "shared_signed_support": support,
        "weighted_gram_identity": gram,
        "amicability_identity": amicable,
        "valid": support and gram and amicable,
    }


def lift_unreal_pair(
    first: Sequence[Sequence[int]], second: Sequence[Sequence[int]], weight: int
) -> IntegerMatrix:
    """Construct a real W(4n,4w), rejecting every unverified component pair."""
    a = _square_integer_matrix(first)
    b = _square_integer_matrix(second)
    if not unreal_pair_conditions(a, b, weight)["valid"]:
        raise ValueError("the component pair does not satisfy the exact unreal Gram identities")
    return tuple(
        tuple(
            a[i][j]
            if inside_row == inside_col
            else b[i][j] * SKEW_CONFERENCE_4[inside_row][inside_col]
            for j in range(len(a))
            for inside_col in range(4)
        )
        for i in range(len(a))
        for inside_row in range(4)
    )


def lift_unreal_cgw(matrix: ExponentMatrix, weight: int, phase: int = 6) -> IntegerMatrix:
    """Lift a verified unreal CGW(n,w;6) to a real W(4n,4w)."""
    return lift_unreal_pair(*unreal_components(matrix, phase), weight)


def is_real_weighing_matrix(matrix: Sequence[Sequence[int]], weight: int) -> bool:
    """Independently check the signed alphabet and both Euclidean Gram matrices."""
    rows = _square_integer_matrix(matrix)
    _parameters(len(rows), weight)
    if any(entry not in (-1, 0, 1) for row in rows for entry in row):
        return False
    columns = tuple(zip(*rows, strict=True))
    return all(
        sum(x * y for x, y in zip(first, second, strict=True)) == (weight if i == j else 0)
        for vectors in (rows, columns)
        for i, first in enumerate(vectors)
        for j, second in enumerate(vectors)
    )


def complete_half_weighing(
    first: Sequence[Sequence[int]], second: Sequence[Sequence[int]]
) -> IntegerMatrix:
    """Verify a complementary weighing pair and return its dense real Hadamard sum.

    Each matrix must have weight N/2, their supports must partition every
    entry, and their mixed Gram terms must cancel.  Complementary supports
    alone are insufficient.
    """
    a = _square_integer_matrix(first)
    b = _square_integer_matrix(second)
    order = len(a)
    if len(b) != order or order % 2:
        raise ValueError("expected equal even matrix orders")
    if not (is_real_weighing_matrix(a, order // 2) and is_real_weighing_matrix(b, order // 2)):
        raise ValueError("both matrices must be weighing matrices of weight N/2")
    if any(
        (x == 0) == (y == 0)
        for ar, br in zip(a, b, strict=True)
        for x, y in zip(ar, br, strict=True)
    ):
        raise ValueError("the supports do not partition the full matrix")
    if any(
        sum(x * v + u * y for x, y, u, v in zip(ar, ac, br, bc, strict=True)) != 0
        for ar, br in zip(a, b, strict=True)
        for ac, bc in zip(a, b, strict=True)
    ):
        raise ValueError("the mixed Gram terms do not cancel")
    return tuple(
        tuple(x + y for x, y in zip(ar, br, strict=True)) for ar, br in zip(a, b, strict=True)
    )


def eisenstein_inert_prime_obstructions(order: int, weight: int) -> tuple[int, ...]:
    """Return primes ruling out *any* CGW(order,weight;6), before gauge search.

    The determinant is an Eisenstein integer of norm ``weight**order``.
    Primes congruent to two modulo three are inert and must have even norm
    valuation.  An empty tuple is only a passed necessary condition.
    """
    _parameters(order, weight)
    if order % 2 == 0:
        return ()
    obstructions = []
    for prime in prime_divisors(weight):
        if prime % 3 != 2:
            continue
        valuation = 0
        remaining = weight
        while remaining % prime == 0:
            valuation += 1
            remaining //= prime
        if valuation % 2:
            obstructions.append(prime)
    return tuple(obstructions)

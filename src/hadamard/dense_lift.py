"""Exact cyclotomic transfers whose outputs have only +1 and -1 entries.

The tenth-root transfer is the partial morphism of Egan and O Cathain,
Morphisms of Butson classes, Example 12, followed by Turyn's transfer.
Its source must have no real entries. No zero-filling step is involved.
"""

from __future__ import annotations

from collections.abc import Sequence

from .cgw import ExponentMatrix, is_cgw_exponent_matrix

SignMatrix = tuple[tuple[int, ...], ...]

# Exponents of i in 2M, 2M^2, 2M^3, 2M^4, where
# M=H_4 diag(-1,1,i,i)/2.  M is unitary and Phi_5(M)=0.
FIFTH_POWER_QUATERNARY_BLOCKS: tuple[ExponentMatrix, ...] = (
    ((2, 0, 1, 1), (2, 2, 1, 3), (2, 0, 3, 3), (2, 2, 3, 1)),
    ((3, 2, 0, 3), (0, 1, 3, 0), (1, 2, 2, 3), (0, 3, 3, 2)),
    ((1, 0, 3, 0), (2, 3, 2, 1), (0, 1, 2, 1), (1, 0, 1, 2)),
    ((2, 2, 2, 2), (0, 2, 0, 2), (3, 3, 1, 1), (3, 1, 1, 3)),
)

# rho(i)=[[0,-1],[1,0]], frame H_2: each rho(i)^e H_2 is a sign block.
QUATERNARY_REAL_BLOCKS: tuple[SignMatrix, ...] = (
    ((1, 1), (1, -1)),
    ((-1, 1), (1, 1)),
    ((-1, -1), (-1, 1)),
    ((1, -1), (-1, -1)),
)

# (1+i) T^e for T=i(I-S)/2, S=[[i,1+i],[-1+i,-i]].
# T is unitary with Phi_12(T)=0. The eight displayed powers are flat;
# the four omitted powers are not in the domain of this partial transfer.
TWELFTH_POWER_QUATERNARY_BLOCKS: dict[int, ExponentMatrix] = {
    1: ((1, 0), (1, 2)),
    2: ((1, 1), (2, 0)),
    4: ((2, 1), (2, 3)),
    5: ((2, 2), (3, 1)),
    7: ((3, 2), (3, 0)),
    8: ((3, 3), (0, 2)),
    10: ((0, 3), (0, 1)),
    11: ((0, 0), (1, 3)),
}


def _dense_exponents(matrix: Sequence[Sequence[int]], phase: int) -> ExponentMatrix:
    rows = tuple(tuple(row) for row in matrix)
    if not rows or any(len(row) != len(rows) for row in rows):
        raise ValueError("expected a nonempty square exponent matrix")
    if any(type(x) is not int for row in rows for x in row):
        raise ValueError("dense Hadamard inputs require integer exponents and no zeros")
    return tuple(tuple(x % phase for x in row) for row in rows)


def lift_quaternary_hadamard(matrix: Sequence[Sequence[int]]) -> SignMatrix:
    """Return a real Hadamard H(2n) from an exactly verified BH(n,4)."""
    source = _dense_exponents(matrix, 4)
    if not is_cgw_exponent_matrix(source, 4, len(source)):
        raise ValueError("the source is not an exact quaternary Hadamard matrix")
    return tuple(
        tuple(QUATERNARY_REAL_BLOCKS[e][r][c] for e in row for c in range(2))
        for row in source
        for r in range(2)
    )


def unreal_tenth_to_quaternary(matrix: Sequence[Sequence[int]], phase: int = 10) -> ExponentMatrix:
    """Return BH(4n,4) from BH(n,10) without +/-1, or BH(n,5) without 1.

    Entry zeta_5^k maps to 2M^k. For tenth roots use
    zeta_10^e=(-1)^e zeta_5^(3e). All permitted scalar blocks are dense.
    """
    if type(phase) is not int or phase not in (5, 10):
        raise ValueError("the source alphabet must be fifth or tenth roots")
    source = _dense_exponents(matrix, phase)
    if any(e % 5 == 0 for row in source for e in row):
        raise ValueError("the source contains a forbidden real entry; the transfer is partial")
    if not is_cgw_exponent_matrix(source, phase, len(source)):
        raise ValueError("the source is not an exact fifth- or tenth-root Hadamard matrix")

    def block(e: int) -> ExponentMatrix:
        power = e if phase == 5 else 3 * e % 5
        sign = 0 if phase == 5 else 2 * (e % 2)
        return tuple(
            tuple((x + sign) % 4 for x in row) for row in FIFTH_POWER_QUATERNARY_BLOCKS[power - 1]
        )

    blocks = {e: block(e) for row in source for e in row}
    return tuple(
        tuple(blocks[e][r][c] for e in row for c in range(4)) for row in source for r in range(4)
    )


def lift_unreal_tenth_hadamard(matrix: Sequence[Sequence[int]], phase: int = 10) -> SignMatrix:
    """Return a dense real H(8n) from a verified unreal fifth/tenth-root input.

    Every output entry is +/-1. Orthogonality follows from two exact
    star-representation identities, with scale factors four and two.
    """
    return lift_quaternary_hadamard(unreal_tenth_to_quaternary(matrix, phase))


def unreal_twelfth_to_quaternary(
    matrix: Sequence[Sequence[int]], phase: int = 12
) -> ExponentMatrix:
    """Return BH(2n,4) from a twelfth-root source avoiding the fourth roots.

    Sixth- or third-root inputs are embedded in the twelfth roots first.
    This keeps Q(i) as the intermediate coefficient field instead of
    turning the first-stage matrix into rational coefficients.
    """
    if type(phase) is not int or phase not in (3, 6, 12):
        raise ValueError("the source alphabet must be third, sixth or twelfth roots")
    source = _dense_exponents(matrix, phase)
    adjusted = tuple(tuple((12 // phase) * e for e in row) for row in source)
    if any(e % 3 == 0 for row in adjusted for e in row):
        raise ValueError("the source contains a forbidden fourth root; the transfer is partial")
    if not is_cgw_exponent_matrix(source, phase, len(source)):
        raise ValueError("the source is not an exact Hadamard matrix")
    return tuple(
        tuple(TWELFTH_POWER_QUATERNARY_BLOCKS[e][r][c] for e in row for c in range(2))
        for row in adjusted
        for r in range(2)
    )


def lift_unreal_twelfth_hadamard(matrix: Sequence[Sequence[int]], phase: int = 12) -> SignMatrix:
    """Compose the exact twelfth-to-fourth-root and fourth-to-sign transfers."""
    return lift_quaternary_hadamard(unreal_twelfth_to_quaternary(matrix, phase))

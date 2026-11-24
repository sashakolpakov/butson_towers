"""Exhaust the integral symmetric order-four coefficients for the 668 template.

This is a bounded exact classification, not a search for a Hadamard matrix.
Every B with B^T=B, B^2=501I and B=I modulo two is conjugate by a signed
permutation to one of the returned matrices.  The proof of the normalization
and bounds accompanies the family obstruction in the documentation.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations_with_replacement, product
from math import isqrt

from sympy import Matrix


def normalized_coefficient_forms():
    """Enumerate every normalized B using four linear equations and two entries.

    The polynomial x^2-501 is irreducible, so trace(B)=0. Choose a positive
    diagonal entry first, make its off-diagonal entries nonnegative by sign
    conjugation, and sort them by permuting the last three coordinates.
    All entries have magnitude at most floor(sqrt(501))=22.
    """
    norm = 501
    bound = isqrt(norm)
    first_rows = [
        (a, *row)
        for a in range(1, bound + 1, 2)
        for row in combinations_with_replacement(range(0, bound + 1, 2), 3)
        if a * a + sum(x * x for x in row) == norm
    ]
    matrices = []
    trials = 0
    for a, u, v, w in first_rows:
        # B=[[a,r^T],[r,S]]. Unknowns of S: three diagonal entries b,c,d,
        # then entries x,y,z at positions (0,1),(0,2),(1,2).
        equations = Matrix(
            [
                [u, 0, 0, v, w, 0, -a * u],
                [0, v, 0, u, 0, w, -a * v],
                [0, 0, w, 0, u, v, -a * w],
                [1, 1, 1, 0, 0, 0, -a],
            ]
        )
        reduced, pivots = equations.rref()
        if len(pivots) != 4 or any(j >= 6 for j in pivots):
            raise ValueError("the normalized first-row equations do not have rank four")
        free = [j for j in range(6) if j not in pivots]
        coefficients = [[Fraction(reduced[i, j]) for j in range(7)] for i in range(4)]
        ranges = [
            range(-bound + 1, bound, 2) if j < 3 else range(-bound, bound + 1, 2)
            for j in free
        ]
        for values in product(*ranges):
            trials += 1
            entries = [None] * 6
            for j, value in zip(free, values, strict=True):
                entries[j] = value
            for i, j in enumerate(pivots):
                value = coefficients[i][6] - sum(coefficients[i][k] * entries[k] for k in free)
                if value.denominator != 1 or abs(value) > bound or value.numerator % 2 != (j < 3):
                    break
                entries[j] = int(value)
            else:
                b, c, d, x, y, z = entries
                matrix = [[a, u, v, w], [u, b, x, y], [v, x, c, z], [w, y, z, d]]
                if all(
                    sum(matrix[i][k] * matrix[k][j] for k in range(4))
                    == (norm if i == j else 0)
                    for i in range(4)
                    for j in range(4)
                ):
                    matrices.append(matrix)
    return {
        "norm": norm,
        "entry_bound": bound,
        "normalization": "positive first diagonal, nonnegative sorted first-row off-diagonals",
        "equivalence": "simultaneous conjugation by a signed permutation",
        "normalized_first_rows": [list(row) for row in first_rows],
        "bounded_free_entry_assignments_checked": trials,
        "matrices": matrices,
        "scope": "All integral symmetric B of order four with B^2=501I and B=I mod2",
    }

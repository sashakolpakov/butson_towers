"""Exact integer-matrix operations with optional FLINT acceleration."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp

from .types import IntegerMatrix, gcd_many


def _as_matrix(rows: Sequence[Sequence[int]]) -> sp.Matrix:
    return sp.Matrix([[int(value) for value in row] for row in rows])


def _as_integer_tuple(matrix: sp.Matrix) -> IntegerMatrix:
    out = []
    for row in matrix.tolist():
        converted = []
        for value in row:
            rational = sp.Rational(value)
            if rational.q != 1:
                raise ValueError("matrix is not integral")
            converted.append(int(rational.p))
        out.append(tuple(converted))
    return tuple(out)


def matrix_rank(rows: Sequence[Sequence[int]], width: int | None = None) -> int:
    if not rows:
        return 0
    values = [tuple(int(x) for x in row) for row in rows]
    if width is not None:
        values = [row + (0,) * (width - len(row)) for row in values]
    try:
        from flint import fmpz_mat
    except ImportError:
        pass
    else:
        return int(fmpz_mat(values).rank())
    return int(_as_matrix(values).rank())


def determinant(matrix: Sequence[Sequence[int]]) -> int:
    return int(_as_matrix(matrix).det())


def inverse_unimodular(matrix: Sequence[Sequence[int]]) -> IntegerMatrix:
    value = _as_matrix(matrix)
    determinant_value = int(value.det())
    if abs(determinant_value) != 1:
        raise ValueError("matrix is not unimodular")
    return _as_integer_tuple(value.inv())


def transpose(matrix: Sequence[Sequence[int]]) -> IntegerMatrix:
    return _as_integer_tuple(_as_matrix(matrix).T)


def matrix_multiply(left: Sequence[Sequence[int]], right: Sequence[Sequence[int]]) -> IntegerMatrix:
    return _as_integer_tuple(_as_matrix(left) * _as_matrix(right))


def column_to_unimodular(vector: Sequence[int], *, return_inverse: bool = False):
    """Complete a primitive integer column to a unimodular matrix.

    This is the Euclidean row-reduction construction used by the KKPR backend,
    implemented with plain integer operations.
    """
    v = tuple(int(value) for value in vector)
    if not v:
        raise ValueError("vector must be nonempty")
    if gcd_many(v) != 1:
        raise ValueError("vector is not primitive")
    n = len(v)
    U = [[1 if i == j else 0 for j in range(n)] for i in range(n)]

    def multiply_vector() -> list[int]:
        return [sum(U[i][j] * v[j] for j in range(n)) for i in range(n)]

    while True:
        w = multiply_vector()
        candidates = [i for i, value in enumerate(w) if value]
        if not candidates:
            raise ArithmeticError("primitive vector unexpectedly reduced to zero")
        pivot = min(candidates, key=lambda i: abs(w[i]))
        if pivot:
            U[0], U[pivot] = U[pivot], U[0]
            w = multiply_vector()
        if w[0] < 0:
            U[0] = [-value for value in U[0]]
            w = multiply_vector()
        for row in range(1, n):
            quotient = w[row] // w[0]
            if quotient:
                U[row] = [U[row][column] - quotient * U[0][column] for column in range(n)]
        w = multiply_vector()
        if not any(w[row] for row in range(1, n)):
            break
    reducer = tuple(tuple(row) for row in U)
    completion = inverse_unimodular(reducer)
    if return_inverse:
        return completion, reducer
    return completion


def row_to_unimodular(vector: Sequence[int], *, return_inverse: bool = False):
    if return_inverse:
        completion, reducer = column_to_unimodular(vector, return_inverse=True)
        return transpose(completion), transpose(reducer)
    return transpose(column_to_unimodular(vector))


def hermite_form(rows: Sequence[Sequence[int]]) -> IntegerMatrix:
    if not rows:
        return tuple()
    matrix = _as_matrix(rows)
    # SymPy computes HNF for matrices with rows <= columns most naturally after transpose.
    hnf = hermite_normal_form(matrix.T).T
    nonzero = [tuple(int(value) for value in row) for row in hnf.tolist() if any(row)]
    return tuple(nonzero)


def hermite_form_with_transform(
    rows: Sequence[Sequence[int]],
) -> tuple[IntegerMatrix, IntegerMatrix]:
    """Return row HNF and the integral row operation producing it.

    For a full-row-rank matrix ``A``, the result is ``(H, U)`` with
    ``H = U A``. The transformation lets callers transport character values
    carried by the rows along with the exponent lattice.
    """

    if not rows:
        return tuple(), tuple()
    values = tuple(tuple(int(value) for value in row) for row in rows)
    if not values[0] or any(len(row) != len(values[0]) for row in values):
        raise ValueError("HNF requires a nonempty rectangular matrix")
    try:
        from flint import fmpz_mat
    except ImportError:
        matrix = _as_matrix(values)
        hnf_matrix = hermite_normal_form(matrix.T).T
        nonzero_indices = [
            index
            for index, row in enumerate(hnf_matrix.tolist())
            if any(row)
        ]
        if len(nonzero_indices) != len(values):
            raise ValueError("HNF transformation requires full row rank") from None
        hnf_matrix = hnf_matrix[nonzero_indices, :]
        transform = hnf_matrix * matrix.T * (matrix * matrix.T).inv()
        hnf = _as_integer_tuple(hnf_matrix)
        integral_transform = _as_integer_tuple(transform)
        if transform * matrix != hnf_matrix:
            raise ArithmeticError(
                "failed to reconstruct the HNF row transform"
            ) from None
        return hnf, integral_transform

    # FLINT's row-HNF convention pivots from the left, whereas the historical
    # project key is the transpose of SymPy's column HNF. Reversing ambient
    # columns before row HNF, then reversing rows and columns of the result,
    # gives exactly that established canonical form.
    matrix = fmpz_mat(
        tuple(tuple(reversed(row)) for row in values)
    )
    hnf_matrix, transform_matrix = matrix.hnf(transform=True)
    nonzero_rows = [
        row
        for row in range(hnf_matrix.nrows())
        if any(
            hnf_matrix[row, column]
            for column in range(hnf_matrix.ncols())
        )
    ]
    hnf = tuple(
        tuple(
            int(hnf_matrix[row, column])
            for column in reversed(range(hnf_matrix.ncols()))
        )
        for row in reversed(nonzero_rows)
    )
    if len(hnf) != len(values):
        raise ValueError("HNF transformation requires full row rank")
    transform = tuple(
        tuple(
            int(transform_matrix[row, column])
            for column in range(transform_matrix.ncols())
        )
        for row in reversed(nonzero_rows)
    )
    return hnf, transform


def smith_diagonal(rows: Sequence[Sequence[int]]) -> tuple[int, ...]:
    if not rows:
        return tuple()
    diagonal = smith_normal_form(_as_matrix(rows), domain=sp.ZZ)
    return tuple(
        abs(int(diagonal[i, i]))
        for i in range(min(diagonal.rows, diagonal.cols))
        if diagonal[i, i]
    )


@dataclass(frozen=True, slots=True)
class SmithDecomposition:
    """Smith form together with its two unimodular basis changes.

    The matrices satisfy ``normal_form = left * source * right``.  Keeping
    both transforms is essential for translated lattices: the left transform
    acts on character constants, while the inverse right transform gives a
    primitive basis for the saturated exponent lattice.
    """

    normal_form: IntegerMatrix
    left_transform: IntegerMatrix
    right_transform: IntegerMatrix

    @property
    def invariant_factors(self) -> tuple[int, ...]:
        if not self.normal_form:
            return tuple()
        return tuple(
            abs(self.normal_form[index][index])
            for index in range(
                min(len(self.normal_form), len(self.normal_form[0]))
            )
            if self.normal_form[index][index]
        )


def smith_decomposition(
    rows: Sequence[Sequence[int]],
) -> SmithDecomposition:
    """Return an exact Smith decomposition with transformation matrices."""

    if not rows:
        raise ValueError("Smith decomposition requires at least one row")
    width = len(rows[0])
    if width == 0 or any(len(row) != width for row in rows):
        raise ValueError("Smith decomposition requires a nonempty rectangular matrix")
    values = [[ZZ(int(value)) for value in row] for row in rows]
    source = DomainMatrix(values, (len(values), width), ZZ)
    normal, left, right = smith_normal_decomp(source)
    return SmithDecomposition(
        _as_integer_tuple(normal.to_Matrix()),
        _as_integer_tuple(left.to_Matrix()),
        _as_integer_tuple(right.to_Matrix()),
    )


@dataclass(frozen=True, slots=True)
class LatticeSummary:
    rank: int
    hermite_rows: IntegerMatrix
    smith_diagonal: tuple[int, ...]


def summarize_lattice(rows: Sequence[Sequence[int]]) -> LatticeSummary:
    return LatticeSummary(matrix_rank(rows), hermite_form(rows), smith_diagonal(rows))

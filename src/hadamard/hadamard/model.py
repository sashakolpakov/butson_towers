"""The common dephased Hadamard Laurent variety."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass

from ..coefficients.cyclotomic import CyclotomicField
from ..laurent.ideal import LaurentIdeal
from ..laurent.polynomial import LaurentPolynomial, LaurentRing


@dataclass(frozen=True, slots=True)
class EquationTag:
    kind: str
    first: int
    second: int
    power: int = 1


class HadamardTorus:
    """Laurent coordinates for a dephased order-``n`` complex Hadamard matrix.

    The first row and first column are constants equal to one.  The remaining
    ``(n-1)^2`` entries are independent torus coordinates before orthogonality
    equations are imposed.
    """

    def __init__(
        self,
        n: int,
        *,
        field: CyclotomicField | None = None,
        term_order: str = "degrevlex",
        include_column_equations: bool = True,
        include_reverse_equations: bool = True,
    ) -> None:
        if int(n) < 2:
            raise ValueError("Hadamard order must be at least 2")
        self.n = int(n)
        self.m = self.n - 1
        self.field = field or CyclotomicField.rational()
        self.term_order = str(term_order)
        self.include_column_equations = bool(include_column_equations)
        self.include_reverse_equations = bool(include_reverse_equations)
        names = tuple(f"x_{r}_{c}" for r in range(1, self.n) for c in range(1, self.n))
        self.ring = LaurentRing(self.field, names, self.term_order)
        self.variables = self.ring.variables()
        self._tagged_generator_cache: dict[
            int,
            tuple[tuple[EquationTag, LaurentPolynomial], ...],
        ] = {}
        self._index = {
            (r, c): (r - 1) * self.m + (c - 1)
            for r in range(1, self.n)
            for c in range(1, self.n)
        }

    @property
    def ambient_dimension(self) -> int:
        return self.m * self.m

    def variable_index(self, row: int, column: int) -> int:
        if row == 0 or column == 0:
            raise ValueError("the dephased border consists of constants")
        try:
            return self._index[(int(row), int(column))]
        except KeyError as exc:
            raise IndexError((row, column)) from exc

    def h(self, row: int, column: int) -> LaurentPolynomial:
        if not 0 <= int(row) < self.n or not 0 <= int(column) < self.n:
            raise IndexError((row, column))
        if row == 0 or column == 0:
            return self.ring.one()
        return self.variables[self.variable_index(row, column)]

    def ratio(self, row: int, other_row: int, column: int, *, power: int = 1) -> LaurentPolynomial:
        if row == other_row:
            raise ValueError("row indices must be distinct")
        return (self.h(row, column) / self.h(other_row, column)) ** int(power)

    def column_ratio(self, column: int, other_column: int, row: int, *, power: int = 1) -> LaurentPolynomial:
        if column == other_column:
            raise ValueError("column indices must be distinct")
        return (self.h(row, column) / self.h(row, other_column)) ** int(power)

    def row_equation(
        self,
        row: int,
        other_row: int,
        *,
        power: int = 1,
        columns: Iterable[int] | None = None,
    ) -> LaurentPolynomial:
        indices = range(self.n) if columns is None else tuple(int(c) for c in columns)
        return sum((self.ratio(row, other_row, c, power=power) for c in indices), self.ring.zero())

    def column_equation(
        self,
        column: int,
        other_column: int,
        *,
        power: int = 1,
        rows: Iterable[int] | None = None,
    ) -> LaurentPolynomial:
        indices = range(self.n) if rows is None else tuple(int(r) for r in rows)
        return sum(
            (self.column_ratio(column, other_column, r, power=power) for r in indices),
            self.ring.zero(),
        )

    def ordered_pairs(self) -> Iterable[tuple[int, int]]:
        if self.include_reverse_equations:
            for first in range(self.n):
                for second in range(self.n):
                    if first != second:
                        yield first, second
        else:
            for first in range(self.n):
                for second in range(first + 1, self.n):
                    yield first, second

    # Compatibility alias for the original prototype.
    _ordered_pairs = ordered_pairs

    def tagged_generators(self, *, power: int = 1) -> tuple[tuple[EquationTag, LaurentPolynomial], ...]:
        power = int(power)
        cached = self._tagged_generator_cache.get(power)
        if cached is not None:
            return cached
        out: list[tuple[EquationTag, LaurentPolynomial]] = []
        for row, other_row in self.ordered_pairs():
            out.append((EquationTag("row", row, other_row, power), self.row_equation(row, other_row, power=power)))
        if self.include_column_equations:
            for column, other_column in self.ordered_pairs():
                out.append(
                    (
                        EquationTag("column", column, other_column, power),
                        self.column_equation(column, other_column, power=power),
                    )
                )
        result = tuple(out)
        self._tagged_generator_cache[power] = result
        return result

    def generators(self, *, power: int = 1) -> tuple[LaurentPolynomial, ...]:
        return tuple(polynomial for _, polynomial in self.tagged_generators(power=power))

    def ideal(self, *, power: int = 1) -> LaurentIdeal:
        return LaurentIdeal(self.ring, self.generators(power=power), hint=f"hadamard-order-{self.n}")

    def square_cover_generators(self, sign_pattern, *, include_original: bool = True) -> tuple[LaurentPolynomial, ...]:
        """Construct ``I(x) + I(epsilon*x^2)`` for one global sign pattern."""
        if sign_pattern.n != self.n:
            raise ValueError("sign pattern has the wrong order")
        out = list(self.generators()) if include_original else []
        for tag, _ in self.tagged_generators(power=2):
            if tag.kind == "row":
                terms = []
                for column in range(self.n):
                    sign = -1 if sign_pattern.bit(tag.first, column) ^ sign_pattern.bit(tag.second, column) else 1
                    terms.append(sign * self.ratio(tag.first, tag.second, column, power=2))
            else:
                terms = []
                for row in range(self.n):
                    sign = -1 if sign_pattern.bit(row, tag.first) ^ sign_pattern.bit(row, tag.second) else 1
                    terms.append(sign * self.column_ratio(tag.first, tag.second, row, power=2))
            out.append(sum(terms, self.ring.zero()))
        return tuple(out)

    @staticmethod
    def _validate_permutation(permutation: Sequence[int], n: int) -> tuple[int, ...]:
        value = tuple(int(x) for x in permutation)
        if len(value) != n or set(value) != set(range(n)):
            raise ValueError(f"not a permutation of range({n})")
        return value

    def variable_substitution_for_dephased_symmetry(
        self,
        row_permutation: Sequence[int],
        column_permutation: Sequence[int],
        *,
        transpose: bool = False,
        invert: bool = False,
    ) -> tuple[LaurentPolynomial, ...]:
        """Monomial substitution induced by equivalence followed by re-dephasing."""
        rows = self._validate_permutation(row_permutation, self.n)
        columns = self._validate_permutation(column_permutation, self.n)

        def source(a: int, b: int) -> LaurentPolynomial:
            value = self.h(b, a) if transpose else self.h(a, b)
            return value.inverse_monomial() if invert else value

        row_zero, column_zero = rows[0], columns[0]
        substitutions: list[LaurentPolynomial] = []
        for row in range(1, self.n):
            for column in range(1, self.n):
                numerator = source(rows[row], columns[column]) * source(row_zero, column_zero)
                denominator = source(rows[row], column_zero) * source(row_zero, columns[column])
                substitutions.append(numerator / denominator)
        return tuple(substitutions)

    def apply_dephased_symmetry(
        self,
        polynomial: LaurentPolynomial,
        row_permutation: Sequence[int],
        column_permutation: Sequence[int],
        *,
        transpose: bool = False,
        invert: bool = False,
    ) -> LaurentPolynomial:
        substitutions = self.variable_substitution_for_dephased_symmetry(
            row_permutation,
            column_permutation,
            transpose=transpose,
            invert=invert,
        )
        return polynomial.substitute_variables(substitutions)

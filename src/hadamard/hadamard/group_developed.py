"""Laurent models for finite-group-developed Hadamard matrices.

A function ``a: G -> C^*`` develops the matrix

    H[x, y] = a(x^{-1} y).

After the harmless normalization ``a(e) = 1``, only ``|G|-1`` torus
coordinates remain.  All row orthogonality conditions are translates of

    sum(a(hg) / a(g) for g in G) = 0,    h != e.

Consequently the complete group-developed Hadamard ideal has ``|G|-1``
generators of weight ``|G|``.  The implementation below deliberately stops at
constructing that Laurent ideal; it does not choose covers or invoke the
torsion-closure orchestration.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Any

from ..coefficients.cyclotomic import CyclotomicElement, CyclotomicField
from ..laurent.ideal import LaurentIdeal, SupportProfile, support_profile
from ..laurent.polynomial import LaurentPolynomial, LaurentRing


@dataclass(frozen=True, slots=True)
class FiniteGroupTable:
    """A finite group presented by a checked multiplication table.

    Elements are the integers ``0, ..., order-1``.  The identity is detected
    from the table, so a caller is not required to number it by zero.
    Associativity is checked exhaustively: these tables are intended for the
    small symmetry groups for which group-developed classification is useful.
    """

    table: tuple[tuple[int, ...], ...]
    name: str | None = None
    identity: int = field(init=False)
    _inverses: tuple[int, ...] = field(init=False, repr=False)

    def __post_init__(self) -> None:
        table = tuple(tuple(int(value) for value in row) for row in self.table)
        order = len(table)
        if order == 0:
            raise ValueError("a finite-group table must be nonempty")
        if any(len(row) != order for row in table):
            raise ValueError("a finite-group multiplication table must be square")
        if any(not 0 <= value < order for row in table for value in row):
            raise ValueError("group-table entries must index elements of the table")

        elements = range(order)
        identities = tuple(
            candidate
            for candidate in elements
            if all(
                table[candidate][element] == element
                and table[element][candidate] == element
                for element in elements
            )
        )
        if len(identities) != 1:
            raise ValueError("multiplication table does not have a unique two-sided identity")
        identity = identities[0]

        for left in elements:
            for middle in elements:
                left_middle = table[left][middle]
                for right in elements:
                    if table[left_middle][right] != table[left][table[middle][right]]:
                        raise ValueError(
                            "multiplication table is not associative at "
                            f"({left}, {middle}, {right})"
                        )

        inverses: list[int] = []
        for element in elements:
            candidates = tuple(
                candidate
                for candidate in elements
                if table[element][candidate] == identity
                and table[candidate][element] == identity
            )
            if len(candidates) != 1:
                raise ValueError(f"element {element} does not have a unique two-sided inverse")
            inverses.append(candidates[0])

        object.__setattr__(self, "table", table)
        object.__setattr__(self, "identity", identity)
        object.__setattr__(self, "_inverses", tuple(inverses))
        object.__setattr__(self, "name", str(self.name or f"G{order}"))

    @classmethod
    def cyclic(cls, order: int) -> FiniteGroupTable:
        """Return the additive cyclic group ``C_order``."""

        order = int(order)
        if order < 1:
            raise ValueError("cyclic-group order must be positive")
        return cls(
            tuple(
                tuple((left + right) % order for right in range(order))
                for left in range(order)
            ),
            name=f"C{order}",
        )

    @property
    def order(self) -> int:
        return len(self.table)

    @property
    def elements(self) -> tuple[int, ...]:
        return tuple(range(self.order))

    def multiply(self, left: int, right: int) -> int:
        left = int(left)
        right = int(right)
        if not 0 <= left < self.order or not 0 <= right < self.order:
            raise IndexError((left, right))
        return self.table[left][right]

    def inverse(self, element: int) -> int:
        element = int(element)
        if not 0 <= element < self.order:
            raise IndexError(element)
        return self._inverses[element]


def cyclic_group(order: int) -> FiniteGroupTable:
    """Convenience constructor for ``C_order``."""

    return FiniteGroupTable.cyclic(order)


class GroupDevelopedTorus:
    """The normalized Laurent torus of matrices developed over a finite group."""

    def __init__(
        self,
        group: FiniteGroupTable,
        *,
        field: CyclotomicField | None = None,
        term_order: str = "degrevlex",
    ) -> None:
        if not isinstance(group, FiniteGroupTable):
            raise TypeError("group must be a FiniteGroupTable")
        if group.order < 2:
            raise ValueError("a group-developed Hadamard model must have order at least 2")
        self.group = group
        self.n = group.order
        self.field = field or CyclotomicField.rational()
        self.term_order = str(term_order)
        self.coordinate_elements = tuple(
            element for element in group.elements if element != group.identity
        )
        self.ring = LaurentRing(
            self.field,
            tuple(f"a_{element}" for element in self.coordinate_elements),
            self.term_order,
        )
        self.variables = self.ring.variables()
        self._index = {
            element: index for index, element in enumerate(self.coordinate_elements)
        }
        self._generator_cache: dict[int, tuple[LaurentPolynomial, ...]] = {}

    @classmethod
    def cyclic(
        cls,
        order: int,
        *,
        field: CyclotomicField | None = None,
        term_order: str = "degrevlex",
    ) -> GroupDevelopedTorus:
        """Construct the group-developed torus for ``C_order``."""

        return cls(
            FiniteGroupTable.cyclic(order),
            field=field,
            term_order=term_order,
        )

    @property
    def ambient_dimension(self) -> int:
        return self.n - 1

    @property
    def generator_shifts(self) -> tuple[int, ...]:
        return self.coordinate_elements

    def variable_index(self, element: int) -> int:
        element = int(element)
        if element == self.group.identity:
            raise ValueError("the normalized identity coordinate is constant")
        try:
            return self._index[element]
        except KeyError as exc:
            raise IndexError(element) from exc

    def phase(self, element: int) -> LaurentPolynomial:
        """Return ``a(element)``, with ``a(identity) = 1``."""

        element = int(element)
        if not 0 <= element < self.n:
            raise IndexError(element)
        if element == self.group.identity:
            return self.ring.one()
        return self.variables[self.variable_index(element)]

    # The short notation mirrors the mathematical function a: G -> C^*.
    a = phase

    def matrix_entry(self, row: int, column: int) -> LaurentPolynomial:
        """Return ``H[row,column] = a(row^-1 * column)``."""

        quotient = self.group.multiply(self.group.inverse(row), column)
        return self.phase(quotient)

    def row_inner_product(
        self,
        first: int,
        second: int,
        *,
        power: int = 1,
    ) -> LaurentPolynomial:
        """Return the formal row quotient sum for two distinct rows."""

        first = int(first)
        second = int(second)
        if first == second:
            raise ValueError("row indices must be distinct")
        return sum(
            (
                (self.matrix_entry(first, column) / self.matrix_entry(second, column))
                ** int(power)
                for column in self.group.elements
            ),
            self.ring.zero(),
        )

    def correlation(self, shift: int, *, power: int = 1) -> LaurentPolynomial:
        """Return ``sum_g (a(shift*g)/a(g))**power``.

        For a nonidentity ``shift`` this is the row inner product of the
        identity row with the ``shift`` row.  Those ``|G|-1`` correlations
        therefore impose all row orthogonality conditions.
        """

        shift = int(shift)
        if not 0 <= shift < self.n:
            raise IndexError(shift)
        return sum(
            (
                (
                    self.phase(self.group.multiply(shift, element))
                    / self.phase(element)
                )
                ** int(power)
                for element in self.group.elements
            ),
            self.ring.zero(),
        )

    def tagged_generators(
        self,
        *,
        power: int = 1,
    ) -> tuple[tuple[int, LaurentPolynomial], ...]:
        return tuple(
            (shift, polynomial)
            for shift, polynomial in zip(
                self.generator_shifts,
                self.generators(power=power),
                strict=True,
            )
        )

    def generators(self, *, power: int = 1) -> tuple[LaurentPolynomial, ...]:
        power = int(power)
        cached = self._generator_cache.get(power)
        if cached is None:
            cached = tuple(
                self.correlation(shift, power=power)
                for shift in self.generator_shifts
            )
            self._generator_cache[power] = cached
        return cached

    def ideal(self, *, power: int = 1) -> LaurentIdeal:
        return LaurentIdeal(
            self.ring,
            self.generators(power=power),
            hint=f"group-developed-{self.group.name}",
        )

    def support_profile(self, *, power: int = 1) -> SupportProfile:
        return support_profile(self.generators(power=power), self.ambient_dimension)

    @property
    def profile(self) -> SupportProfile:
        return self.support_profile()

    def evaluate_correlations(
        self,
        values: Sequence[Any],
        *,
        power: int = 1,
    ) -> tuple[CyclotomicElement, ...]:
        """Evaluate every nontrivial correlation using exact field arithmetic."""

        return tuple(
            polynomial.evaluate(values) for polynomial in self.generators(power=power)
        )

    def is_solution(self, values: Sequence[Any], *, power: int = 1) -> bool:
        return all(
            value.is_zero
            for value in self.evaluate_correlations(values, power=power)
        )


__all__ = ["FiniteGroupTable", "GroupDevelopedTorus", "cyclic_group"]

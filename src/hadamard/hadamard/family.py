"""Machine-readable monomial families of complex Hadamard matrices.

A family entry has the form

    omega[r,c] * t_1**a[r,c,1] * ... * t_d**a[r,c,d],

where ``omega[r,c]`` is a root of unity in an explicitly presented
cyclotomic field.  Orthogonality is verified as an exact Laurent-polynomial
identity, so no numerical sampling is involved.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from operator import index
from pathlib import Path
from typing import Any

from ..coefficients.cyclotomic import CyclotomicElement, CyclotomicField
from ..laurent.polynomial import LaurentPolynomial, LaurentRing
from ..provenance import write_json


def _matrix_shape(matrix: Sequence[Sequence[Any]]) -> tuple[int, int]:
    rows = len(matrix)
    columns = len(matrix[0]) if rows else 0
    if any(len(row) != columns for row in matrix):
        raise ValueError("matrix rows have unequal length")
    return rows, columns


@dataclass(frozen=True, slots=True)
class FamilyVerification:
    valid: bool
    dephased: bool
    row_failures: tuple[tuple[int, int, LaurentPolynomial], ...]
    column_failures: tuple[tuple[int, int, LaurentPolynomial], ...]

    def to_json(self) -> dict[str, Any]:
        return {
            "valid": self.valid,
            "dephased": self.dephased,
            "row_failures": [
                {"first": first, "second": second, "polynomial": polynomial.to_json()}
                for first, second, polynomial in self.row_failures
            ],
            "column_failures": [
                {"first": first, "second": second, "polynomial": polynomial.to_json()}
                for first, second, polynomial in self.column_failures
            ],
        }


@dataclass(frozen=True, slots=True)
class MonomialHadamardFamily:
    """A dephased square matrix with monomial torus parameters."""

    name: str
    field: CyclotomicField
    parameter_names: tuple[str, ...]
    coefficients: tuple[tuple[CyclotomicElement, ...], ...]
    exponents: tuple[tuple[tuple[int, ...], ...], ...]
    provenance: Mapping[str, Any] | None = None

    def __post_init__(self) -> None:
        n, columns = _matrix_shape(self.coefficients)
        if n == 0 or n != columns:
            raise ValueError("Hadamard family must be a nonempty square matrix")
        exponent_rows, exponent_columns = _matrix_shape(self.exponents)
        if (exponent_rows, exponent_columns) != (n, n):
            raise ValueError("coefficient and exponent matrices have different shapes")
        parameter_count = len(self.parameter_names)
        for row in self.coefficients:
            for coefficient in row:
                if coefficient.field != self.field:
                    raise ValueError("all coefficients must use the declared cyclotomic field")
                if coefficient.multiplicative_order() is None:
                    raise ValueError("family coefficients must be roots of unity")
        normalized_exponents = []
        for row in self.exponents:
            normalized_row = []
            for exponent in row:
                if len(exponent) != parameter_count:
                    raise ValueError("entry exponent has the wrong parameter dimension")
                normalized_exponent = []
                for value in exponent:
                    if isinstance(value, bool):
                        raise ValueError("family exponents must be exact integers")
                    try:
                        normalized_exponent.append(index(value))
                    except TypeError as exc:
                        raise ValueError("family exponents must be exact integers") from exc
                normalized_row.append(tuple(normalized_exponent))
            normalized_exponents.append(tuple(normalized_row))
        object.__setattr__(self, "exponents", tuple(normalized_exponents))

    @property
    def order(self) -> int:
        return len(self.coefficients)

    @property
    def parameter_dimension(self) -> int:
        return len(self.parameter_names)

    @property
    def ring(self) -> LaurentRing:
        return LaurentRing(self.field, self.parameter_names)

    def entry(self, row: int, column: int) -> LaurentPolynomial:
        return self.ring.monomial(self.exponents[row][column], self.coefficients[row][column])

    def row_inner_product(self, first: int, second: int) -> LaurentPolynomial:
        ring = self.ring
        return sum(
            (
                ring.monomial(
                    tuple(
                        a - b
                        for a, b in zip(
                            self.exponents[first][column],
                            self.exponents[second][column],
                            strict=True,
                        )
                    ),
                    self.coefficients[first][column] / self.coefficients[second][column],
                )
                for column in range(self.order)
            ),
            ring.zero(),
        )

    def column_inner_product(self, first: int, second: int) -> LaurentPolynomial:
        ring = self.ring
        return sum(
            (
                ring.monomial(
                    tuple(
                        a - b
                        for a, b in zip(
                            self.exponents[row][first],
                            self.exponents[row][second],
                            strict=True,
                        )
                    ),
                    self.coefficients[row][first] / self.coefficients[row][second],
                )
                for row in range(self.order)
            ),
            ring.zero(),
        )

    def is_dephased(self) -> bool:
        zero = (0,) * self.parameter_dimension
        one = self.field.one()
        return all(
            self.coefficients[0][index] == one
            and self.coefficients[index][0] == one
            and self.exponents[0][index] == zero
            and self.exponents[index][0] == zero
            for index in range(self.order)
        )

    def verify(self, *, check_columns: bool = True) -> FamilyVerification:
        row_failures = []
        for first in range(self.order):
            for second in range(first + 1, self.order):
                polynomial = self.row_inner_product(first, second)
                if polynomial:
                    row_failures.append((first, second, polynomial))
        column_failures = []
        if check_columns:
            for first in range(self.order):
                for second in range(first + 1, self.order):
                    polynomial = self.column_inner_product(first, second)
                    if polynomial:
                        column_failures.append((first, second, polynomial))
        dephased = self.is_dephased()
        return FamilyVerification(
            valid=dephased and not row_failures and not column_failures,
            dephased=dephased,
            row_failures=tuple(row_failures),
            column_failures=tuple(column_failures),
        )

    def to_json(self) -> dict[str, Any]:
        return {
            "schema": "hadamard.monomial-family.v1",
            "name": self.name,
            "field_conductor": self.field.conductor,
            "parameters": list(self.parameter_names),
            "coefficients": [
                [coefficient.to_json() for coefficient in row]
                for row in self.coefficients
            ],
            "exponents": [
                [[int(value) for value in exponent] for exponent in row]
                for row in self.exponents
            ],
            "provenance": None if self.provenance is None else dict(self.provenance),
        }

    @classmethod
    def from_json(cls, payload: Mapping[str, Any]) -> MonomialHadamardFamily:
        schema = payload.get("schema")
        if schema not in {None, "hadamard.monomial-family.v1"}:
            raise ValueError(f"unsupported family schema: {schema!r}")
        field = CyclotomicField(int(payload["field_conductor"]))
        coefficients = tuple(
            tuple(CyclotomicElement.from_json(value) for value in row)
            for row in payload["coefficients"]
        )
        exponents = tuple(
            tuple(tuple(exponent) for exponent in row)
            for row in payload["exponents"]
        )
        return cls(
            name=str(payload.get("name", "unnamed")),
            field=field,
            parameter_names=tuple(str(value) for value in payload.get("parameters", ())),
            coefficients=coefficients,
            exponents=exponents,
            provenance=payload.get("provenance"),
        )

    @classmethod
    def read(cls, path: str | Path) -> MonomialHadamardFamily:
        import json

        with Path(path).open("r", encoding="utf-8") as handle:
            return cls.from_json(json.load(handle))

    def write(self, path: str | Path) -> None:
        write_json(path, self.to_json())

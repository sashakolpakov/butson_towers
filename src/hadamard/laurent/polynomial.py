"""Sparse Laurent polynomials over explicit cyclotomic fields."""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any

from ..coefficients.cyclotomic import CyclotomicElement, CyclotomicField
from ..types import Exponent


def _validate_exponent(exponent: Sequence[int], dimension: int) -> Exponent:
    value = tuple(int(x) for x in exponent)
    if len(value) != dimension:
        raise ValueError(f"expected exponent length {dimension}, got {len(value)}")
    return value


def _zero_exponent(dimension: int) -> Exponent:
    return (0,) * dimension


@dataclass(frozen=True, slots=True)
class LaurentRing:
    """A Laurent polynomial ring with named variables."""

    field: CyclotomicField
    names: tuple[str, ...]
    order: str = "degrevlex"

    def __post_init__(self) -> None:
        names = tuple(str(name) for name in self.names)
        if len(set(names)) != len(names):
            raise ValueError("variable names must be distinct")
        object.__setattr__(self, "names", names)

    @classmethod
    def rational(cls, names: Sequence[str], *, order: str = "degrevlex") -> LaurentRing:
        return cls(CyclotomicField.rational(), tuple(names), order)

    @property
    def dimension(self) -> int:
        return len(self.names)

    def zero(self) -> LaurentPolynomial:
        return LaurentPolynomial(self, ())

    def one(self) -> LaurentPolynomial:
        return self.constant(1)

    def constant(self, value: Any) -> LaurentPolynomial:
        coefficient = self.field.coerce(value)
        if not coefficient:
            return self.zero()
        return LaurentPolynomial(self, ((_zero_exponent(self.dimension), coefficient),))

    def monomial(self, exponent: Sequence[int], coefficient: Any = 1) -> LaurentPolynomial:
        exp = _validate_exponent(exponent, self.dimension)
        coeff = self.field.coerce(coefficient)
        if not coeff:
            return self.zero()
        return LaurentPolynomial(self, ((exp, coeff),))

    def variable(self, index: int) -> LaurentPolynomial:
        if not 0 <= int(index) < self.dimension:
            raise IndexError(index)
        exponent = [0] * self.dimension
        exponent[int(index)] = 1
        return self.monomial(exponent)

    def variables(self) -> tuple[LaurentPolynomial, ...]:
        return tuple(self.variable(i) for i in range(self.dimension))

    def polynomial(self, terms: Mapping[Sequence[int], Any] | Iterable[tuple[Sequence[int], Any]]) -> LaurentPolynomial:
        items = terms.items() if isinstance(terms, Mapping) else terms
        accumulator: dict[Exponent, CyclotomicElement] = {}
        for exponent, coefficient in items:
            exp = _validate_exponent(exponent, self.dimension)
            coeff = self.field.coerce(coefficient)
            if not coeff:
                continue
            accumulator[exp] = accumulator.get(exp, self.field.zero()) + coeff
            if not accumulator[exp]:
                del accumulator[exp]
        return LaurentPolynomial(self, tuple(sorted(accumulator.items())))

    def change_field(self, field: CyclotomicField) -> LaurentRing:
        return LaurentRing(field, self.names, self.order)

    def __str__(self) -> str:
        return f"{self.field}[{', '.join(name + '^±1' for name in self.names)}]"


@dataclass(frozen=True, slots=True)
class LaurentPolynomial:
    """An immutable sparse Laurent polynomial.

    ``items`` is a sorted tuple of ``(exponent, coefficient)`` pairs.  Public
    construction should normally go through :class:`LaurentRing`.
    """

    ring: LaurentRing
    items: tuple[tuple[Exponent, CyclotomicElement], ...]

    def __post_init__(self) -> None:
        accumulator: dict[Exponent, CyclotomicElement] = {}
        for exponent, coefficient in self.items:
            exp = _validate_exponent(exponent, self.ring.dimension)
            coeff = self.ring.field.coerce(coefficient)
            if not coeff:
                continue
            accumulator[exp] = accumulator.get(exp, self.ring.field.zero()) + coeff
            if not accumulator[exp]:
                del accumulator[exp]
        object.__setattr__(self, "items", tuple(sorted(accumulator.items())))

    @property
    def terms(self) -> dict[Exponent, CyclotomicElement]:
        return dict(self.items)

    @property
    def support(self) -> tuple[Exponent, ...]:
        return tuple(exponent for exponent, _ in self.items)

    @property
    def coefficients(self) -> tuple[CyclotomicElement, ...]:
        return tuple(coefficient for _, coefficient in self.items)

    @property
    def number_of_terms(self) -> int:
        return len(self.items)

    @property
    def is_zero(self) -> bool:
        return not self.items

    @property
    def is_monomial(self) -> bool:
        return len(self.items) == 1

    @property
    def is_unit(self) -> bool:
        return self.is_monomial and bool(self.items[0][1])

    @property
    def leading_term(self) -> tuple[Exponent, CyclotomicElement] | None:
        return self.items[-1] if self.items else None

    def coefficient(self, exponent: Sequence[int]) -> CyclotomicElement:
        return self.terms.get(_validate_exponent(exponent, self.ring.dimension), self.ring.field.zero())

    def _coerce_polynomial(self, other: Any) -> LaurentPolynomial:
        if isinstance(other, LaurentPolynomial):
            if other.ring == self.ring:
                return other
            if other.ring.names != self.ring.names:
                raise TypeError("cannot mix Laurent rings with different variables")
            return other.change_field(self.ring.field)
        return self.ring.constant(other)

    def __add__(self, other: Any) -> LaurentPolynomial:
        other_poly = self._coerce_polynomial(other)
        accumulator = self.terms
        for exponent, coefficient in other_poly.items:
            accumulator[exponent] = accumulator.get(exponent, self.ring.field.zero()) + coefficient
            if not accumulator[exponent]:
                del accumulator[exponent]
        return LaurentPolynomial(self.ring, tuple(sorted(accumulator.items())))

    def __radd__(self, other: Any) -> LaurentPolynomial:
        return self + other

    def __neg__(self) -> LaurentPolynomial:
        return LaurentPolynomial(self.ring, tuple((e, -c) for e, c in self.items))

    def __sub__(self, other: Any) -> LaurentPolynomial:
        return self + (-self._coerce_polynomial(other))

    def __rsub__(self, other: Any) -> LaurentPolynomial:
        return self._coerce_polynomial(other) - self

    def __mul__(self, other: Any) -> LaurentPolynomial:
        other_poly = self._coerce_polynomial(other)
        if self.is_zero or other_poly.is_zero:
            return self.ring.zero()
        accumulator: dict[Exponent, CyclotomicElement] = {}
        for exponent_left, coefficient_left in self.items:
            for exponent_right, coefficient_right in other_poly.items:
                exponent = tuple(a + b for a, b in zip(exponent_left, exponent_right, strict=True))
                coefficient = coefficient_left * coefficient_right
                accumulator[exponent] = accumulator.get(exponent, self.ring.field.zero()) + coefficient
                if not accumulator[exponent]:
                    del accumulator[exponent]
        return LaurentPolynomial(self.ring, tuple(sorted(accumulator.items())))

    def __rmul__(self, other: Any) -> LaurentPolynomial:
        return self * other

    def inverse_monomial(self) -> LaurentPolynomial:
        if not self.is_monomial:
            raise ValueError("only Laurent monomials are invertible as polynomials")
        exponent, coefficient = self.items[0]
        return self.ring.monomial(tuple(-x for x in exponent), coefficient.inverse())

    def __truediv__(self, other: Any) -> LaurentPolynomial:
        if isinstance(other, LaurentPolynomial):
            return self * other.inverse_monomial()
        return self * self.ring.field.coerce(other).inverse()

    def __pow__(self, exponent: int) -> LaurentPolynomial:
        exponent = int(exponent)
        if exponent < 0:
            return self.inverse_monomial() ** (-exponent)
        result = self.ring.one()
        base = self
        while exponent:
            if exponent & 1:
                result *= base
            base *= base
            exponent >>= 1
        return result

    def change_field(self, field: CyclotomicField) -> LaurentPolynomial:
        if field == self.ring.field:
            return self
        target_ring = self.ring.change_field(field)
        return LaurentPolynomial(
            target_ring,
            tuple((exponent, coefficient.embed_into(field)) for exponent, coefficient in self.items),
        )

    def map_coefficients(self, function) -> LaurentPolynomial:
        mapped = tuple((exponent, function(coefficient)) for exponent, coefficient in self.items)
        if not mapped:
            return self
        target_field = mapped[0][1].field
        target_ring = self.ring.change_field(target_field)
        return LaurentPolynomial(target_ring, mapped)

    def apply_galois(self, exponent: int) -> LaurentPolynomial:
        return self.map_coefficients(lambda coefficient: coefficient.galois(exponent))

    def rescale_variables(self, scalars: Mapping[int, Any]) -> LaurentPolynomial:
        normalized = {int(i): self.ring.field.coerce(value) for i, value in scalars.items()}
        items: list[tuple[Exponent, CyclotomicElement]] = []
        for exponent, coefficient in self.items:
            factor = self.ring.field.one()
            for index, scalar in normalized.items():
                if not 0 <= index < self.ring.dimension:
                    raise IndexError(index)
                if exponent[index]:
                    factor *= scalar ** exponent[index]
            items.append((exponent, coefficient * factor))
        return LaurentPolynomial(self.ring, tuple(items))

    # Compatibility spelling used in the previous Sage frontend.
    rescale_vars = rescale_variables

    def toric_coordinate_change(self, matrix: Sequence[Sequence[int]]) -> LaurentPolynomial:
        rows = tuple(tuple(int(value) for value in row) for row in matrix)
        if len(rows) != self.ring.dimension or any(len(row) != self.ring.dimension for row in rows):
            raise ValueError("toric coordinate matrix must be square of ambient dimension")
        transformed: list[tuple[Exponent, CyclotomicElement]] = []
        for exponent, coefficient in self.items:
            new_exponent = tuple(
                sum(exponent[i] * rows[i][j] for i in range(self.ring.dimension))
                for j in range(self.ring.dimension)
            )
            transformed.append((new_exponent, coefficient))
        return LaurentPolynomial(self.ring, tuple(transformed))

    def power_map(self, multiplier: int) -> LaurentPolynomial:
        multiplier = int(multiplier)
        return LaurentPolynomial(
            self.ring,
            tuple((tuple(multiplier * value for value in exponent), coefficient) for exponent, coefficient in self.items),
        )

    def invert_variables(self) -> LaurentPolynomial:
        return LaurentPolynomial(
            self.ring,
            tuple((tuple(-value for value in exponent), coefficient) for exponent, coefficient in self.items),
        )

    def decimate(self, variable_index: int, modulus: int, *, require_divisible: bool = False) -> LaurentPolynomial:
        variable_index = int(variable_index)
        modulus = int(modulus)
        if modulus <= 0:
            raise ValueError("decimation modulus must be positive")
        transformed = []
        for exponent, coefficient in self.items:
            if require_divisible and exponent[variable_index] % modulus:
                raise ValueError("exponent is not divisible by the decimation modulus")
            values = list(exponent)
            values[variable_index] //= modulus
            transformed.append((tuple(values), coefficient))
        return LaurentPolynomial(self.ring, tuple(transformed))

    def clear_negative_exponents(self) -> tuple[Exponent, LaurentPolynomial]:
        if self.is_zero:
            return _zero_exponent(self.ring.dimension), self
        shift = tuple(-min(0, min(exponent[i] for exponent in self.support)) for i in range(self.ring.dimension))
        polynomial = self * self.ring.monomial(shift)
        return shift, polynomial

    def monomial_content(self) -> Exponent:
        if self.is_zero:
            return _zero_exponent(self.ring.dimension)
        return tuple(min(exponent[i] for exponent in self.support) for i in range(self.ring.dimension))

    def normalized(self) -> LaurentPolynomial:
        """Normalize up to multiplication by a Laurent monomial and scalar."""
        if self.is_zero:
            return self
        base_exponent, base_coefficient = self.items[0]
        inverse = base_coefficient.inverse()
        return LaurentPolynomial(
            self.ring,
            tuple(
                (
                    tuple(
                        value - base
                        for value, base in zip(exponent, base_exponent, strict=True)
                    ),
                    coefficient * inverse,
                )
                for exponent, coefficient in self.items
            ),
        )

    def substitute_variables(self, images: Sequence[LaurentPolynomial]) -> LaurentPolynomial:
        if len(images) != self.ring.dimension:
            raise ValueError("wrong number of variable images")
        if any(image.ring != self.ring for image in images):
            raise ValueError("all substitution images must lie in the same ring")
        result = self.ring.zero()
        for exponent, coefficient in self.items:
            term = self.ring.constant(coefficient)
            for image, power in zip(images, exponent, strict=True):
                if power:
                    term *= image**power
            result += term
        return result

    def evaluate(self, values: Sequence[Any]) -> CyclotomicElement:
        if len(values) != self.ring.dimension:
            raise ValueError("wrong number of values")
        coerced = tuple(self.ring.field.coerce(value) for value in values)
        result = self.ring.field.zero()
        for exponent, coefficient in self.items:
            term = coefficient
            for value, power in zip(coerced, exponent, strict=True):
                if power:
                    term *= value**power
            result += term
        return result

    def to_json(self) -> dict[str, Any]:
        return {
            "field": self.ring.field.conductor,
            "variables": list(self.ring.names),
            "order": self.ring.order,
            "terms": [
                {"exponent": list(exponent), "coefficient": coefficient.to_json()}
                for exponent, coefficient in self.items
            ],
        }

    @classmethod
    def from_json(cls, payload: Mapping[str, Any]) -> LaurentPolynomial:
        field = CyclotomicField(int(payload["field"]))
        ring = LaurentRing(field, tuple(payload["variables"]), str(payload.get("order", "degrevlex")))
        return ring.polynomial(
            (
                tuple(item["exponent"]),
                CyclotomicElement.from_json(item["coefficient"]),
            )
            for item in payload["terms"]
        )

    def _monomial_string(self, exponent: Exponent) -> str:
        factors = []
        for name, power in zip(self.ring.names, exponent, strict=True):
            if power == 0:
                continue
            factors.append(name if power == 1 else f"{name}^{power}")
        return "*".join(factors) or "1"

    def __bool__(self) -> bool:
        return not self.is_zero

    def __str__(self) -> str:
        if self.is_zero:
            return "0"
        chunks: list[str] = []
        for exponent, coefficient in self.items:
            monomial = self._monomial_string(exponent)
            if monomial == "1":
                body = str(abs(coefficient.rational_value)) if coefficient.is_rational else str(coefficient)
            elif coefficient == self.ring.field.one() or coefficient == -self.ring.field.one():
                body = monomial
            else:
                body = f"({coefficient})*{monomial}"
            negative = coefficient.is_rational and coefficient.rational_value < 0 or coefficient == -self.ring.field.one()
            if not chunks:
                chunks.append(f"-{body}" if negative and not body.startswith("-") else body)
            else:
                chunks.append((" - " if negative else " + ") + (body.lstrip("-") if negative else body))
        return "".join(chunks)

    def __repr__(self) -> str:
        return f"LaurentPolynomial({self})"

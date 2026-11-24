r"""Exact elements of :math:`\mathbf Q(\zeta_N)` without SageMath.

The representation is deliberately small and auditable: an element is a rational
polynomial in a distinguished primitive ``N``-th root, reduced modulo the
cyclotomic polynomial ``Phi_N``.  SymPy is used only for univariate exact
polynomial arithmetic; multivariate ideal algebra is delegated to Singular.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from fractions import Fraction
from functools import cache, lru_cache
from math import gcd
from typing import Any

import sympy as sp

from ..types import lcm

_X = sp.Symbol("x")


def _as_fraction(value: Any) -> Fraction:
    if isinstance(value, Fraction):
        return value
    if isinstance(value, int):
        return Fraction(value)
    if isinstance(value, sp.Integer):
        return Fraction(int(value))
    if isinstance(value, sp.Rational):
        return Fraction(int(value.p), int(value.q))
    raise TypeError(f"cannot coerce {type(value).__name__} to an exact rational")


def _trim(values: list[Fraction]) -> list[Fraction]:
    while values and values[-1] == 0:
        values.pop()
    return values


@cache
def _modulus_ascending(conductor: int) -> tuple[Fraction, ...]:
    polynomial = sp.Poly(sp.cyclotomic_poly(conductor, _X), _X, domain=sp.QQ)
    return tuple(_as_fraction(c) for c in reversed(polynomial.all_coeffs()))


@cache
def _degree(conductor: int) -> int:
    return len(_modulus_ascending(conductor)) - 1


def _to_sympy_poly(coefficients: Iterable[Fraction]) -> sp.Poly:
    expression = sum(
        sp.Rational(c.numerator, c.denominator) * _X**i for i, c in enumerate(coefficients)
    )
    return sp.Poly(expression, _X, domain=sp.QQ)


def _from_sympy_poly(polynomial: sp.Poly, degree: int) -> tuple[Fraction, ...]:
    out = [Fraction(0)] * degree
    for (exponent,), coefficient in polynomial.terms():
        if exponent >= degree:
            raise ValueError("polynomial was not reduced modulo the cyclotomic polynomial")
        out[exponent] = _as_fraction(coefficient)
    return tuple(out)


def _reduce(conductor: int, coefficients: Iterable[Fraction]) -> tuple[Fraction, ...]:
    """Reduce coefficient lists modulo the monic cyclotomic polynomial.

    This is the coefficient-arithmetic hot path.  Direct fraction long
    division avoids constructing two SymPy ``Poly`` objects for every field
    multiplication while retaining exactly the same normal form.
    """
    degree = _degree(conductor)
    # A coefficient vector of exactly the field degree already represents a
    # polynomial of degree < deg(Phi_N), hence is in canonical normal form.
    # Arithmetic methods overwhelmingly pass tuples of Fractions in precisely
    # this form.  Returning them directly avoids a second normalization pass
    # after multiplication and all normalization work after addition, without
    # changing the represented element.
    if (
        isinstance(coefficients, tuple)
        and len(coefficients) == degree
        and all(isinstance(value, Fraction) for value in coefficients)
    ):
        return coefficients
    modulus = _modulus_ascending(conductor)
    values = _trim([_as_fraction(value) for value in coefficients])
    while len(values) > degree:
        leading = values[-1] / modulus[-1]
        shift = len(values) - len(modulus)
        for index, coefficient in enumerate(modulus):
            values[shift + index] -= leading * coefficient
        _trim(values)
    values.extend([Fraction(0)] * (degree - len(values)))
    return tuple(values)


@dataclass(frozen=True, slots=True)
class CyclotomicField:
    """The explicitly presented field ``Q(zeta_N)``.

    Conductors are not silently minimized.  This keeps Galois actions and
    embeddings transparent.  ``Q`` is represented by conductor 1.
    """

    conductor: int = 1

    def __post_init__(self) -> None:
        if int(self.conductor) < 1:
            raise ValueError("cyclotomic conductor must be positive")
        object.__setattr__(self, "conductor", int(self.conductor))

    @classmethod
    def rational(cls) -> CyclotomicField:
        return cls(1)

    @property
    def degree(self) -> int:
        return _degree(self.conductor)

    @property
    def modulus(self) -> tuple[Fraction, ...]:
        """Cyclotomic polynomial coefficients in ascending order."""
        return _modulus_ascending(self.conductor)

    @property
    def root_group_order(self) -> int:
        """Order of the full root-of-unity group in this field."""
        if self.conductor <= 2:
            return 2
        return self.conductor if self.conductor % 2 == 0 else 2 * self.conductor

    def zero(self) -> CyclotomicElement:
        return _cached_field_zero(self.conductor)

    def one(self) -> CyclotomicElement:
        return _cached_field_one(self.conductor)

    def zeta(self) -> CyclotomicElement:
        return _cached_field_zeta(self.conductor)

    def full_root_generator(self) -> CyclotomicElement:
        return _cached_full_root_generator(self.conductor)

    def element(self, coefficients: Iterable[int | Fraction]) -> CyclotomicElement:
        return CyclotomicElement(
            self, _reduce(self.conductor, tuple(_as_fraction(c) for c in coefficients))
        )

    def coerce(self, value: Any) -> CyclotomicElement:
        if isinstance(value, CyclotomicElement):
            if value.field == self:
                return value
            return value.embed_into(self)
        rational = _as_fraction(value)
        coefficients = [Fraction(0)] * self.degree
        coefficients[0] = rational
        return CyclotomicElement(self, tuple(coefficients))

    def common_field(self, other: CyclotomicField) -> CyclotomicField:
        if self == other:
            return self
        return CyclotomicField(lcm(self.conductor, other.conductor))

    def defining_polynomial_string(self, symbol: str = "a", power_operator: str = "^") -> str:
        terms: list[str] = []
        for exponent, coefficient in enumerate(self.modulus):
            if coefficient == 0:
                continue
            abs_coefficient = abs(coefficient)
            if exponent == 0:
                monomial = ""
            elif exponent == 1:
                monomial = symbol
            else:
                monomial = f"{symbol}{power_operator}{exponent}"
            if exponent and abs_coefficient == 1:
                body = monomial
            elif monomial:
                body = f"{abs_coefficient}*{monomial}"
            else:
                body = str(abs_coefficient)
            if not terms:
                terms.append(body if coefficient > 0 else f"-{body}")
            else:
                terms.append(("+" if coefficient > 0 else "-") + body)
        return "".join(terms) or "0"

    def __str__(self) -> str:
        return "Q" if self.conductor == 1 else f"Q(zeta_{self.conductor})"


@dataclass(frozen=True, slots=True)
class CyclotomicElement:
    """An exact element in a specified :class:`CyclotomicField`."""

    field: CyclotomicField
    coefficients: tuple[Fraction, ...]

    def __post_init__(self) -> None:
        reduced = _reduce(self.field.conductor, tuple(_as_fraction(c) for c in self.coefficients))
        object.__setattr__(self, "coefficients", reduced)

    @property
    def is_zero(self) -> bool:
        return not any(self.coefficients)

    @property
    def is_one(self) -> bool:
        return self == self.field.one()

    @property
    def is_rational(self) -> bool:
        return not any(self.coefficients[1:])

    @property
    def rational_value(self) -> Fraction:
        if not self.is_rational:
            raise ValueError("element is not rational")
        return self.coefficients[0]

    @lru_cache(maxsize=32_768)  # noqa: B019 - bounded immutable-value memoization
    def embed_into(self, target: CyclotomicField) -> CyclotomicElement:
        if target == self.field:
            return self
        if target.conductor % self.field.conductor != 0:
            raise ValueError(
                f"only explicit conductor-divisibility embeddings are supported: "
                f"{self.field.conductor} does not divide {target.conductor}"
            )
        image = target.zeta() ** (target.conductor // self.field.conductor)
        result = target.zero()
        power = target.one()
        for coefficient in self.coefficients:
            result += target.coerce(coefficient) * power
            power *= image
        return result

    def _common(self, other: Any) -> tuple[CyclotomicElement, CyclotomicElement]:
        if isinstance(other, CyclotomicElement):
            target = self.field.common_field(other.field)
            return self.embed_into(target), other.embed_into(target)
        return self, self.field.coerce(other)

    def __add__(self, other: Any) -> CyclotomicElement:
        try:
            left, right = self._common(other)
        except TypeError:
            return NotImplemented
        if left.is_zero:
            return right
        if right.is_zero:
            return left
        return CyclotomicElement(
            left.field,
            tuple(a + b for a, b in zip(left.coefficients, right.coefficients, strict=True)),
        )

    def __radd__(self, other: Any) -> CyclotomicElement:
        return self + other

    def __neg__(self) -> CyclotomicElement:
        if self.is_zero:
            return self
        return CyclotomicElement(self.field, tuple(-c for c in self.coefficients))

    def __sub__(self, other: Any) -> CyclotomicElement:
        try:
            _left, right = self._common(other)
        except TypeError:
            return NotImplemented
        return self + (-right)

    def __rsub__(self, other: Any) -> CyclotomicElement:
        return (-self) + other

    def __mul__(self, other: Any) -> CyclotomicElement:
        try:
            left, right = self._common(other)
        except TypeError:
            return NotImplemented
        if left.is_rational:
            scalar = left.rational_value
            if scalar == 0:
                return left.field.zero()
            if scalar == 1:
                return right
            return CyclotomicElement(
                left.field,
                tuple(scalar * coefficient for coefficient in right.coefficients),
            )
        if right.is_rational:
            scalar = right.rational_value
            if scalar == 0:
                return left.field.zero()
            if scalar == 1:
                return left
            return CyclotomicElement(
                left.field,
                tuple(scalar * coefficient for coefficient in left.coefficients),
            )
        degree = left.field.degree
        raw = [Fraction(0)] * (2 * degree - 1)
        for i, a in enumerate(left.coefficients):
            if a == 0:
                continue
            for j, b in enumerate(right.coefficients):
                if b:
                    raw[i + j] += a * b
        return CyclotomicElement(left.field, _reduce(left.field.conductor, raw))

    def __rmul__(self, other: Any) -> CyclotomicElement:
        return self * other

    @lru_cache(maxsize=8_192)  # noqa: B019 - bounded immutable-value memoization
    def inverse(self) -> CyclotomicElement:
        if self.is_zero:
            raise ZeroDivisionError("division by zero in a cyclotomic field")
        if self.is_rational:
            return self.field.coerce(1 / self.rational_value)
        polynomial = _to_sympy_poly(self.coefficients)
        modulus = _to_sympy_poly(self.field.modulus)
        inverse = sp.invert(polynomial, modulus)
        if not isinstance(inverse, sp.Poly):
            inverse = sp.Poly(inverse, _X, domain=sp.QQ)
        return CyclotomicElement(
            self.field, _from_sympy_poly(inverse.rem(modulus), self.field.degree)
        )

    def __truediv__(self, other: Any) -> CyclotomicElement:
        try:
            left, right = self._common(other)
        except TypeError:
            return NotImplemented
        return left * right.inverse()

    def __rtruediv__(self, other: Any) -> CyclotomicElement:
        return self.field.coerce(other) / self

    def __pow__(self, exponent: int) -> CyclotomicElement:
        exponent = int(exponent)
        if exponent < 0:
            return self.inverse() ** (-exponent)
        result = self.field.one()
        base = self
        while exponent:
            if exponent & 1:
                result *= base
            base *= base
            exponent >>= 1
        return result

    @lru_cache(maxsize=32_768)  # noqa: B019 - bounded immutable-value memoization
    def galois(self, exponent: int) -> CyclotomicElement:
        """Apply ``zeta_N -> zeta_N^exponent``.

        The exponent must be prime to ``N`` unless the field is rational.
        """
        if self.field.conductor > 2 and gcd(int(exponent), self.field.conductor) != 1:
            raise ValueError("Galois exponent must be coprime to the conductor")
        image = self.field.zeta() ** int(exponent)
        result = self.field.zero()
        power = self.field.one()
        for coefficient in self.coefficients:
            result += self.field.coerce(coefficient) * power
            power *= image
        return result

    def conjugate(self) -> CyclotomicElement:
        if self.field.conductor <= 2:
            return self
        return self.galois(-1)

    @lru_cache(maxsize=32_768)  # noqa: B019 - bounded immutable-value memoization
    def multiplicative_order(self) -> int | None:
        if self.is_zero:
            return None
        order = self.field.root_group_order
        if self**order != self.field.one():
            return None
        for prime in sp.factorint(order):
            while order % prime == 0 and self ** (order // prime) == self.field.one():
                order //= prime
        return order

    @lru_cache(maxsize=32_768)  # noqa: B019 - bounded immutable-value memoization
    def root_of_unity_exponent(self) -> int | None:
        """Return ``k`` with ``self = omega^k`` for the full root generator."""
        order = self.multiplicative_order()
        if order is None:
            return None
        generator = self.field.full_root_generator()
        value = self.field.one()
        for exponent in range(self.field.root_group_order):
            if value == self:
                return exponent
            value *= generator
        return None

    def valuation_at_two_if_rational(self) -> int:
        """The normalized 2-adic valuation for a rational element."""
        value = self.rational_value
        if value == 0:
            raise ValueError("valuation of zero is infinite")

        def v2(integer: int) -> int:
            integer = abs(integer)
            count = 0
            while integer % 2 == 0:
                integer //= 2
                count += 1
            return count

        return v2(value.numerator) - v2(value.denominator)

    def to_sympy(self, symbol: sp.Symbol | None = None) -> sp.Expr:
        symbol = symbol or sp.Symbol(f"zeta{self.field.conductor}")
        return sp.expand(
            sum(
                sp.Rational(c.numerator, c.denominator) * symbol**i
                for i, c in enumerate(self.coefficients)
            )
        )

    def to_json(self) -> dict[str, Any]:
        return {
            "conductor": self.field.conductor,
            "coefficients": [[c.numerator, c.denominator] for c in self.coefficients],
        }

    @classmethod
    def from_json(cls, payload: dict[str, Any]) -> CyclotomicElement:
        field = CyclotomicField(int(payload["conductor"]))
        coefficients = tuple(Fraction(int(a), int(b)) for a, b in payload["coefficients"])
        return cls(field, coefficients)

    def singular_string(self, symbol: str = "a") -> str:
        terms: list[str] = []
        for exponent, coefficient in enumerate(self.coefficients):
            if coefficient == 0:
                continue
            abs_coefficient = abs(coefficient)
            if exponent == 0:
                monomial = ""
            elif exponent == 1:
                monomial = symbol
            else:
                monomial = f"{symbol}^{exponent}"
            if exponent and abs_coefficient == 1:
                body = monomial
            elif monomial:
                body = f"({abs_coefficient})*{monomial}"
            else:
                body = str(abs_coefficient)
            if not terms:
                terms.append(body if coefficient > 0 else f"-{body}")
            else:
                terms.append(("+" if coefficient > 0 else "-") + body)
        return "".join(terms) or "0"

    def __bool__(self) -> bool:
        return not self.is_zero

    def __str__(self) -> str:
        if self.is_zero:
            return "0"
        if self.is_rational:
            return str(self.rational_value)
        symbol = f"zeta_{self.field.conductor}"
        terms: list[str] = []
        for exponent, coefficient in enumerate(self.coefficients):
            if coefficient == 0:
                continue
            abs_coefficient = abs(coefficient)
            monomial = "" if exponent == 0 else symbol if exponent == 1 else f"{symbol}^{exponent}"
            if exponent and abs_coefficient == 1:
                body = monomial
            elif monomial:
                body = f"{abs_coefficient}*{monomial}"
            else:
                body = str(abs_coefficient)
            if not terms:
                terms.append(body if coefficient > 0 else f"-{body}")
            else:
                terms.append((" + " if coefficient > 0 else " - ") + body)
        return "".join(terms)


@lru_cache(maxsize=128)
def _cached_field_zero(conductor: int) -> CyclotomicElement:
    field = CyclotomicField(conductor)
    return CyclotomicElement(field, (Fraction(0),) * field.degree)


@lru_cache(maxsize=128)
def _cached_field_one(conductor: int) -> CyclotomicElement:
    field = CyclotomicField(conductor)
    coefficients = [Fraction(0)] * field.degree
    coefficients[0] = Fraction(1)
    return CyclotomicElement(field, tuple(coefficients))


@lru_cache(maxsize=128)
def _cached_field_zeta(conductor: int) -> CyclotomicElement:
    field = CyclotomicField(conductor)
    if conductor == 1:
        return _cached_field_one(conductor)
    coefficients = [Fraction(0)] * max(2, field.degree)
    coefficients[1] = Fraction(1)
    return CyclotomicElement(field, _reduce(conductor, coefficients))


@lru_cache(maxsize=128)
def _cached_full_root_generator(conductor: int) -> CyclotomicElement:
    if conductor <= 2:
        return -_cached_field_one(conductor)
    zeta = _cached_field_zeta(conductor)
    return zeta if conductor % 2 == 0 else -zeta

"""Laurent ideals and backend-independent structural measurements."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from typing import Any

from ..integer_lattice import matrix_rank
from ..types import Exponent
from .polynomial import LaurentPolynomial, LaurentRing


@dataclass(frozen=True, slots=True)
class SupportProfile:
    generator_count: int
    minimum_support: int
    maximum_support: int
    total_support: int
    active_exponent_rank: int


def exponent_difference_rows(polynomials: Iterable[LaurentPolynomial]) -> list[Exponent]:
    rows: list[Exponent] = []
    for polynomial in polynomials:
        support = polynomial.support
        if not support:
            continue
        base = support[0]
        for exponent in support[1:]:
            row = tuple(a - b for a, b in zip(exponent, base, strict=True))
            if any(row):
                rows.append(row)
    return rows


def active_exponent_rank(polynomials: Iterable[LaurentPolynomial], ambient_dimension: int | None = None) -> int:
    rows = exponent_difference_rows(polynomials)
    if not rows:
        return 0
    width = ambient_dimension if ambient_dimension is not None else len(rows[0])
    return matrix_rank(rows, width)


def support_profile(polynomials: Iterable[LaurentPolynomial], ambient_dimension: int | None = None) -> SupportProfile:
    values = tuple(polynomials)
    sizes = tuple(polynomial.number_of_terms for polynomial in values if polynomial)
    return SupportProfile(
        generator_count=len(values),
        minimum_support=min(sizes) if sizes else 0,
        maximum_support=max(sizes) if sizes else 0,
        total_support=sum(sizes),
        active_exponent_rank=active_exponent_rank(values, ambient_dimension),
    )


@dataclass(frozen=True, slots=True)
class LaurentIdeal:
    ring: LaurentRing
    generators: tuple[LaurentPolynomial, ...]
    hint: str | None = None
    _normalized_key: tuple[str, ...] | None = field(
        default=None,
        init=False,
        repr=False,
        compare=False,
        hash=False,
    )
    _normalized_set: frozenset[LaurentPolynomial] | None = field(
        default=None,
        init=False,
        repr=False,
        compare=False,
        hash=False,
    )

    def __post_init__(self) -> None:
        normalized: list[LaurentPolynomial] = []
        seen: set[LaurentPolynomial] = set()
        has_unit = False
        for generator in self.generators:
            if generator.ring != self.ring:
                raise ValueError("ideal generator belongs to a different ring")
            if not generator:
                continue
            has_unit = has_unit or generator.is_unit
            if generator not in seen:
                seen.add(generator)
                normalized.append(generator)
        if has_unit:
            # Every nonzero Laurent monomial is invertible, so the presence
            # of any such generator makes the whole ideal equal to (1).
            # Canonicalizing here avoids retaining and checkpointing all of
            # the other generators of an already empty branch.
            normalized = [self.ring.one()]
        object.__setattr__(self, "generators", tuple(normalized))

    @classmethod
    def create(cls, ring: LaurentRing, generators: Iterable[LaurentPolynomial], *, hint: str | None = None) -> LaurentIdeal:
        return cls(ring, tuple(generators), hint)

    @property
    def is_zero(self) -> bool:
        return not self.generators

    @property
    def is_unit(self) -> bool:
        return any(generator.is_unit for generator in self.generators)

    @property
    def active_rank(self) -> int:
        return active_exponent_rank(self.generators, self.ring.dimension)

    @property
    def profile(self) -> SupportProfile:
        return support_profile(self.generators, self.ring.dimension)

    def add(self, *generators: LaurentPolynomial) -> LaurentIdeal:
        return LaurentIdeal(self.ring, self.generators + tuple(generators), self.hint)

    def __add__(self, other: LaurentIdeal) -> LaurentIdeal:
        if self.ring != other.ring:
            raise ValueError("cannot add ideals from different rings")
        return LaurentIdeal(self.ring, self.generators + other.generators, self.hint or other.hint)

    def map(self, function, *, target_ring: LaurentRing | None = None) -> LaurentIdeal:
        mapped = tuple(function(generator) for generator in self.generators)
        ring = target_ring or (mapped[0].ring if mapped else self.ring)
        return LaurentIdeal(ring, mapped, self.hint)

    def change_field(self, field) -> LaurentIdeal:
        target_ring = self.ring.change_field(field)
        return LaurentIdeal(target_ring, tuple(generator.change_field(field) for generator in self.generators), self.hint)

    def normalized_generator_key(self) -> tuple[str, ...]:
        cached = self._normalized_key
        if cached is None:
            cached = tuple(
                sorted(
                    str(generator.normalized())
                    for generator in self.generators
                )
            )
            object.__setattr__(self, "_normalized_key", cached)
        return cached

    def normalized_generator_set(self) -> frozenset[LaurentPolynomial]:
        """Return an exact order-free key without formatting polynomials."""

        cached = self._normalized_set
        if cached is None:
            cached = frozenset(
                generator.normalized() for generator in self.generators
            )
            object.__setattr__(self, "_normalized_set", cached)
        return cached

    def to_json(self) -> dict[str, Any]:
        return {
            "ring": {
                "field": self.ring.field.conductor,
                "variables": list(self.ring.names),
                "order": self.ring.order,
            },
            "hint": self.hint,
            "generators": [generator.to_json() for generator in self.generators],
        }

    @classmethod
    def from_json(cls, payload: Mapping[str, Any]) -> LaurentIdeal:
        ring_payload = payload["ring"]
        from ..coefficients.cyclotomic import CyclotomicField

        ring = LaurentRing(
            CyclotomicField(int(ring_payload["field"])),
            tuple(ring_payload["variables"]),
            str(ring_payload.get("order", "degrevlex")),
        )
        generators = []
        for generator_payload in payload["generators"]:
            terms = []
            from ..coefficients.cyclotomic import CyclotomicElement

            for item in generator_payload["terms"]:
                terms.append((tuple(item["exponent"]), CyclotomicElement.from_json(item["coefficient"])))
            generators.append(ring.polynomial(terms))
        return cls(ring, tuple(generators), payload.get("hint"))

    def __str__(self) -> str:
        if self.is_zero:
            return f"(0) in {self.ring}"
        return "<" + ", ".join(str(generator) for generator in self.generators) + ">"

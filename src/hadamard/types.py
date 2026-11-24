"""Shared immutable type aliases and helpers."""

from __future__ import annotations

from collections.abc import Iterable
from fractions import Fraction
from math import gcd
from typing import TypeAlias

Exponent: TypeAlias = tuple[int, ...]
IntegerMatrix: TypeAlias = tuple[tuple[int, ...], ...]
RationalLike: TypeAlias = int | Fraction


def gcd_many(values: Iterable[int]) -> int:
    result = 0
    for value in values:
        result = gcd(result, abs(int(value)))
    return result


def lcm(a: int, b: int) -> int:
    if a == 0 or b == 0:
        return 0
    return abs(a // gcd(a, b) * b)


def lcm_many(values: Iterable[int]) -> int:
    result = 1
    for value in values:
        result = lcm(result, int(value))
    return result

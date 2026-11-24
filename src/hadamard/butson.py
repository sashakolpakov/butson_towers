"""General exact utilities for Butson--Hadamard parameter checks."""

from __future__ import annotations


def prime_divisors(value: int) -> tuple[int, ...]:
    """Return the distinct prime divisors of a positive integer."""
    if value < 1:
        raise ValueError("value must be positive")
    remaining = value
    divisors: list[int] = []
    candidate = 2
    while candidate * candidate <= remaining:
        if remaining % candidate == 0:
            divisors.append(candidate)
            while remaining % candidate == 0:
                remaining //= candidate
        candidate = 3 if candidate == 2 else candidate + 2
    if remaining > 1:
        divisors.append(remaining)
    return tuple(divisors)


def vanishing_sum_weight_possible(weight: int, q: int) -> bool:
    """Apply the Lam--Leung characterization for a ``q``th-root sum.

    A vanishing sum of ``weight`` many ``q``th roots exists exactly when the
    weight is a nonnegative-integer combination of the prime divisors of
    ``q``.
    """
    if weight < 0:
        raise ValueError("weight must be nonnegative")
    primes = prime_divisors(q)
    reachable = [False] * (weight + 1)
    reachable[0] = True
    for current in range(weight + 1):
        if not reachable[current]:
            continue
        for prime in primes:
            target = current + prime
            if target <= weight:
                reachable[target] = True
    return reachable[weight]


def lam_leung_obstructs_butson(order: int, q: int) -> bool:
    """Return whether pairwise row orthogonality is Lam--Leung obstructed."""
    if order < 1:
        raise ValueError("order must be positive")
    return not vanishing_sum_weight_possible(order, q)

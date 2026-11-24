"""Bounded-space exact monomial equivalence of finite exponent matrices.

The canonical key uses the same column-major convention as the older
``torsion.equivalence.canonical_butson_key`` on dephased input.  At each
column prefix, rows tied on that prefix can still be permuted freely.  Sort
the next column inside those blocks and retain every column attaining the
least possible next signature.  All other choices are lexicographically
larger and cannot contribute to the canonical form.  This avoids eagerly
traversing every row permutation and retains only one recursion path and
its tied choices; there is no unbounded orbit or memoization table.

Only independent row/column permutations and phase multipliers are included.
Neither transposition nor a coefficient Galois automorphism is added.
"""

from __future__ import annotations

import hashlib
import operator
from collections.abc import Iterable
from math import gcd

ExponentMatrix = tuple[tuple[int, ...], ...]
FiniteButsonKey = tuple[int, int, tuple[int, ...]]


def _residual_fingerprint(
    interior: ExponentMatrix,
    columns: tuple[int, ...],
    blocks: tuple[tuple[int, ...], ...],
) -> tuple[ExponentMatrix, ...]:
    """A sufficient, not necessary, test of equivalent continuation states.

    Every operation here is a permitted permutation of the remaining columns
    or of rows inside one equal-prefix block.  Therefore equal fingerprints
    imply identical sets of continuations. Unequal fingerprints merely leave
    redundant branches; they never reject a possible canonical form.
    """
    rows = blocks
    for _ in range(2):
        columns = tuple(sorted(
            columns,
            key=lambda c: tuple(interior[r][c] for block in rows for r in block),
        ))
        rows = tuple(tuple(sorted(
            block, key=lambda r: tuple(interior[r][c] for c in columns),
        )) for block in rows)
    return tuple(
        tuple(tuple(interior[r][c] for c in columns) for r in block)
        for block in rows
    )


def _exact_integer(value: object, name: str) -> int:
    if isinstance(value, bool):
        raise TypeError(f"{name} must be an exact integer, not bool")
    try:
        return operator.index(value)
    except TypeError as exc:
        raise TypeError(f"{name} must be an exact integer") from exc


def _normalized_exponents(
    matrix: Iterable[Iterable[int]], modulus: int,
) -> tuple[int, ExponentMatrix]:
    modulus = _exact_integer(modulus, "modulus")
    if modulus < 1:
        raise ValueError("modulus must be positive")
    value = tuple(
        tuple(_exact_integer(entry, "exponent") % modulus for entry in row)
        for row in matrix
    )
    if not value or any(len(row) != len(value) for row in value):
        raise ValueError("an exponent matrix must be nonempty and square")
    # Dephase before reducing the alphabet: arbitrary input row/column phases
    # must not change the minimal alphabet or the monomial class key.
    value = tuple(
        tuple((entry - row[0] - value[0][c] + value[0][0]) % modulus
              for c, entry in enumerate(row))
        for row in value
    )
    divisor = modulus
    for row in value:
        for entry in row:
            divisor = gcd(divisor, entry)
    return modulus // divisor, tuple(
        tuple(entry // divisor for entry in row) for row in value
    )


def canonical_finite_butson_key(
    matrix: Iterable[Iterable[int]], modulus: int,
) -> FiniteButsonKey:
    """Return the exact full monomial key, without transpose or Galois action.

    The result is ``(order, reduced_modulus, column_major_interior)``.  Input
    need not be dephased, and integral exponents may lie outside the given
    residue range.  Nonintegral values are rejected rather than truncated.
    The method also accepts nonorthogonal square matrices, which is useful
    for independent equivalence tests and adversarial regressions.
    """
    modulus, value = _normalized_exponents(matrix, modulus)
    n = len(value)
    width = n - 1
    if modulus == 1:
        return n, 1, (0,) * (width * width)
    best: tuple[int, ...] | None = None
    seen_pivots: set[tuple[ExponentMatrix, ...]] = set()

    for first_row in range(n):
        row_indices = tuple(r for r in range(n) if r != first_row)
        for first_column in range(n):
            column_indices = tuple(c for c in range(n) if c != first_column)
            interior = tuple(
                tuple(
                    (value[r][c] - value[r][first_column]
                     - value[first_row][c] + value[first_row][first_column]) % modulus
                    for c in column_indices
                )
                for r in row_indices
            )
            initial_columns = tuple(range(width))
            initial_blocks = (tuple(range(width)),)
            pivot_fingerprint = _residual_fingerprint(
                interior, initial_columns, initial_blocks,
            )
            if pivot_fingerprint in seen_pivots:
                continue
            seen_pivots.add(pivot_fingerprint)

            def visit(
                columns: tuple[int, ...],
                blocks: tuple[tuple[int, ...], ...],
                prefix: tuple[int, ...],
                interior: ExponentMatrix = interior,
            ) -> None:
                nonlocal best
                least: tuple[int, ...] | None = None
                choices: list[tuple[int, tuple[tuple[int, ...], ...]]] = []
                seen_columns: set[tuple[int, ...]] = set()
                for column in columns:
                    # Identical columns have identical continuations after
                    # exchanging their labels. This matters for degenerate
                    # inputs, but never identifies distinct phase values.
                    column_values = tuple(row[column] for row in interior)
                    if column_values in seen_columns:
                        continue
                    seen_columns.add(column_values)
                    entries: list[int] = []
                    refined: list[tuple[int, ...]] = []
                    for block in blocks:
                        groups: dict[int, list[int]] = {}
                        for row in block:
                            groups.setdefault(interior[row][column], []).append(row)
                        for exponent in sorted(groups):
                            members = tuple(groups[exponent])
                            entries.extend((exponent,) * len(members))
                            refined.append(members)
                    signature = tuple(entries)
                    if least is None or signature < least:
                        least = signature
                        choices = [(column, tuple(refined))]
                    elif signature == least:
                        choices.append((column, tuple(refined)))
                assert least is not None
                candidate = prefix + least
                if best is not None and candidate > best[:len(candidate)]:
                    return
                if len(columns) == 1:
                    if best is None or candidate < best:
                        best = candidate
                    return
                seen_continuations: set[tuple[ExponentMatrix, ...]] = set()
                for column, refined in choices:
                    remaining = tuple(c for c in columns if c != column)
                    if len(choices) > 1:
                        fingerprint = _residual_fingerprint(interior, remaining, refined)
                        if fingerprint in seen_continuations:
                            continue
                        seen_continuations.add(fingerprint)
                    visit(remaining, refined, candidate)

            visit(initial_columns, initial_blocks, ())
    assert best is not None
    return n, modulus, best


def finite_butson_class_digest(key: FiniteButsonKey) -> str:
    """Hash a canonical key with the existing point-class digest convention."""
    return hashlib.sha256(repr(key).encode()).hexdigest()


def canonical_finite_butson_matrix(
    matrix: Iterable[Iterable[int]], modulus: int,
) -> tuple[int, ExponentMatrix]:
    """Return the reduced modulus and dephased canonical exponent matrix."""
    n, reduced_modulus, interior = canonical_finite_butson_key(matrix, modulus)
    width = n - 1
    return reduced_modulus, ((0,) * n,) + tuple(
        (0,) + tuple(interior[c * width + r] for c in range(width))
        for r in range(width)
    )


__all__ = [
    "FiniteButsonKey",
    "canonical_finite_butson_key",
    "canonical_finite_butson_matrix",
    "finite_butson_class_digest",
]

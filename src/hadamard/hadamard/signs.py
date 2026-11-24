"""Pure-Python combinatorics for KKPR-style 2-torsion covers.

This module deliberately has no SageMath imports.  It is used both by the
Sage frontend and by the ordinary-Python unit tests shipped with the bundle.

Convention
----------
A dephased sign pattern of order ``n`` is an ``n x n`` matrix over F_2 whose
row 0 and column 0 are zero.  The stored tuple contains the lower-right
``(n-1) x (n-1)`` block in row-major order.  A bit 1 means multiplication by
-1.  Permutations are written in *pullback* convention: ``p[new]`` is the old
index used at the new position.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass
from itertools import product

Permutation = tuple[int, ...]


def identity_permutation(n: int) -> Permutation:
    return tuple(range(n))


def adjacent_transposition(n: int, i: int) -> Permutation:
    """Return the pullback permutation swapping ``i`` and ``i+1``."""
    if not 0 <= i < n - 1:
        raise ValueError("adjacent-transposition index out of range")
    p = list(range(n))
    p[i], p[i + 1] = p[i + 1], p[i]
    return tuple(p)


def validate_permutation(p: Sequence[int], n: int) -> Permutation:
    p = tuple(int(x) for x in p)
    if len(p) != n or set(p) != set(range(n)):
        raise ValueError(f"not a permutation of range({n}): {p!r}")
    return p


@dataclass(frozen=True, slots=True)
class SignPattern:
    """A dephased ``n x n`` sign matrix over F_2."""

    n: int
    bits: tuple[int, ...]

    def __post_init__(self) -> None:
        if self.n < 2:
            raise ValueError("matrix order n must be at least 2")
        expected = (self.n - 1) ** 2
        if len(self.bits) != expected:
            raise ValueError(f"expected {expected} bits, got {len(self.bits)}")
        if any(bit not in (0, 1, False, True) for bit in self.bits):
            raise ValueError("sign bits must lie in {0,1}")
        object.__setattr__(self, "bits", tuple(int(bit) for bit in self.bits))

    @classmethod
    def zero(cls, n: int) -> SignPattern:
        return cls(n, (0,) * ((n - 1) ** 2))

    @classmethod
    def from_rows(cls, rows: Sequence[Sequence[int]]) -> SignPattern:
        """Build from a full dephased square bit matrix."""
        n = len(rows)
        if n < 2 or any(len(row) != n for row in rows):
            raise ValueError("rows must form an n x n matrix")
        if any(int(rows[0][c]) != 0 for c in range(n)):
            raise ValueError("first row must be zero")
        if any(int(rows[r][0]) != 0 for r in range(n)):
            raise ValueError("first column must be zero")
        bits = tuple(int(rows[r][c]) & 1 for r in range(1, n) for c in range(1, n))
        return cls(n, bits)

    def _index(self, r: int, c: int) -> int:
        return (r - 1) * (self.n - 1) + (c - 1)

    def bit(self, r: int, c: int) -> int:
        if not 0 <= r < self.n or not 0 <= c < self.n:
            raise IndexError((r, c))
        if r == 0 or c == 0:
            return 0
        return self.bits[self._index(r, c)]

    def full_rows(self) -> tuple[tuple[int, ...], ...]:
        return tuple(tuple(self.bit(r, c) for c in range(self.n)) for r in range(self.n))

    def variable_row_ints(self) -> tuple[int, ...]:
        """Encode non-first rows as integers, least-significant column first."""
        out = []
        for r in range(1, self.n):
            value = 0
            for c in range(1, self.n):
                value |= self.bit(r, c) << (c - 1)
            out.append(value)
        return tuple(out)

    def is_zero(self) -> bool:
        return not any(self.bits)

    def transpose(self) -> SignPattern:
        return SignPattern(
            self.n,
            tuple(self.bit(c, r) for r in range(1, self.n) for c in range(1, self.n)),
        )

    def row_difference(self, r: int, s: int) -> tuple[int, ...]:
        return tuple(self.bit(r, c) ^ self.bit(s, c) for c in range(self.n))

    def column_difference(self, c: int, d: int) -> tuple[int, ...]:
        return tuple(self.bit(r, c) ^ self.bit(r, d) for r in range(self.n))

    def row_blocks(self, r: int, s: int) -> tuple[tuple[int, ...], tuple[int, ...]]:
        """Return (plus-columns, minus-columns) for the row-pair cover split."""
        delta = self.row_difference(r, s)
        plus = tuple(c for c, bit in enumerate(delta) if bit == 0)
        minus = tuple(c for c, bit in enumerate(delta) if bit == 1)
        return plus, minus

    def column_blocks(self, c: int, d: int) -> tuple[tuple[int, ...], tuple[int, ...]]:
        """Return (plus-rows, minus-rows) for the column-pair cover split."""
        delta = self.column_difference(c, d)
        plus = tuple(r for r, bit in enumerate(delta) if bit == 0)
        minus = tuple(r for r, bit in enumerate(delta) if bit == 1)
        return plus, minus

    def row_distance_spectrum(self) -> tuple[int, ...]:
        return tuple(
            sum(self.row_difference(r, s))
            for r in range(self.n)
            for s in range(r + 1, self.n)
        )

    def column_distance_spectrum(self) -> tuple[int, ...]:
        return tuple(
            sum(self.column_difference(c, d))
            for c in range(self.n)
            for d in range(c + 1, self.n)
        )

    def singleton_obstructions(self) -> tuple[tuple[str, int, int, int], ...]:
        """List row/column pairs for which one cover block has size one.

        A root-of-unity sum with one term cannot vanish.  Hence every listed
        obstruction is a rigorous reason to discard the sign cover before any
        ideal is constructed.
        """
        bad: list[tuple[str, int, int, int]] = []
        for r in range(self.n):
            for s in range(r + 1, self.n):
                d = sum(self.row_difference(r, s))
                if d in (1, self.n - 1):
                    bad.append(("row", r, s, d))
        for c in range(self.n):
            for dcol in range(c + 1, self.n):
                d = sum(self.column_difference(c, dcol))
                if d in (1, self.n - 1):
                    bad.append(("column", c, dcol, d))
        return tuple(bad)

    def passes_singleton_test(self) -> bool:
        return not self.singleton_obstructions()

    def cheap_invariant(self) -> tuple[object, ...]:
        """A cheap invariant of dephased row/column permutations and transpose.

        Re-dephasing can complement a row- or column-difference word, so the
        invariant datum is ``min(d, n-d)`` rather than the raw Hamming distance.
        Transpose exchanges the row and column spectra, so the two spectra are
        stored as an unordered pair.  This key is suitable for bucketing only;
        equality still requires an exact orbit or ideal test.
        """
        row = tuple(sorted(min(d, self.n - d) for d in self.row_distance_spectrum()))
        column = tuple(sorted(min(d, self.n - d) for d in self.column_distance_spectrum()))
        return self.n, tuple(sorted((row, column)))

    def dephased_permute(
        self,
        row_permutation: Sequence[int] | None = None,
        column_permutation: Sequence[int] | None = None,
    ) -> SignPattern:
        """Apply arbitrary row/column permutations and re-dephase.

        If ``p`` and ``q`` are pullback permutations, the new bit is

        ``e[p(r),q(c)] + e[p(r),q(0)] + e[p(0),q(c)] + e[p(0),q(0)]`` mod 2.

        This is the F_2 shadow of the standard monomial re-dephasing action on
        a complex Hadamard matrix.
        """
        p = validate_permutation(
            row_permutation if row_permutation is not None else identity_permutation(self.n),
            self.n,
        )
        q = validate_permutation(
            column_permutation if column_permutation is not None else identity_permutation(self.n),
            self.n,
        )
        p0, q0 = p[0], q[0]
        bits = []
        for r in range(1, self.n):
            for c in range(1, self.n):
                value = (
                    self.bit(p[r], q[c])
                    ^ self.bit(p[r], q0)
                    ^ self.bit(p0, q[c])
                    ^ self.bit(p0, q0)
                )
                bits.append(value)
        return SignPattern(self.n, tuple(bits))

    def apply_generator(self, generator: SignSymmetryGenerator) -> SignPattern:
        pattern = self.transpose() if generator.transpose else self
        return pattern.dephased_permute(generator.rows, generator.columns)


@dataclass(frozen=True, slots=True)
class SignSymmetryGenerator:
    n: int
    rows: Permutation
    columns: Permutation
    transpose: bool = False
    name: str = ""

    def __post_init__(self) -> None:
        object.__setattr__(self, "rows", validate_permutation(self.rows, self.n))
        object.__setattr__(self, "columns", validate_permutation(self.columns, self.n))


def full_dephased_symmetry_generators(n: int, include_transpose: bool = True) -> tuple[SignSymmetryGenerator, ...]:
    """Generators for ``(S_n x S_n) semidirect C_2`` after re-dephasing."""
    ident = identity_permutation(n)
    gens: list[SignSymmetryGenerator] = []
    for i in range(n - 1):
        p = adjacent_transposition(n, i)
        gens.append(SignSymmetryGenerator(n, p, ident, False, f"row_s{i}"))
        gens.append(SignSymmetryGenerator(n, ident, p, False, f"col_s{i}"))
    if include_transpose:
        gens.append(SignSymmetryGenerator(n, ident, ident, True, "transpose"))
    return tuple(gens)


def exact_orbit(
    pattern: SignPattern,
    generators: Sequence[SignSymmetryGenerator] | None = None,
    *,
    max_nodes: int | None = None,
) -> tuple[SignPattern, ...] | None:
    """Return the exact generator orbit, or ``None`` if the safety cap is hit.

    Returning ``None`` rather than a partial orbit is important: callers must
    not use an incomplete orbit to discard branches.
    """
    gens = tuple(generators or full_dephased_symmetry_generators(pattern.n))
    seen = {pattern}
    queue = [pattern]
    cursor = 0
    while cursor < len(queue):
        current = queue[cursor]
        cursor += 1
        for generator in gens:
            nxt = current.apply_generator(generator)
            if nxt not in seen:
                seen.add(nxt)
                queue.append(nxt)
                if max_nodes is not None and len(seen) > max_nodes:
                    return None
    return tuple(sorted(seen, key=lambda p: p.bits))


def exact_canonical_key(
    pattern: SignPattern,
    generators: Sequence[SignSymmetryGenerator] | None = None,
    *,
    max_nodes: int | None = None,
) -> tuple[int, ...] | None:
    orbit = exact_orbit(pattern, generators, max_nodes=max_nodes)
    if orbit is None:
        return None
    return min(p.bits for p in orbit)


def _row_candidates(n: int) -> tuple[int, ...]:
    """All non-first sign rows that pass the row-vs-zero singleton test."""
    width = n - 1
    return tuple(
        word
        for word in range(1 << width)
        if word.bit_count() not in (1, n - 1)
    )


def _word_distance(a: int, b: int) -> int:
    return (a ^ b).bit_count()


def _bits_from_row_words(n: int, words: Sequence[int]) -> tuple[int, ...]:
    return tuple((words[r - 1] >> (c - 1)) & 1 for r in range(1, n) for c in range(1, n))


def _partial_column_feasible(n: int, rows: Sequence[int]) -> bool:
    """Sound partial pruning for future column-distance singleton tests."""
    width = n - 1
    used = len(rows)
    remaining = width - used
    # Column 0 is the all-zero column.  Other columns have a partial word.
    partial_columns = [0]
    for c in range(width):
        value = 0
        for r, word in enumerate(rows):
            value |= ((word >> c) & 1) << r
        partial_columns.append(value)
    for c in range(n):
        for d in range(c + 1, n):
            current = _word_distance(partial_columns[c], partial_columns[d])
            possible = range(current, current + remaining + 1)
            if all(value in (1, n - 1) for value in possible):
                return False
    return True


def iter_singleton_admissible_patterns(
    n: int,
    *,
    nondecreasing_rows: bool = True,
    exclude_zero: bool = True,
    limit: int | None = None,
) -> Iterator[SignPattern]:
    """Lazily enumerate sign patterns surviving all singleton tests.

    ``nondecreasing_rows=True`` quotients the evident permutation action on
    non-first rows by canonical augmentation.  This is exact for that subgroup,
    not a heuristic.  Full row/column/transposition orbit reduction is a
    separate optional step.
    """
    if n < 2:
        raise ValueError("n must be at least 2")
    candidates = _row_candidates(n)
    rows: list[int] = []
    yielded = 0

    def recurse(start_index: int) -> Iterator[SignPattern]:
        nonlocal yielded
        if limit is not None and yielded >= limit:
            return
        if len(rows) == n - 1:
            pattern = SignPattern(n, _bits_from_row_words(n, rows))
            if exclude_zero and pattern.is_zero():
                return
            if pattern.passes_singleton_test():
                yielded += 1
                yield pattern
            return

        iterable = range(start_index, len(candidates)) if nondecreasing_rows else range(len(candidates))
        for index in iterable:
            word = candidates[index]
            if any(_word_distance(word, previous) in (1, n - 1) for previous in rows):
                continue
            rows.append(word)
            if _partial_column_feasible(n, rows):
                next_start = index if nondecreasing_rows else 0
                yield from recurse(next_start)
            rows.pop()
            if limit is not None and yielded >= limit:
                return

    yield from recurse(0)


def iter_orbit_representatives(
    patterns: Iterable[SignPattern],
    *,
    generators: Sequence[SignSymmetryGenerator] | None = None,
    max_orbit_nodes: int | None = None,
) -> Iterator[SignPattern]:
    """Deduplicate by exact orbits whenever the orbit cap permits.

    If an orbit exceeds ``max_orbit_nodes``, the pattern is yielded and *not*
    entered into the seen set.  This can create duplicates but can never delete
    a legitimate branch.
    """
    seen: set[tuple[int, ...]] = set()
    for pattern in patterns:
        if pattern.bits in seen:
            continue
        orbit = exact_orbit(pattern, generators, max_nodes=max_orbit_nodes)
        if orbit is None:
            yield pattern
            continue
        for member in orbit:
            seen.add(member.bits)
        yield orbit[0]


def exhaustive_patterns(n: int) -> Iterator[SignPattern]:
    """Reference enumerator used only by tests and tiny orders."""
    for bits in product((0, 1), repeat=(n - 1) ** 2):
        yield SignPattern(n, tuple(bits))


def iter_row_sorted_patterns(
    n: int,
    *,
    exclude_zero: bool = False,
    limit: int | None = None,
) -> Iterator[SignPattern]:
    """Enumerate all sign patterns modulo permutations of non-first rows.

    This is intended for the sign-and-squaring/odd-rescaling cover, where the
    singleton test is not valid.  It is exact for the stated row subgroup but
    remains exponential and should therefore be used only with strict caps or
    additional canonical augmentation.
    """
    if n < 2:
        raise ValueError("n must be at least 2")
    width = n - 1
    candidates = tuple(range(1 << width))
    rows: list[int] = []
    yielded = 0

    def recurse(start_index: int) -> Iterator[SignPattern]:
        nonlocal yielded
        if limit is not None and yielded >= limit:
            return
        if len(rows) == width:
            pattern = SignPattern(n, _bits_from_row_words(n, rows))
            if exclude_zero and pattern.is_zero():
                return
            yielded += 1
            yield pattern
            return
        for index in range(start_index, len(candidates)):
            rows.append(candidates[index])
            yield from recurse(index)
            rows.pop()
            if limit is not None and yielded >= limit:
                return

    yield from recurse(0)

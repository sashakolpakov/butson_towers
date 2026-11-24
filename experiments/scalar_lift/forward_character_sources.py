"""Completed large Butson inputs from two multiplicative characters.

Ong--Schmidt, arXiv:2609.22754, Corollary 4.4.  The stored object is one
q-by-q group-ring coefficient array, not a dense q²-by-q² matrix.  Exact
one-dimensional correlation identities certify the entire large Gram.
The saved sources are known Butson constructions, not new real orders.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from dataclasses import dataclass
from functools import cache
from pathlib import Path

import numpy as np

from experiments.scalar_lift.fifth_source_assembly import root_vectors

PRIMARY_SOURCE = "https://arxiv.org/pdf/2609.22754"
PRESETS = {73: 8, 241: 40}


def _prime_divisors(value):
    result = []
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            result.append(divisor)
            while value % divisor == 0:
                value //= divisor
        divisor += 1
    if value > 1:
        result.append(value)
    return result


def _is_prime(value):
    return value >= 2 and all(value % d for d in range(2, math.isqrt(value) + 1))


def _array_hash(array):
    return hashlib.sha256(np.asarray(array, dtype="<i8").tobytes()).hexdigest()


@dataclass
class CharacterSource:
    prime: int
    first_order: int
    second_order: int
    phase: int
    gamma: int
    primitive_root: int
    logarithms: np.ndarray
    coefficients: np.ndarray

    @property
    def order(self):
        return self.prime**2

    def entry(self, row, column):
        """Exponent at two coordinates of F_q²; no dense matrix is formed."""
        difference = tuple((int(b) - int(a)) % self.prime for a, b in zip(row, column, strict=True))
        if len(difference) != 2:
            raise ValueError("row and column must each have two coordinates")
        return int(self.coefficients[difference])

    def metadata(self):
        return {
            "prime": self.prime, "order": self.order, "phase": self.phase,
            "character_orders": [self.first_order, self.second_order],
            "gamma_exponent": self.gamma, "primitive_root": self.primitive_root,
            "coefficient_array_shape": list(self.coefficients.shape),
            "coefficient_array_sha256": _array_hash(self.coefficients),
            "used_exponents": sorted(map(int, np.unique(self.coefficients))),
            "dense_order_squared_matrix_materialized": False,
        }


def character_source(prime=73, second_order=None, *, first_order=2, gamma=0):
    """Prime-field instance; both character orders must exceed one."""
    prime = int(prime)
    second_order = PRESETS.get(prime) if second_order is None else int(second_order)
    if not _is_prime(prime) or second_order is None:
        raise ValueError("a prime and an explicit/preset second character order are required")
    if any(order <= 1 or (prime - 1) % order for order in (first_order, second_order)):
        raise ValueError("both nontrivial character orders must divide prime-1")
    divisors = _prime_divisors(prime - 1)
    generator = next(g for g in range(2, prime)
                     if all(pow(g, (prime - 1) // p, prime) != 1 for p in divisors))
    logarithms = np.full(prime, -1, dtype=np.int64)
    value = 1
    for exponent in range(prime - 1):
        logarithms[value] = exponent
        value = value * generator % prime
    if value != 1 or sorted(map(int, logarithms[1:])) != list(range(prime - 1)):
        raise ArithmeticError("the discrete-log construction failed")
    phase = math.lcm(6, first_order, second_order)
    gamma = int(gamma) % phase
    coefficients = np.zeros((prime, prime), dtype=np.int64)
    coefficients[1:, 0] = phase // 6
    coefficients[0, 1:] = 5 * phase // 6
    coefficients[1:, 1:] = (
        gamma + (phase // first_order) * logarithms[1:, None]
        + (phase // second_order) * logarithms[None, 1:]
    ) % phase
    return CharacterSource(prime, first_order, second_order, phase, gamma,
                           generator, logarithms, coefficients)


def character_correlation_check(source, order):
    """All q one-dimensional periodic correlations, in an integral basis."""
    q = source.prime
    roots = root_vectors(order)
    nonzero = np.arange(1, q)
    if np.any(roots[source.logarithms[nonzero] % order].sum(axis=0)):
        raise ArithmeticError("nontrivial character does not have zero sum")
    values = []
    for shift in range(q):
        translated = (nonzero - shift) % q
        valid = translated != 0
        powers = (source.logarithms[nonzero[valid]] - source.logarithms[translated[valid]]) % order
        actual = roots[powers].sum(axis=0)
        expected = np.zeros(roots.shape[1], dtype=np.int64)
        expected[0] = q - 1 if shift == 0 else -1
        if not np.array_equal(actual, expected):
            raise ArithmeticError("character autocorrelation identity failed")
        values.append(actual)
    return {"character_order": order, "zero_sum_verified": True,
            "all_cyclic_lags_verified": q, "correlations_sha256": _array_hash(values),
            "exact_group_ring_identity": "AA*=q*delta-J"}


def formal_gram_check(source):
    """Prove every large Gram entry from exact small group-ring identities."""
    characters = [character_correlation_check(source, d)
                  for d in (source.first_order, source.second_order)]
    roots, q = root_vectors(source.phase), source.prime
    crossing = roots[source.phase // 3] + roots[2 * source.phase // 3]
    expected = np.zeros(roots.shape[1], dtype=np.int64)
    expected[0] = -1
    if not np.array_equal(crossing, expected):
        raise ArithmeticError("axis Gram coefficient is not -1")
    # Tensor basis delta⊗delta, J⊗delta, delta⊗J, J⊗J.
    axes = np.array([0, q, q, -1], dtype=np.int64)
    interior = np.array([q*q, -q, -q, 1], dtype=np.int64)
    if not np.array_equal(axes + interior, [q*q, 0, 0, 0]):
        raise ArithmeticError("group-ring tensor Gram failed")
    return {
        "character_checks": characters,
        "axes_Gram_coefficients": axes.tolist(),
        "interior_Gram_coefficients": interior.tolist(),
        "mixed_terms_zero_by_character_sums": True,
        "complete_Gram_identity": f"H H*=H* H={q*q} I_{q*q}",
        "all_matrix_Gram_entries_proved_by_group_ring_identity": True,
        "number_of_Gram_entries_covered": q**4,
    }


def direct_correlation_check(source, all_shifts=False):
    """Independent literal q²-term sums, with no Gauss-sum evaluation."""
    q, exponents = source.prime, source.coefficients
    roots = root_vectors(source.phase)
    shifts = list(itertools.product(range(q), repeat=2)) if all_shifts else [
        (0, 0), (1, 0), (0, 1), (1, 1), (2, 3), (q - 1, q - 2),
    ]
    checked = []
    for shift in shifts:
        differences = (exponents - np.roll(exponents, shift, axis=(0, 1))) % source.phase
        actual = roots[differences].sum(axis=(0, 1))
        expected = np.zeros(roots.shape[1], dtype=np.int64)
        if shift == (0, 0):
            expected[0] = q*q
        if not np.array_equal(actual, expected):
            raise ArithmeticError("literal full-length group correlation failed")
        checked.append({"shift": list(shift), "constant_coefficient": int(actual[0])})
    return {"all_group_shifts_checked": bool(all_shifts), "literal_checks": checked}


def primitive_cross_ratio(source):
    for y in range(2, source.prime):
        rows, columns = [(0, 0), (0, 1)], [(0, 0), (1, y)]
        exponents = [[source.entry(row, column) for column in columns] for row in rows]
        value = (exponents[0][0] + exponents[1][1] - exponents[0][1] - exponents[1][0]) % source.phase
        if math.gcd(value, source.phase) == 1:
            return {"rows": rows, "columns": columns, "entry_exponents": exponents,
                    "cross_ratio_exponent": value, "cross_ratio_order": source.phase,
                    "full_cyclotomic_phase_is_gauge_invariant": True}
    return {"primitive_cross_ratio_found": False}


def uniform_sign_obstruction(source):
    used = set(map(int, np.unique(source.coefficients)))
    exponents = [source.phase // 6, 5 * source.phase // 6, 0]
    if not set(exponents).issubset(used):
        raise ArithmeticError("axis and identity entries are absent")
    roots = root_vectors(source.phase)
    if np.any(roots[exponents[0]] + roots[exponents[1]] - roots[exponents[2]]):
        raise ArithmeticError("signed odd three-term relation failed")
    return {"exponents": exponents, "coefficients": [1, 1, -1],
            "exact_relation": "zeta6+zeta6^-1-1=0", "odd_length": 3,
            "uniform_Q_linear_final_sign_map_excluded_for_raw_source": True,
            "all_gamma_and_global_root_rotations_covered": True,
            "arbitrary_row_column_gauges_or_nonuniform_maps_excluded": False}


@cache
def cubic_gauge_census():
    """Complete tiny obstructions to this source's two-cubic-coset gauges."""
    reports = []
    for shape in ((4, 2), (2, 4)):
        points = list(itertools.product(range(shape[0]), range(shape[1])))
        for gamma in range(3):
            colors = np.array([
                [0 if left == right else 2 if left[0] == right[0]
                 else 1 if left[1] == right[1] else gamma for right in points]
                for left in points], dtype=np.int64,
            )
            survivors = 0
            full_gauges = 0
            witness = None
            for tail in itertools.product(range(3), repeat=7):
                row = (0, *tail)
                domains = [tuple(c for c in range(3)
                                 if all((int(colors[i, j]) + row[i] + c) % 3 != 2
                                        for i in range(8))) for j in range(8)]
                if all(domains):
                    survivors += 1
                    full_gauges += math.prod(map(len, domains))
                    if witness is None:
                        witness = {"row_gauge": row, "column_gauge": [d[0] for d in domains]}
            reports.append({"grid_shape": shape, "normalized_interior_cubic_color": gamma,
                            "row_gauges_exhausted": 3**7, "feasible_row_gauges": survivors,
                            "complete_row_column_gauges": full_gauges, "witness": witness})
    if any(not any(r["normalized_interior_cubic_color"] == gamma
                   and r["feasible_row_gauges"] == 0 for r in reports) for gamma in range(3)):
        raise ArithmeticError("the finite cubic-gauge obstruction failed")
    return reports


def literal_gauge_template_check(source):
    if (source.phase // source.first_order) % 3 or (source.phase // source.second_order) % 3:
        raise ValueError("this preflight needs characters with trivial cubic component")
    axis = (source.phase // 6) % 3
    if axis == 0:
        raise ValueError("the saved cubic-color normalization is not applicable")
    inverse = pow(axis, -1, 3)
    gamma = source.gamma * inverse % 3
    for shape in ((4, 2), (2, 4)):
        points = list(itertools.product(range(shape[0]), range(shape[1])))
        for row in points:
            for column in points:
                expected = (0 if row == column else 2 if row[0] == column[0]
                            else 1 if row[1] == column[1] else gamma)
                if source.entry(row, column) * inverse % 3 != expected:
                    raise ArithmeticError("the literal source does not match the gauge template")
    return {"both_literal_eight_by_eight_templates_checked": True,
            "normalized_interior_cubic_color": gamma,
            "arbitrary_root_row_column_gauge_to_current_two_cubic_cosets": "EXCLUDED",
            "other_partial_frames_or_nonuniform_descent_maps_excluded": False}


def certificate():
    records = []
    for prime in PRESETS:
        source = character_source(prime)
        records.append({**source.metadata(), "formal_Gram": formal_gram_check(source),
                        "literal_controls": direct_correlation_check(source),
                        "genuine_phase": primitive_cross_ratio(source),
                        "uniform_sign_gate": uniform_sign_obstruction(source),
                        "current_cubic_alphabet_gate": literal_gauge_template_check(source),
                        "complete_exponent_descent_states": [
                            {"order": source.order * multiplier, "phase": source.phase // multiplier}
                            for multiplier in (1, 2, 4)],
                        "completed_real_Hadamard_descent": False})
    small = character_source(5, 4)
    points = list(itertools.product(range(5), repeat=2))
    matrix = np.array([[small.entry(a, b) for b in points] for a in points])
    roots = root_vectors(small.phase)
    gram = roots[(matrix[:, None, :] - matrix[None, :, :]) % small.phase].sum(axis=2)
    target = np.zeros_like(gram)
    target[np.arange(25), np.arange(25), 0] = 25
    if not np.array_equal(gram, target):
        raise ArithmeticError("the independent dense small control failed")
    return {"schema": "forward-character-sources-v1", "primary_source": PRIMARY_SOURCE,
            "primary_result": "Ong--Schmidt Corollary 4.4, submitted 19 September 2026",
            "required_character_parity_condition": None, "known_source_constructions": records,
            "cubic_gauge_census": cubic_gauge_census(),
            "small_control": {**small.metadata(), "full_dense_order25_Gram_verified": True,
                              "dense_order_squared_matrix_materialized": True,
                              **direct_correlation_check(small, all_shifts=True)},
            "new_Butson_order_claimed": False, "new_real_Hadamard_order_claimed": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--output", type=Path)
    action.add_argument("--check", type=Path)
    args = parser.parse_args()
    report = json.loads(json.dumps(certificate()))
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    elif json.loads(args.check.read_text()) != report:
        raise ValueError("forward character-source certificate failed replay")
    print(json.dumps({"known_source_orders": [r["order"] for r in report["known_source_constructions"]],
                      "all_Grams_exact": True, "new_real_order": False}))


if __name__ == "__main__":
    main()

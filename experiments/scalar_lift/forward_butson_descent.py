"""Forward restriction of genuine Butson sources through repeated-prime fields.

BH(n,p*q) -> BH(p*n,q) when p is prime and p divides q. Pure dyadic
towers terminate in actual sign matrices. Odd factors need a separately
proved terminal lift; this module does not assume one exists.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from functools import cache
from pathlib import Path

import numpy as np
from sympy import Poly, Symbol, cyclotomic_poly, factorint, isprime, totient

from experiments.scalar_lift.fifth_source_assembly import exact_column_gram, root_vectors


def edge_parameters(upper_phase, prime, generator_mode="primary"):
    upper, p = int(upper_phase), int(prime)
    if upper != upper_phase or upper < 4 or p != prime or not isprime(p) or upper % (p * p):
        raise ValueError("the prime must occur at least twice in the upper root order")
    if generator_mode not in ("primary", "primitive"):
        raise ValueError("generator mode must be primary or primitive")
    lower = upper // p
    coprime = upper
    while coprime % p == 0:
        coprime //= p
    generator = coprime if generator_mode == "primary" else 1
    if int(totient(upper)) != p * int(totient(lower)):
        raise ValueError("the proposed cyclotomic field extension has the wrong degree")
    return {"upper_phase": upper, "lower_phase": lower, "prime": p,
            "relative_field_degree": p, "order_multiplier": p,
            "generator_mode": generator_mode, "generator_upper_exponent": generator,
            "generator_order": upper // generator, "xi_lower_exponent": generator,
            "coprime_part": coprime,
            "lower_field_conductor": lower // 2 if lower % 4 == 2 else lower,
            "lower_root_alphabet_retained": lower}


def available_steps(phase):
    return [int(p) for p, exponent in factorint(int(phase)).items() if exponent >= 2]


def block_entry(exponent, row, column, upper_phase, prime, generator_mode="primary"):
    edge = edge_parameters(upper_phase, prime, generator_mode)
    p, q = edge["prime"], edge["lower_phase"]
    if not 0 <= row < p or not 0 <= column < p:
        raise ValueError("scalar block indices lie outside the relative field degree")
    e = int(exponent) % upper_phase
    generator = edge["generator_upper_exponent"]
    remainder = e * pow(generator, -1, p) % p
    quotient = (e - generator * remainder) // p
    return (quotient + generator * int(row < remainder) + (q // p) * ((row - remainder) % p) * column) % q


@cache
def scalar_blocks(upper_phase, prime, generator_mode="primary"):
    edge = edge_parameters(upper_phase, prime, generator_mode)
    p = edge["prime"]
    return np.asarray([[[block_entry(e, r, s, upper_phase, p, generator_mode) for s in range(p)]
                        for r in range(p)] for e in range(upper_phase)], dtype=np.int64)


def _monomial_power(exponent, upper_phase, prime, generator_mode="primary"):
    q = upper_phase // prime
    generator = edge_parameters(upper_phase, prime, generator_mode)["generator_upper_exponent"]
    e = int(exponent) % upper_phase
    remainder = e * pow(generator, -1, prime) % prime
    quotient = (e - generator * remainder) // prime
    return [(r, (r - remainder) % prime, (quotient + generator * int(r < remainder)) % q) for r in range(prime)]


@cache
def scalar_certificate(upper_phase, prime, generator_mode="primary"):
    edge = edge_parameters(upper_phase, prime, generator_mode)
    p, q = edge["prime"], edge["lower_phase"]
    roots = root_vectors(q)
    blocks = scalar_blocks(upper_phase, p, generator_mode)
    x = Symbol("x")
    polynomial_checks = []
    for order, multiplier in ((upper_phase, 1), (edge["generator_order"], edge["generator_upper_exponent"])):
        coefficients = [int(c) for c in reversed(Poly(cyclotomic_poly(order, x), x).all_coeffs())]
        polynomial = np.zeros((p, p, roots.shape[1]), dtype=np.int64)
        for e, coefficient in enumerate(coefficients):
            for r, s, value in _monomial_power(multiplier * e, upper_phase, p, generator_mode):
                polynomial[r, s] += coefficient * roots[value]
        if np.any(polynomial):
            raise ValueError("a cyclic action does not annihilate its cyclotomic polynomial")
        polynomial_checks.append({"order": order, "scalar_upper_exponent": multiplier,
                                  "polynomial_coefficients_ascending": coefficients, "exact_zero": True})
    for e in range(upper_phase):
        for f in range(upper_phase):
            actual = np.array([[roots[(blocks[e, r] - blocks[f, s]) % q].sum(axis=0)
                                for s in range(p)] for r in range(p)])
            expected = np.zeros_like(actual)
            for r, s, value in _monomial_power(e - f, upper_phase, p, generator_mode):
                expected[r, s] = p * roots[value]
            if not np.array_equal(actual, expected):
                raise ValueError("a complete scalar cross-Gram identity failed")
    return {**edge, "chosen_generator_minimal_polynomial_over_lower_field": f"X^p-zeta_q^{edge['xi_lower_exponent']}",
            "cyclotomic_polynomial_checks": polynomial_checks,
            "Phi_upper_of_represented_primitive_root_equals_zero_exactly": True,
            "all_scalar_cross_Grams_checked": upper_phase**2,
            "cross_Gram_identity": "B_e B_f* = p rho(zeta_upper^(e-f))",
            "all_blocks_have_only_lower_roots": True,
            "scalar_blocks": blocks.tolist(), "terminal_sign_blocks": q == 2,
            "arbitrary_partial_alphabet_preservation_claimed": False}


def exponent_hash(source):
    return hashlib.sha256(np.asarray(source, dtype="<i8").tobytes()).hexdigest()


def verify_source(source, phase):
    source = np.asarray(source)
    if source.dtype.kind not in "iu":
        raise ValueError("source exponents must be integers")
    source = source.astype(np.int64) % phase
    if source.ndim != 2 or source.shape[0] != source.shape[1] or not len(source):
        raise ValueError("a nonempty square Butson source is required")
    gram = exact_column_gram(source, phase)
    expected = np.zeros_like(gram)
    expected[np.arange(len(source)), np.arange(len(source)), 0] = len(source)
    if not np.array_equal(gram, expected):
        raise ValueError("the input is not a Butson matrix")
    return {"order": len(source), "phase": int(phase), "exponent_sha256": exponent_hash(source),
            "full_Gram_verified": True,
            "verification": "complete exact column Gram; square invertibility gives the row Gram"}


class ForwardTower:
    """A literal nested matrix, with entries available without dense expansion.

    The constructor validates dimensions/alphabet. A Butson claim additionally
    requires the source Gram proof passed to ``proof``; it is never inferred
    merely from a root exponent array.
    """

    def __init__(self, source, phase, primes=None, generator_mode="primary"):
        source = np.asarray(source)
        if source.dtype.kind not in "iu" or source.ndim != 2 or source.shape[0] != source.shape[1] or not len(source):
            raise ValueError("a nonempty square integer-exponent source is required")
        if int(phase) != phase or phase < 2:
            raise ValueError("the root alphabet order must be an integer at least two")
        self.source = source.astype(np.int64) % int(phase)
        self.phase = int(phase)
        self.edges = []
        self.orders = [len(source)]
        self.phases = [self.phase]
        if primes is None:
            primes = []
            current = self.phase
            while available_steps(current):
                p = min(available_steps(current))
                primes.append(p)
                current //= p
        for p in primes:
            edge = edge_parameters(self.phases[-1], p, generator_mode)
            self.edges.append(edge)
            self.orders.append(self.orders[-1] * edge["prime"])
            self.phases.append(edge["lower_phase"])

    def entry(self, row, column, stage=None):
        stage = len(self.edges) if stage is None else int(stage)
        if not 0 <= stage <= len(self.edges) or not 0 <= row < self.orders[stage] or not 0 <= column < self.orders[stage]:
            raise ValueError("matrix entry lies outside the chosen tower stage")
        digits = []
        for level in range(stage - 1, -1, -1):
            p = self.edges[level]["prime"]
            row, r = divmod(row, p)
            column, s = divmod(column, p)
            digits.append((level, r, s))
        exponent = int(self.source[row, column])
        for level, r, s in reversed(digits):
            edge = self.edges[level]
            exponent = block_entry(exponent, r, s, edge["upper_phase"], edge["prime"], edge["generator_mode"])
        return exponent

    def sign_entry(self, row, column):
        if self.phases[-1] != 2:
            raise ValueError("the tower has not reached the real sign alphabet")
        return 1 - 2 * self.entry(row, column)

    def materialize(self, stage=None, max_entries=1_000_000):
        stage = len(self.edges) if stage is None else int(stage)
        if not 0 <= stage <= len(self.edges) or self.orders[stage]**2 > max_entries:
            raise ValueError("the requested dense stage exceeds the entry budget")
        result = self.source.copy()
        for edge in self.edges[:stage]:
            p, old = edge["prime"], len(result)
            blocks = scalar_blocks(edge["upper_phase"], p, edge["generator_mode"])
            result = blocks[result].transpose(0, 2, 1, 3).reshape(p * old, p * old)
        return result

    def proof(self, source_proof):
        if (source_proof.get("order") != len(self.source) or source_proof.get("phase") != self.phase
                or source_proof.get("exponent_sha256") != exponent_hash(self.source)
                or source_proof.get("full_Gram_verified") is not True):
            raise ValueError("the supplied source proof does not match this exact source")
        return {"source_proof": source_proof,
                "stages": [{"order": n, "phase": q} for n, q in zip(self.orders, self.phases, strict=True)],
                "edges": [scalar_certificate(edge["upper_phase"], edge["prime"], edge["generator_mode"]) for edge in self.edges],
                "all_stage_Grams_follow_from_exact_scalar_cross_Grams": True,
                "terminal_is_real_Hadamard": self.phases[-1] == 2,
                "terminal_order": self.orders[-1], "terminal_phase": self.phases[-1],
                "entry_accessor_is_literal_matrix_not_a_missing_completion": True,
                "new_Hadamard_order_claimed": False}


def fourier_control(phase):
    source = np.outer(np.arange(phase), np.arange(phase)) % phase
    tower = ForwardTower(source, phase)
    report = tower.proof(verify_source(source, phase))
    checked = []
    for stage in range(1, len(tower.edges) + 1):
        matrix = tower.materialize(stage)
        checked.append(verify_source(matrix, tower.phases[stage]))
        for row, column in ((0, 0), (1, 2), (len(matrix) - 1, len(matrix) - 2)):
            if tower.entry(row, column, stage) != int(matrix[row, column]):
                raise ValueError("the lazy nested entry disagrees with dense restriction")
    report["independent_complete_stage_Grams"] = checked
    if tower.phases[-1] == 2:
        signs = 1 - 2 * tower.materialize()
        if not np.array_equal(signs @ signs.T, len(signs) * np.eye(len(signs), dtype=np.int64)):
            raise ValueError("the terminal literal sign Gram failed")
        report["terminal_full_sign_Gram_verified"] = True
    return report


def high_conference_control():
    from experiments.scalar_lift.even_petrescu_conference_family import (
        formal_family_certificate,
        source_exponents,
    )
    from experiments.scalar_lift.even_petrescu_lattice import lattice_data

    # lattice_data replays the exact order-500 real/skew identities.
    _, _, skew_plus_identity, _, _ = lattice_data()
    identity = np.eye(500, dtype=np.int64)
    conference = skew_plus_identity - identity
    if not np.array_equal(conference.T, -conference):
        raise ValueError("the actual conference is not skew")
    formal = formal_family_certificate(8)
    source = source_exponents(1000, 8)
    prior_path = Path(__file__).resolve().parents[2] / "certificates/scalar-lift/even-petrescu-conference-family.json"
    prior = json.loads(prior_path.read_text())
    prior_source = next(row for row in prior["complete_known_sources"] if row["order"] == 1000 and row["phase"] == 8)
    if hashlib.sha256(source.astype(np.uint8).tobytes()).hexdigest() != prior_source["normalized_exponent_sha256"]:
        raise ValueError("the actual large source differs from its existing certificate")
    if source[1, 2] != 1:
        raise ValueError("the genuine primitive eighth-root cross ratio disappeared")
    proof = {"order": 1000, "phase": 8, "exponent_sha256": exponent_hash(source),
             "full_Gram_verified": True,
             "verification": "exact actual skew conference Gram and all four symbolic order-1000 Gram blocks",
             "formal_Gram": formal, "prior_source_certificate": str(prior_path.relative_to(prior_path.parents[2])),
             "prior_source_certificate_sha256": hashlib.sha256(prior_path.read_bytes()).hexdigest(),
             "primitive_cross_ratio_exponent": 1, "primitive_cross_ratio_order": 8}
    tower = ForwardTower(source, 8)
    report = tower.proof(proof)
    report["terminal_dense_4000_matrix_materialized"] = False
    report["literal_terminal_sign_samples"] = [
        {"row": i, "column": j, "value": tower.sign_entry(i, j)}
        for i, j in ((0, 0), (1, 2), (1999, 2000), (3999, 3998))]
    return report


def certificate():
    return {"schema": "forward-butson-descent-v1", "direction": "verified source to explicit descendants",
            "general_edge": "BH(n,p*q) -> BH(p*n,q), p prime and p divides q",
            "mixed_phase_generator_controls": [scalar_certificate(12, 2, mode) for mode in ("primary", "primitive")],
            "fourier_eight": fourier_control(8), "fourier_twenty_seven": fourier_control(27),
            "actual_order1000_eighth_root_source": high_conference_control(),
            "odd_full_alphabet_obstruction": "an odd prime polygon cannot map to an entrywise odd sum of signs equal to zero",
            "odd_alphabets_need_a_separately_verified_terminal_lift": True,
            "new_Hadamard_order_claimed": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--output", type=Path)
    action.add_argument("--check", type=Path)
    args = parser.parse_args()
    report = certificate()
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    elif json.loads(args.check.read_text()) != report:
        raise ValueError("the forward descent certificate failed replay")
    print(json.dumps({"real_orders": [32, 4000], "odd_terminal": {"order": 243, "phase": 3},
                      "new_Hadamard_order_claimed": False}))


if __name__ == "__main__":
    main()

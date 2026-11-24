"""Exact relative cyclotomic restriction, retaining its Hermitian metric.

No frame or flat-entry condition is imposed. The lower field and its
declared root label are kept separate, including the rational labels 1,2.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass
from functools import cache

import numpy as np
from sympy import ImmutableMatrix, Matrix, Poly, Rational, Symbol, cyclotomic_poly, totient

from experiments.scalar_lift.vanishing_collision_search import integral_roots


@dataclass(frozen=True)
class RationalArray:
    numerators: np.ndarray
    denominator: int = 1

    @classmethod
    def exact(cls, values, shape=None):
        values = np.asarray(values, dtype=object)
        shape = values.shape if shape is None else shape
        rationals = [Rational(value) for value in values.flat]
        denominator = math.lcm(*(int(value.q) for value in rationals)) if rationals else 1
        numbers = [int(value * denominator) for value in rationals]
        divisor = math.gcd(denominator, *numbers)
        denominator //= divisor
        numbers = [value // divisor for value in numbers]
        if any(abs(value) >= 2**63 for value in numbers):
            raise OverflowError("exact restriction coefficients exceed the int64 storage range")
        array = np.asarray(numbers, dtype=np.int64).reshape(shape)
        array.setflags(write=False)
        return cls(array, denominator)

    def __getitem__(self, index):
        return RationalArray(self.numerators[index], self.denominator)

    def to_json(self):
        return {"numerators": self.numerators.tolist(), "denominator": self.denominator}

    def to_sympy(self):
        values = np.asarray([Rational(int(value), self.denominator) for value in self.numerators.flat], dtype=object)
        return values.reshape(self.numerators.shape)


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


@cache
def _lower_tables(phase):
    roots = np.asarray(integral_roots(phase), dtype=np.int64)
    h = roots.shape[1]
    multiplication = np.asarray([np.array([roots[(a+b) % phase] for b in range(h)]).T for a in range(h)])
    conjugation = np.array([roots[-a % phase] for a in range(h)]).T
    return roots, multiplication, conjugation


def regular_numerators(value, lower_phase):
    """Expand field coefficients without introducing approximate arithmetic."""
    _, multiplication, _ = _lower_tables(lower_phase)
    numbers = value.numerators
    if numbers.ndim != 3 or numbers.shape[2] != len(multiplication):
        raise ValueError("expected [rows,columns,phi(lower)] field coefficients")
    rows, columns, h = numbers.shape
    bound = h * int(np.max(np.abs(numbers))) * int(np.max(np.abs(multiplication)))
    if bound >= 2**62:
        raise OverflowError("regular coefficient multiplication exceeds the exact int64 range")
    return np.einsum("rcu,uab->racb", numbers, multiplication).reshape(rows*h, columns*h)


def regular_matrix(value, lower_phase):
    """Expand a K-valued matrix in the rational power basis of K."""
    expanded = regular_numerators(value, lower_phase)
    return ImmutableMatrix(expanded.tolist()) / value.denominator


def conjugate_transpose(value, lower_phase):
    _, _, conjugation = _lower_tables(lower_phase)
    if value.numerators.ndim != 3:
        raise ValueError("expected a field-valued matrix")
    bound = len(conjugation) * int(np.max(np.abs(value.numerators))) * int(np.max(np.abs(conjugation)))
    if bound >= 2**62:
        raise OverflowError("field conjugation exceeds the exact int64 range")
    return RationalArray(np.tensordot(value.numerators.transpose(1, 0, 2), conjugation.T, axes=(2, 0)),
                         value.denominator)


class RelativeRestriction:
    def __init__(self, upper_phase, lower_phase, basis_exponents=None):
        Q, q = int(upper_phase), int(lower_phase)
        if upper_phase != Q or lower_phase != q or Q < 1 or q < 1 or Q % q:
            raise ValueError("positive root labels with lower dividing upper are required")
        self.upper_phase, self.lower_phase = Q, q
        self.upper_degree, self.lower_degree = int(totient(Q)), int(totient(q))
        if self.upper_degree % self.lower_degree:
            raise ValueError("the relative field degree is not integral")
        self.width = self.upper_degree // self.lower_degree
        d, h, D = self.width, self.lower_degree, self.upper_degree
        if basis_exponents is not None and any(int(a) != a for a in basis_exponents):
            raise ValueError("relative basis exponents must be exact integers")
        self.basis_exponents = tuple(range(d)) if basis_exponents is None else tuple(int(a) % Q for a in basis_exponents)
        if len(self.basis_exponents) != d or self.basis_exponents[0] != 0:
            raise ValueError("a relative basis of the required degree starting with 1 is required")
        self.combined_exponents = tuple(a + (Q // q) * b for a in self.basis_exponents for b in range(h))
        roots = np.asarray(integral_roots(Q), dtype=np.int64)
        self.basis_matrix = ImmutableMatrix(np.array([roots[e % Q] for e in self.combined_exponents]).T.tolist())
        if self.basis_matrix.det() == 0:
            raise ValueError("the proposed relative power basis is singular")
        self.inverse_basis_matrix = ImmutableMatrix(self.basis_matrix.inv())
        self.basis_matrix_data = RationalArray.exact(self.basis_matrix.tolist())
        self.inverse_basis_matrix_data = RationalArray.exact(self.inverse_basis_matrix.tolist())
        all_coordinates = self.inverse_basis_matrix * Matrix(roots.T.tolist())
        self._coordinates = RationalArray.exact(np.asarray(all_coordinates.T.tolist(), dtype=object).reshape(Q, d, h))
        actions = np.empty((Q, d, d, h), dtype=object)
        coordinates = self._coordinates.to_sympy()
        for e in range(Q):
            for j in range(d):
                actions[e, :, j, :] = coordinates[(e+self.basis_exponents[j]) % Q]
        self.actions = RationalArray.exact(actions)
        self.trace_automorphisms = tuple(k for k in range(Q) if math.gcd(k, Q) == 1 and (k-1) % q == 0)
        if len(self.trace_automorphisms) != d:
            raise ValueError("the relative trace stabilizer has the wrong order")
        traces = []
        for e in range(Q):
            value = roots[[(k*e) % Q for k in self.trace_automorphisms]].sum(axis=0)
            relative = self.inverse_basis_matrix * Matrix(value.tolist())
            if any(relative[i] for i in range(h, D)):
                raise ValueError("the relative trace did not lie in the declared lower field")
            traces.append(list(relative[:h]))
        self._traces = RationalArray.exact(traces)
        trace_values = self._traces.to_sympy()
        self.metric = RationalArray.exact([[trace_values[(b-a) % Q] for b in self.basis_exponents]
                                          for a in self.basis_exponents])
        self.gram = self.metric
        metric_regular = regular_matrix(self.metric, q)
        inverse_regular = metric_regular.inv()
        inverse = [[[inverse_regular[a*h+b, c*h] for b in range(h)] for c in range(d)] for a in range(d)]
        self.inverse_metric = RationalArray.exact(inverse)
        self.inverse_gram = self.inverse_metric
        if regular_matrix(self.inverse_metric, q) != inverse_regular:
            raise ValueError("the rational inverse is not the advertised lower-field metric inverse")
        if regular_matrix(conjugate_transpose(self.metric, q), q) != metric_regular:
            raise ValueError("the relative trace metric is not Hermitian")
        # All scalar powers are literal coordinate transforms of multiplication
        # in the original rational cyclotomic field. Check the complete root
        # table against those transforms, independent of lower-field arithmetic.
        upper_multiplications = roots[(np.arange(Q)[:, None] + np.arange(D)[None, :]) % Q].transpose(0, 2, 1)
        inverse_numbers = self.inverse_basis_matrix_data.numerators
        basis_numbers = self.basis_matrix_data.numerators
        bound = (D*D * int(np.max(np.abs(inverse_numbers))) * int(np.max(np.abs(roots)))
                 * int(np.max(np.abs(basis_numbers))) * self.actions.denominator)
        if bound >= 2**62:
            raise OverflowError("basis conjugacy exceeds the exact int64 range")
        expected = inverse_numbers @ upper_multiplications @ basis_numbers
        _, lower_multiplications, _ = _lower_tables(q)
        bound = (h * int(np.max(np.abs(self.actions.numerators))) * int(np.max(np.abs(lower_multiplications)))
                 * self.inverse_basis_matrix_data.denominator)
        if bound >= 2**62:
            raise OverflowError("all-root regular coefficients exceed the exact int64 range")
        actual = np.einsum("ercu,uab->eracb", self.actions.numerators, lower_multiplications).reshape(Q, D, D)
        if not np.array_equal(actual * self.inverse_basis_matrix_data.denominator,
                              expected * self.actions.denominator):
            raise ValueError("a relative scalar action disagrees with exact multiplication")
        generator = regular_matrix(self.actions[1 % Q], q)
        conjugate_generator = regular_matrix(self.actions[-1 % Q], q)
        generator_adjoint = regular_matrix(conjugate_transpose(self.actions[1 % Q], q), q)
        if generator_adjoint * metric_regular != metric_regular * conjugate_generator:
            raise ValueError("the scalar action fails its exact star-metric identity")
        x = Symbol("x")
        polynomial = list(reversed([int(c) for c in Poly(cyclotomic_poly(Q, x), x).all_coeffs()]))
        if sum(abs(c) for c in polynomial) * int(np.max(np.abs(actual))) >= 2**62:
            raise OverflowError("cyclotomic evaluation exceeds the exact int64 range")
        total = sum((coefficient * actual[e % Q] for e, coefficient in enumerate(polynomial)), np.zeros((D, D), dtype=np.int64))
        if np.any(total) or generator * regular_matrix(self.actions[-1 % Q], q) != Matrix.eye(D):
            raise ValueError("the relative action fails cyclotomic or inverse identities")
        if generator * inverse_regular * generator_adjoint != inverse_regular:
            raise ValueError("the weighted inverse metric is not preserved")
        self._proof = {
            "all_upper_scalar_powers_verified_by_exact_basis_conjugacy": Q,
            "upper_cyclotomic_polynomial_coefficients_ascending": polynomial,
            "Phi_upper_of_generator_exact_zero": True,
            "metric_equals_unnormalized_relative_trace": True,
            "metric_Hermitian_and_inverse_verified": True,
            "star_metric_identity": "rho(z)^* G = G rho(conjugate(z))",
            "generator_star_metric_verified_all_powers_follow_by_multiplication": True,
            "weighted_unitarity_identity": "rho(zeta) G^-1 rho(zeta)^* = G^-1",
            "source_Gram_implication": "rho(H)(I_n tensor G^-1)rho(H)^* = n(I_n tensor G^-1)",
            "metric_positive_definite": "relative trace is the sum of squared complex embedding magnitudes",
            "ordinary_identity_metric_or_flat_entries_claimed": False,
        }

    def coordinates(self, exponent):
        return self._coordinates[int(exponent) % self.upper_phase]

    def trace_exponent(self, exponent):
        return self._traces[int(exponent) % self.upper_phase]

    def action(self, exponent):
        return self.actions[int(exponent) % self.upper_phase]

    def entry(self, exponent, row, column):
        if not 0 <= row < self.width or not 0 <= column < self.width:
            raise ValueError("relative scalar entry outside the block")
        return self.action(exponent)[row, column].to_json()

    def proof(self):
        return dict(self._proof)

    def metadata(self):
        data = {"schema": "relative-scalar-restriction-v1", "upper_phase": self.upper_phase,
                "lower_phase": self.lower_phase, "upper_degree": self.upper_degree,
                "lower_degree": self.lower_degree, "relative_degree": self.width,
                "lower_field_conductor": self.lower_phase // 2 if self.lower_phase % 4 == 2 else self.lower_phase,
                "declared_lower_root_label_retained": self.lower_phase,
                "relative_basis_exponents": list(self.basis_exponents),
                "combined_rational_basis_exponents": list(self.combined_exponents),
                "basis_matrix": self.basis_matrix_data.to_json(),
                "inverse_basis_matrix": self.inverse_basis_matrix_data.to_json(),
                "actions": self.actions.to_json(), "metric": self.metric.to_json(),
                "inverse_metric": self.inverse_metric.to_json(),
                "trace_automorphisms": list(self.trace_automorphisms), "proof": self.proof()}
        data["exact_definition_sha256"] = _digest(data)
        return data


@cache
def edge(upper_phase, lower_phase, basis_exponents=None):
    return RelativeRestriction(upper_phase, lower_phase, basis_exponents)


def restriction_in_basis(upper_phase, lower_phase, basis_exponents):
    return edge(upper_phase, lower_phase, tuple(basis_exponents))

"""Necessary arithmetic filters for target-first cyclotomic scalar towers.

Passing these tests proves neither a star representation in the standard
Euclidean metric nor a flat frame nor existence of the source Butson matrix.
The accompanying note proves the determinant and local-at-two filters.
"""

from __future__ import annotations

from math import gcd, lcm

from sympy import Matrix, divisors, eye, factorint, isprime, ones, primitive_root, totient

from hadamard.butson import lam_leung_obstructs_butson


def cyclotomic_orders_with_degree_dividing(d):
    """Invert phi(k)|d exactly, without an arbitrary conductor cutoff."""
    if type(d) is not int or d < 1:
        raise ValueError("the representation dimension must be a positive integer")
    states = {(1, 1)}
    for p in sorted(int(t) + 1 for t in divisors(d) if isprime(t + 1)):
        options = [(1, 1)]
        power, degree = p, p - 1
        while d % degree == 0:
            options.append((power, degree))
            power *= p
            degree *= p
        states = {
            (n * power, phi * degree)
            for n, phi in states
            for power, degree in options
            if d % (phi * degree) == 0
        }
    return sorted(n for n, _ in states)


def minus_one_power(a, modulus):
    """Return j>=1 with a**j=-1 mod modulus, or exhaust the finite orbit."""
    if modulus < 2 or gcd(a, modulus) != 1:
        return None
    value = a % modulus
    exponent = 1
    while True:
        if value == modulus - 1:
            return exponent
        if value == 1:
            return None
        value = value * a % modulus
        exponent += 1


def representation_obstructions(block_order, phase):
    """Necessary conditions for Q(zeta_phase) -> (M_d(Q), transpose)."""
    degree = int(totient(phase))
    if block_order % degree:
        return [{"kind": "field_degree", "degree": degree}]
    if phase <= 2:
        return []
    failures = []
    multiplicity = block_order // degree
    # |disc Q(zeta_k)| = product_{p|k} p**(v_p(k)*phi(k)-phi(k)/(p-1)).
    discriminant_parities = {
        int(p): int((int(a) * degree - degree // (int(p) - 1)) % 2)
        for p, a in factorint(phase).items()
    }
    nonsquare_primes = sorted(p for p, parity in discriminant_parities.items() if parity)
    if multiplicity % 2 and nonsquare_primes:
        failures.append(
            {
                "kind": "hermitian_trace_determinant",
                "multiplicity": multiplicity,
                "discriminant_squareclass_primes": nonsquare_primes,
            }
        )
    odd_part = phase
    while odd_part % 2 == 0:
        odd_part //= 2
    if odd_part > 1 and minus_one_power(2, odd_part) is None and block_order % 8:
        failures.append(
            {
                "kind": "dyadic_split_conjugation",
                "odd_conductor_part": odd_part,
                "minus_one_not_in_subgroup_generated_by_two": True,
                "required_dimension_divisor": 8,
            }
        )
    return failures


def source_obstructions(order, phase):
    """Vanishing-sum, real-order and determinant-norm obstructions only."""
    failures = []
    if phase == 2 and order not in (1, 2) and order % 4:
        failures.append({"kind": "real_hadamard_order"})
    if order > 1 and lam_leung_obstructs_butson(order, phase):
        failures.append({"kind": "vanishing_sum_length"})
    for p, valuation in factorint(order).items():
        p, valuation = int(p), int(valuation)
        if phase % p == 0 or (order * valuation) % 2 == 0:
            continue
        exponent = minus_one_power(p, phase)
        if exponent is not None:
            failures.append(
                {
                    "kind": "source_determinant_norm",
                    "prime": p,
                    "valuation_of_source_order": valuation,
                    "self_conjugate_power": exponent,
                }
            )
    return failures


def gaussian_terminal_degree_test(block_order, phase):
    """Necessary degree budget after adjoining i for a Gaussian last stage."""
    compositum_phase = lcm(phase, 4)
    absolute_degree = int(totient(compositum_phase))
    relative_degree = absolute_degree // 2
    first_stage_degree = block_order // 2
    return {
        "compositum_phase": compositum_phase,
        "compositum_absolute_degree": absolute_degree,
        "minimum_relative_degree_over_Q_i": relative_degree,
        "allocated_first_stage_degree": first_stage_degree,
        "passes_necessary_degree_test": block_order % 2 == 0
        and first_stage_degree % relative_degree == 0,
    }


def dyadic_flat_source_reduction(block_order, source_order, phase):
    """The fixed-frame flat-alphabet obstruction when d=4 mod 8 and n is odd.

    Flat Hadamard blocks have one common row-product sign each. Equal
    signs make their quotient integral in Z_(2). A rational orthogonal
    field representation cannot carry nontrivial odd torsion integrally
    at two: reduction modulo two fixes the all-ones vector, whereas the
    determinant of U-I is an odd cyclotomic norm. Thus all flat phases
    occupy at most two cosets of a 2-power root subgroup. Odd source
    order then forces, up to a common phase, the sixth/twelfth alphabets.

    This classifies the used entry alphabet. Unused auxiliary field
    extensions, including proper cyclotomic subfields, remain possible.
    """
    applicable = block_order % 8 == 4 and source_order > 1 and source_order % 2 == 1
    if not applicable:
        return {"applicable": False}
    full_torsion = phase if phase % 2 == 0 else 2 * phase
    available = [k for k in (6, 12) if full_torsion % k == 0]
    reduced = [
        {"phase": k, "source_obstructions": source_obstructions(source_order, k)} for k in available
    ]
    admitted = [r["phase"] for r in reduced if not r["source_obstructions"]]
    return {
        "applicable": True,
        "scope": "uniform faithful field representation with one common scalar frame",
        "global_phase_change_allowed": True,
        "required_used_source_alphabets": [6, 12],
        "available_reduced_alphabets": available,
        "reduced_source_tests": reduced,
        "surviving_reduced_alphabets": admitted,
        "excluded": not admitted,
        "does_not_exclude_unused_auxiliary_extensions": True,
    }


def target_factor_audit(target):
    """Enumerate every possible nontrivial sign-block factor and coefficient field.

    The even phases are the full torsion orders of the cyclotomic fields;
    odd conductors give duplicate fields and are retained in the degree
    inventory but omitted from the detailed source tests. A small-block
    classification already proved in HADAMARD_SCALAR_TRANSFER.md restricts
    the block-order-two case to Gaussian scalars.
    """
    branches = []
    for divisor in divisors(target):
        d = int(divisor)
        n = target // d
        if n <= 1 or (d != 2 and d % 4):
            continue
        phases = cyclotomic_orders_with_degree_dividing(d)
        candidates = []
        for k in phases:
            if k % 2:
                continue
            failures = representation_obstructions(d, k) + source_obstructions(n, k)
            flat_reduction = dyadic_flat_source_reduction(d, n, k)
            if flat_reduction["applicable"] and flat_reduction["excluded"]:
                failures.append(
                    {
                        "kind": "dyadic_flat_source_alphabet",
                        "reduced_alphabets": flat_reduction["available_reduced_alphabets"],
                        "reduced_source_tests": flat_reduction["reduced_source_tests"],
                    }
                )
            if d == 2 and k != 4:
                failures.append({"kind": "order_two_flat_block_classification"})
            candidates.append(
                {
                    "phase": k,
                    "field_degree": int(totient(k)),
                    "field_multiplicity": d // int(totient(k)),
                    "gaussian_terminal_degree_test": gaussian_terminal_degree_test(d, k),
                    "dyadic_flat_source_reduction": flat_reduction,
                    "obstructions": failures,
                    "survives_necessary_tests": not failures,
                }
            )
        branches.append(
            {
                "source_order": n,
                "block_order": d,
                "block_prime_factorization": {str(p): int(a) for p, a in factorint(d).items()},
                "all_root_orders_passing_degree_test": phases,
                "coefficient_fields": candidates,
                "surviving_phases": [c["phase"] for c in candidates if not c["obstructions"]],
            }
        )
    return {
        "target_order": target,
        "target_prime_factorization": {str(p): int(a) for p, a in factorint(target).items()},
        "branches": branches,
        "necessary_only": True,
        "source_or_lift_existence_not_asserted": True,
    }


def auxiliary_degree_167_subfield():
    """Explicit real Gaussian-period field of degree 167 and its trace basis.

    The field is a proper subfield of Q(zeta_2339), so representing this
    field does not require representing a primitive 2339th root. This is
    a coefficient-field and trace-metric witness, not a dense-frame lift.
    """
    p, q, h = 167, 2339, 14
    if not (isprime(p) and isprime(q) and q - 1 == p * h):
        raise ValueError("invalid auxiliary field parameters")
    g = int(primitive_root(q))
    subgroup = {pow(g, p * j, q) for j in range(h)}
    cosets = [{pow(g, i, q) * x % q for x in subgroup} for i in range(p)]
    if (
        len(subgroup) != h
        or q - 1 not in subgroup
        or any(len(c) != h for c in cosets)
        or set().union(*cosets) != set(range(1, q))
        or sum(map(len, cosets)) != q - 1
    ):
        raise ValueError("the real period field cosets do not form the required partition")
    # eta_j=sum_{a in C_j} zeta_q^a. Each eta_j has trace -1.
    # Tr(eta_i eta_j)=(q*#{a+b=0}-h*h)/h = q delta_ij-h.
    if any(
        len(cosets[0].intersection({(-x) % q for x in c})) != (h if j == 0 else 0)
        for j, c in enumerate(cosets)
    ):
        raise ValueError("the Gaussian-period trace calculation failed")
    jacobsthal = Matrix(
        p,
        p,
        lambda i, j: 0 if i == j else (1 if pow((j - i) % p, (p - 1) // 2, p) == 1 else -1),
    )
    unit, all_ones = eye(p), ones(p)
    if not (
        jacobsthal + jacobsthal.T == Matrix.zeros(p)
        and jacobsthal * jacobsthal.T == p * unit - all_ones
        and jacobsthal * all_ones == Matrix.zeros(p)
    ):
        raise ValueError("the rational trace-basis construction failed its character identity")
    gram = q * unit - h * all_ones
    numerator = 171 * unit - 7 * jacobsthal + 55 * all_ones
    # S=numerator/(4q): S^T G S=I. Verify with integers throughout.
    if numerator.T * gram * numerator != 16 * q * q * unit:
        raise ValueError("the auxiliary period basis is not exactly trace orthonormal")
    return {
        "prime_conductor": q,
        "primitive_root": g,
        "subgroup_order": h,
        "subgroup_elements": sorted(subgroup),
        "subgroup_contains_minus_one": True,
        "period_field_degree": p,
        "period_field_is_totally_real": True,
        "period_definition": "eta_j=sum(zeta_2339**a for a in g**j H)",
        "trace_gram": "2339 I_167 - 14 J_167",
        "orthonormal_basis_change": "(171 I_167 - 7 Q_167 + 55 J_167)/(4*2339)",
        "jacobsthal_definition": "Q_ij=quadratic_character_167(j-i)",
        "basis_trace_gram_exactly_verified": True,
        "compositum_with_Q_zeta12_degree": 668,
        "relative_degree_over_Q_i": 334,
        "full_ambient_cyclotomic_field_not_represented": True,
        "full_ambient_conductor": 12 * q,
        "full_ambient_degree": 4 * (q - 1),
        "displayed_separable_basis_is_not_a_dense_lift": True,
        "dense_sign_frame_constructed": False,
    }

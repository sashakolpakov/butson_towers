"""Replay the finite necessary-arithmetic scan for odd C3 multiplier charts.

Let ``m=3*k`` be odd, let ``N=(p-1)/m``, and normalize the value at zero
to ``+e0``.  If the nonzero coordinate-zero orbits have signed count ``a``
and the other orbits have signed count ``b``, the four coordinate sums are
``(1+m*a, k*b, k*b, k*b)``.  Their squared norm must equal ``p``.  Hence

    N = 2*a + k*(3*a*a+b*b),    p = 1+3*k*N.

We normalize ``b>0``; then ``a`` and ``b`` have the same parity.  Put
``q=3*a*a+b*b``.  These conditions imply ``q>=4`` and ``q+2*a>=2``, so
``N-(4*k-2) >= q+2*a-2 >= 0``.  Consequently a bound on nonzero orbits
gives the unconditional arithmetic bound ``p <= 1+3*N*(N+2)/4`` within
this chart.  The prime-first enumeration uses this bound, independently
of the parameter-first enumeration.

This module constructs no T-function, parity map, frame or Hadamard matrix.
The sign dimensions are explicitly conditional estimates, and the Turyn
tag tests only its prime-power arithmetic hypothesis.  Neither a survivor
nor failure of that tag establishes an existence or novelty result.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import time

from sympy import divisors, factorint, isprime, primerange
from sympy.ntheory import n_order


ROOT = Path(__file__).resolve().parents[2]
WILLIAMSON_INPUT = Path("experiments/scalar_lift/williamson_inputs.json")
SCHEMA = "finite-odd-C3-multiplier-target-scan-v1"


def prime_bound(max_orbits: int) -> int:
    """Return the proved prime cutoff for at most this many nonzero orbits."""
    if isinstance(max_orbits, bool) or not isinstance(max_orbits, int) or max_orbits < 2:
        raise ValueError("max_orbits must be an integer at least 2")
    return 1 + 3 * max_orbits * (max_orbits + 2) // 4


def parameter_candidates(max_orbits: int) -> set[tuple[int, int, int, int]]:
    """Enumerate ``(p,m,a,b)`` from the expanded integer norm equation."""
    cutoff = prime_bound(max_orbits)
    result = set()
    for a in range(-max_orbits, max_orbits + 1):
        b_squared_bound = max_orbits - 2 * a - 3 * a * a
        if b_squared_bound < 1:
            continue
        for b in range(1, math.isqrt(b_squared_bound) + 1):
            if (a - b) % 2:
                continue
            q = 3 * a * a + b * b
            for k in range(1, (max_orbits - 2 * a) // q + 1, 2):
                N = 2 * a + k * q
                p, m = 1 + 3 * k * N, 3 * k
                if not 0 < N <= max_orbits or not isprime(p):
                    continue
                if N % 2 or p > cutoff or p != (1 + m * a) ** 2 + 3 * (k * b) ** 2:
                    raise ArithmeticError("parameter enumeration violated the norm or prime bound")
                if abs(a) + b > N:
                    raise ArithmeticError("parameter enumeration violated the count bound")
                result.add((int(p), int(m), a, b))
    return result


def prime_candidates(max_orbits: int) -> tuple[set[tuple[int, int, int, int]], int]:
    """Independently enumerate norm representations ``p=x²+3y²`` first.

    Choose ``y>0`` and the sign of ``x`` with ``x=1 mod3``.  Every allowed
    ``k`` is an odd positive divisor of ``gcd((x-1)/3,y)``.  This routine
    does not assume uniqueness of the norm representation.
    """
    result = set()
    representations = 0
    for p in primerange(7, prime_bound(max_orbits) + 1):
        if p % 3 != 1:
            continue
        for y in range(1, math.isqrt(p // 3) + 1):
            square = int(p) - 3 * y * y
            x = math.isqrt(square)
            if x * x != square:
                continue
            representations += 1
            x = x if x % 3 == 1 else -x
            common_divisor = math.gcd((x - 1) // 3, y)
            for k in divisors(common_divisor):
                k = int(k)
                if k % 2 == 0:
                    continue
                m = 3 * k
                if (p - 1) % m:
                    raise ArithmeticError("norm representation did not give an allowed multiplier")
                N = int(p - 1) // m
                if N > max_orbits:
                    continue
                a, b = (x - 1) // m, y // k
                if N % 2 or b <= 0 or (a - b) % 2:
                    raise ArithmeticError("prime enumeration violated parity or positivity")
                result.add((int(p), m, a, b))
    return result, representations


def _candidate_row(parameter: tuple[int, int, int, int]) -> dict:
    p, m, a, b = parameter
    N = (p - 1) // m
    L = N // 2
    order2 = int(n_order(2, p))
    order4 = order2 // math.gcd(order2, 2)
    frequency_subgroup_order = math.lcm(order4, m)
    extension_index = frequency_subgroup_order // m
    half_frequency_orbits = (p - 1) // (2 * frequency_subgroup_order)
    if N != 2 * half_frequency_orbits * extension_index:
        raise ArithmeticError("frequency subgroup indices disagree")
    prime_factors = {str(q): int(power) for q, power in factorint(2 * p - 1).items()}
    return {
        "p": p,
        "m": m,
        "a": a,
        "b": b,
        "N": N,
        "n": N + 1,
        "correlation_classes": L,
        "coordinate_sums": [1 + m * a] + [(m // 3) * b] * 3,
        "admissible_head_counts": {"minimum": abs(a), "maximum": N - b, "step": 2},
        "headless_allowed": a == 0,
        "conditional_generic_free_bits_estimate": max(0, L - 2),
        "conditional_extreme_head_free_bits_estimate": max(0, L - max(2, abs(a))) if a else None,
        "ord_p_2": order2,
        "ord_p_4": order4,
        "frequency_extension_index": extension_index,
        "half_frequency_orbit_count": half_frequency_orbits,
        "frequency_inverse_structure": (
            "paired distinct H-orbits" if extension_index % 2 else "self-inverse H-orbits"
        ),
        "odd_single_factor_size_compatible": (
            half_frequency_orbits == 1 and extension_index % 2 == 1
        ),
        "turyn_arithmetic_gate": {
            "q": 2 * p - 1,
            "q_factorization": prime_factors,
            "q_is_prime_power": len(prime_factors) == 1,
        },
    }


def scan(max_orbits: int = 160) -> dict:
    """Return the complete replay report after both finite enumerations agree."""
    start = time.monotonic()
    parameter_first = parameter_candidates(max_orbits)
    after_parameters = time.monotonic()
    prime_first, representations = prime_candidates(max_orbits)
    after_primes = time.monotonic()
    if parameter_first != prime_first:
        raise ArithmeticError(
            f"independent scans disagree: parameter-only={sorted(parameter_first-prime_first)}, "
            f"prime-only={sorted(prime_first-parameter_first)}"
        )
    parameters = sorted(parameter_first, key=lambda z: ((z[0] - 1) // z[1], z[0], -z[1], z[2], z[3]))
    rows = [_candidate_row(parameter) for parameter in parameters]
    source_bytes = (ROOT / WILLIAMSON_INPUT).read_bytes()
    literal_source = json.loads(source_bytes)
    widths = [entry["order"] for entry in literal_source["entries"]]
    if widths != sorted(set(widths)) or any(type(width) is not int or width <= 0 for width in widths):
        raise ValueError("literal Williamson source widths must be distinct increasing positive integers")
    parameter_digest = hashlib.sha256(
        json.dumps([list(parameter) for parameter in parameters], separators=(",", ":")).encode()
    ).hexdigest()
    return {
        "schema": SCHEMA,
        "scope": {
            "object": "Necessary integer row-sum arithmetic for prime periodic T-functions with an odd C3 multiplier chart.",
            "normalization": "f(0)=+e0; the multiplier fixes coordinate0 and cycles coordinates1,2,3; b>0.",
            "construction_or_search_performed": False,
            "existence_nonexistence_or_novelty_claim": False,
            "unlisted_orders_assumed_unknown": False,
        },
        "limits": {"max_nonorigin_orbits": max_orbits, "derived_prime_max": prime_bound(max_orbits)},
        "parameter_definitions": {
            "m": "3*k, with k positive and odd; m divides p-1.",
            "N": "(p-1)/m, the number of nonzero multiplication-by-mu orbits; N is even.",
            "n": "N+1, including the representative0.",
            "a": "Signed count of the N orbits assigned to coordinate0.",
            "b": "Positive signed count of all other orbits; a and b have the same parity.",
            "norm_equation": "p=(1+m*a)^2+3*((m/3)*b)^2; equivalently N=2*a+k*(3*a^2+b^2).",
            "admissible_head_counts": "h=abs(a),abs(a)+2,...,N-b; signed counts then determine all four sign totals.",
        },
        "prime_bound_proof": {
            "formula": "p <= 1+3*N*(N+2)/4",
            "steps": [
                "For b>0 and a=b mod2, q=3*a^2+b^2 satisfies q>=4 and q+2*a>=2.",
                "N-(4*k-2)=k*(q-4)+2*a+2 >= q+2*a-2 >= 0.",
                "Therefore k <= (N+2)/4 and p=1+3*k*N <= 1+3*N*(N+2)/4.",
            ],
            "even_multiplier_exclusion": (
                "For even m divisible by3, mu^(m/2)=-1 and the coordinate permutation to that power is identity. "
                "Thus f(-x)=f(x); pairing the summands of C(s) by x -> -x-s leaves one summand1, so C(s) is odd."
            ),
        },
        "independent_enumeration": {
            "parameter_first": "Enumerate a,b,k from N=2*a+k*(3*a^2+b^2), with b>0, a=b mod2, k positive odd, and p prime.",
            "prime_first": "Enumerate prime p=x^2+3*y^2 within the proved bound; y>0, x=1 mod3; k runs over odd divisors of gcd((x-1)/3,y).",
            "agree": True,
            "prime_first_norm_representations_checked": representations,
            "candidate_parameter_sha256": parameter_digest,
        },
        "sign_dimension_estimate_scope": {
            "generic_formula": "max(0,L-2), where L=N/2.",
            "assumptions": "The full correlation rank and count-parity dependency pattern of the validated determinant charts, with no additional forced signs; actual axis words require exact elimination.",
            "extreme_head_formula": "max(0,L-max(2,abs(a))) when a is nonzero and h=abs(a), under the corresponding additional fixed-sign rank assumptions.",
            "headless_caveat": "a=0,h=0 requires a tail gauge; no extreme-head estimate is supplied.",
            "unconditional_kernel_dimension_claim": False,
        },
        "frequency_index_scope": {
            "definitions": "H=<4,mu> in F_p^*; |H|=lcm(ord_p(4),m); e=|H|/m; r=(p-1)/(2*|H|); N=2*r*e.",
            "inverse_structure": "If e is odd, inverse H-orbits pair distinctly; if e is even, -1 lies in H and H-orbits are self-inverse.",
            "size_compatibility_only": "r=1 and e odd is only an arithmetic triage label, not a determinant-map validation.",
        },
        "turyn_tag_scope": {
            "test": "q=2*p-1 is a prime power; q=1 mod4 automatically for odd p.",
            "interpretation": "Only this arithmetic hypothesis is tagged; failure does not establish novelty or rule out other constructions.",
            "retained_theorem_reference": "docs/HADAMARD_FRAME_TARGET_SCAN.md",
        },
        "conditional_endpoint_orders": {
            "formula": "4*p*w",
            "hypotheses": "A completed exact periodic T_p witness and a verified retained Williamson source of width w are required.",
            "williamson_input": WILLIAMSON_INPUT.as_posix(),
            "williamson_input_sha256": hashlib.sha256(source_bytes).hexdigest(),
            "literal_widths": widths,
            "endpoint_constructed_by_this_scan": False,
        },
        "summary": {
            "candidate_rows": len(rows),
            "distinct_primes": len({row["p"] for row in rows}),
            "maximum_p_found": max(row["p"] for row in rows),
        },
        "candidates": rows,
        "runtime_seconds": {
            "parameter_first": after_parameters - start,
            "prime_first": after_primes - after_parameters,
            "total": time.monotonic() - start,
        },
    }


def stable_report(report: dict) -> dict:
    """Copy a report without the sole nonmathematical timing field."""
    return {key: value for key, value in report.items() if key != "runtime_seconds"}


def main(argv: list[str] | None = None) -> int:
    """Write a report and/or compare its stable content with a saved report."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-orbits", type=int, default=160, help="maximum number of nonzero multiplier orbits (default:160)")
    parser.add_argument("--output", type=Path, help="write the complete JSON report to this path")
    parser.add_argument("--check", type=Path, help="require exact agreement with this JSON report, ignoring runtime_seconds only")
    args = parser.parse_args(argv)
    try:
        result = scan(args.max_orbits)
        if args.check:
            expected = json.loads(args.check.read_text())
            if stable_report(result) != stable_report(expected):
                raise ValueError(f"stable mathematical report differs from {args.check}")
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(json.dumps(result, indent=2) + "\n")
    except (ArithmeticError, OSError, ValueError) as exc:
        parser.exit(1, f"multiplier scan failed: {exc}\n")
    if not args.output and not args.check:
        print(json.dumps(result, indent=2))
    else:
        print(json.dumps({"summary": result["summary"], "limits": result["limits"], "check_passed": bool(args.check)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

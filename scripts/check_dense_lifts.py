#!/usr/bin/env python3
"""Replay dense sign-valued lifts and their exact field and Gram identities."""

from __future__ import annotations

import argparse
import hashlib
import json
from itertools import product
from pathlib import Path

from sympy import I, Matrix, eye, kronecker_product, simplify

from hadamard.cgw import is_cgw_exponent_matrix
from hadamard.dense_lift import (
    FIFTH_POWER_QUATERNARY_BLOCKS,
    lift_unreal_tenth_hadamard,
)

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REPORT = ROOT / "certificates/scalar-lift/dense-lifts.json"


def _clean(matrix):
    return matrix.applyfunc(lambda x: simplify(x.expand()))


def norm_f9_hadamard() -> Matrix:
    """Paley H20 from N_(F9/F3)(a+bu)=a^2+b^2, u^2=-1."""
    points = tuple(product(range(3), repeat=2))

    def character(a, b):
        norm = (a * a + b * b) % 3
        return 0 if norm == 0 else (1 if norm == 1 else -1)

    q = Matrix(
        9, 9, lambda i, j: character(points[j][0] - points[i][0], points[j][1] - points[i][1])
    )
    if not (q == q.T and q * q.T == 9 * eye(9) - Matrix.ones(9)):
        raise ValueError("the finite-field norm character identity failed")
    c = Matrix.zeros(10)
    c[0, 1:] = Matrix.ones(1, 9)
    c[1:, 0] = Matrix.ones(9, 1)
    c[1:, 1:] = q
    if c * c.T != 9 * eye(10):
        raise ValueError("the norm-character conference identity failed")
    a, b = c + eye(10), c - eye(10)
    return a.row_join(b).col_join(b.row_join(-a))


def _record(name, h):
    n = h.rows
    if not (h.cols == n and set(h) <= {-1, 1} and h * h.T == h.T * h == n * eye(n)):
        raise ValueError(f"not a dense real Hadamard matrix: {name}")
    quotient_counts = [
        sum(h[i, k] == h[j, k] for k in range(n)) for i in range(n) for j in range(i)
    ]
    if quotient_counts != [n // 2] * (n * (n - 1) // 2):
        raise ValueError("integer C2-character balance failed")
    rows = [[int(x) for x in row] for row in h.tolist()]
    digest = hashlib.sha256(json.dumps(rows, separators=(",", ":")).encode()).hexdigest()
    return {
        "name": name,
        "real_order": n,
        "entry_alphabet": [-1, 1],
        "zero_count": 0,
        "both_integer_grams_valid": True,
        "every_row_quotient_has_equal_sign_counts": True,
        "matrix_sha256": digest,
        "known_order_control": True,
    }


def build_report():
    blocks = [Matrix([[I**e for e in row] for row in b]) for b in FIFTH_POWER_QUATERNARY_BLOCKS]
    m = blocks[0] / 2
    if not (
        _clean(m * m.conjugate().T) == eye(4)
        and _clean(sum((m**e for e in range(5)), Matrix.zeros(4))) == Matrix.zeros(4)
        and all(_clean(2 * m**e) == blocks[e - 1] for e in range(1, 5))
        and all(set(b) <= {1, -1, I, -I} for b in blocks)
    ):
        raise ValueError("the fifth-root star representation or flat powers failed")

    p = Matrix.hstack(*(m**j * eye(4)[:, 0] for j in range(4)))
    trace = Matrix(4, 4, lambda i, j: _clean(m ** (i + j)).trace())
    moments = Matrix([_clean(m**j)[0, 0] for j in range(4)])
    coefficients = _clean(trace.inv() * moments)
    theta = _clean(sum((coefficients[j] * m**j for j in range(4)), Matrix.zeros(4)))
    gram = Matrix(4, 4, lambda i, j: _clean(theta * m ** (j - i)).trace())
    basis = _clean(p.inv())
    if not (
        theta == _clean(theta.conjugate().T)
        and all(simplify(theta[:j, :j].det()) > 0 for j in range(1, 5))
        and _clean(basis.conjugate().T * gram * basis) == eye(4)
    ):
        raise ValueError("the weighted relative trace basis failed")

    five = tuple(tuple((3 * (i + j) ** 2 + 1) % 5 for j in range(5)) for i in range(5))
    ten = tuple(
        tuple((2 * five[i // 2][j // 2] + 5 * (i % 2) * (j % 2)) % 10 for j in range(10))
        for i in range(10)
    )
    h40 = Matrix(lift_unreal_tenth_hadamard(five, 5))
    h80 = Matrix(lift_unreal_tenth_hadamard(ten, 10))

    h4 = kronecker_product(Matrix([[1, 1], [1, -1]]), Matrix([[1, 1], [1, -1]]))
    projected_source = tuple(
        tuple(
            (
                5 * int(bool(h4[i % 4, j % 4] < 0))
                + (2 * (j % 4) if i >= 4 else 0)
                + (5 if i >= 4 and j >= 4 else 0)
            )
            % 10
            for j in range(8)
        )
        for i in range(8)
    )
    if not is_cgw_exponent_matrix(projected_source, 10, 8):
        raise ValueError("the projection control source is not Hadamard")
    projected = Matrix([[(-1) ** e for e in row] for row in projected_source])

    return {
        "scope": "Exact dense-output identities; no weighing-matrix completion or repair",
        "sign_group": "F3^times = C2, interpreted as integer +1 and -1",
        "tenth_root_transfer": {
            "source": "BH(n,10) with no entries +/-1",
            "intermediate": "BH(4n,4)",
            "output": "real Hadamard of order 8n, all entries +/-1",
            "source_field": "Q(zeta_5)",
            "extension_field": "Q(i,zeta_5)=Q(zeta_20)",
            "extension_degree": 2,
            "relative_trace_base": "Q(i)",
            "theta_coefficients_in_1_t_t2_t3": [str(v) for v in coefficients],
            "orthonormal_basis_power_coordinates_by_column": [
                [str(v) for v in row] for row in basis.tolist()
            ],
            "all_four_nonidentity_powers_quaternary_flat": True,
            "reference": "https://arxiv.org/html/1707.08815v3#S3.SS1",
        },
        "sources": {"unreal_fifth_order_five": five, "unreal_tenth_order_ten": ten},
        "known_order_controls": [
            _record("Finite-field norm F9 to F3, Paley", norm_f9_hadamard()),
            _record("Cyclotomic fifth-root extension", h40),
            _record("Cyclotomic tenth-root extension", h80),
            _record("Multiplicative projection BH(8,10) to signs", projected),
        ],
        "new_real_hadamard_order_claimed": False,
        "novelty_claimed_for_transfer": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", type=Path)
    args = parser.parse_args()
    if args.output is not None and args.check is not None:
        parser.error("use either --output or --check")
    encoded = json.dumps(build_report(), indent=2) + "\n"
    if args.check is not None:
        if json.loads(encoded) != json.loads(args.check.read_text()):
            raise SystemExit("report replay differs")
        print(f"Verified {args.check}: four dense controls and exact field identities.")
    else:
        destination = args.output or DEFAULT_REPORT
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(encoded)
        print(f"Wrote {destination}; all outputs are dense signs; orders already known.")


if __name__ == "__main__":
    main()

"""Verified small-source inputs for a subsequently excluded bordered assembly.

The odd-border obstruction now excludes the stated two-coset completions.
The seed matrices remain valid. See docs/BORDERED_ODD_BORDER_OBSTRUCTION.md.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from sympy import isprime

from hadamard.cgw import is_cgw_exponent_matrix

ROOT = Path(__file__).resolve().parents[2]


def paley_748():
    p = 373
    assert isprime(p) and p % 4 == 1
    character = np.zeros(p, dtype=np.int64)
    character[(np.arange(1, p, dtype=np.int64) ** 2) % p] = 1
    character[(character == 0) & (np.arange(p) != 0)] = -1
    conference = np.ones((p + 1, p + 1), dtype=np.int64)
    conference[0, 0] = 0
    conference[1:, 1:] = character[(np.arange(p)[:, None] - np.arange(p)) % p]
    eye = np.eye(p + 1, dtype=np.int64)
    if not np.array_equal(conference, conference.T) or not np.array_equal(
        conference @ conference.T, p * eye
    ):
        raise ValueError("the Paley conference identity failed")
    return np.block([[conference + eye, conference - eye],
                     [conference - eye, -conference - eye]])


def check_sign_hadamard(matrix):
    n = len(matrix)
    if matrix.shape != (n, n) or not np.all(np.abs(matrix) == 1):
        raise ValueError("a square all-sign matrix is required")
    words = [sum(1 << j for j, value in enumerate(row) if value == 1) for row in matrix]
    for i, word in enumerate(words):
        if any((word ^ words[j]).bit_count() * 2 != n for j in range(i)):
            raise ValueError("an exact row inner product is nonzero")
    return {"order": n, "all_row_pairs_checked": n * (n - 1) // 2,
            "binary_sha256": hashlib.sha256(((matrix + 1) // 2).astype(np.uint8).tobytes()).hexdigest()}


def border_source(order, phase):
    if order == 3:
        cube = next(e for e in (phase // 3, 2 * phase // 3) if e % 3 == 1)
        source = np.array([[cube if i == j else 0 for j in range(3)] for i in range(3)])
    elif order == 7:
        records = json.loads((ROOT / "certificates/scalar-lift/large-joint-targets.json").read_text())
        seed = next(s["exponents"] for s in records["verified_sources"] if s["order"] == 7)
        source = (np.array(seed, dtype=np.int64) * (phase // 6) + 2) % phase
        support = (source % 3 == 1).astype(np.int64)
        if int(support[0].sum()) % 2 == 0:
            cube = next(e for e in (phase // 3, 2 * phase // 3) if e % 3 == 1)
            source = (cube - source) % phase
    else:
        raise ValueError("the certified border orders are three and seven")
    if not is_cgw_exponent_matrix(source.tolist(), phase, order):
        raise ValueError("the border source failed its exact cyclotomic Gram identity")
    support = (source % 3 == 1).astype(np.int64)
    if np.any(source % 3 == 2) or not np.array_equal(
        (support @ support.T) % 2, np.eye(order, dtype=np.int64)
    ):
        raise ValueError("the border source failed its required coset gauge")
    return source


def assemble_exponents(x, y, t, d, phase):
    x, y, t, d = (np.asarray(item, dtype=np.int64) for item in (x, y, t, d))
    s, border = len(x), len(d)
    if x.shape != (s, s) or y.shape != (s, s) or t.shape != (s, border) or d.shape != (border, border):
        raise ValueError("incompatible bordered block sizes")
    return np.block([[x, y, t], [y, x, t], [-t.T, -t.T, d]]) % phase


def build_report():
    base = paley_748()
    doubled = np.block([[base, base], [base, -base]])
    inputs = {748: base, 1496: doubled}
    full_inputs = [check_sign_hadamard(inputs[s]) for s in (748, 1496)]
    routes = []
    for s, border, phase in ((748, 3, 96), (1496, 7, 48), (1496, 7, 96)):
        t = inputs[s][:, :border]
        if not np.array_equal(t.T @ t, s * np.eye(border, dtype=np.int64)):
            raise ValueError("the partial Butson columns are not orthogonal")
        d = border_source(border, phase)
        routes.append({
            "source_order": 2 * s + border, "phase": phase,
            "real_target_order": (2 * s + border) * (phase // 3),
            "square_seed_order": s, "border_source_order": border,
            "T_signs": t.tolist(), "D_exponents": d.tolist(),
            "T_star_T_equals_sI": True, "D_D_star_equals_tI": True,
            "D_binary_coset_support_orthogonal": True,
            "X_and_Y_constructed": False,
            "specified_adjoint_border_completion_excluded": True,
            "exclusion_holds_throughout_dyadic_extensions": True,
        })
    return {"scope": "Verified seed inputs; specified adjoint-border source completions are excluded",
            "Paley_prime": 373, "full_real_seed_checks": full_inputs, "routes": routes,
            "exclusion_proof": "docs/BORDERED_ODD_BORDER_OBSTRUCTION.md",
            "general_Butson_sources_excluded": False,
            "new_Butson_sources_constructed": False, "new_real_Hadamard_order_constructed": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", type=Path)
    args = parser.parse_args()
    report = build_report()
    if args.check and json.loads(args.check.read_text()) != report:
        raise SystemExit("bordered seed certificate differs from its replay")
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"verified_real_seed_orders": [748, 1496],
                      "excluded_adjoint_border_source_orders": [1499, 2999],
                      "source_completions_constructed": False}))


if __name__ == "__main__":
    main()

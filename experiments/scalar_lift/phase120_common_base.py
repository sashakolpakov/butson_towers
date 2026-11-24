"""Exact phase-120 partial sign lift of width 32 over a shared eighth-root field.

The BH(15,120) control is complete. The order-1499 source needed for the
conditional real order 47968 is not supplied by this construction.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from functools import cache
from math import gcd
from pathlib import Path

import numpy as np
from sympy import Poly, Symbol, cyclotomic_poly

from experiments.scalar_lift.coupled_seven_relative_extension import extension_control
from experiments.scalar_lift.dyadic_partial_lift import local_frame
from experiments.scalar_lift.fifth_source_assembly import exact_column_gram
from hadamard.dense_lift import FIFTH_POWER_QUATERNARY_BLOCKS


@cache
def common_base():
    """Return the commuting rational actions, and raw and sign frames.

    Coordinate ordering is (fifth index, cubic index, eighth-root index).
    Each base-field block has size four, so this is a relative tensor
    product of ranks four and two over Q(zeta_8), of rational degree 32.
    """
    identity4 = np.eye(4, dtype=np.int64)
    r = np.zeros((4, 4), dtype=np.int64)
    r[:3, 1:] = np.eye(3, dtype=np.int64)
    r[3, 0] = -1
    a = identity4 + r + np.linalg.matrix_power(r, 3)
    cubic8 = -np.block([[identity4, a], [-a.T, identity4]])
    r8 = np.kron(np.eye(2, dtype=np.int64), r)
    j8 = r8 @ r8
    fifth = [np.block([[np.linalg.matrix_power(j8, int(e)) for e in row]
                       for row in exponents])
             for exponents in FIFTH_POWER_QUATERNARY_BLOCKS]
    cubic = np.kron(np.eye(4, dtype=np.int64), cubic8)
    dyadic = np.kron(np.eye(4, dtype=np.int64), r8)
    raw = 2 * np.kron(np.eye(4, dtype=np.int64), local_frame())
    frame_numerator = fifth[0] @ raw
    if np.any(frame_numerator % 2):
        raise ValueError("the common-base sign frame is not integral")
    frame = frame_numerator // 2
    return {
        "base_R": r, "norm_three_A": a,
        "cubic_numerator": cubic, "fifth_numerators": fifth,
        "dyadic_action": dyadic, "raw_frame": raw, "sign_frame": frame,
    }


def action_power(exponent):
    """Return an exact integer numerator and denominator for S^exponent."""
    if not isinstance(exponent, (int, np.integer)) or isinstance(exponent, (bool, np.bool_)):
        raise ValueError("an exact integer exponent is required")
    e = int(exponent) % 120
    data = common_base()
    identity = np.eye(32, dtype=np.int64)
    cp = e % 3
    fp = 4 * e % 5
    c = identity if cp == 0 else (data["cubic_numerator"] if cp == 1
                                  else data["cubic_numerator"].T)
    f = identity if fp == 0 else data["fifth_numerators"][fp - 1]
    denominator = (1 if cp == 0 else 2) * (1 if fp == 0 else 2)
    numerator = np.linalg.matrix_power(data["dyadic_action"], 7 * e % 8) @ c @ f
    return numerator, denominator


@cache
def phase_blocks(shifted_flat_frame=True):
    """All 120 integral phase blocks, including those outside the sign alphabet."""
    key = "sign_frame" if shifted_flat_frame else "raw_frame"
    frame = common_base()[key]
    blocks = []
    for e in range(120):
        numerator, denominator = action_power(e)
        product = numerator @ frame
        if np.any(product % denominator):
            raise ValueError("a phase block unexpectedly lost integrality")
        blocks.append(product // denominator)
    return blocks


def admitted_exponents(shifted_flat_frame=True):
    return [e for e in range(120)
            if e % 3 in (0, 1) and e % 5 != (1 if shifted_flat_frame else 0)]


def lift_source(exponents, shifted_flat_frame=True):
    """Verify an exact restricted BH(n,120) and return its sign lift."""
    rows = [list(row) for row in exponents]
    if not rows or any(len(row) != len(rows) for row in rows):
        raise ValueError("a nonempty square source is required")
    if any(not isinstance(e, (int, np.integer)) or isinstance(e, (bool, np.bool_))
           for row in rows for e in row):
        raise ValueError("source exponents must be exact integers")
    source = [[int(e) % 120 for e in row] for row in rows]
    allowed = set(admitted_exponents(shifted_flat_frame))
    if any(e not in allowed for row in source for e in row):
        raise ValueError("the source uses an excluded phase")
    gram = exact_column_gram(source, 120)
    expected = np.zeros_like(gram)
    expected[np.arange(len(source)), np.arange(len(source)), 0] = len(source)
    if not np.array_equal(gram, expected):
        raise ValueError("the source fails its exact cyclotomic Gram")
    blocks = phase_blocks(shifted_flat_frame)
    return np.block([[blocks[e] for e in row] for row in source])


def source_fifteen():
    """A genuine-phase-120 Diţă control in the raw frame's alphabet.

    The outer source is the cubic diagonal-omega matrix of order three;
    the inner source is the admitted chirped Fourier matrix of order five.
    The eighth-root diagonal twist is indexed by inner row and outer column.
    """
    return [[(40 * (i == j) + 15 * a * j
              + 24 * ((3 * (a + b) ** 2 + 1) % 5)) % 120
             for j in range(3) for b in range(5)]
            for i in range(3) for a in range(5)]


def source_thirty_five():
    """Tensor the genuinely coupled BH(7,24) with a fifth-root chirp.

    The field embedding zeta_24 -> zeta_120^25 is an automorphic embedding,
    not the conventional exponent-five embedding. It preserves the admitted
    cubic residue classes. The fifth shift matches the sign-frame alphabet.
    """
    seven = extension_control(24, 0)["source_exponents"]
    return [[(25 * seven[i][j] + 24 * ((3 * (a + b) ** 2 + 3) % 5)) % 120
             for j in range(7) for b in range(5)]
            for i in range(7) for a in range(5)]


def _digest(matrix):
    return hashlib.sha256(np.asarray(matrix, dtype="<i8").tobytes()).hexdigest()


@cache
def build_report():
    data = common_base()
    identity = np.eye(32, dtype=np.int64)
    r, a = data["base_R"], data["norm_three_A"]
    c, fifth, rd = data["cubic_numerator"], data["fifth_numerators"], data["dyadic_action"]
    if not np.array_equal(a @ a.T, 3 * np.eye(4, dtype=np.int64)):
        raise ValueError("the base-field norm-three identity failed")
    if not np.array_equal(c @ c.T, 4 * identity) or not np.array_equal(c @ c, 2 * c.T):
        raise ValueError("the cubic action failed its exact identities")
    if np.any(2 * identity + c + c.T):
        raise ValueError("Phi3 failed")
    for j in range(1, 5):
        if not np.array_equal(fifth[j - 1] @ fifth[j - 1].T, 4 * identity):
            raise ValueError("a fifth action is not orthogonal")
        for k in range(1, 5):
            p = (j + k) % 5
            expected = 4 * identity if p == 0 else 2 * fifth[p - 1]
            if not np.array_equal(fifth[j - 1] @ fifth[k - 1], expected):
                raise ValueError("the fifth power table failed")
    if np.any(2 * identity + sum(fifth)):
        raise ValueError("Phi5 failed")
    for left, right in ((c, fifth[0]), (c, rd), (fifth[0], rd)):
        if not np.array_equal(left @ right, right @ left):
            raise ValueError("the common-base actions do not commute")
    powers = [action_power(e) for e in range(120)]
    generator, generator_denominator = powers[1]
    for e, (numerator, denominator) in enumerate(powers):
        if not np.array_equal(numerator @ numerator.T, denominator**2 * identity):
            raise ValueError("a scalar power is not rational orthogonal")
        following, following_denominator = powers[(e + 1) % 120]
        if not np.array_equal(following_denominator * generator @ numerator,
                              generator_denominator * denominator * following):
            raise ValueError("the exact scalar recurrence failed")
    for exponent, matrix, denominator in ((15, rd, 1), (40, c, 2), (24, fifth[0], 2)):
        numerator, actual_denominator = powers[exponent]
        if not np.array_equal(denominator * numerator, actual_denominator * matrix):
            raise ValueError("a relative tower restriction failed")
    polynomial = Poly(cyclotomic_poly(120, Symbol("x")), Symbol("x"))
    residual = np.zeros((32, 32), dtype=np.int64)
    for (e,), coefficient in polynomial.terms():
        numerator, denominator = powers[e]
        residual += int(coefficient) * (4 // denominator) * numerator
    if np.any(residual) or polynomial.degree() != 32:
        raise ValueError("the degree-32 cyclotomic polynomial failed")
    frame_reports = []
    for shifted in (False, True):
        blocks = phase_blocks(shifted)
        actual = [e for e, block in enumerate(blocks) if np.all(np.abs(block) == 1)]
        if actual != admitted_exponents(shifted):
            raise ValueError("the complete flat alphabet differs from its formula")
        for block in blocks:
            if not np.array_equal(block @ block.T, 32 * identity):
                raise ValueError("a phase frame has the wrong Gram")
        # A direct exact check of every cross-Gram, independent of entry flatness.
        for e in range(120):
            for f in range(120):
                numerator, denominator = powers[(e - f) % 120]
                if not np.array_equal(denominator * (blocks[e] @ blocks[f].T),
                                      32 * numerator):
                    raise ValueError("a phase cross-Gram failed")
        frame_reports.append({"uses_shifted_sign_frame": shifted,
                              "allowed_exponents": actual,
                              "all_120_integer_frame_Grams_verified": True,
                              "all_14400_cross_Grams_verified": True,
                              "frame_sha256": _digest(blocks[0])})
    source = source_fifteen()
    sign = lift_source(source, False)
    shifted_source = [[(e - 24) % 120 for e in row] for row in source]
    if not np.array_equal(sign, lift_source(shifted_source, True)):
        raise ValueError("the source/frame phase-shift equivalence failed")
    if not np.all(np.abs(sign) == 1) or not np.array_equal(
        sign @ sign.T, 480 * np.eye(480, dtype=np.int64)
    ):
        raise ValueError("the complete H480 control failed")
    cross_ratio = (source[0][0] + source[6][6] - source[0][6] - source[6][0]) % 120
    if cross_ratio != 119 or gcd(cross_ratio, 120) != 1:
        raise ValueError("the control lost its primitive gauge-invariant cross-ratio")
    source35 = source_thirty_five()
    sign35 = lift_source(source35)
    if not np.all(np.abs(sign35) == 1) or not np.array_equal(
        sign35 @ sign35.T, 1120 * np.eye(1120, dtype=np.int64)
    ):
        raise ValueError("the coupled order-35 source-to-H1120 control failed")
    return {
        "phase": 120, "dimension": 32, "base_field": "Q(zeta_8)",
        "relative_degrees": [2, 4], "ordinary_rational_tensor_dimension_avoided": 64,
        "coordinate_order": ["fifth index", "cubic index", "base-field index"],
        "base_R": r.tolist(), "norm_three_A": a.tolist(),
        "generator_numerator": generator.tolist(), "generator_denominator": generator_denominator,
        "cyclotomic_polynomial": str(polynomial.as_expr()),
        "Phi120_annihilation_verified": True,
        "minimal_and_characteristic_polynomial_equal_Phi120": True,
        "faithful_order_120_orthogonal_action": True,
        "relative_restrictions": {"S^15": "R", "S^40": "C", "S^24": "F"},
        "frames": frame_reports,
        "sign_frame": data["sign_frame"].tolist(),
        "control": {"source_order": 15, "source_phase": 120,
                    "source_exponents_raw_frame": source,
                    "exact_cyclotomic_Gram_verified": True,
                    "primitive_cross_ratio_exponent": cross_ratio,
                    "cross_ratio_rows_and_columns": [0, 6],
                    "real_order": 480, "full_integer_Gram_verified": True,
                    "sign_matrix_sha256": _digest(sign),
                    "previously_unknown_real_order_claimed": False},
        "coupled_control": {"source_order": 35, "source_phase": 120,
                            "source_exponents_sign_frame": source35,
                            "order_seven_branch": 0,
                            "order_seven_relative_AB_star_nonzero": True,
                            "order_seven_field_embedding_exponent": 25,
                            "exact_cyclotomic_Gram_verified": True,
                            "real_order": 1120, "full_integer_Gram_verified": True,
                            "sign_matrix_sha256": _digest(sign35),
                            "previously_unknown_real_order_claimed": False},
        "conditional_high_route": {"source_order": 1499, "source_phase": 120,
                                   "real_order": 47968, "source_constructed": False,
                                   "source_alphabet_raw_frame": "e mod 3 in {0,1}; e mod 5 != 0",
                                   "new_real_Hadamard_claimed": False},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--verify", type=Path)
    args = parser.parse_args()
    report = build_report()
    if args.verify and json.loads(args.verify.read_text()) != report:
        raise SystemExit("certificate replay differs")
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"phase": 120, "dimension": 32, "admitted_phases": 64,
                      "control_real_orders": [480, 1120], "high_source_constructed": False}))


if __name__ == "__main__":
    main()

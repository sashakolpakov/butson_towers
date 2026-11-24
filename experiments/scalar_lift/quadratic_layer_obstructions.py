"""Exact controls for two necessary quadratic-source layer constraints.

The proofs are in docs/QUADRATIC_LAYER_OBSTRUCTIONS.md.  The controls are
known even-order sources; no missing odd-order source is constructed.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from experiments.scalar_lift.fifth_source_assembly import root_vectors
from experiments.scalar_lift.fifth_source_descent import (
    check_relative_equations,
    check_source,
    coefficient_product,
    split_coefficients,
)
from experiments.scalar_lift.relative_source_tower import sparse_row_product


def binary_rank(matrix):
    work = np.array(matrix, dtype=np.uint8, copy=True) % 2
    rank = 0
    for column in range(work.shape[1]):
        candidates = np.flatnonzero(work[rank:, column])
        if not len(candidates):
            continue
        pivot = rank + int(candidates[0])
        work[[rank, pivot]] = work[[pivot, rank]]
        for row in range(rank + 1, len(work)):
            if work[row, column]:
                work[row] ^= work[rank]
        rank += 1
        if rank == len(work):
            break
    return rank


def parity_dephase(source, upper_phase):
    """Use only dyadic row/column phases, preserving every fifth color."""
    source = np.asarray(source, dtype=np.int64) % upper_phase
    row_bits = source[:, 0] % 2
    column_bits = (source[0] + source[0, 0]) % 2
    result = (source - 5 * row_bits[:, None] - 5 * column_bits[None, :]) % upper_phase
    expected = (source + source[:, :1] + source[:1] + source[0, 0]) % 2
    if not np.array_equal(result % 2, expected):
        raise ValueError("dyadic gauges failed to dephase the parity matrix")
    if np.any((result - source) % 5):
        raise ValueError("dyadic gauges changed a fifth color")
    return result, expected


def collapse_rank_one(source, upper_phase):
    """Certify the conditional lower-field collapse when parity rank <=1."""
    check_source(source, upper_phase)
    gauged, parity = parity_dephase(source, upper_phase)
    rank = binary_rank(parity)
    if rank > 1:
        raise ValueError("this lemma does not cover parity rank greater than one")
    a, b = split_coefficients(gauged, upper_phase)
    lower = upper_phase // 2
    if np.any(coefficient_product(a, b, lower)):
        raise ValueError("the rank-one rectangular support did not decouple")
    collapsed = np.maximum(a, b)
    check_source(collapsed, lower)
    return collapsed, {"dephased_parity_rank": rank,
                       "dyadic_gauge_preserves_fifth_colors": True,
                       "ABstar_zero_after_parity_dephasing": True,
                       "lower_source_complete_Gram_verified": True}


def fourier_control(upper_phase):
    order = upper_phase // 5
    exponents = (upper_phase // 5 + 5 * np.outer(np.arange(order), np.arange(order))) % upper_phase
    gauged = exponents.copy()
    gauged[1] = (gauged[1] + 5) % upper_phase
    reports = []
    for name, source in (("canonical_Fourier", exponents), ("one_odd_row_multiplied_by_eta", gauged)):
        check_source(source, upper_phase)
        relative = check_relative_equations(source, upper_phase)
        collapsed, collapse = collapse_rank_one(source, upper_phase)
        a, b = split_coefficients(source, upper_phase)
        cross = coefficient_product(a, b, upper_phase // 2)
        positions = np.argwhere(np.any(cross != 0, axis=2))
        reports.append({"name": name, "source_exponents": source.tolist(),
                        "relative_equations": relative, "rank_one_collapse": collapse,
                        "collapsed_lower_exponents": collapsed.tolist(),
                        "nonzero_ABstar_entries": [{"row": int(i), "column": int(j),
                                                     "coefficients": cross[i, j].tolist()}
                                                    for i, j in positions]})
    if reports[0]["relative_equations"]["ABstar_is_nonzero"]:
        raise ValueError("the Fourier split unexpectedly has a cross term")
    if not reports[1]["relative_equations"]["ABstar_is_nonzero"]:
        raise ValueError("the gauge control lost its nonzero cross term")
    _, first = parity_dephase(exponents, upper_phase)
    _, second = parity_dephase(gauged, upper_phase)
    if not np.array_equal(first, second):
        raise ValueError("dephased parity changed under a phase gauge")
    return {"upper_phase": upper_phase, "order": order,
            "same_dephased_parity_under_dyadic_gauge": True,
            "ABstar_nonzero_is_not_gauge_invariant": True, "controls": reports}


def even_weighing_control(upper_phase):
    # H=c[[1,eta],[eta,-xi]], where eta^2=xi.  Each layer has weight1.
    color = upper_phase // 5
    source = (color + np.array([[0, 5], [5, 10 + upper_phase // 2]])) % upper_phase
    check_source(source, upper_phase)
    relative = check_relative_equations(source, upper_phase)
    a, b = split_coefficients(source, upper_phase)
    lower = upper_phase // 2
    aa = coefficient_product(a, a, lower)
    bb = coefficient_product(b, b, lower)
    expected = np.zeros_like(aa)
    expected[np.arange(2), np.arange(2), 0] = 1
    if not np.array_equal(aa, expected) or not np.array_equal(bb, expected):
        raise ValueError("the even control lost its two weight-one layers")
    return {"order": 2, "upper_phase": upper_phase,
            "source_exponents": source.tolist(), "relative_equations": relative,
            "both_layer_Grams_equal_identity": True,
            "shows_odd_order_hypothesis_is_necessary": True}


def rank_two_block_identity(upper_phase):
    """Audit the proposed block equations on non-Hadamard test data.

    A nonzero residual is intentional: this checks every coefficient,
    including the xi^-1 cross terms, without claiming a completion.
    """
    lower = upper_phase // 2
    sizes = (2, 2, 2, 2)
    cuts = np.cumsum((0, *sizes))
    rng = np.random.default_rng(721 + upper_phase)
    admitted = np.array([e for e in range(lower) if e % 5])
    blocks = [[rng.choice(admitted, (sizes[i], sizes[j])) for j in range(4)] for i in range(4)]
    parity = np.ones((4, 4), dtype=np.int64) - np.eye(4, dtype=np.int64)
    parity[0] = parity[:, 0] = 0
    upper = np.block([[(2 * blocks[i][j] + 5 * parity[i, j]) % upper_phase
                       for j in range(4)] for i in range(4)])
    a, b = split_coefficients(upper, upper_phase)
    c0 = coefficient_product(a, a, lower) + coefficient_product(b, b, lower)
    c1 = coefficient_product(b, a, lower) + coefficient_product(a, b, lower, shift=-5)
    direct0, direct1 = np.zeros_like(c0), np.zeros_like(c1)

    def product(i, j, k, column, shift=0):
        value = sparse_row_product(blocks[i][j], blocks[k][column], lower)
        roots = root_vectors(lower)
        return value @ roots[(np.arange(roots.shape[1]) + shift) % lower]

    for i in range(4):
        for k in range(4):
            region = np.s_[cuts[i]:cuts[i + 1], cuts[k]:cuts[k + 1]]
            if i == k == 0:
                direct0[region] = sum(product(0, j, 0, j) for j in range(4))
            elif i == k:
                direct0[region] = sum(product(i, j, i, j) for j in range(4))
            elif i == 0:
                direct0[region] = product(0, 0, k, 0) + product(0, k, k, k)
                direct1[region] = sum(product(0, j, k, j, -5) for j in range(1, 4) if j != k)
            elif k == 0:
                direct0[region] = product(i, 0, 0, 0) + product(i, i, 0, i)
                direct1[region] = sum(product(i, j, 0, j) for j in range(1, 4) if j != i)
            else:
                third = 6 - i - k
                direct0[region] = product(i, 0, k, 0) + product(i, third, k, third)
                direct1[region] = product(i, k, k, k) + product(i, i, k, i, -5)
    if not np.array_equal(c0, direct0) or not np.array_equal(c1, direct1):
        raise ValueError("the three-way block equations disagree with the coefficient Gram")
    roots = root_vectors(upper_phase)
    degree = root_vectors(lower).shape[1]
    reconstructed = c0 @ roots[2 * np.arange(degree)] + c1 @ roots[(2 * np.arange(degree) + 5) % upper_phase]
    if not np.array_equal(reconstructed, sparse_row_product(upper, upper, upper_phase)):
        raise ValueError("the block coefficients disagree with the independent upper Gram")
    if not np.any(c1) or binary_rank(parity) != 2:
        raise ValueError("the test data fail to exercise the rank-two cross coupling")
    return {"upper_phase": upper_phase, "test_order": len(upper),
            "block_sizes": list(sizes), "source_exponents": upper.tolist(),
            "parity_type_matrix": parity.tolist(), "dephased_parity_rank": 2,
            "all_block_Gram_formulas_match_full_exact_Gram": True,
            "nonzero_cross_residual_present": True,
            "test_data_are_Not_a_Hadamard_source": True}


def certificate():
    return {"schema": "quadratic-layer-obstructions-v1",
            "coordinate": "eta=zeta_(2d), xi=eta^2=zeta_d; both coefficient channels preserve fifth color",
            "proved_constraints": {
                "parity_rank_at_most_one": "After dyadic gauges the B support is a rectangle; AB*=BA*=0 and A+B is a restricted lower-field Butson.",
                "nontrivial_scalar_Gram_layer_at_odd_order": "Impossible: U=A*B/k is scalar-unitary, a scalar multiple of a Hermitian matrix, and trace zero, forcing even order.",
                "block_diagonal_smaller_Butsons": "If A is block diagonal with complete smaller Butsons and B occupies its complement, B commutes with the diagonal block-size matrix; all block sizes must agree, then the odd-order scalar-Gram obstruction applies.",
                "scope": "Necessary conditions on this quadratic split; not nonexistence of the high sources."},
            "essential_target_steps": [
                {"source_order": n, "upper_phase": q, "excluded_lower_phase": q // 2,
                 "dephased_dyadic_parity_rank_must_be_at_least_two": True,
                 "neither_nontrivial_layer_can_have_scalar_Gram": True,
                 "source_constructed": False} for n, q in ((1499, 40), (2999, 80))],
            "Fourier_gauge_controls": [fourier_control(q) for q in (40, 80)],
            "even_weighing_controls": [even_weighing_control(q) for q in (40, 80)],
            "rank_two_block_formula_audits": [rank_two_block_identity(q) for q in (40, 80)],
            "rank_two_additive_charts": [
                {"source_order": 3 * s + 2, "upper_phase": q, "border_size": 2,
                 "three_core_sizes": [s, s, s],
                 "each_off_diagonal_core_block_rank_at_least": (s * s + 3 * s + 1) // (3 * s + 2),
                 "two_off_diagonal_block_ranks_in_each_column_sum_at_most": s + 2,
                 "source_constructed_or_feasibility_asserted": False}
                for s, q in ((499, 40), (999, 80))],
            "new_unknown_source_constructed": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--output", type=Path)
    group.add_argument("--check", type=Path)
    args = parser.parse_args()
    report = certificate()
    if args.check:
        if json.loads(args.check.read_text()) != report:
            raise ValueError("the quadratic-layer certificate differs from replay")
    else:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"Fourier_control_orders": [8, 16], "even_weighing_control_order": 2,
                      "exact_controls_verified": True, "new_unknown_source_constructed": False}))


if __name__ == "__main__":
    main()

"""Exact variable Butson boundaries and a noncollapsing second-bit residue.

Physical real-Hadamard row pairs receive independent cyclotomic parameters.
The current modular construction uses one parameter for each ten-row copy.
No characteristic-zero square source is claimed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from itertools import product
from math import gcd
from pathlib import Path

import numpy as np

from experiments.scalar_lift.cubic_core_support_lift import (
    boundary_data,
    phase_coordinates,
    template,
)
from experiments.scalar_lift.fifth_source_assembly import exact_column_gram, root_vectors
from experiments.scalar_lift.mixed_cubic_boundaries import doubled_real_frame, real_frame
from experiments.scalar_lift.mixed_cubic_second_ramification import polynomial_gram_from_roots


def primitive_frame_cross_ratio(u, phase):
    ratios = (u-u[:, :1]-u[:1]+u[0, 0]) % phase
    for i, row in enumerate(ratios):
        for j, value in enumerate(row):
            if gcd(int(value), phase) == 1:
                return {"rows": [0, i], "columns": [0, j], "exponent": int(value), "order": phase}
    raise ValueError("the variable frame has no certified primitive cross ratio")


def variable_problem(order, phase, copy_second_bits=None):
    """Reconstruct actual frames, with three physical pairs in each copy."""
    boundary = boundary_data(order, phase)
    s = len(boundary["active"])
    pairs = real_frame()[3] if order == 1499 else [[i, i + 748] for i in range(372)]
    pairs = [list(pair) for pair in pairs if boundary["active"][pair[0]]]
    on = np.asarray(pairs, dtype=np.int64).ravel()
    off = np.flatnonzero(1 - boundary["active"])
    copies = len(on) // 6
    if (len(set(on.tolist())) != len(on) or len(on) != int(boundary["active"].sum())
            or any(not boundary["active"][j] for j in on) or len(on) % 6):
        raise ValueError("the physical row pairs do not partition the active rows")
    if copy_second_bits is None:
        copy_second_bits = [int(i >= 31) for i in range(copies)]
    bits = np.asarray(copy_second_bits, dtype=np.int64)
    if bits.shape != (copies,) or np.any((bits != 0) & (bits != 1)):
        raise ValueError("one binary phase choice is required per ten-row copy")
    small = template()
    correction = np.eye(s, dtype=np.int64)
    b, w = np.zeros((s, s), dtype=np.int64), np.zeros((s, s), dtype=np.int64)
    row_bits = np.zeros(s, dtype=np.int64)
    groups = []
    for index in range(copies):
        group = np.concatenate((on[6*index:6*index+6], off[4*index:4*index+4]))
        where = np.ix_(group, group)
        correction[where], b[where], w[where] = small["T"], small["B"], small["W"]
        row_bits[group] = bits[index]
        groups.append(group.tolist())
    u = (boundary["U_exponents"] + 6*row_bits[:, None]*boundary["A"]) % phase
    for first, second in pairs:
        if row_bits[first] != row_bits[second]:
            raise ValueError("a phase parameter differs within its physical unitary pair")
    t = len(boundary["border_C"])
    expected = np.zeros((t, t, root_vectors(phase).shape[1]), dtype=np.int64)
    expected[np.arange(t), np.arange(t), 0] = s
    if not np.array_equal(exact_column_gram(u, phase), expected):
        raise ValueError("the variable actual frame lost its exact characteristic-zero Gram")
    colors, dyadic = phase_coordinates(u, phase)
    if (not np.array_equal(colors, boundary["A"])
            or not np.array_equal(dyadic % 4, boundary["A"]*(1+2*row_bits[:, None]))):
        raise ValueError("the changed actual frame has the wrong low phase bits")
    return {**boundary, "source_order": order, "phase": phase,
            "U_exponents": u, "C": 1 ^ correction, "B": b, "W": w, "T": correction,
            "primitive_cross_ratio": primitive_frame_cross_ratio(u, phase),
            "template_groups": groups, "inactive_fillers": off[4*copies:].tolist(),
            "actual_row_pairs": pairs, "copy_second_bits": bits,
            "row_copy_bits": row_bits, "core_base_shift": 2*int(order == 2999),
            "boundary_next_bits": (dyadic // 2) % 2,
            "actual_variable_boundary_Gram_verified": True}


def representative_root_rows(problem):
    groups, bits = problem["template_groups"], problem["copy_second_bits"]
    chosen = []
    for value in (0, 1):
        indices = np.flatnonzero(bits == value)
        if len(indices) < 2:
            raise ValueError("the orbit check needs two copies of each phase type")
        chosen.extend(int(index) for index in indices[:2])
    fillers = problem["inactive_fillers"]
    if len(fillers) < 2:
        raise ValueError("the orbit check needs two inactive fillers")
    rows = np.array([i for index in chosen for i in groups[index]] + fillers[:2], dtype=np.int64)
    c, w = problem["C"][rows], problem["W"][rows]
    uc, ub = phase_coordinates(problem["U_exponents"], problem["phase"])
    dc, db = phase_coordinates(problem["D_exponents"], problem["phase"])
    colors = np.block([[np.zeros_like(c), c, uc[rows]], [c, np.zeros_like(c), uc[rows]],
                       [uc.T, uc.T, dc]])
    shift = problem["core_base_shift"]
    x = np.full_like(c, shift)
    y = (shift + w*(1+2*problem["row_copy_bits"][rows, None])) % 4
    exponents = np.block([[x, y, ub[rows] % 4], [y, x, ub[rows] % 4],
                          [ub.T % 4, ub.T % 4, db % 4]])
    s = len(problem["C"])
    return colors, exponents, rows.tolist() + (s+rows).tolist() + list(range(2*s, problem["source_order"]))


def full_parameter_problem(order, phase, copy_parameters=None):
    """Actual pair phases representing every odd class in the mod-two ring.

The root source is specified by X exponent alpha and Y exponents
alpha+W*lambda_i, with the actual changed U and fixed D. Its modular
Gram is established separately by the projector certificate.
"""
    p = variable_problem(order, phase)
    precision = phase // 6
    copies = len(p["template_groups"])
    if copy_parameters is None:
        copy_parameters = [2*(i % (precision//2))+1 for i in range(copies)]
    parameters = np.asarray(copy_parameters, dtype=np.int64)
    if (parameters.shape != (copies,) or np.any(parameters < 0)
            or np.any(parameters >= precision) or np.any(parameters % 2 != 1)):
        raise ValueError("each copy requires an odd phase representative modulo d/2")
    row_parameters = np.ones(len(p["C"]), dtype=np.int64)
    for parameter, group in zip(parameters, p["template_groups"], strict=True):
        row_parameters[group] = parameter
    previous = 1+2*p["row_copy_bits"]
    u = (p["U_exponents"] + 3*(row_parameters-previous)[:, None]*p["A"]) % phase
    uc, ub = phase_coordinates(u, phase)
    _, db = phase_coordinates(p["D_exponents"], phase)
    t = len(db)
    expected = np.zeros((t, t, root_vectors(phase).shape[1]), dtype=np.int64)
    expected[np.arange(t), np.arange(t), 0] = len(u)
    if not np.array_equal(exact_column_gram(u, phase), expected):
        raise ValueError("the full-parameter actual boundary lost its exact Gram")
    beta = int(db[0, 0] % precision)
    if (np.any(db % precision != beta)
            or not np.array_equal(uc, p["A"])
            or not np.array_equal(ub % precision, p["A"]*row_parameters[:, None])):
        raise ValueError("the actual phases do not have the required normalized form")
    alpha = (-beta) % precision
    if alpha % 2:
        raise ValueError("the core/border phase gauge has no integer half shift")
    return {**p, "U_exponents": u, "copy_parameters": parameters,
            "primitive_cross_ratio": primitive_frame_cross_ratio(u, phase),
            "row_parameters": row_parameters, "mod_two_precision": precision,
            "copy_second_bits": (parameters//2) % 2,
            "row_copy_bits": (row_parameters//2) % 2,
            "boundary_next_bits": (ub//2) % 2,
            "core_base_shift": alpha, "border_base_shift": beta,
            "full_parameter_source_Gram_claimed_here": False}


def orthogonal_parameter_problem(order, phase):
    """Align distinct phase classes with orthogonal real boundary patterns.

The resulting actual boundary Gram is checked exactly. The separate
projector and compatibility audits must verify the modular source and
its necessary characteristic-zero conditions.
"""
    if (order, phase) not in ((1499, 24), (2999, 48)):
        raise ValueError("the orthogonal-pattern construction covers1499/24 and2999/48")
    p = full_parameter_problem(order, phase)
    v = p["boundary_vector"]
    columns = np.flatnonzero(v)
    if order == 1499:
        _, signs, support, raw_pairs, _ = real_frame()
        raw_pairs = [pair for pair in raw_pairs if np.any(support[pair[0]])]
    else:
        _, signs, _, raw_pairs, _, _, _, _, _ = doubled_real_frame()
    buckets = {}
    for pair in raw_pairs:
        pattern = signs[pair[0], columns]
        pattern = pattern*pattern[0]
        if order == 2999 and int(np.prod(pattern)) != 1:
            continue
        buckets.setdefault(tuple(int(value) for value in pattern), []).append(list(pair))
    patterns = sorted(buckets)
    if len(patterns) != phase//12:
        raise ValueError("the required number of orthogonal pattern classes is missing")
    pattern_array = np.asarray(patterns, dtype=np.int64)
    if not np.array_equal(pattern_array @ pattern_array.T, len(columns)*np.eye(len(patterns), dtype=np.int64)):
        raise ValueError("the selected boundary sign classes are not mutually orthogonal")
    selected_pairs, group_parameters, pattern_reports = [], [], []
    for index, pattern in enumerate(patterns):
        if len(buckets[pattern]) < 93:
            raise ValueError("the chosen sign class does not supply93 physical row pairs")
        selected_pairs.extend(buckets[pattern][:93])
        group_parameters.extend([2*index+1]*31)
        pattern_reports.append({"projective_sign_pattern": list(pattern),
                                "available_physical_pairs": len(buckets[pattern]),
                                "selected_physical_pairs": 93, "template_copies": 31,
                                "odd_phase_parameter": 2*index+1})
    s = len(signs)
    on = np.asarray(selected_pairs, dtype=np.int64).ravel()
    active = np.zeros(s, dtype=np.int64)
    active[on] = 1
    off = np.flatnonzero(1-active)
    small = template()
    correction = np.eye(s, dtype=np.int64)
    b, w = np.zeros((s, s), dtype=np.int64), np.zeros((s, s), dtype=np.int64)
    row_parameters = np.ones(s, dtype=np.int64)
    groups = []
    for index, parameter in enumerate(group_parameters):
        group = np.concatenate((on[6*index:6*index+6], off[4*index:4*index+4]))
        where = np.ix_(group, group)
        correction[where], b[where], w[where] = small["T"], small["B"], small["W"]
        row_parameters[group] = parameter
        groups.append(group.tolist())
    a = active[:, None]*v[None, :]
    cube = next(e for e in (phase//3, 2*phase//3) if e % 3 == 1)
    u = ((signs == -1)*(phase//2)+(cube+3*row_parameters[:, None])*a) % phase
    expected = np.zeros((len(v), len(v), root_vectors(phase).shape[1]), dtype=np.int64)
    expected[np.arange(len(v)), np.arange(len(v)), 0] = s
    if not np.array_equal(exact_column_gram(u, phase), expected):
        raise ValueError("the orthogonal-pattern actual boundary lost its exact Gram")
    uc, ub = phase_coordinates(u, phase)
    if not np.array_equal(uc, a) or not np.array_equal(ub % (phase//6), a*row_parameters[:, None]):
        raise ValueError("the orthogonal-pattern root phases differ from the declared parameters")
    return {**p, "U_exponents": u, "A": a, "active": active,
            "C": 1 ^ correction, "B": b, "W": w, "T": correction,
            "template_groups": groups, "inactive_fillers": off[4*len(groups):].tolist(),
            "actual_row_pairs": selected_pairs,
            "copy_parameters": np.array(group_parameters, dtype=np.int64),
            "row_parameters": row_parameters,
            "copy_second_bits": (np.array(group_parameters)//2) % 2,
            "row_copy_bits": (row_parameters//2) % 2,
            "boundary_next_bits": (ub//2) % 2,
            "primitive_cross_ratio": primitive_frame_cross_ratio(u, phase),
            "orthogonal_pattern_classes": pattern_reports,
            "boundary_pattern_Gram": (pattern_array @ pattern_array.T).tolist()}


def full_parameter_row_data(problem, core_rows):
    """Actual columns for selected rows in both core halves and the border."""
    rows = np.asarray(core_rows, dtype=np.int64)
    c, w = problem["C"][rows], problem["W"][rows]
    uc, ub = phase_coordinates(problem["U_exponents"], problem["phase"])
    dc, db = phase_coordinates(problem["D_exponents"], problem["phase"])
    colors = np.block([[np.zeros_like(c), c, uc[rows]], [c, np.zeros_like(c), uc[rows]],
                       [uc.T, uc.T, dc]])
    precision, shift = problem["mod_two_precision"], problem["core_base_shift"]
    x = np.full_like(c, shift)
    y = shift+w*problem["row_parameters"][rows, None]
    exponents = np.block([[x, y, ub[rows]], [y, x, ub[rows]], [ub.T, ub.T, db]]) % precision
    s = len(problem["C"])
    full_rows = rows.tolist() + (s+rows).tolist() + list(range(2*s, problem["source_order"]))
    return colors, exponents, full_rows


def full_parameter_minor(problem, modulus):
    """Extract a star minor containing every odd relative exponent."""
    if modulus < 4 or problem["mod_two_precision"] % modulus:
        raise ValueError("the minor modulus must divide the mod-two precision")
    indices = []
    for value in range(1, modulus, 2):
        found = np.flatnonzero(problem["copy_parameters"] % modulus == value)
        if not len(found):
            raise ValueError("the chosen parameters omit a required odd relative phase")
        indices.append(problem["template_groups"][int(found[0])][0])
    border = int(np.flatnonzero(problem["boundary_vector"])[0])
    _, ub = phase_coordinates(problem["U_exponents"], problem["phase"])
    _, db = phase_coordinates(problem["D_exponents"], problem["phase"])
    size, shift = len(indices)+1, problem["core_base_shift"]
    raw = np.full((size, size), shift, dtype=np.int64)
    raw[:-1, -1] = ub[indices, border]
    raw[-1, :-1] = ub[indices, border]
    raw[-1, -1] = db[border, border]
    gauge = np.array([-shift//2]*(size-1)+[shift//2])
    normalized = (raw+gauge[:, None]+gauge[None, :]) % modulus
    expected = np.zeros((size, size), dtype=np.int64)
    expected[:-1, -1] = expected[-1, :-1] = np.arange(1, modulus, 2)
    if not np.array_equal(normalized, expected):
        raise ValueError("the actual root phases do not produce the claimed star minor")
    return {"source_order": problem["source_order"], "phase": problem["phase"],
            "relative_modulus": modulus,
            "actual_source_row_and_column_indices": indices+[2*len(problem["C"])+border],
            "raw_dyadic_exponents_modulo": (raw % modulus).tolist(),
            "same_row_and_column_gauge": gauge.tolist(),
            "normalized_star_minor": normalized.tolist()}


def adjacent_gauge_count(matrix, modulus=4):
    """Exhaust all normalized row/column gauges for a small minor.

The first row gauge can be fixed to zero. Each column then has exactly
two possible gauges, determined by its first-row output in {0,1}.
"""
    matrix = np.asarray(matrix, dtype=np.int64) % modulus
    count = checked = 0
    for other_rows in product(range(modulus), repeat=len(matrix)-1):
        row = np.array((0, *other_rows), dtype=np.int64)[:, None]
        for first_outputs in product((0, 1), repeat=matrix.shape[1]):
            column = np.array(first_outputs)[None, :] - matrix[:1]
            checked += 1
            if np.all(((matrix+row+column) % modulus) <= 1):
                count += 1
    return {"normalized_gauges_checked": checked, "admissible_gauges": count}


def verify_variable_residue(problem):
    colors, exponents, rows = representative_root_rows(problem)
    grams = polynomial_gram_from_roots(colors, exponents)
    if (not np.array_equal(grams[0], np.eye(len(rows), dtype=np.int64))
            or any(np.any(value) for value in grams[1:])):
        raise ValueError("the actual variable-boundary source fails modulo pi fourth")
    # The two selected active rows and a boundary coordinate produce a
    # small gauge obstruction, independent of any copy symmetry of gauges.
    bits = problem["copy_second_bits"]
    first = problem["template_groups"][int(np.flatnonzero(bits == 0)[0])][0]
    second = problem["template_groups"][int(np.flatnonzero(bits == 1)[0])][0]
    border = int(np.flatnonzero(problem["boundary_vector"])[0])
    s = len(problem["C"])
    actual_indices = [first, second, 2*s+border]
    row_positions = [rows.index(index) for index in actual_indices]
    minor = exponents[np.ix_(row_positions, actual_indices)]
    gauge = problem["core_base_shift"] // 2
    shifts = np.array([-gauge, -gauge, gauge])
    gauged = (minor+shifts[:, None]+shifts[None, :]) % 4
    expected = np.array([[0, 0, 1], [0, 0, 3], [1, 3, 0]])
    if not np.array_equal(gauged, expected):
        raise ValueError("the actual source does not contain the claimed gauge obstruction")
    obstruction = adjacent_gauge_count(gauged)
    if obstruction["admissible_gauges"]:
        raise ValueError("the minor still permits two adjacent phase channels")
    return {"all_actual_source_columns_used": colors.shape[1],
            "representative_source_rows": rows,
            "ordered_representative_row_pairs_checked": len(rows)**2,
            "all_four_Gram_coefficients_verified": True,
            "row_pair_orbits_covered_by_permutations_within_each_copy_type_and_fillers": True,
            "source_row_pairs_covered": problem["source_order"]*(problem["source_order"]+1)//2,
            "phase_rows_sha256": hashlib.sha256(colors.astype(np.uint8).tobytes()
                                                + exponents.astype(np.uint8).tobytes()).hexdigest(),
            "obstruction_actual_row_and_column_indices": actual_indices,
            "obstruction_dyadic_minor_mod4": minor.tolist(),
            "normalized_obstruction_minor_mod4": gauged.tolist(),
            "adjacent_channel_gauge_obstruction": obstruction,
            "characteristic_zero_source_constructed": False}


def certificate():
    reports = []
    for order, phase in ((1499, 96), (2999, 48), (2999, 96)):
        p = variable_problem(order, phase)
        reports.append({"source_order": order, "phase": phase,
                        "conditional_real_order": order*(phase//3),
                        "copy_second_bits": p["copy_second_bits"].tolist(),
                        "copy_type_counts": np.bincount(p["copy_second_bits"], minlength=2).tolist(),
                        "template_groups": p["template_groups"],
                        "actual_row_pairs": p["actual_row_pairs"],
                        "U_exponents": p["U_exponents"].tolist(),
                        "D_exponents": p["D_exponents"].tolist(),
                        "actual_variable_boundary_Gram_verified": True,
                        "full_modular_source_check": verify_variable_residue(p)})
    return {"schema": "variable-cubic-frame-bits-v1",
            "scope": "Exact variable Butson frames and noncollapsing pi-fourth root residues",
            "pair_parameter": "omega*eta^(1+2h), h constant on each physical row pair",
            "core_second_bit_rule": "X: source_order==2999; Y: X xor (W*copy_h)",
            "routes": reports, "no_general_or_higher_precision_feasibility_inferred": True,
            "new_real_Hadamard_order_constructed": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--output", type=Path)
    action.add_argument("--check", type=Path)
    args = parser.parse_args()
    report = certificate()
    if args.output:
        args.output.write_text(json.dumps(report, indent=2)+"\n")
    elif json.loads(args.check.read_text()) != report:
        raise ValueError("the variable-frame certificate differs from replay")
    print(json.dumps({"verified_routes": [[r["source_order"], r["phase"]] for r in report["routes"]],
                      "all_routes_escape_adjacent_channel_collapse_mod4": True,
                      "characteristic_zero_source_constructed": False}))


if __name__ == "__main__":
    main()

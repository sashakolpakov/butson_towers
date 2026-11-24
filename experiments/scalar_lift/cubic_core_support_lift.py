"""Constructive full-source first ramification lifts at1499 and2999.

The sources are unitary modulo pi squared, with actual root phase bits.
No characteristic-zero completion or real Hadamard matrix is asserted.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from math import gcd
from pathlib import Path

import numpy as np

from experiments.scalar_lift.fifth_source_assembly import exact_column_gram, root_vectors
from experiments.scalar_lift.mixed_cubic_boundaries import mixed_frames, mixed_frames_2999
from experiments.scalar_lift.mixed_cubic_ramification import CONJUGATE, derivative, product

C0 = np.array([
    [0,1,0,0,1,1,0,0,1,1],
    [0,1,1,1,0,0,1,0,1,0],
    [1,1,0,0,0,1,0,0,1,1],
    [1,0,0,1,1,0,0,1,1,0],
    [1,0,1,0,1,0,0,1,1,0],
    [0,0,1,1,0,1,1,0,1,0],
    [0,0,1,1,0,0,1,0,1,1],
    [1,1,1,1,1,1,1,1,0,1],
    [1,0,0,0,1,0,1,1,1,0],
    [0,1,0,0,0,1,0,1,1,1],
], dtype=np.int64)
B0 = np.array([
    [0,0,0,0,0,1,0,0,0,1],
    [0,0,0,1,0,0,1,0,0,0],
    [0,1,0,0,0,0,0,0,0,1],
    [0,0,0,0,1,0,0,1,0,0],
    [1,0,0,0,0,0,0,1,0,0],
    [0,0,1,0,0,0,1,0,0,0],
    [0,0,1,1,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0],
    [1,0,0,0,1,0,0,0,0,0],
    [0,1,0,0,0,1,0,0,0,0],
], dtype=np.int64)


def template():
    c, b = C0.copy(), B0.copy()
    active = np.array([1]*6+[0]*4, dtype=np.int64)
    eye = np.eye(10, dtype=np.int64)
    d = c @ b.T @ c % 2
    if (not np.array_equal(c @ c.T % 2, eye) or not np.array_equal(c @ active % 2, active)
            or not np.array_equal(c.T @ active % 2, active)
            or np.any(b*(1-c)) or np.any(c*d)
            or np.any(b.sum(axis=0) % 2) or np.any(b.sum(axis=1) % 2)
            or not np.array_equal(b @ active % 2, active)
            or not np.array_equal(b.T @ active % 2, active)):
        raise ValueError("the constructive ten-row support template failed")
    return {"C": c, "B": b, "a": active, "D": d, "W": b ^ d, "T": 1 ^ c}


def phase_coordinates(exponents, phase):
    cube = next(e for e in (phase//3, 2*phase//3) if e % 3 == 1)
    exponents = np.asarray(exponents, dtype=np.int64)
    color = exponents % 3
    dyadic = ((exponents-cube*color)//3) % (phase//3)
    if np.any((3*dyadic+cube*color-exponents) % phase):
        raise ValueError("the cubic/dyadic phase decomposition failed")
    return color, dyadic


def boundary_data(source_order, phase):
    if source_order == 1499 and phase in (24, 48, 96):
        old = mixed_frames(phase)
        u = np.array(old["U_exponents"], dtype=np.int64)
    elif source_order == 2999 and phase in (48, 96):
        old = mixed_frames_2999(phase)
        signs = np.array(old["T_signs"], dtype=np.int64)
        support = np.zeros_like(signs)
        #372 active pairs give744 active rows; the earlier chart used374.
        for first in range(372):
            second = first+748
            support[first] = support[second] = signs[first] != signs[second]
        eta = old["scalar_pair_action"]["eta_exponent"]
        u = ((signs == -1)*(phase//2)+eta*support) % phase
    else:
        raise ValueError("the constructed routes are1499/24,48,96 and2999/48,96")
    d = np.array(old["D_exponents"], dtype=np.int64)
    a, dyadic = phase_coordinates(u, phase)
    c, border_dyadic = phase_coordinates(d, phase)
    if np.any(a > 1) or np.any(c > 1):
        raise ValueError("a boundary left the admitted cubic alphabet")
    active = np.any(a, axis=1).astype(np.int64)
    vector = a[np.flatnonzero(active)[0]]
    if (not np.array_equal(a, active[:, None]*vector[None, :])
            or active.sum() % 2 or vector.sum() % 2
            or not np.array_equal(dyadic % 2, a) or np.any((dyadic//2) % 2)
            or np.any(border_dyadic % 2)):
        raise ValueError("the fixed boundary phase bits are not the intended ones")
    expected = np.zeros((len(vector), len(vector), root_vectors(phase).shape[1]), dtype=np.int64)
    expected[np.arange(len(vector)), np.arange(len(vector)), 0] = len(u)
    if not np.array_equal(exact_column_gram(u, phase), expected):
        raise ValueError("the modified actual boundary lost its characteristic-zero Gram")
    ratios = (u-u[:, :1]-u[:1]+u[0, 0]) % phase
    position = next(((i, j) for i, row in enumerate(ratios) for j, value in enumerate(row)
                     if gcd(int(value), phase) == 1), None)
    if position is None:
        raise ValueError("the modified boundary lost its primitive phase cross ratio")
    i, j = position
    return {"U_exponents": u, "D_exponents": d, "A": a, "border_C": c,
            "active": active, "boundary_vector": vector,
            "primitive_cross_ratio": {"rows": [0, i], "columns": [0, j],
                                      "exponent": int(ratios[i, j]), "order": phase},
            "boundary_next_bits": (dyadic//2) % 2,
            "border_next_bits": (border_dyadic//2) % 2}


def lifted_core(active):
    active = np.asarray(active, dtype=np.int64)
    size, weight = len(active), int(active.sum())
    if weight % 6 or size-weight < 4*(weight//6):
        raise ValueError("the active/inactive counts cannot fit the ten-row template")
    small = template()
    on, off = np.flatnonzero(active), np.flatnonzero(1-active)
    copies = weight//6
    t = np.eye(size, dtype=np.int64)
    b, w = np.zeros((size, size), dtype=np.int64), np.zeros((size, size), dtype=np.int64)
    groups = []
    for index in range(copies):
        coordinates = np.concatenate((on[6*index:6*index+6], off[4*index:4*index+4]))
        where = np.ix_(coordinates, coordinates)
        t[where], b[where], w[where] = small["T"], small["B"], small["W"]
        groups.append(coordinates.tolist())
    c = 1 ^ t
    if np.any(b*(1-c)) or np.any(c*(b ^ w)):
        raise ValueError("the replicated cubic supports overlap incorrectly")
    return {"C": c, "B": b, "W": w, "T": t,
            "template_groups": groups, "inactive_fillers": off[4*copies:].tolist()}


def large_problem(source_order, phase):
    boundary = boundary_data(source_order, phase)
    core = lifted_core(boundary["active"])
    return {"source_order": source_order, "phase": phase, **boundary, **core}


def source_layers(problem):
    c, w, a, border = (problem[key] for key in ("C", "W", "A", "border_C"))
    ones, zeros = np.ones_like(c), np.zeros_like(c)
    y0, u0, d0 = ones ^ (3*c), np.ones_like(a) ^ (3*a), np.ones_like(border) ^ (3*border)
    y1, u1 = y0*w, 2*a
    h0 = np.block([[ones, y0, u0], [y0, ones, u0], [u0.T, u0.T, d0]])
    h1 = np.block([[zeros, y1, u1], [y1, zeros, u1], [u1.T, u1.T, np.zeros_like(border)]])
    return h0, h1


def root_candidate_exponents(problem):
    """Actual admitted roots whose Gram is correct only modulo pi squared."""
    phase = problem["phase"]
    cube = next(e for e in (phase//3, 2*phase//3) if e % 3 == 1)
    x = np.zeros_like(problem["C"])
    y = (cube*problem["C"]+3*problem["W"]) % phase
    u, d = problem["U_exponents"], problem["D_exponents"]
    candidate = np.block([[x, y, u], [y, x, u], [u.T, u.T, d]])
    if np.any(candidate % 3 == 2):
        raise ValueError("the constructed root candidate left the admitted alphabet")
    return candidate


def direct_small_control():
    small = template()
    a = small["a"][:, None]*np.array([[0, 1, 1]])
    problem = {**small, "A": a, "border_C": np.eye(3, dtype=np.int64)}
    h0, h1 = source_layers(problem)
    if not np.array_equal(product(h0, CONJUGATE[h0].T), np.eye(23, dtype=np.int64)):
        raise ValueError("the independent order23 F4 Gram failed")
    if np.any(derivative(h0, h1)):
        raise ValueError("the independent full order23 first Gram derivative failed")
    return {"source_size": 23, "complete_F4_Gram_and_first_derivative_verified": True,
            "H0": h0.tolist(), "H1": h1.tolist()}


def packed_rows(matrix, bit):
    values = (np.asarray(matrix, dtype=np.uint8) >> bit) & 1
    return [int.from_bytes(row.tobytes(), "little") for row in np.packbits(values, axis=1, bitorder="little")]


def full_bitset_check(h0, h1):
    """Check every row pair of the full source, without cubic multiplication."""
    n = len(h0)
    zero0, zero1 = packed_rows(h0, 0), packed_rows(h0, 1)
    first0, first1 = packed_rows(h1, 0), packed_rows(h1, 1)
    zero_bar = [a ^ b for a, b in zip(zero0, zero1, strict=True)]
    first_bar = [a ^ b for a, b in zip(first0, first1, strict=True)]
    # h0=J+omega²P; zero1 is exactly the binary support P.
    if any(row.bit_count() % 2 != 1 for row in zero1):
        raise ValueError("the full cubic support has an even row")
    count = 0
    for i in range(n):
        a, b, x, y = zero0[i], zero1[i], first0[i], first1[i]
        for j in range(i+1):
            if ((b & zero1[j]).bit_count() & 1) != int(i == j):
                raise ValueError("the full binary support Gram is not I")
            constant = ((x & zero_bar[j]).bit_count() ^ (y & zero1[j]).bit_count()
                        ^ (a & first_bar[j]).bit_count() ^ (b & first1[j]).bit_count()) & 1
            omega = ((x & zero1[j]).bit_count() ^ (y & zero0[j]).bit_count()
                     ^ (a & first1[j]).bit_count() ^ (b & first0[j]).bit_count()) & 1
            if constant or omega:
                raise ValueError("the full source first Gram derivative is nonzero")
            count += 1
    return {"source_order": n, "all_row_pairs_including_diagonal_checked": count,
            "full_binary_support_Gram_equals_I": True,
            "full_F4_Gram_equals_I_from_exact_support_identity": True,
            "full_first_Gram_derivative_zero": True,
            "H0_uint8_sha256": hashlib.sha256(h0.astype(np.uint8).tobytes()).hexdigest(),
            "H1_uint8_sha256": hashlib.sha256(h1.astype(np.uint8).tobytes()).hexdigest()}


def build_report():
    start = time.monotonic()
    small = template()
    reports, seen = [], {}
    for order, phase in ((1499, 96), (2999, 48), (2999, 96)):
        problem = large_problem(order, phase)
        h0, h1 = source_layers(problem)
        digest = hashlib.sha256(h0.astype(np.uint8).tobytes()+h1.astype(np.uint8).tobytes()).hexdigest()
        reused = digest in seen
        check = seen[digest] if reused else full_bitset_check(h0, h1)
        seen[digest] = check
        reports.append({"source_order": order, "phase": phase,
                        "conditional_real_order": order*(phase//3),
                        "core_order": len(problem["C"]), "active_rows": int(problem["active"].sum()),
                        "template_copies": len(problem["template_groups"]),
                        "template_groups": problem["template_groups"],
                        "inactive_fillers": problem["inactive_fillers"],
                        "U_exponents": problem["U_exponents"].tolist(),
                        "D_exponents": problem["D_exponents"].tolist(),
                        "boundary_next_bits": problem["boundary_next_bits"].tolist(),
                        "border_next_bits": problem["border_next_bits"].tolist(),
                        "actual_boundary_Gram_reverified_exactly": True,
                        "boundary_primitive_cross_ratio": problem["primitive_cross_ratio"],
                        "cubic_support_connected_through_all_cross_group_edges": True,
                        "full_source_check": check, "same_verified_residue_reused": reused,
                        "root_phases_unitary_mod_pi_squared": True,
                        "unitary_mod_pi_fourth_claimed": False,
                        "characteristic_zero_source_constructed": False})
    return {"schema": "cubic-core-support-lift-v1", "template": {key: value.tolist() for key, value in small.items()},
            "independent_direct_control": direct_small_control(), "routes": reports,
            "large_source_residue_lift_constructed": True,
            "new_real_Hadamard_constructed": False, "wall_seconds": round(time.monotonic()-start, 6)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = build_report()
    args.output.write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps({"routes": len(report["routes"]), "wall_seconds": report["wall_seconds"]}))


if __name__ == "__main__":
    main()

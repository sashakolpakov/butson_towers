"""Exact modular audit of every Q=F diag(epsilon) F^T/d basis at once.

These rational orthogonal bases can have the required odd denominators.
The audit uses finite-field kernels, not enumeration or sampled signs.
"""
from __future__ import annotations

import argparse
import json
from itertools import product
from pathlib import Path

import numpy as np
from flint import nmod_mat
from sympy import isprime

from experiments.scalar_lift.alpoge_mixed_basis_probe import binary_rank, canonical_fraction
from experiments.scalar_lift.broad_phase_basis_probe import FRAME_SPECS, array_hash, load_frame
from experiments.scalar_lift.q12_basis_probe import block_statistics, twice_phase_action


def native(array, prime):
    return nmod_mat(np.asarray(array).tolist(), prime)


def field_array(matrix, rows=None):
    rows = matrix.nrows() if rows is None else rows
    columns = matrix.ncols()
    return np.fromiter((int(matrix[i, j]) for i in range(rows) for j in range(columns)),
                       dtype=np.uint16, count=rows*columns).reshape(rows, columns)


def support_graph(rref, pivots):
    """All diagonal +/-1 stabilizers of ker[I A], encoded by components."""
    rows, columns = rref.shape
    pivots = list(map(int, pivots))
    free = sorted(set(range(columns))-set(pivots))
    if len(pivots) != rows or not np.array_equal(rref[:, pivots], np.eye(rows, dtype=rref.dtype)):
        raise ArithmeticError("the displayed RREF pivots do not form an identity")
    parent = list(range(columns))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    spanning, edge_count = [], 0
    for i, pivot in enumerate(pivots):
        occupied = np.flatnonzero(rref[i, free])
        edge_count += len(occupied)
        for j in occupied:
            column = free[int(j)]
            a, b = find(pivot), find(column)
            if a != b:
                parent[b] = a
                spanning.append({"pivot_coordinate": pivot, "free_coordinate": column,
                                 "rref_row": i, "coefficient": int(rref[i, column])})
    roots = [find(i) for i in range(columns)]
    labels = {root: i for i, root in enumerate(sorted(set(roots)))}
    components = [labels[root] for root in roots]
    return {"pivot_coordinates": pivots, "free_coordinates": free,
            "support_edges": int(edge_count), "component_count": len(labels),
            "component_labels": components, "spanning_forest_nonzero_edges": spanning,
            "all_diagonal_sign_stabilizers": "one independent sign on each connected component",
            "stabilizer_count_power_of_two": len(labels),
            "constant_signs_are_only_stabilizers": len(labels) == 1}


def column_factorization(frame, prime):
    """Verify F=V*C modulo p, with V independent columns and C=[I A]."""
    matrix = native(frame, prime)
    reduced, rank = matrix.rref()
    top = field_array(reduced, rank)
    del reduced
    pivots = [int(np.flatnonzero(row)[0]) for row in top]
    if len(set(pivots)) != rank:
        raise ArithmeticError("the native RREF has inconsistent pivots")
    columns = np.asarray(frame)[:, pivots]
    V = native(columns, prime)
    C = native(top, prime)
    if matrix != V*C:
        raise ArithmeticError("the complete exact column factorization failed")
    del C, matrix
    return columns, V, top, pivots


def kernel_sign_points(reduced, pivots, prime, maximum_nullity=12):
    """Enumerate only the 2^nullity free-coordinate sign assignments."""
    width = reduced.shape[1]
    free = sorted(set(range(width))-set(pivots))
    if len(free) > maximum_nullity:
        return free, None
    points = []
    for signs in product((-1, 1), repeat=len(free)):
        values = np.zeros(width, dtype=np.int64)
        values[free] = signs
        values[pivots] = (-reduced[:, free].astype(np.int64)@values[free]) % prime
        values %= prime
        if np.all((values == 1) | (values == prime-1)):
            points.append((list(signs), np.where(values == 1, 1, -1).astype(np.int64)))
    return free, points


def sign_point_rank_gate(frame, exponent, points):
    numerator, denominator = canonical_fraction(twice_phase_action(frame, exponent), 2)
    if denominator != 1:
        raise ArithmeticError("the quaternion baseline should be integer on a sign frame")
    rank = binary_rank(numerator-1)
    records = []
    for signs, epsilon in points:
        global_sign = 1 if np.count_nonzero(epsilon < 0) <= len(frame)//2 else -1
        minority = np.flatnonzero(epsilon != global_sign)
        records.append({"free_coordinate_signs": signs, "global_sign": global_sign,
                        "minority_coordinates": minority.tolist(),
                        "sign_vector_sha256": array_hash(epsilon),
                        "rank_Q_minus_global_sign_I": len(minority),
                        "baseline_binary_rank": rank,
                        "rank_gate_excludes": rank > 2*len(minority)})
    return rank, records


def small_linear_gate(frame, exponent):
    """Every F diag(eps) K equation is linear in eps; use three F rows."""
    width, prime = len(frame), len(frame)//4
    F = native(frame, prime)
    K = F.transpose()*native(twice_phase_action(frame, exponent), prime)
    values = field_array(K).astype(np.int64)
    coefficients = np.vstack([(values.T*frame[i]) % prime for i in (0, 1, 2)])
    if np.any(coefficients.sum(axis=1) % prime):
        raise ArithmeticError("constant signs do not satisfy the necessary equations")
    identity = array_hash(coefficients)
    reduced, rank = native(coefficients, prime).rref()
    del coefficients, values, F, K
    top = field_array(reduced, rank)
    del reduced
    pivots = [int(np.flatnonzero(row)[0]) for row in top]
    free, points = kernel_sign_points(top, pivots, prime)
    rank_gate, survivors = (None, []) if points is None else sign_point_rank_gate(frame, exponent, points)
    return {"selected_frame_rows": [0, 1, 2], "modulus": prime,
            "constraint_matrix_sha256": identity, "constraint_rank": rank,
            "constraint_nullity": width-rank, "free_coordinates": free,
            "rref_sha256": array_hash(top),
            "free_sign_assignments_checked": 0 if points is None else 2**len(free),
            "all_kernel_sign_points_exhausted": points is not None,
            "sign_points": survivors, "baseline_binary_rank": rank_gate,
            "rank_gate": "rank_F2(rho0(z)F-J) <= 2 rank_Q(Q-global_sign*I)",
            "all_sign_points_excluded": points is not None and all(p["rank_gate_excludes"] for p in survivors),
            "scope": "Necessary three-row modular equations, followed by exhaustive sign points in their small kernel and the exact low-rank gate"}


def audit_frame(identifier):
    frame, proof = load_frame(identifier)
    width, prime = len(frame), len(frame)//4
    if width % 4 or not isprime(prime) or prime == 2:
        raise ValueError("the exact odd-prime frame model is required")
    baselines = {}
    for exponent in (1, 2):
        numerator, denominator = canonical_fraction(twice_phase_action(frame, exponent), 2)
        baselines[str(exponent)] = block_statistics(numerator, denominator)
    columns, V, reduced, pivots = column_factorization(frame, prime)
    rank = len(pivots)
    graph = support_graph(reduced, pivots)
    roots = {}
    for exponent in (1, 2):
        transformed = twice_phase_action(columns, exponent)
        T_V = native(transformed, prime)
        compressed = V.transpose()*T_V
        del transformed, T_V
        compressed_rank = compressed.rank()
        compressed_hash = array_hash(field_array(compressed))
        del compressed
        full_image = rank == width//2 and compressed_rank == width//2
        connected = graph["constant_signs_are_only_stabilizers"]
        closed = (compressed_rank == 0 or
                  (full_image and connected and not baselines[str(exponent)]["all_entries_signs"]))
        linear = None
        if not closed:
            linear = small_linear_gate(frame, exponent)
            closed = linear["all_sign_points_excluded"]
        roots[str(exponent)] = {
            "root_order": 12//int(np.gcd(12, exponent)),
            "compressed_K_rank": compressed_rank, "full_K_rank": compressed_rank,
            "compressed_K_sha256": compressed_hash,
            "rank_identity": "F=V*C, K=F^T*T*F=C^T*(V^T*T*V)*C; C has full row rank",
            "column_K_equals_kernel_F": full_image,
            "necessary_modular_condition": "F diag(epsilon) K=0 mod r",
            "diagonal_signs_must_stabilize_kernel_F": full_image,
            "rank_zero_dyadic_right_action_gate": compressed_rank == 0,
            "small_linear_gate": linear,
            "family_closed": closed,
            "status": ("EXCLUDED_all_2^d_sign_choices" if closed else
                       "OPEN_rank_deficient_image" if not full_image else "OPEN_multiple_support_components"),
            "baseline_statistics": baselines[str(exponent)]}
    related = {}
    for exponent, base in ((4, 2), (5, 1)):
        parent = roots[str(base)]
        inherited = parent["column_K_equals_kernel_F"] or parent["rank_zero_dyadic_right_action_gate"]
        if inherited:
            related[str(exponent)] = {"family_closed": parent["family_closed"],
                "exact_transfer": "K4=K2 mod r" if exponent == 4 else "K5=-K1^T has full kernel image or zero image",
                "baseline_exclusion": "the primitive odd-part dyadic baseline cannot map an order4r Hadamard to signs"}
        else:
            linear = small_linear_gate(frame, exponent)
            related[str(exponent)] = {"family_closed": linear["all_sign_points_excluded"],
                                      "small_linear_gate": linear}
    return {"frame_id": identifier, "width": width, "odd_prime": prime,
            "published_frame_proof": proof, "frame_sha256": array_hash(frame),
            "modular_F_rank": rank, "required_half_rank": width//2,
            "complete_column_factorization_verified": True,
            "rref_nonzero_rows_sha256": array_hash(reduced), "kernel_support_graph": graph,
            "root_cases": roots,
            "related_root_cases": related,
            "negative_roots": "7,8,10,11 are negative scalar actions; integrality and sign conditions unchanged",
            "all_non_two_power_root_families_closed": all(row["family_closed"] for row in (*roots.values(), *related.values())),
            "no_random_epsilon_or_basis_search_performed": True}


def exhaustive_small_control():
    """Check the stabilizer criterion against all16 signs in a nontrivial code."""
    prime = 3
    frame = np.array([[1, 0, 1, 1], [0, 1, 1, -1],
                      [1, 1, -1, 0], [1, -1, 0, -1]], dtype=np.int64)
    columns, V, reduced, pivots = column_factorization(frame, prime)
    graph = support_graph(reduced, pivots)
    free = graph["free_coordinates"]
    kernel = np.empty((4, 2), dtype=np.int64)
    kernel[pivots] = -reduced[:, free].astype(np.int64)
    kernel[free] = np.eye(2, dtype=np.int64)
    T = np.diag([1, 1, -1, -1])
    K = frame.T@T@frame
    if not np.all((frame@kernel) % prime == 0) or native(K, prime).rank() != 2:
        raise ArithmeticError("the small full-image control is invalid")
    stabilizers = []
    for signs in product((-1, 1), repeat=4):
        E = np.diag(signs)
        direct = not np.any((frame@E@kernel) % prime)
        predicted = all(signs[i] == signs[j] for i in range(4) for j in range(4)
                        if graph["component_labels"][i] == graph["component_labels"][j])
        image_test = not np.any((frame@E@K) % prime)
        if direct != predicted or direct != image_test:
            raise ArithmeticError("the complete small diagonal-stabilizer enumeration disagrees")
        if direct:
            stabilizers.append(list(signs))
    if len(stabilizers) != 2:
        raise ArithmeticError("the connected small support graph should have only constant signs")
    return {"modulus": 3, "dimension": 4, "sign_choices_checked": 16,
            "stabilizers": stabilizers, "compressed_rank": (V.transpose()*native(T@columns, prime)).rank(),
            "all_three_tests_agree": True}


def report(progress=False):
    controls = exhaustive_small_control()
    results = []
    for identifier, _width in FRAME_SPECS:
        result = audit_frame(identifier)
        results.append(result)
        if progress:
            print(json.dumps({"frame": identifier, "rank_F": result["modular_F_rank"],
                              "rank_K": {e: r["full_K_rank"] for e, r in result["root_cases"].items()},
                              "support_components": result["kernel_support_graph"]["component_count"],
                              "closed": result["all_non_two_power_root_families_closed"]}), flush=True)
    return {"schema": "odd-denominator-frame-diagonal-sign-kernel-audit-v1",
            "basis_family": "Q=F diag(epsilon) F^T/d for ALL epsilon in {+1,-1}^d",
            "orthogonality": "Q Q^T=I from F F^T=F^T F=dI",
            "new_root_block": "G=F diag(epsilon) K diag(epsilon)/(2d), K=F^T*(2rho0(z))*F",
            "integer_entry_gate": "F diag(epsilon) K=0 mod r, where d=4r",
            "kernel_gate": "If rankF=rankK=d/2 and FF^T=0 mod r, imK=kerF; diagonal signs must stabilize that kernel",
            "support_graph_gate": "For RREF(F)=[I A], epsilon_pivot_i=epsilon_free_j whenever Aij!=0",
            "connected_case": "Only epsilon=+1 or-1 globally; then Q=+I or-I and the failed baseline returns",
            "rank_zero_case": "K divisible by r makes C=K/(2d) dyadic; G=F diag(eps) C diag(eps) is excluded by the transposed dyadic primitive-root gate",
            "deficient_case": "Stack K^T diag(F_i) for i0,1,2. Enumerate all2^nullity free-coordinate sign assignments if nullity<=12; prune their exact low-rank basis cases",
            "compressed_native_computation": "Verify F=V*C exactly, then rankK=rank(V^T*T*V) using FLINT nmod_mat",
            "small_exhaustive_control": controls, "frames": results,
            "frame_count": len(results),
            "families_closed": sum(r["all_non_two_power_root_families_closed"] for r in results),
            "open_frames": [r["frame_id"] for r in results if not r["all_non_two_power_root_families_closed"]],
            "pending_declared_audits": 0, "new_Hadamard_found": False,
            "scope": "These entire frame-diagonal sign-conjugacy families for all four non-two-power root representatives1,2,4,5; not all rational bases with the needed odd denominators"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    choice = parser.add_mutually_exclusive_group(required=True)
    choice.add_argument("--output", type=Path)
    choice.add_argument("--check", type=Path)
    parser.add_argument("--progress", action="store_true")
    args = parser.parse_args()
    result = report(args.progress)
    if args.output:
        args.output.write_text(json.dumps(result, separators=(",", ":"))+"\n")
    elif json.loads(args.check.read_text()) != result:
        raise ValueError("the odd-denominator exact kernel audit differs from replay")
    print(json.dumps({k: result[k] for k in ("frame_count", "families_closed", "open_frames", "pending_declared_audits")}))


if __name__ == "__main__":
    main()

"""Exact controls for using a real GS design as a scalar-restriction chart.

This replays an existing order-twelve control. It runs no solver, creates no
source-search tree edge, and does not claim a new Hadamard order.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from sympy import I, Matrix

from experiments.scalar_lift.ambient_field_candidates import load_ambient_source
from experiments.scalar_lift.joint_gs_223 import goethals_seidel_array


def _hash(matrix):
    return hashlib.sha256(np.asarray(matrix, dtype="<i8").tobytes()).hexdigest()


def report():
    # Existing literal control from test_goethals_seidel_gram_and_design_to_pair_identities.
    a = np.array([[0, 1, 1], [-1, 0, 1], [-1, -1, 0], [-1, 1, -1]], dtype=np.int64)
    b = np.eye(4, 3, dtype=np.int64)
    U, V = goethals_seidel_array(a), goethals_seidel_array(b)
    identity = np.eye(12, dtype=np.int64)
    if (not np.array_equal(U @ U.T, 9*identity)
            or not np.array_equal(V @ V.T, 3*identity)
            or np.any(U @ V.T+V @ U.T)
            or not np.all(np.abs(U)+np.abs(V) == 1)):
        raise ArithmeticError("the known formal real GS control failed")
    F, G = U-V, U+V
    if (not np.all(np.abs(F) == 1) or not np.all(np.abs(G) == 1)
            or not np.array_equal(F @ F.T, 12*identity)
            or not np.array_equal(G @ G.T, 12*identity)
            or not np.array_equal(G @ F.T+F @ G.T, 12*identity)):
        raise ArithmeticError("the exact sign-pair control failed")
    R = Matrix((G @ F.T).tolist())/12
    if Matrix.eye(12) != R*R.T or R*R-R+Matrix.eye(12) != Matrix.zeros(12):
        raise ArithmeticError("the recovered sextic action failed")
    blocks = (F, G, G-F, -F, -G, F-G)
    controls = []
    for source_id in ("unreal3-q6", "catalogue-p7-q6", "unreal13-q6"):
        source = load_ambient_source(source_id)
        restricted = np.block([[blocks[source.entry(i, j)] for j in range(source.order)]
                               for i in range(source.order)])
        order = 12*source.order
        if not np.all(np.abs(restricted) == 1) or not np.array_equal(restricted @ restricted.T, order*np.eye(order, dtype=np.int64)):
            raise ArithmeticError("the complete known scalar-restriction control failed")
        controls.append({"source_id": source_id, "source_hash": source.source_hash,
                         "source_Gram_replayed": source.proof["full_Gram_verified"],
                         "matrix_order": order, "full_sign_Gram_checked": True,
                         "matrix_sha256": _hash(restricted), "new_order_claimed": False})
    # The ordinary real GS layout is not automatically a complex Hadamard layout.
    complex_array = Matrix([[1, I, 1, 1], [-I, 1, -1, 1],
                            [-1, 1, 1, -I], [-1, -1, I, 1]])
    gram = complex_array*complex_array.H
    if any(gram[j, j] != 4 for j in range(4)) or gram[0, 1] != 2*I:
        raise ArithmeticError("the exact complex negative control changed")
    return {"schema": "real-GS-scalar-chart-audit-v1",
            "existing_control_source": "tests/test_joint_gs_223.py:test_goethals_seidel_gram_and_design_to_pair_identities",
            "a_sequences": a.tolist(), "b_sequences": b.tolist(),
            "formal_GS_design": "D(x,y)D(x,y)^T=(9x^2+3y^2)I12 for commuting involution-fixed x,y",
            "all_three_real_polynomial_Gram_coefficients_checked": True,
            "complementary_supports_checked": True,
            "F_sha256": _hash(F), "G_sha256": _hash(G),
            "recovered_R_definition": "R=G F^T/12",
            "R_orthogonal_and_R_squared_minus_R_plus_I_checked": True,
            "complete_scalar_restriction_controls": controls,
            "naive_complex_GS_counterexample": {
                "four_length_one_inputs": ["1", "i", "1", "1"],
                "entries_are_fourth_roots": True, "all_diagonal_Gram_entries": 4,
                "nonzero_offdiagonal_Gram_witness": {"row": 0, "column": 1, "value": "2i"}},
            "solver_or_source_search_run": False, "new_tree_edge_added": False,
            "new_Hadamard_order_claimed": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    choices = parser.add_mutually_exclusive_group(required=True)
    choices.add_argument("--output", type=Path)
    choices.add_argument("--check", type=Path)
    args = parser.parse_args()
    result = report()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+"\n")
    elif json.loads(args.check.read_text()) != result:
        raise ValueError("the GS scalar-chart audit differs from exact replay")
    print(json.dumps({"known_control_orders": [c["matrix_order"] for c in result["complete_scalar_restriction_controls"]],
                      "complex_negative_control_checked": True, "new_Hadamard": False}))


if __name__ == "__main__":
    main()

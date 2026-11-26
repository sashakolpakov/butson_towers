"""Exact T631 endpoints from published symmetric amicable sign quadruples.

The first circulant in a near-Williamson input is multiplied by reversal.
It then commutes with the three symmetric circulants.  These matrices
generate a commutative reduced rational algebra, represented on Q^w with
the identity metric.  They need not be elements of the cyclic group algebra.
"""
from __future__ import annotations

import argparse
from functools import cache
import hashlib
import json
from pathlib import Path

import numpy as np

from experiments.scalar_lift.williamson_scalar_endpoints import (
    ROOT, ScalarEndpoint, array_hash, require, row_controls,
    verified_t_design,
)

INPUTS = Path(__file__).with_name("near_williamson_inputs.json")
STATUSES = Path(__file__).with_name("near_williamson_order_status.json")
ALGEBRAS = ROOT / "certificates/scalar-lift/near-williamson-scalar-algebras.json"


def literal_matrices(record):
    """Reconstruct the exact four matrices, including the recorded reversal."""
    w = record["width"]
    words = record["first_rows_pm"]
    require(type(w) is int and w > 0 and len(words) == 4,
            "invalid source width or number of rows")
    require(all(isinstance(s, str) and len(s) == w and set(s) <= {"+", "-"}
                for s in words), "invalid literal sign row")
    require(type(record["first_matrix_reversed"]) is bool,
            "the first-matrix reversal must be specified")
    rows = np.array([[1 if c == "+" else -1 for c in s] for s in words],
                    dtype=np.int64)
    index = (np.arange(w)[None, :] - np.arange(w)[:, None]) % w
    matrices = rows[:, index].copy()
    if record["first_matrix_reversed"]:
        matrices[0] = matrices[0, :, ::-1]
    return rows, matrices


def verify_quadruple(record):
    """Check every entry of the norm, symmetry and all pairwise identities."""
    rows, matrices = literal_matrices(record)
    w = record["width"]
    require(np.all(abs(matrices) == 1), "a source entry is not a sign")
    require(np.array_equal(matrices, matrices.transpose(0, 2, 1)),
            "a source matrix is not symmetric")
    for i in range(4):
        for j in range(i):
            require(np.array_equal(matrices[i] @ matrices[j],
                                   matrices[j] @ matrices[i]),
                    "source matrices do not commute")
            require(np.array_equal(matrices[i] @ matrices[j].T,
                                   matrices[j] @ matrices[i].T),
                    "source matrices are not amicable")
    gram = sum(a @ a.T for a in matrices)
    require(np.array_equal(gram, 4*w*np.eye(w, dtype=np.int64)),
            "complete summed source Gram failed")
    expected = record.get("verification", {})
    require(array_hash(matrices) == expected["WT_matrices_int8_row_major_sha256"],
            "reconstructed matrices differ from the acquired literal data")
    proof = {
        "width": w,
        "literal_first_rows_pm": record["first_rows_pm"],
        "first_matrix_reversed": record["first_matrix_reversed"],
        "circulant_entry": "C_j[r,s]=row_j[(s-r) mod w]",
        "first_matrix_entry": (
            "A_0[r,s]=row_0[(w-1-s-r) mod w]"
            if record["first_matrix_reversed"] else
            "A_0[r,s]=row_0[(s-r) mod w]"
        ),
        "all_four_sign_and_symmetric": True,
        "all_six_complete_commutation_and_amicability_identities_checked": True,
        "all_summed_integer_Gram_entries_checked": w*w,
        "source_norm": f"sum_j A_j A_j^T={4*w}I_{w}",
        "literal_rows_sha256": array_hash(rows),
        "matrices_sha256": array_hash(matrices),
        "primary_literal_source": record["source"],
    }
    matrices = matrices.astype(np.int8)
    matrices.setflags(write=False)
    return matrices, proof


@cache
def verified_source(w):
    require(type(w) is int, "the width must be an integer")
    data = json.loads(INPUTS.read_text())
    record = next((r for r in data["entries"] if r["width"] == w), None)
    require(record is not None, "no retained compatible input at this width")
    return verify_quadruple(record)


def endpoint(w, p=631):
    """Return a literal entry/row recipe; never allocate the dense endpoint."""
    x, _ = verified_t_design(p)
    a, _ = verified_source(w)
    return ScalarEndpoint(x, a)


def certificate(widths=None):
    data = json.loads(INPUTS.read_text())
    if widths is None:
        widths = [r["width"] for r in data["entries"]]
    statuses = json.loads(STATUSES.read_text()) if STATUSES.exists() else {}
    statuses = statuses.get("entries", statuses)
    x, t_proof = verified_t_design(631)
    outputs = []
    for w in widths:
        a, source = verified_source(w)
        literal = ScalarEndpoint(x, a)
        used = [{"source_index": j, "sign": sign,
                 "operator_int8_sha256": array_hash(sign*a[j])}
                for j in range(4)
                for sign in sorted(set(map(int, np.unique(x[j]))) - {0})]
        outputs.append({
            "width": w,
            "output_order": literal.order,
            "source": source,
            "scalar_algebra": "Q[A_0,A_1,A_2,A_3] inside M_w(Q)",
            "algebra_is_commutative_and_reduced": True,
            "reason_reduced": "Its generators are commuting real symmetric matrices, hence simultaneously diagonalizable.",
            "scalar_action": "the explicitly reconstructed matrices on V=Q^w",
            "conjugation": "identity on the scalar algebra; its represented matrices are symmetric",
            "module_metric": f"I_{w} in the standard basis of V",
            "weighted_source_order": 2524,
            "weighted_source_Gram_scale": literal.order,
            "weighted_source_is_Butson": False,
            "actual_source_alphabet": "the eight recorded operators +A_j,-A_j, with coincident operators identified",
            "actual_used_operators": used,
            "distinct_used_operator_count": len({r["operator_int8_sha256"] for r in used}),
            "terminal_frame": f"I_{w}",
            "terminal_frame_Gram_scale": 1,
            "rational_order_after_scalar_action": literal.order,
            "endpoint_Gram_scale": literal.order,
            "endpoint_entry_alphabet": [-1, 1],
            "endpoint_recipe": "sum_j X_j tensor A_j, with the certified T631 signed partition X",
            "complete_Gram_proof": "sum_i X_i X_i^T tensor A_i A_i^T =631 I tensor sum_i A_i A_i^T; cross terms cancel by the checked amicability identities",
            "flatness_proof": "Each X position selects exactly one signed A_j block; every entry of each block is a sign.",
            "complete_algebraic_endpoint_Gram_verified": True,
            "dense_endpoint_allocated": False,
            "literal_entry_and_row_recipe_available": True,
            "literal_row_controls": row_controls(literal),
            "order_status": statuses.get(str(w)),
        })
    return {
        "schema": "complete-commutative-operator-scalar-endpoints-v1",
        "T_design": t_proof,
        "input_sha256": hashlib.sha256(INPUTS.read_bytes()).hexdigest(),
        "implementation_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "outputs": outputs,
        "all_outputs_are_completed_literal_sign_matrix_recipes": True,
        "scalar_factor_audit": str(ALGEBRAS.relative_to(ROOT)),
        "scope": "The represented scalar algebra acts on the stated rational module. Its dimension is not assumed equal to the module dimension. Field factors, multiplicities and rational reconstruction are recorded in the separate exact algebra audit.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--widths", type=int, nargs="+")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = certificate(args.widths)
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"completed_orders": [x["output_order"] for x in report["outputs"]]}))


if __name__ == "__main__":
    main()

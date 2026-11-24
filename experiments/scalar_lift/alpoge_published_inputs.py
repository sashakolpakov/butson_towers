"""Replay fixed published real Hadamard inputs, without searching for sources.

The Table A2 order-716 input is decoded by the paper's width-four recipe.
It is an existing real matrix, not a cubic scalar action or a new order.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from experiments.scalar_lift.h668_width12 import decode_signs, goethals_seidel
from hadamard.inverse_lift import is_sign_hadamard

SOURCE_URL = "https://hexagonmath.org/pdf/2611.00637v1"
SOURCE_TEX_URL = "https://hexagonmath.org/api/v1/works/2611.00637v1/source?path=hadamard-revised.tex"
H716_HEX = (
    "6370A5D5F65EAF7C823604CE47AD9504D157BE511B024",
    "1C3826C46D2311AC3EADBA71E3EC1DC96988D61F56DD0",
    "7F530ED2DD2826F0F1CCBAC05678969162137878E65D4",
    "FFF4473346A7786A13B7920005DC665CAABC3509DBC90",
)


def known_h716():
    """Return the exact published matrix and independently checked metadata."""
    m, order = 178, 716
    sequences = np.array([decode_signs(word, m) for word in H716_HEX], dtype=np.int64)
    if sequences.sum(axis=1).tolist() != [2, 0, 0, 0]:
        raise ValueError("the published input has incorrect row sums")
    for shift in range(1, m):
        if int(np.sum(sequences*np.roll(sequences, shift, axis=1))) != -4:
            raise ValueError("the exact Legendre-quadruple condition failed")
    corner = np.array([[-1, 1, 1, -1], [1, -1, 1, -1],
                       [1, 1, -1, -1], [-1, -1, -1, -1]], dtype=np.int64)
    left = np.array([[1, 1, 1, -1], [-1, -1, 1, -1],
                     [-1, 1, -1, -1], [1, -1, -1, -1]], dtype=np.int64)
    numerator = -(corner @ left.T)
    if np.any(numerator % 2):
        raise ArithmeticError("the fixed top border is nonintegral")
    top = numerator//2
    core = goethals_seidel(sequences, (m,))
    matrix = np.block([[corner, np.kron(top, np.ones((1, m), dtype=np.int64))],
                       [np.kron(left, np.ones((m, 1), dtype=np.int64)), core]])
    if matrix.shape != (order, order) or not np.all(np.abs(matrix) == 1):
        raise ValueError("the fixed input has incorrect dimensions or entries")
    if not np.array_equal(matrix @ matrix.T, order*np.eye(order, dtype=np.int64)):
        raise ValueError("the full exact integer Gram check failed")
    if not is_sign_hadamard(matrix.tolist()):
        raise ValueError("the independent bit-distance Gram check failed")
    report = {
        "schema": "published-order716-input-v1",
        "source_url": SOURCE_URL,
        "source_tex_url": SOURCE_TEX_URL,
        "source_identifier": "hexagon:2611.00637v1",
        "source_license": "CC0-1.0",
        "source_location": "Table A2, order716; Proposition2; section6.1 width4 recipe",
        "source_author_as_printed": "Claude",
        "announced_by": ["Levent Alpöge", "Philippe Voinov", "Saul Reynolds-Haertle"],
        "order": order, "phase": 2, "core_group_moduli": [m], "border_width": 4,
        "literal_core_hex": list(H716_HEX), "core_row_sums": [2, 0, 0, 0],
        "full_integer_Gram_verified": True,
        "independent_bit_distance_Gram_verified": True,
        "all_nonzero_periodic_autocorrelation_sums_equal_minus4": True,
        "matrix_binary_sha256": hashlib.sha256(((matrix+1)//2).astype(np.uint8).tobytes()).hexdigest(),
        "source_construction_search_performed": False,
        "tree_edge_added": False, "new_order_claimed": False,
        "compatible_nonreal_scalar_action_supplied": False,
    }
    return matrix, report


def natural_cubic_reuse_audit():
    """Test the specified four-row quaternion chart, with an exact census.

    The closed-quadruple argument in ALPOGE_REUSE_AUDIT.md excludes its
    entire signed row-permutation conjugacy class for order4r, r odd >1.
    It does not exclude a dense rational scalar operator.
    """
    matrix, source = known_h716()
    r, m = 179, 178
    permutation = np.concatenate([np.r_[f, np.arange(4+m*f, 4+m*(f+1))]
                                  for f in range(4)])
    frame = matrix[permutation]
    twice_action = np.array([[1, -1, -1, -1], [1, 1, -1, 1],
                             [1, 1, 1, -1], [1, -1, 1, 1]], dtype=np.int64)
    unit = np.eye(4, dtype=np.int64)
    if not np.array_equal(twice_action @ twice_action.T, 4*unit):
        raise ArithmeticError("the quaternion action is not orthogonal")
    if np.any(twice_action @ twice_action-2*twice_action+4*unit):
        raise ArithmeticError("the exact sextic scalar polynomial failed")
    numerator = np.einsum("ab,brc->arc", twice_action,
                          frame.reshape(4, r, 4*r)).reshape(4*r, 4*r)
    if np.any(numerator % 2):
        raise ArithmeticError("the specified scalar images are nonintegral")
    image = numerator//2
    if not np.array_equal(image @ image.T, 716*np.eye(716, dtype=np.int64)):
        raise ArithmeticError("the normalized scalar-image Gram failed")
    values, counts = np.unique(image, return_counts=True)
    parities, parity_counts = np.unique(np.prod(frame.reshape(4, r, 4*r), axis=0), return_counts=True)
    return {
        "schema": "published-h716-natural-cubic-chart-audit-v1",
        "source_matrix_binary_sha256": source["matrix_binary_sha256"],
        "row_partition": "channel f contains border row f followed by its178 core rows",
        "scalar_action": "R=(I+L_i+L_j+L_k)/2 tensor I_179",
        "twice_local_action": twice_action.tolist(),
        "scalar_action_is_orthogonal": True,
        "scalar_polynomial_R_squared_minus_R_plus_I": True,
        "complete_image_entry_counts": {str(int(v)): int(c) for v, c in zip(values, counts, strict=True)},
        "complete_quartet_product_counts": {str(int(v)): int(c) for v, c in zip(parities, parity_counts, strict=True)},
        "image_is_integer": True, "image_Gram_scale": 716,
        "image_all_signs": bool(np.all(np.abs(image) == 1)),
        "all_signed_row_quartetings_excluded_by_closed_quadruple_parity": True,
        "arbitrary_dense_rational_scalar_operator_excluded": False,
        "new_source_construction_or_tree_edge_performed": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--check", type=Path)
    parser.add_argument("--bundle", type=Path)
    parser.add_argument("--reuse-audit-report", type=Path)
    parser.add_argument("--check-reuse-audit", type=Path)
    args = parser.parse_args()
    matrix, report = known_h716()
    if args.check is not None and json.loads(args.check.read_text()) != report:
        raise ValueError("the published-input certificate differs from replay")
    if args.report is not None:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, indent=2)+"\n")
    if args.bundle is not None:
        args.bundle.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(args.bundle, matrix=matrix.astype(np.int8),
                            metadata_json=np.array(json.dumps(report, sort_keys=True)))
    if args.reuse_audit_report is not None:
        args.reuse_audit_report.write_text(json.dumps(natural_cubic_reuse_audit(), indent=2)+"\n")
    if (args.check_reuse_audit is not None
            and json.loads(args.check_reuse_audit.read_text()) != natural_cubic_reuse_audit()):
        raise ValueError("the natural cubic-reuse audit differs from replay")
    print(json.dumps({key: report[key] for key in
                      ("order", "matrix_binary_sha256", "full_integer_Gram_verified",
                       "independent_bit_distance_Gram_verified", "new_order_claimed")}))


if __name__ == "__main__":
    main()

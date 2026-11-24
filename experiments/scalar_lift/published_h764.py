"""Reconstruct the two fixed published order764 Hadamard inputs.

The literal coset data are from Djokovic, arXiv:math/0703312v1, pp2–3.
No difference sets, signs, sources, or frames are searched.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from experiments.scalar_lift.h668_width12 import goethals_seidel
from hadamard.inverse_lift import is_sign_hadamard

SOURCE = "https://arxiv.org/pdf/math/0703312v1"
REPRESENTATIVES = (1, 2, 3, 4, 6, 8, 9, 11, 12, 13, 16, 17, 18, 19, 22, 32, 36, 38, 41)
SUBGROUP = (1, 39, 184, 109, 49)
INDEX_SETS = (
    ((1, 7, 9, 10, 11, 13, 17, 18, 25, 26, 30, 31, 33, 34, 35, 36, 37),
     (1, 4, 7, 9, 11, 12, 13, 14, 19, 21, 22, 23, 24, 25, 26, 29, 36, 37),
     (0, 3, 4, 5, 7, 8, 9, 16, 17, 19, 24, 25, 29, 30, 31, 33, 35, 37),
     (1, 3, 4, 5, 8, 11, 14, 18, 19, 20, 21, 23, 24, 25, 28, 29, 30, 32, 34, 35)),
    ((0, 1, 6, 8, 9, 11, 12, 16, 18, 20, 21, 23, 28, 31, 33, 36, 37),
     (0, 1, 3, 4, 10, 12, 13, 17, 20, 22, 24, 31, 32, 33, 34, 35, 36, 37),
     (4, 8, 9, 10, 12, 13, 14, 16, 17, 20, 21, 24, 26, 27, 29, 31, 32, 34),
     (1, 7, 9, 10, 11, 12, 14, 15, 16, 17, 20, 22, 23, 25, 28, 29, 32, 33, 34, 37)),
)


def known_h764(index=0):
    if index not in (0, 1):
        raise ValueError("the printed paper contains exactly two inputs")
    cosets = []
    for representative in REPRESENTATIVES:
        coset = {representative*x % 191 for x in SUBGROUP}
        cosets.extend((coset, {-x % 191 for x in coset}))
    if len(set().union(*cosets)) != 190 or any(len(c) != 5 for c in cosets):
        raise ArithmeticError("the printed subgroup cosets do not partition the nonzero residues")
    negatives = [set().union(*(cosets[j] for j in indices)) for indices in INDEX_SETS[index]]
    if list(map(len, negatives)) != [85, 90, 90, 100]:
        raise ArithmeticError("the printed SDS cardinalities changed")
    for shift in range(1, 191):
        differences = sum(sum((x+shift) % 191 in block for x in block) for block in negatives)
        if differences != 174:
            raise ArithmeticError("the exact supplementary difference condition failed")
    sequences = np.array([[-1 if j in block else 1 for j in range(191)]
                          for block in negatives], dtype=np.int64)
    matrix = goethals_seidel(sequences, (191,))
    if not np.all(np.abs(matrix) == 1) or not np.array_equal(
            matrix @ matrix.T, 764*np.eye(764, dtype=np.int64)):
        raise ArithmeticError("the full integer sign Gram failed")
    if not is_sign_hadamard(matrix.tolist()):
        raise ArithmeticError("the independent bit-distance Gram failed")
    report = {
        "schema": "published-h764-input-v1", "source_url": SOURCE,
        "author": "Dragomir Z. Djokovic", "published_seed_index_zero_based": index,
        "source_location": "pp2–3, the two explicitly printed SDS index quadruples",
        "subgroup": list(SUBGROUP), "positive_coset_representatives": list(REPRESENTATIVES),
        "SDS_index_sets": [list(indices) for indices in INDEX_SETS[index]],
        "SDS_cardinalities": [85, 90, 90, 100], "SDS_lambda": 174,
        "all190_difference_conditions_exactly_verified": True,
        "core_row_sums": sequences.sum(axis=1).tolist(),
        "order": 764, "phase": 2,
        "full_integer_sign_Gram_verified": True,
        "independent_bit_distance_Gram_verified": True,
        "matrix_binary_sha256": hashlib.sha256(((matrix+1)//2).astype(np.uint8).tobytes()).hexdigest(),
        "new_source_or_frame_search_performed": False, "tree_edge_added": False,
        "compatible_nonreal_scalar_action_supplied": False,
        "new_Hadamard_order_claimed": False,
    }
    return matrix, report


def build_report(bundle_directory=None):
    results = []
    for index in (0, 1):
        matrix, report = known_h764(index)
        results.append(report)
        if bundle_directory is not None:
            directory = Path(bundle_directory)
            directory.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(directory/f"h764-{index}.npz", matrix=matrix.astype(np.int8),
                                metadata_json=np.array(json.dumps(report, sort_keys=True)))
    return {"schema": "published-h764-two-input-acquisition-v1", "inputs": results,
            "known_H764_frames_acquired": 2, "new_order6876_constructed": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--check", type=Path)
    parser.add_argument("--bundle-directory", type=Path)
    args = parser.parse_args()
    report = build_report(args.bundle_directory)
    if args.report is not None:
        args.report.write_text(json.dumps(report, indent=2)+"\n")
    if args.check is not None and json.loads(args.check.read_text()) != report:
        raise ValueError("the published order764 certificate differs from exact replay")
    print(json.dumps({"verified_published_H764_inputs": 2,
                      "matrix_hashes": [r["matrix_binary_sha256"] for r in report["inputs"]]}))


if __name__ == "__main__":
    main()

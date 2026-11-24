"""Replay all five unbordered published Alpoge/Claude GS frame seeds.

Only fixed Table A1 data are decoded.  No source or frame is searched.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from experiments.scalar_lift.h668_width12 import goethals_seidel
from hadamard.inverse_lift import is_sign_hadamard

DATA = Path(__file__).with_name("broad_catalogue_inputs")/"alpoge-table-a1.json"
PARAMETERS = {892: (223, 3, 3), 1132: (283, 3, 3), 1244: (311, 5, 17),
              1948: (487, 9, 3), 1964: (491, 7, 2)}
TEX_SHA256 = "b40d2420405ce96ac56f154cb36adb071c7f67fa2f6b7bd7a23e5ae3b485e15f"
DATA_SHA256 = "8d1967f7585091fa5e266722bc7e57ca5d9ccaf1146c8472c472f061957daee2"


def known_gs_frame(order):
    if order not in PARAMETERS:
        raise ValueError("order is not in the published unbordered TableA1 list")
    encoded = DATA.read_bytes()
    if hashlib.sha256(encoded).hexdigest() != DATA_SHA256:
        raise ValueError("the immutable extracted published seed changed")
    data = json.loads(encoded)
    if data["source_file_sha256"] != TEX_SHA256:
        raise ValueError("the recorded primary TeX source identity changed")
    record = next(row for row in data["records"] if row["order"] == order)
    prime, multiplier_size, generator = PARAMETERS[order]
    orbit_count = (prime-1)//multiplier_size
    sequences = []
    for word in record["orbit_signs"]:
        if len(word) != orbit_count+1 or set(word)-set("+-"):
            raise ValueError("the printed orbit word has incorrect length or alphabet")
        sequence = np.empty(prime, dtype=np.int64)
        sequence[0] = 1 if word[0] == "+" else -1
        filled = {0}
        for j in range(orbit_count):
            sign = 1 if word[j+1] == "+" else -1
            for u in range(multiplier_size):
                position = pow(generator, j+orbit_count*u, prime)
                if position in filled:
                    raise ArithmeticError("the printed multiplier orbits overlap")
                filled.add(position)
                sequence[position] = sign
        if len(filled) != prime:
            raise ArithmeticError("the multiplier orbits omit a core coordinate")
        sequences.append(sequence)
    sequences = np.array(sequences)
    for shift in range(1, prime):
        if np.sum(sequences*np.roll(sequences, shift, axis=1)) != 0:
            raise ArithmeticError("the full periodic complementarity condition failed")
    matrix = goethals_seidel(sequences, (prime,))
    if matrix.shape != (order, order) or not np.all(np.abs(matrix) == 1):
        raise ArithmeticError("the published GS output has incorrect signs or dimensions")
    if not np.array_equal(matrix @ matrix.T, order*np.eye(order, dtype=np.int64)):
        raise ArithmeticError("the full exact integer Gram failed")
    if not is_sign_hadamard(matrix.tolist()):
        raise ArithmeticError("the independent bit-distance Gram failed")
    minus = (-np.arange(prime)) % prime
    report = {
        "schema": "published-unbordered-alpoge-gs-input-v1",
        "source_url": "https://hexagonmath.org/pdf/2611.00637v1",
        "source_tex_url": data["source_url"], "source_tex_sha256": TEX_SHA256,
        "extracted_data_sha256": hashlib.sha256(DATA.read_bytes()).hexdigest(),
        "source_location": "TableA1, Table3, section6.1 unbordered recipe",
        "source_license": data["source_license"],
        "order": order, "prime_core_order": prime,
        "multiplier_size": multiplier_size, "primitive_generator": generator,
        "literal_orbit_words": record["orbit_signs"],
        "core_row_sums": sequences.sum(axis=1).tolist(),
        "symmetric_core_blocks": [bool(np.array_equal(row, row[minus])) for row in sequences],
        "all_nonzero_periodic_autocorrelation_sums_zero": True,
        "full_integer_sign_Gram_verified": True,
        "independent_bit_distance_Gram_verified": True,
        "matrix_binary_sha256": hashlib.sha256(((matrix+1)//2).astype(np.uint8).tobytes()).hexdigest(),
        "new_source_or_frame_search_performed": False, "tree_edge_added": False,
        "compatible_nonreal_scalar_action_supplied": False,
        "new_Hadamard_order_claimed": False,
    }
    return matrix, report


def build_report(bundle_directory=None, orders=None):
    result = []
    for order in (PARAMETERS if orders is None else orders):
        matrix, report = known_gs_frame(order)
        result.append(report)
        if bundle_directory is not None:
            directory = Path(bundle_directory)
            directory.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(directory/f"h{order}.npz", matrix=matrix.astype(np.int8),
                                metadata_json=np.array(json.dumps(report, sort_keys=True)))
    return {"schema": "published-unbordered-alpoge-gs-acquisition-v1", "inputs": result,
            "published_frames_acquired": len(result), "new_Hadamard_order_claimed": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--check", type=Path)
    parser.add_argument("--bundle-directory", type=Path)
    parser.add_argument("--order", type=int, action="append")
    args = parser.parse_args()
    report = build_report(args.bundle_directory, args.order)
    if args.report is not None:
        args.report.write_text(json.dumps(report, indent=2)+"\n")
    if args.check is not None and json.loads(args.check.read_text()) != report:
        raise ValueError("the published frame certificate differs from exact replay")
    print(json.dumps({"verified_published_frames": [r["order"] for r in report["inputs"]],
                      "new_Hadamard": False}))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Replay periodic T-matrices of order127; optionally verify/write known H1524.

Only the standard library and NumPy are required. No search or repository
module is imported, and no reference matrix archive is needed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np

HERE = Path(__file__).resolve().parent
P = 127
MU = 22
MULTIPLIER_ORDER = 9
WIDTH = 508
ORDER = 1524
STATES = (0, 0, 4, 0, 1, 6, 2, 3, 6, 1, 0, 3, 5, 4, 3)
D4 = np.array(((0, 1, 1, 1), (-1, 0, 1, -1),
               (-1, -1, 0, 1), (-1, 1, -1, 0)), dtype=np.int64)
P_HASH = "907e3feb8e602a6e11338654e068ded900e03cc065b8675322d3c5809a4a261e"
Q_HASH = "72c2d61022b7cfd28349afa2e89fd915477bd757f941a93bad1ffba68c22c269"
PACKED_HASH = "2f90301a69292ef94dfb8cb5d2bd2329b094eebf882e5b1a4bfebfecf1934ef8"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def expand_sequences():
    record = json.loads((HERE / "t127.json").read_text())
    require((record["p"], record["multiplier"], record["multiplier_order"])
            == (P, MU, MULTIPLIER_ORDER), "wrong T127 parameters")
    states = record["states"]
    require(tuple(states) == STATES and all(type(s) is int and 0 <= s < 8 for s in states),
            "wrong literal orbit states")
    require(pow(MU, 9, P) == 1 and pow(MU, 3, P) != 1,
            "multiplier does not have order9")
    b = np.zeros((4, P), dtype=np.int64)
    b[0, 0] = 1
    seen, representatives, orbit = {0}, [0], 0
    for first in range(1, P):
        if first in seen:
            continue
        orbit += 1
        representatives.append(first)
        axis, negative = divmod(states[orbit], 2)
        x = first
        for k in range(MULTIPLIER_ORDER):
            require(x not in seen, "overlapping multiplier orbits")
            seen.add(x)
            channel = 0 if axis == 0 else 1 + (axis - 1 + k) % 3
            b[channel, x] = 1 - 2 * negative
            x = x * MU % P
        require(x == first, "multiplier orbit did not close")
    require(len(seen) == P and orbit + 1 == len(states), "incomplete orbit cover")
    require(representatives == record["orbit_representatives"], "different orbit ordering")
    require(np.array_equal(b, np.asarray(record["b_sequences"])), "saved sequences differ")
    require(np.all(np.sum(b * b, axis=0) == 1), "supports do not partition the positions")
    require(b.sum(axis=1).tolist() == [10, 3, 3, 3], "wrong T127 row sums")
    for shift in range(P):
        value = int(np.sum(b * np.roll(b, -shift, axis=1)))
        require(value == (P if shift == 0 else 0), f"periodic correlation {shift} failed")
    return b


def t_matrices(b):
    indices = np.arange(P)
    return np.array([row[(indices[None, :] - indices[:, None]) % P] for row in b],
                    dtype=np.int64)


def coefficient_pair(b):
    require(np.array_equal(D4.T, -D4)
            and np.array_equal(D4.T @ D4, 3 * np.eye(4, dtype=np.int64)),
            "skew conference identities failed")
    a = D4 @ b
    record = json.loads((HERE / "t127.json").read_text())
    require(np.array_equal(a, np.asarray(record["a_sequences"])), "saved a differs from D4 b")
    require(np.all(a * b == 0) and np.all(a * a + b * b == 1), "coefficient flatness failed")
    for shift in range(P):
        ar, br = np.roll(a, -shift, axis=1), np.roll(b, -shift, axis=1)
        value = (int(np.sum(a * ar)), int(np.sum(b * br)),
                 int(np.sum(a * br + b * ar)))
        require(value == ((3 * P, P, 0) if shift == 0 else (0, 0, 0)),
                f"coefficient correlation {shift} failed")
    return a


def goethals_seidel(sequences):
    indices = np.arange(P)
    A, B, C, E = [row[(indices[None, :] - indices[:, None]) % P]
                  for row in np.asarray(sequences, dtype=np.int8)]
    rev = np.arange(P - 1, -1, -1)
    return np.block([[A, B[:, rev], C[:, rev], E[:, rev]],
                     [-B[:, rev], A, -E.T[:, rev], C.T[:, rev]],
                     [-C[:, rev], E.T[:, rev], A, -B.T[:, rev]],
                     [-E[:, rev], -C.T[:, rev], B.T[:, rev], A]])


def packed_rows(matrix):
    require(matrix.ndim == 2 and matrix.shape[0] == matrix.shape[1]
            and np.all(np.abs(matrix) == 1), "invalid square sign matrix")
    packed = np.packbits(matrix < 0, axis=1, bitorder="little")
    return packed, [int.from_bytes(row.tobytes(), "little") for row in packed]


def full_gram(rows, order):
    require(len(rows) == order, "wrong row count")
    for i, row in enumerate(rows):
        for j in range(i):
            require((row ^ rows[j]).bit_count() == order // 2,
                    f"integer row Gram failed at {i},{j}")


def source_gram():
    source = [[(-1, 0) if i == j else (0, 1) for j in range(3)] for i in range(3)]
    for i in range(3):
        for j in range(3):
            first = second = 0
            for k in range(3):
                a, b = source[i][k]
                c, d = source[j][k]
                c, d = c + d, -d  # conjugate in Z[zeta6], zeta6^2=zeta6-1
                first += a * c - b * d
                second += a * d + b * c + b * d
            require((first, second) == ((3, 0) if i == j else (0, 0)),
                    "exact BH(3,6) source Gram failed")


def hadamard_control(a, b, write=None, compare=None):
    source_gram()
    frame, image = goethals_seidel(a - b), goethals_seidel(a + b)
    pp, pr = packed_rows(frame)
    qp, qr = packed_rows(image)
    require(digest(pp.tobytes()) == P_HASH and digest(qp.tobytes()) == Q_HASH,
            "frame differs from the independent reference audit")
    full_gram(pr, WIDTH)
    full_gram(qr, WIDTH)
    numerators = set()
    for i in range(WIDTH):
        diagonal = WIDTH - 2 * (qr[i] ^ pr[i]).bit_count()
        require(diagonal == WIDTH // 2, "cross diagonal failed")
        numerators.add(diagonal)
        for j in range(i):
            ij = WIDTH - 2 * (qr[i] ^ pr[j]).bit_count()
            ji = WIDTH - 2 * (qr[j] ^ pr[i]).bit_count()
            require(ij + ji == 0, "symmetric cross identity failed")
            numerators.update((ij, ji))
    matrix = np.block([[-frame, image, image], [image, -frame, image],
                       [image, image, -frame]])
    packed, rows = packed_rows(matrix)
    require(matrix.shape == (ORDER, ORDER) and digest(packed.tobytes()) == PACKED_HASH,
            "H1524 differs from the independent reference audit")
    require(packed.shape == (ORDER, 191) and np.all((packed[:, -1] & 0xf0) == 0),
            "incorrect packed shape or padding")
    full_gram(rows, ORDER)
    if compare is not None:
        with np.load(compare, allow_pickle=False) as saved:
            require(set(saved.files) == {"negative_entry_bits", "order", "bitorder"},
                    "unexpected matrix archive fields")
            require(saved["order"].shape == () and saved["order"].item() == ORDER
                    and saved["bitorder"].shape == () and saved["bitorder"].item() == "little",
                    "incorrect archive parameters")
            require(saved["negative_entry_bits"].dtype == np.uint8
                    and np.array_equal(saved["negative_entry_bits"], packed),
                    "archive differs from the reconstructed matrix")
    if write is not None:
        require(write.suffix == ".npz" and not write.exists(), "output must be a new NPZ file")
        np.savez_compressed(write, negative_entry_bits=packed, order=np.array(ORDER),
                            bitorder=np.array("little"))
    return {"known_order_control": ORDER, "actual_alphabet": [-1, 1],
            "all_P_Q_and_symmetric_cross_integer_Grams_verified": True,
            "full_H1524_integer_Gram_verified": True,
            "distinct_Hadamard_row_pairs_checked": ORDER * (ORDER - 1) // 2,
            "packed_sha256": PACKED_HASH, "field": "Q(zeta_6)/Q", "field_degree": 2,
            "module_rank": WIDTH // 2, "rational_dimension": WIDTH,
            "source_gram_scale": 3, "frame_gram_scale": WIDTH,
            "endpoint_gram_scale": ORDER, "actual_R_entry_numerators": sorted(numerators),
            "actual_R_entry_denominator": WIDTH, "new_real_order_claimed": False}


def verify(hadamard=False, write=None, compare=None):
    start = time.monotonic()
    b = expand_sequences()
    t = t_matrices(b)
    require(np.all(np.sum(np.abs(t), axis=0) == 1), "T-matrix supports are not disjoint")
    require(np.array_equal(sum(matrix @ matrix.T for matrix in t),
                           P * np.eye(P, dtype=np.int64)), "complete T-matrix Gram failed")
    a = coefficient_pair(b)
    result = {"periodic_T_matrix_order": P, "multiplier": MU,
              "multiplier_order": MULTIPLIER_ORDER,
              "actual_T_entry_alphabet": [-1, 0, 1],
              "b_row_sums": b.sum(axis=1).tolist(),
              "b_support_sizes": np.sum(b * b, axis=1).tolist(),
              "all127_periodic_correlations_verified": True,
              "full_T_matrix_disjoint_support_and_integer_Gram_verified": True,
              "T_matrix_Gram_scale": P,
              "all127_coefficient_pair_correlations_verified": True,
              "aperiodic_T_sequence_claimed": False}
    if hadamard or write is not None or compare is not None:
        result["scalar_Hadamard_control"] = hadamard_control(a, b, write, compare)
    result["elapsed_seconds"] = time.monotonic() - start
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hadamard", action="store_true", help="also check every H1524 row pair")
    parser.add_argument("--write-hadamard", type=Path, help="also write a new verified H1524 NPZ")
    parser.add_argument("--compare-npz", type=Path, help="also compare a saved H1524 NPZ")
    args = parser.parse_args()
    print(json.dumps(verify(args.hadamard, args.write_hadamard, args.compare_npz), indent=2))

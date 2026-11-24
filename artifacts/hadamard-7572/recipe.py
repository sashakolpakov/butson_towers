#!/usr/bin/env python3
"""Reconstruct H7572 and verify every integer row inner product.

Only Python's standard library and NumPy are required. No repository imports
or search programs are used. Run from any directory: python /path/recipe.py.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np


HERE = Path(__file__).resolve().parent
P = 631
MU = 458
MULTIPLIER_ORDER = 21
D = 4 * P
ORDER = 3 * D
D4 = np.array(((0, 1, 1, 1), (-1, 0, 1, -1),
               (-1, -1, 0, 1), (-1, 1, -1, 0)), dtype=np.int64)
# Byte hash of the original reference archive, not a requirement on new ZIP files.
NPZ_SHA256 = "cb0541a9019643d1e76698405c260731458627c01c8516bcb6d4d13c626d6364"
PACKED_SHA256 = "bc3cb6fca9243a6ef01c82de2a446080c8866e0f30d76c0f23ef0508ae5e2e9b"
P_SHA256 = "6797e67ca52723c3a0d787ab10099742dbe765f806b937d69b01e656f6a42110"
Q_SHA256 = "3f0b69cc4a2ab3b5f1a0d2f4d13e6e5f51247ad1741f59be3156e092d988ef36"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def verify_source():
    """Exact arithmetic in Z[z]/(z^2-z+1), with z conjugate = 1-z."""
    def multiply(x, y):
        a, b = x
        c, d = y
        return a * c - b * d, a * d + b * c + b * d

    source = [[(-1, 0) if i == j else (0, 1) for j in range(3)]
              for i in range(3)]
    for i in range(3):
        for j in range(3):
            total = [0, 0]
            for k in range(3):
                a, b = source[j][k]
                term = multiply(source[i][k], (a + b, -b))
                total = [total[t] + term[t] for t in range(2)]
            require(total == ([3, 0] if i == j else [0, 0]),
                    "BH(3,6) source Gram failed")


def expand_sequences():
    data = json.loads((HERE / "t631.json").read_text())
    require((data["p"], data["multiplier"], data["multiplier_order"])
            == (P, MU, MULTIPLIER_ORDER), "wrong T631 parameters")
    states = data["states"]
    require(len(states) == 31 and all(type(s) is int and 0 <= s < 8 for s in states)
            and states[:2] == [0, 0], "invalid orbit states")
    require(pow(MU, 21, P) == 1 and pow(MU, 7, P) != 1
            and pow(MU, 3, P) != 1, "multiplier does not have order 21")
    seen = {0}
    orbits = [[(0, 0)]]
    for representative in range(1, P):
        if representative in seen:
            continue
        orbit = []
        x = representative
        for k in range(MULTIPLIER_ORDER):
            require(x not in seen, "overlapping multiplier orbits")
            seen.add(x)
            orbit.append((x, k % 3))
            x = MU * x % P
        require(x == representative, "multiplier orbit did not close")
        orbits.append(orbit)
    require(len(seen) == P and len(orbits) == len(states), "incomplete orbit cover")
    require([o[0][0] for o in orbits] == data["orbit_representatives"],
            "orbit ordering differs from recipe")
    b = np.zeros((4, P), dtype=np.int64)
    for state, orbit in zip(states, orbits):
        axis, negative = divmod(state, 2)
        for x, phase in orbit:
            row = 0 if axis == 0 else 1 + (axis - 1 + phase) % 3
            b[row, x] = -1 if negative else 1
    require(np.array_equal(b, np.asarray(data["b_sequences"])),
            "saved b sequences differ from orbit expansion")
    require(b.sum(axis=1).tolist() == [22, 7, 7, 7], "wrong b row sums")
    require(np.all(np.sum(b * b, axis=0) == 1), "b is not one-hot")
    for shift in range(P):
        value = int(np.sum(b * np.roll(b, -shift, axis=1)))
        require(value == (P if shift == 0 else 0), f"T631 correlation {shift} failed")
    require(np.array_equal(D4.T, -D4)
            and np.array_equal(D4 @ D4.T, 3 * np.eye(4, dtype=np.int64)),
            "skew conference identity failed")
    a = D4 @ b
    require(np.array_equal(a, np.asarray(data["a_sequences"])),
            "saved a sequences differ from D4 b")
    require(np.all(a * b == 0) and np.all(a * a + b * b == 1),
            "coefficient entry flatness failed")
    for shift in range(P):
        ar = np.roll(a, -shift, axis=1)
        br = np.roll(b, -shift, axis=1)
        observed = (int(np.sum(a * ar)), int(np.sum(b * br)),
                    int(np.sum(a * br + b * ar)))
        require(observed == ((3 * P, P, 0) if shift == 0 else (0, 0, 0)),
                f"coefficient correlation {shift} failed")
    return a, b


def goethals_seidel(sequences):
    """Literal circulants and reversal; this fixes the full matrix convention."""
    indices = np.arange(P)
    A, B, C, E = [row[(indices[None, :] - indices[:, None]) % P]
                  for row in np.asarray(sequences, dtype=np.int8)]
    rev = np.arange(P - 1, -1, -1)
    return np.block([[A, B[:, rev], C[:, rev], E[:, rev]],
                     [-B[:, rev], A, -E.T[:, rev], C.T[:, rev]],
                     [-C[:, rev], E.T[:, rev], A, -B.T[:, rev]],
                     [-E[:, rev], -C.T[:, rev], B.T[:, rev], A]])


def pack_signs(matrix):
    require(matrix.ndim == 2 and matrix.shape[0] == matrix.shape[1]
            and np.all((matrix == -1) | (matrix == 1)), "matrix is not square/sign")
    packed = np.packbits(matrix < 0, axis=1, bitorder="little")
    integers = [int.from_bytes(row.tobytes(), "little") for row in packed]
    return packed, integers


def verify_gram(rows, order):
    """For signs, <r_i,r_j> = order - 2*HammingDistance(r_i,r_j)."""
    require(len(rows) == order, "wrong number of rows")
    for i, row in enumerate(rows):
        for j in range(i):
            require((row ^ rows[j]).bit_count() == order // 2,
                    f"integer Gram failed at rows {i},{j}")


def verify(write=None):
    started = time.monotonic()
    verify_source()
    a, b = expand_sequences()
    frame, image = goethals_seidel(a - b), goethals_seidel(a + b)
    packed_p, rows_p = pack_signs(frame)
    packed_q, rows_q = pack_signs(image)
    require(sha256(packed_p.tobytes()) == P_SHA256, "P hash differs")
    require(sha256(packed_q.tobytes()) == Q_SHA256, "Q hash differs")
    verify_gram(rows_p, D)
    verify_gram(rows_q, D)
    r_numerators = {D // 2}
    for i in range(D):
        require((rows_p[i] ^ rows_q[i]).bit_count() == D // 4,
                "symmetric cross diagonal failed")
        for j in range(i):
            distance_1 = (rows_p[i] ^ rows_q[j]).bit_count()
            distance_2 = (rows_q[i] ^ rows_p[j]).bit_count()
            require(distance_1 + distance_2 == D,
                    f"symmetric cross identity failed at {i},{j}")
            r_numerators.update((D - 2 * distance_1, D - 2 * distance_2))
    # Source entries: diagonal -1 -> -P, off-diagonal zeta_6 -> Q.
    matrix = np.block([[-frame, image, image], [image, -frame, image],
                       [image, image, -frame]])
    packed, rows = pack_signs(matrix)
    require(matrix.shape == (ORDER, ORDER), "wrong endpoint dimension")
    require(sha256(packed.tobytes()) == PACKED_SHA256, "endpoint packed hash differs")
    verify_gram(rows, ORDER)
    artifact = HERE / "h7572.npz"
    archive_checked = artifact.exists()
    if archive_checked:
        require(sha256(artifact.read_bytes()) == NPZ_SHA256, "immutable NPZ hash differs")
        with np.load(artifact, allow_pickle=False) as saved:
            require(set(saved.files) == {"negative_entry_bits", "order", "bitorder"},
                    "unexpected matrix archive fields")
            require(saved["order"].shape == () and int(saved["order"]) == ORDER
                    and saved["bitorder"].shape == () and str(saved["bitorder"]) == "little",
                    "wrong matrix archive parameters")
            bits = saved["negative_entry_bits"]
            require(bits.dtype == np.uint8 and bits.shape == packed.shape
                    and np.array_equal(bits, packed), "archive differs from reconstruction")
            require(np.all((bits[:, -1] >> (ORDER % 8)) == 0), "nonzero padding bits")
    else:
        require(write is not None, "reference archive missing; pass --write PATH.npz to reconstruct")
    if write is not None:
        require(write.suffix == ".npz", "output name must end in .npz")
        require(not write.exists(), "refusing to overwrite an existing output")
        np.savez_compressed(write, negative_entry_bits=packed,
                            order=np.array(ORDER), bitorder=np.array("little"))
    return {"order": ORDER, "actual_alphabet": [-1, 1],
            "all_expanded_coefficient_correlations_verified": True,
            "full_P_and_Q_integer_Grams_verified": True,
            "full_symmetric_cross_identity_verified": True,
            "full_endpoint_integer_Gram_verified": True,
            "endpoint_row_pairs_checked": ORDER * (ORDER - 1) // 2,
            "matrix_npz_sha256": NPZ_SHA256, "packed_sign_bits_sha256": PACKED_SHA256,
            "reference_archive_checked": archive_checked,
            "scalar_action_entry_numerators": sorted(r_numerators),
            "scalar_action_entry_denominator": D,
            "novelty_audit": "novelty.json",
            "elapsed_seconds": time.monotonic() - started}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", type=Path, help="also write a new reconstructed NPZ")
    arguments = parser.parse_args()
    print(json.dumps(verify(arguments.write), indent=2))

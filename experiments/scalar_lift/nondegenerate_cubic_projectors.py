"""Exact binary first-lift templates without zero W rows or isolated T components.

Projector idempotence is recorded separately and is not required for a
first-lift witness. Verification does not rerun the bounded search.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import random
import time
from pathlib import Path

import numpy as np

from experiments.scalar_lift.cubic_core_support_lift import source_layers
from experiments.scalar_lift.mixed_cubic_boundaries import seven_border
from experiments.scalar_lift.mixed_cubic_ramification import CONJUGATE, derivative, product

DEFAULT_CERTIFICATE = Path(__file__).resolve().parents[2] / "certificates/scalar-lift/nondegenerate-cubic-first-lifts.json"


def bipartite_components(matrix):
    n = len(matrix)
    unseen = set(range(2 * n))
    components = []
    while unseen:
        stack = [min(unseen)]
        component = []
        while stack:
            node = stack.pop()
            if node not in unseen:
                continue
            unseen.remove(node)
            component.append(node)
            neighbors = (n + np.flatnonzero(matrix[node]) if node < n else
                         np.flatnonzero(matrix[:, node - n]))
            stack.extend(int(value) for value in neighbors)
        components.append(sorted(component))
    return components


def binary_affine_solutions(c, active):
    rows, columns = np.nonzero(c)
    count = len(rows)
    equations = (c[rows[:, None], columns[None, :]] & c[rows[None, :], columns[:, None]])
    boundary = []
    rhs = [0] * count
    for i in range(len(c)):
        boundary.extend((rows == i, columns == i,
                         (rows == i) & (active[columns] != 0),
                         (columns == i) & (active[rows] != 0)))
        rhs.extend((0, 0, int(active[i]), int(active[i])))
    matrix = np.vstack((equations, np.array(boundary, dtype=np.uint8)))
    packed = np.packbits(matrix, axis=1, bitorder="little")
    bitrows = [int.from_bytes(row.tobytes(), "little") | (value << count)
               for row, value in zip(packed, rhs, strict=True)]
    mask, pivots = (1 << count) - 1, {}
    for row in bitrows:
        while row & mask:
            pivot = ((row & mask) & -(row & mask)).bit_length() - 1
            if pivot in pivots:
                row ^= pivots[pivot]
            else:
                pivots[pivot] = row
                break
        if row & mask == 0 and row >> count:
            return None
    descending = sorted(pivots.items(), reverse=True)
    solution = 0
    for pivot, row in descending:
        if ((row & solution).bit_count() ^ (row >> count)) & 1:
            solution ^= 1 << pivot
    basis = []
    for free in range(count):
        if free in pivots:
            continue
        vector = 1 << free
        for pivot, row in descending:
            if (row & vector).bit_count() & 1:
                vector ^= 1 << pivot
        basis.append(vector)
    return rows, columns, solution, basis


def template_checks(c, b, active, *, require_projector=True):
    n = len(c)
    d = c @ b.T @ c % 2
    w = b ^ d
    identity = np.eye(n, dtype=np.int64)
    linear = (np.array_equal(c @ c.T % 2, identity)
              and np.array_equal(c @ active % 2, active)
              and np.array_equal(c.T @ active % 2, active)
              and not np.any(b * (1 - c)) and not np.any(c * d)
              and not np.any(b.sum(axis=0) % 2) and not np.any(b.sum(axis=1) % 2)
              and np.array_equal(b @ active % 2, active)
              and np.array_equal(b.T @ active % 2, active))
    if not linear:
        raise ValueError("a sampled affine solution fails the exact original equations")
    ell = c @ b.T % 2
    f0, f1 = ell.T, ell ^ ell.T
    square0 = (f0 @ f0 + f1 @ f1) % 2
    square1 = (f0 @ f1 + f1 @ f0 + f1 @ f1) % 2
    projector = np.array_equal(square0, f0) and np.array_equal(square1, f1)
    defect = (square0 ^ f0) ^ (2 * (square1 ^ f1))
    nonzero = np.argwhere(defect != 0)
    failure = None if len(nonzero) == 0 else [int(nonzero[0, 0]), int(nonzero[0, 1]), int(defect[tuple(nonzero[0])])]
    dense = bool(np.all(w.any(axis=0)) and np.all(w.any(axis=1)))
    t = 1 ^ c
    no_singleton = bool(np.min(t.sum(axis=0)) > 1 and np.min(t.sum(axis=1)) > 1)
    components = bipartite_components(t)
    return {"C": c.tolist(), "B": b.tolist(), "a": active.tolist(), "D": d.tolist(),
            "W": w.tolist(), "T": t.tolist(), "linear_criterion_verified": True,
            "projector_F_squared_equals_F": projector,
            "projector_failure_entry_over_F4": failure,
            "W_has_no_zero_rows_or_columns": dense,
            "T_has_no_singleton_rows_or_columns": no_singleton,
            "T_bipartite_components": components,
            "T_bipartite_connected": len(components) == 1,
            "W_row_weights": w.sum(axis=1).tolist(), "W_column_weights": w.sum(axis=0).tolist(),
            "first_lift_acceptable": dense and no_singleton and len(components) == 1,
            "projector_family_acceptable": dense and no_singleton and (projector or not require_projector)}


def search(n=12, weight=6, trials=5000, seed=20261122, seconds=50):
    rng = random.Random(seed)
    active = np.array([1] * weight + [0] * (n - weight), dtype=np.int64)
    allowed = [np.array([(word >> j) & 1 for j in range(n)], dtype=np.int64)
               for word in range(1, 1 << n) if word.bit_count() % 2 == 0
               and (word & ((1 << weight) - 1)).bit_count() % 2 == 0]
    c = np.eye(n, dtype=np.int64)
    consistent = samples = dense_count = projector_count = 0
    retained = {}
    start = time.monotonic()
    for trial in range(trials):
        for _ in range(30 if trial == 0 else 5):
            vector = rng.choice(allowed)
            c ^= vector[:, None] * (vector @ c % 2)[None, :]
        if np.any((1 ^ c).sum(axis=0) == 1) or np.any((1 ^ c).sum(axis=1) == 1):
            continue
        solved = binary_affine_solutions(c, active)
        if solved is None:
            if time.monotonic() - start > seconds:
                break
            continue
        consistent += 1
        rows, columns, particular, basis = solved
        combinations = (range(1 << len(basis)) if len(basis) <= 8 else
                        [0] + [rng.getrandbits(len(basis)) for _ in range(127)])
        for choice in combinations:
            vector = particular
            for j, direction in enumerate(basis):
                if (choice >> j) & 1:
                    vector ^= direction
            b = np.zeros_like(c)
            b[rows, columns] = [(vector >> j) & 1 for j in range(len(rows))]
            result = template_checks(c, b, active)
            samples += 1
            dense_count += result["W_has_no_zero_rows_or_columns"]
            projector_count += result["projector_F_squared_equals_F"]
            if result["W_has_no_zero_rows_or_columns"] and result["T_bipartite_connected"]:
                key = hashlib.sha256(bytes(c.ravel()) + bytes(b.ravel())).hexdigest()
                result["template_sha256"] = key
                retained[key] = result
            if result["projector_family_acceptable"]:
                return {"found": True, "template": result, "trial": trial + 1,
                        "consistent_C": consistent, "B_samples": samples,
                        "dense_W_samples": dense_count, "projector_samples": projector_count,
                        "affine_nullity": len(basis), "seed": seed,
                        "retained_first_lifts": list(retained.values()),
                        "elapsed_seconds": time.monotonic() - start}
        if time.monotonic() - start > seconds:
            break
    return {"found": False, "trial": trial + 1, "consistent_C": consistent,
            "B_samples": samples, "dense_W_samples": dense_count,
            "projector_samples": projector_count, "seed": seed,
            "retained_first_lifts": list(retained.values()),
            "elapsed_seconds": time.monotonic() - start,
            "general_nonexistence_claimed": False}


def retained_templates(path=DEFAULT_CERTIFICATE):
    """Literal independently verifiable first lifts; no search on import."""
    report = json.loads(Path(path).read_text())
    return [{**{key: np.asarray(record[key], dtype=np.int64)
                for key in ("C", "B", "a", "D", "W", "T")},
             "template_sha256": record["template_sha256"]}
            for record in report["retained_first_lifts"]]


def independent_first_layer_controls(record):
    template = {key: np.asarray(record[key], dtype=np.int64) for key in ("C", "B", "a", "D", "W", "T")}
    _, border7, vector7, _ = seven_border(48)
    borders = [(np.eye(3, dtype=np.int64), np.array(vector, dtype=np.int64))
               for vector in ((0, 1, 1), (1, 0, 1), (1, 1, 0))]
    borders.append((border7, vector7))
    controls = []
    for border, vector in borders:
        problem = {**template, "A": template["a"][:, None] * vector[None, :], "border_C": border}
        h0, h1 = source_layers(problem)
        if (not np.array_equal(product(h0, CONJUGATE[h0].T), np.eye(len(h0), dtype=np.int64))
                or np.any(derivative(h0, h1))):
            raise ValueError("a literal template fails the full F4 source Gram or first derivative")
        controls.append({"border_order": len(border), "boundary_vector": vector.tolist(),
                         "complete_local_source_order": len(h0),
                         "entire_F4_Gram_equals_I": True,
                         "entire_first_Gram_derivative_zero": True})
    return controls


def certify_retained(path=DEFAULT_CERTIFICATE, *, update=False):
    path = Path(path)
    report = json.loads(path.read_text())
    hashes = set()
    for record in report["retained_first_lifts"]:
        c, b, active = (np.asarray(record[key], dtype=np.int64) for key in ("C", "B", "a"))
        checks = template_checks(c, b, active)
        key = hashlib.sha256(c.astype("<i8").tobytes() + b.astype("<i8").tobytes()).hexdigest()
        if key != record["template_sha256"] or key in hashes or not checks["first_lift_acceptable"]:
            raise ValueError("a retained first-lift witness or its hash is invalid")
        hashes.add(key)
        controls = independent_first_layer_controls(record)
        if update:
            record.pop("acceptable", None)
            record.update(checks)
            record["independent_first_layer_controls"] = controls
        elif any(record.get(name) != value for name, value in checks.items()) or record.get("independent_first_layer_controls") != controls:
            raise ValueError("a saved exact first-lift audit differs from replay")
    if update:
        report["schema"] = "nondegenerate-cubic-first-lifts-v1"
        report["projector_found"] = report.pop("found", False)
        report["first_lift_witness_count"] = len(hashes)
        report["first_lift_active_weight"] = 4
        report["characteristic_zero_source_constructed"] = False
        report["full_modulo_two_projector_family_claimed"] = False
        report["replication_arithmetic"] = [
            {"source_order": 1499, "core_order": 748, "copies": 62, "template_order": 12,
             "active_rows": 248, "remaining_fillers": 4},
            {"source_order": 2999, "core_order": 1496, "copies": 124, "template_order": 12,
             "active_rows": 496, "remaining_fillers": 8}]
        path.write_text(json.dumps(report, indent=2) + "\n")
    return {"verified_first_lifts": len(hashes), "local_full_Gram_and_derivative_checks": 4 * len(hashes),
            "projector_family_claimed": False, "characteristic_zero_source_constructed": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--search-output", type=Path)
    action.add_argument("--check", type=Path)
    action.add_argument("--certify", type=Path)
    parser.add_argument("--weight", type=int, default=6)
    args = parser.parse_args()
    if args.search_output:
        report = search(weight=args.weight)
        args.search_output.write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps({key: value for key, value in report.items() if key not in ("template", "retained_first_lifts")}))
    else:
        print(json.dumps(certify_retained(args.check or args.certify, update=bool(args.certify))))


if __name__ == "__main__":
    main()

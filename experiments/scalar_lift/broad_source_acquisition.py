"""Acquire and audit 524 fixed published Butson inputs at useful orders.

Only monomial gauges of complete published matrices are considered here.
No Butson source is constructed from smaller matrices or searched for.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from functools import cache
from pathlib import Path

import numpy as np

from experiments.scalar_lift.forward_butson_descent import exponent_hash, verify_source
from experiments.scalar_lift.forward_source_registry import VerifiedSource
from experiments.scalar_lift.scalar_restriction_obstructions import alphabet_obstruction

DATA = Path(__file__).with_name("broad_catalogue_inputs")
CATALOGUE = "https://wiki.aalto.fi/spaces/Butson/pages/120482304/Matrices%2Bup%2Bto%2Bmonomial%2Bequivalence"
BASE = "https://wiki.aalto.fi/download/attachments/120482261/"
SPECS = (
    (7, 6, 2, "04aa420d4e3180c5cfe41b6595bccbace269237b3945452ad8a609e47b6aaab6"),
    (7, 12, 4, "7c906b5f2a69ae7b92a7b6684953024be2452a84f26f5ea3a4ff1dbaf8b1e45a"),
    (9, 6, 17, "e6b0c65c17700e3d9a5141b5ba74e025cbbde9d39809a25b89c98ef1083ff5a7"),
    (9, 12, 65, "d284fdd90d50339963261f872b18f6820250bc264be7c8c4f1567f4b70d76b67"),
    (13, 6, 436, "0c4194ba242ee737c9f06c0e0300aa1329d6ccef0beace67a3a7418134192976"),
)


@cache
def dataset(order, phase):
    spec = next(row for row in SPECS if row[:2] == (order, phase))
    path = DATA/f"BH-{order}-{phase}.txt"
    encoded = path.read_bytes()
    if hashlib.sha256(encoded).hexdigest() != spec[3]:
        raise ValueError("the immutable published dataset changed")
    records = encoded.splitlines()
    if len(records) != spec[2]:
        raise ValueError("the published matrix count changed")
    arrays = []
    header = bytes((66, 72, 49, 48+order, 48+order, 48+phase))
    for record in records:
        if record[:6] != header or len(record) != 6+order*order:
            raise ValueError("the published BH1 header or length is malformed")
        values = np.frombuffer(record[6:], dtype=np.uint8).astype(np.int64)-48
        if np.any(values < 0) or np.any(values >= phase):
            raise ValueError("a published exponent is outside its root alphabet")
        matrix = values.reshape(order, order)
        matrix.setflags(write=False)
        arrays.append(matrix)
    return tuple(arrays)


def sources():
    return [{"id": f"aalto-{n}-{q}-{index:03d}", "order": n, "phase": q,
             "family": "published-fixed-aalto-class", "catalogue_class_index": index}
            for n, q, count, _digest in SPECS for index in range(count)]


def root_gauge(matrix, phase):
    """Complete gauge CSP modulo signs for phases6 and12.

    An odd relation exists iff one full triple of antipodal classes occurs.
    At phase6 omit one residue modulo3. At phase12 omit one even and one
    odd residue modulo6. Global phase sets the first missing residue to a
    fixed value; the first row gauge may be fixed using the column gauges.
    """
    if phase not in (6, 12):
        raise ValueError("only the certified sixth/twelfth-root palettes are supported")
    n, modulus = len(matrix), phase//2
    reduced = (np.asarray(matrix, dtype=np.int64) % modulus).tolist()
    patterns = ((2,),) if phase == 6 else ((0, 1), (0, 3), (0, 5))
    nodes = 0
    for forbidden in patterns:
        allowed = [sum(1 << c for c in range(modulus)
                       if (exponent+c) % modulus not in forbidden)
                   for exponent in range(modulus)]

        def visit(rows, masks, allowed=allowed):
            nonlocal nodes
            nodes += 1
            if len(rows) == n:
                return rows, [(mask & -mask).bit_length()-1 for mask in masks]
            i = len(rows)
            for shift in ((0,) if i == 0 else range(modulus)):
                updated = [mask & allowed[(exponent+shift) % modulus]
                           for mask, exponent in zip(masks, reduced[i], strict=True)]
                if 0 in updated:
                    continue
                result = visit(rows+[shift], updated)
                if result is not None:
                    return result
            return None

        result = visit([], [(1 << modulus)-1]*n)
        if result is not None:
            return {"exists": True, "row_shifts": result[0], "column_shifts": result[1],
                    "missing_antipodal_residues": list(forbidden), "visited_row_prefixes": nodes,
                    "complete_over_all_row_column_root_gauges": True,
                    "higher_root_alphabet_gauges_covered": False}
    return {"exists": False, "row_shifts": None, "column_shifts": None,
            "visited_row_prefixes": nodes, "complete_over_all_row_column_root_gauges": True,
            "higher_root_alphabet_gauges_covered": False}


@cache
def load_source(identifier):
    try:
        prefix, n, q, index = identifier.split("-")
        order, phase, index = int(n), int(q), int(index)
        if prefix != "aalto":
            raise ValueError
        raw = dataset(order, phase)[index]
        if index < 0:
            raise ValueError
    except (ValueError, IndexError, StopIteration) as error:
        raise KeyError(identifier) from error
    raw_proof = verify_source(raw, phase)
    gauge = root_gauge(raw, phase)
    matrix = raw.copy()
    if gauge["exists"]:
        matrix = (matrix+np.asarray(gauge["row_shifts"])[:, None]
                  +np.asarray(gauge["column_shifts"])[None, :]) % phase
    global_anchor = int(matrix[0, 0])
    matrix = (matrix-global_anchor) % phase
    proof = verify_source(matrix, phase)
    used = tuple(map(int, np.unique(matrix)))
    gate = alphabet_obstruction(phase, used)
    if not gate["complete_for_odd_integer_relations"] or gauge["exists"] == gate["excluded"]:
        raise ArithmeticError("the independent exact alphabet gate disagrees with the gauge CSP")
    filename = f"BH-{order}-{phase}.txt"
    proof.update({"alphabet_gate": gate, "monomial_gauge_audit": gauge,
                  "acquisition": {"catalogue": CATALOGUE, "dataset_url": BASE+filename,
                      "dataset_sha256": hashlib.sha256((DATA/filename).read_bytes()).hexdigest(),
                      "class_index_zero_based": index,
                      "raw_source_full_Gram": raw_proof,
                      "raw_source_sha256": exponent_hash(raw),
                      "global_anchor_subtracted_after_gauge": global_anchor,
                      "source_construction_search_performed": False,
                      "tree_edge_added": False}})
    matrix = matrix.astype(np.uint16)
    matrix.setflags(write=False)
    return VerifiedSource(identifier, order, phase, exponent_hash(matrix), used, proof,
                          matrix, lambda i, j: int(matrix[i, j]))


def certificate():
    records = [load_source(row["id"]).metadata() for row in sources()]
    retained = [record for record in records if not record["proof"]["alphabet_gate"]["excluded"]]
    counts = Counter((record["order"], record["phase"]) for record in retained)
    return {"schema": "broad-published-source-acquisition-v1",
            "catalogue_url": CATALOGUE, "source_count": len(records),
            "raw_data_bytes": sum((DATA/f"BH-{n}-{q}.txt").stat().st_size for n, q, *_ in SPECS),
            "all_full_Grams_verified": True,
            "retained_fixed_representative_count": len(retained),
            "retained_ids": [record["id"] for record in retained],
            "retained_by_order_phase": [{"order": n, "phase": q, "count": count}
                                        for (n, q), count in sorted(counts.items())],
            "gauge_orbits_excluded_in_their_own_root_alphabet": len(records)-len(retained),
            "gauge_row_prefixes_visited": sum(r["proof"]["monomial_gauge_audit"]["visited_row_prefixes"] for r in records),
            "scope": "All524 published classes and all same-alphabet root monomial gauges; higher-alphabet gauges not covered",
            "sources": records, "new_Hadamard_order_constructed": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--report", type=Path)
    group.add_argument("--check", type=Path)
    args = parser.parse_args()
    report = certificate()
    if args.report is not None:
        args.report.write_text(json.dumps(report, indent=2)+"\n")
    elif json.loads(args.check.read_text()) != report:
        raise ValueError("the broad acquisition certificate differs from exact replay")
    print(json.dumps({key: report[key] for key in
                      ("source_count", "retained_fixed_representative_count",
                       "gauge_orbits_excluded_in_their_own_root_alphabet",
                       "gauge_row_prefixes_visited", "raw_data_bytes")}))


if __name__ == "__main__":
    main()

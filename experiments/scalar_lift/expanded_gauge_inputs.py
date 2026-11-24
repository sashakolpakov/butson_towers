"""Fixed published sixth-root inputs with exhaustive twelfth-root gauges.

This is a separate initial-input acquisition. Existing source registries,
certificates, and scalar-restriction checkpoints are not modified.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from functools import cache
from pathlib import Path

import numpy as np

from experiments.scalar_lift.broad_source_acquisition import (
    BASE,
    CATALOGUE,
    SPECS,
    dataset,
    root_gauge,
)
from experiments.scalar_lift.forward_butson_descent import exponent_hash, verify_source
from experiments.scalar_lift.forward_source_registry import VerifiedSource
from experiments.scalar_lift.scalar_restriction_obstructions import alphabet_obstruction

COUNTS = ((7, 2), (9, 17), (13, 436))


def sources():
    return [{"id": f"expanded-{n}-{index:03d}-q12", "order": n, "phase": 12,
             "catalogue_phase": 6, "catalogue_class_index": index,
             "family": "published-sixth-root-class-with-twelfth-root-initial-gauge"}
            for n, count in COUNTS for index in range(count)]


def projective_root_order(matrix, phase):
    cross = (matrix-matrix[:, :1]-matrix[:1, :]+matrix[0, 0]) % phase
    return phase//math.gcd(phase, *map(int, cross.flat))


@cache
def load_source(identifier):
    try:
        prefix, n, i, suffix = identifier.split("-")
        order, index = int(n), int(i)
        if prefix != "expanded" or suffix != "q12" or index < 0:
            raise ValueError
        raw = dataset(order, 6)[index]
    except (ValueError, IndexError, StopIteration) as error:
        raise KeyError(identifier) from error
    original_gram = verify_source(raw, 6)
    original_gauge = root_gauge(raw, 6)
    gauge = root_gauge(2*raw, 12)
    if original_gauge["exists"] and not gauge["exists"]:
        raise ArithmeticError("enlarging the initial root alphabet lost an existing gauge")
    matrix = 2*raw
    if gauge["exists"]:
        matrix = (matrix+np.asarray(gauge["row_shifts"])[:, None]
                  +np.asarray(gauge["column_shifts"])[None, :]) % 12
    global_anchor = int(matrix[0, 0])
    matrix = (matrix-global_anchor) % 12
    proof = verify_source(matrix, 12)
    used = tuple(map(int, np.unique(matrix)))
    gate = alphabet_obstruction(12, used)
    if not gate["complete_for_odd_integer_relations"] or gate["excluded"] == gauge["exists"]:
        raise ArithmeticError("the independent alphabet relation engine disagrees with the gauge CSP")
    rescued = not original_gauge["exists"] and gauge["exists"]
    effective = 12//math.gcd(12, *used)
    if rescued and effective != 12:
        raise ArithmeticError("a rescue did not actually use the larger literal entry field")
    projective = projective_root_order(matrix, 12)
    if projective != projective_root_order(raw, 6):
        raise ArithmeticError("a row/column root gauge changed the projective invariant")
    proof.update({
        "alphabet_gate": gate,
        "original_sixth_root_gauge_exists": original_gauge["exists"],
        "original_sixth_root_gauge_audit": original_gauge,
        "twelfth_root_gauge_audit": gauge,
        "rescued_by_larger_initial_root_gauge": rescued,
        "effective_literal_root_order": effective,
        "projective_root_order": projective,
        "not_monomially_equivalent_to_a_sixth_root_matrix": False,
        "acquisition": {
            "catalogue_url": CATALOGUE,
            "dataset_url": BASE+f"BH-{order}-6.txt",
            "dataset_sha256": next(spec[3] for spec in SPECS if spec[:2] == (order, 6)),
            "original_source_id": f"aalto-{order}-6-{index:03d}",
            "original_literal_source_sha256": exponent_hash(raw),
            "original_complete_Gram": original_gram,
            "embedding": "original exponent e maps to2e modulo12",
            "global_anchor_subtracted_after_gauge": global_anchor,
            "initial_representative_only": True,
            "source_from_smaller_matrices_constructed": False,
            "scalar_restriction_tree_edge_added": False,
        }})
    matrix = np.asarray(matrix, dtype=np.uint16)
    matrix.setflags(write=False)
    return VerifiedSource(identifier, order, 12, exponent_hash(matrix), used, proof,
                          matrix, lambda row, col: int(matrix[row, col]))


def admissible_sources():
    return [row for row in sources() if not load_source(row["id"]).proof["alphabet_gate"]["excluded"]]


def certificate():
    records = [load_source(row["id"]).metadata() for row in sources()]
    counts = Counter()
    for record in records:
        proof = record["proof"]
        status = ("originally_admissible" if proof["original_sixth_root_gauge_exists"] else
                  "rescued_by_twelfth_root_gauge" if proof["rescued_by_larger_initial_root_gauge"] else
                  "still_excluded_after_all_twelfth_root_gauges")
        counts[record["order"], status] += 1
    rescued = [record for record in records if record["proof"]["rescued_by_larger_initial_root_gauge"]]
    return {
        "schema": "expanded-initial-twelfth-root-gauge-acquisition-v1",
        "catalogue_url": CATALOGUE,
        "original_published_input_count": len(records),
        "originally_admissible_count": sum(count for (_n, status), count in counts.items()
                                           if status == "originally_admissible"),
        "rescued_input_count": len(rescued),
        "still_excluded_count": sum(count for (_n, status), count in counts.items()
                                   if status == "still_excluded_after_all_twelfth_root_gauges"),
        "total_admissible_twelfth_root_representatives": sum(not record["proof"]["alphabet_gate"]["excluded"]
                                                            for record in records),
        "counts_by_order_and_status": [{"order": n, "status": status, "count": count}
                                       for (n, status), count in sorted(counts.items())],
        "complete_CSP_row_prefixes_visited": sum(record["proof"]["twelfth_root_gauge_audit"]["visited_row_prefixes"]
                                                for record in records),
        "all_declared_cases_decided": True, "unknown_cases": 0,
        "rescued_source_ids": [record["id"] for record in rescued],
        "all_rescues_have_literal_entry_root_order12": all(record["proof"]["effective_literal_root_order"] == 12
                                                          for record in rescued),
        "rescue_creates_a_new_projective_Butson_class": False,
        "scope": "All sixth-root classes at orders7,9,13; all twelfth-root row/column gauges; no larger phases",
        "sources": records,
        "existing_registry_or_checkpoint_modified": False,
        "new_Hadamard_order_constructed": False,
    }


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
        raise ValueError("the expanded gauge certificate differs from exact replay")
    print(json.dumps({key: report[key] for key in
                      ("original_published_input_count", "originally_admissible_count", "rescued_input_count",
                       "still_excluded_count", "complete_CSP_row_prefixes_visited", "unknown_cases")}))


if __name__ == "__main__":
    main()

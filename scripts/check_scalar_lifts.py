#!/usr/bin/env python3
"""Replay exact Eisenstein lifts and audit their inputs against a dated order frontier."""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from dataclasses import asdict
from pathlib import Path

from sympy import Matrix, eye, kronecker_product

from hadamard.cgw import is_cgw_exponent_matrix
from hadamard.scalar_lift import (
    eisenstein_inert_prime_obstructions,
    find_unreal_gauge,
    is_real_weighing_matrix,
    lift_eisenstein_weighing,
    lift_unreal_cgw,
    restrict_eisenstein_scalars,
    unreal_components,
    unreal_pair_conditions,
)

ROOT = Path(__file__).resolve().parents[1]
SNAPSHOT = ROOT / "experiments/scalar_lift/frontier_snapshot.json"


def _load(relative):
    return json.loads((ROOT / relative).read_text())


def _example(name, source, phase, weight):
    source = tuple(tuple(row) for row in source)
    if not is_cgw_exponent_matrix(source, phase, weight):
        raise ValueError(f"invalid cyclotomic source: {name}")
    restriction = Matrix(restrict_eisenstein_scalars(source, phase))
    metric = kronecker_product(eye(len(source)), Matrix([[2, 1], [1, 2]]))
    if restriction * metric * restriction.T != weight * metric:
        raise ValueError("restriction failed its exact metric identity")
    gauge = find_unreal_gauge(source, phase)
    sparse = lift_eisenstein_weighing(source, weight, phase)
    if not is_real_weighing_matrix(sparse, 2 * weight):
        raise ValueError("the universal weighing lift failed its Euclidean Gram check")
    result = {
        "name": name,
        "source_phase": phase,
        "source_order": len(source),
        "source_weight": weight,
        "source_exponents": source,
        "source_exact_cyclotomic_gram_valid": True,
        "restriction_dual_trace_gram_valid": True,
        "restriction_euclidean_gram_valid": restriction * restriction.T
        == weight * eye(2 * len(source)),
        "unreal_gauge": asdict(gauge),
        "universal_sparse_lift": {
            "order": len(sparse),
            "weight": 2 * weight,
            "matrix": sparse,
            "both_euclidean_grams_valid": True,
            "dense_hadamard_completion_established": False,
        },
        "new_real_order_claimed": False,
    }
    if gauge.matrix is not None:
        first, second = unreal_components(gauge.matrix)
        lift = lift_unreal_cgw(gauge.matrix, weight)
        direct = Matrix(lift)
        size = 4 * len(source)
        if not (
            is_cgw_exponent_matrix(gauge.matrix, 6, weight)
            and is_real_weighing_matrix(lift, 4 * weight)
            and direct * direct.T == 4 * weight * eye(size)
            and direct.T * direct == 4 * weight * eye(size)
        ):
            raise ValueError("independent exact checks rejected a real lift")
        result.update(
            {
                "component_pair_conditions": unreal_pair_conditions(first, second, weight),
                "real_order": size,
                "real_weight": 4 * weight,
                "real_matrix": lift,
                "real_gram_both_directions_valid": True,
                "is_dense_real_hadamard": weight == len(source),
            }
        )
    return result


def _target_audit():
    snapshot = json.loads(SNAPSHOT.read_text())
    orders = sorted(
        group["multiplier"] * prime
        for group in snapshot["survey"]["groups"]
        for prime in group["primes"]
    )
    if len(orders) != 95 or len(set(orders)) != 95:
        raise ValueError("the historical frontier transcription must have 95 distinct orders")
    reported = snapshot["later_reported_constructions_intersecting_snapshot"]["orders"]
    rows = []
    for order in orders:
        superseding = [
            seed
            for seed in reported
            if order % seed == 0 and (order // seed) & (order // seed - 1) == 0
        ]
        source_order = order // 4
        obstructions = eisenstein_inert_prime_obstructions(source_order, source_order)
        rows.append(
            {
                "target_order": order,
                "target_status": "later_construction_reported"
                if superseding
                else "snapshot_candidate",
                "later_reported_seed_orders": superseding,
                "unreal_sextic_source_order": source_order,
                "unreal_sextic_inert_prime_obstructions": obstructions,
                "unreal_sextic_route": "excluded_by_norm" if obstructions else "input_required",
                "pure_cubic_source_obstructed_by_length": source_order > 1
                and source_order % 3 != 0,
                "pure_cubic_unreal_lift_obstructed_by_weight": source_order > 4,
                "quaternary_source_order": order // 2,
                "quaternary_route": "input_required",
                "source_matrix_constructed": False,
                "new_real_order_certified": False,
            }
        )
    return {
        "snapshot_sha256": hashlib.sha256(SNAPSHOT.read_bytes()).hexdigest(),
        "scope": snapshot["scope"],
        "priority_status": snapshot["priority_status"],
        "historical_order_count": len(rows),
        "later_reported_or_doubled_count": sum(
            bool(row["later_reported_seed_orders"]) for row in rows
        ),
        "sextic_route_counts_over_historical_snapshot": dict(
            Counter(row["unreal_sextic_route"] for row in rows)
        ),
        "sextic_route_counts_after_later_report_suppression": dict(
            Counter(
                row["unreal_sextic_route"]
                for row in rows
                if row["target_status"] == "snapshot_candidate"
            )
        ),
        "orders": rows,
    }


def report():
    cases = [
        _example("Fourier F3", [[i * j % 3 for j in range(3)] for i in range(3)], 3, 3),
        _example(
            "CGW(5,4;3)",
            _load("certificates/cgw5-4-3/classification.json")["representative_exponent_matrix"],
            3,
            4,
        ),
    ]
    for index, item in enumerate(_load("certificates/finite-butson/bh7-18.json")["classes"][:2]):
        matrix = item["representative"]
        if any(e % 3 for row in matrix for e in row):
            raise ValueError("the selected representative does not lie in mu_6")
        cases.append(
            _example(
                f"BH(7,6) representative {index + 1}",
                [[e // 3 for e in row] for row in matrix],
                6,
                7,
            )
        )
    cases.append(
        _example(
            "Fixed CGW(9,7;6) witness",
            _load("certificates/cgw9-7-6/existence.json")["representative_exponent_matrix"],
            6,
            7,
        )
    )
    cases.append(
        _example(
            "Fourier over F3 squared",
            [[((i // 3) * (j // 3) + (i % 3) * (j % 3)) % 3 for j in range(9)] for i in range(9)],
            3,
            9,
        )
    )
    return {
        "schema": "hadamard.eisenstein-real-lift-pilot.v1",
        "scope": "Exact sparse and known partial lifts, local witness gauges, and necessary conditions for historical open orders",
        "novel_construction_claimed": False,
        "new_real_orders_certified": [],
        "examples": cases,
        "frontier": _target_audit(),
        "next_mathematical_requirement": "Construct a source at a surviving gap, or prove a new sparse completion or bordered lift; passing norm filters is insufficient",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", type=Path)
    args = parser.parse_args()
    payload = report()
    encoded = json.dumps(payload, sort_keys=True, indent=2) + "\n"
    if args.check and json.loads(args.check.read_text()) != json.loads(encoded):
        raise ValueError("the saved pilot does not match exact replay")
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded)
    print(
        json.dumps(
            {
                "examples": [
                    {
                        k: item[k]
                        for k in ("name", "unreal_gauge", "real_order", "real_weight")
                        if k in item
                    }
                    for item in payload["examples"]
                ],
                "frontier_counts": {
                    k: v
                    for k, v in payload["frontier"].items()
                    if k.endswith("count") or k.endswith("snapshot")
                },
                "first_targets": [
                    row
                    for row in payload["frontier"]["orders"]
                    if row["target_status"] == "snapshot_candidate"
                ][:4],
                "new_real_orders_certified": [],
                "saved_pilot_replayed": bool(args.check),
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Work backwards from candidate real orders to certified torsion points and lift equations."""

from __future__ import annotations

import argparse
import hashlib
import json
from functools import reduce
from itertools import product
from math import gcd
from pathlib import Path

from sympy import Matrix, eye, isprime, totient

from experiments.scalar_lift.inverse_blocks import (
    BH13_SEED,
    SEED,
    binary_digest,
    hall_variant,
    known_h668,
    known_order_thirteen_source,
    known_order_twelve_pair,
    pair_operator_arithmetic,
    prime_block_operator_requirements,
    triple_unreal_source,
)
from hadamard.cgw import is_cgw_exponent_matrix
from hadamard.inverse_lift import (
    is_sign_hadamard,
    lift_unreal_sixth_from_pair,
    sixth_root_block_pair_conditions,
)
from hadamard.scalar_lift import SKEW_CONFERENCE_4

ROOT = Path(__file__).resolve().parents[1]
SNAPSHOT = ROOT / "experiments/scalar_lift/frontier_snapshot.json"
DEFAULT_REPORT = ROOT / "certificates/scalar-lift/inverse-dense-lifts.json"


def _snapshot_candidates():
    data = json.loads(SNAPSHOT.read_text())
    historical = sorted(g["multiplier"] * p for g in data["survey"]["groups"] for p in g["primes"])
    reported = data["later_reported_constructions_intersecting_snapshot"]["orders"]
    candidates = [
        n
        for n in historical
        if not any(n % s == 0 and (n // s) & (n // s - 1) == 0 for s in reported)
    ]
    return historical, candidates


def _avoid_third_subgroup(source, phase):
    """Exact monomial gauge: e_ij+r_i+c_j != 0 mod 3; no Gram alteration."""
    n = len(source)
    for checked, tail in enumerate(product(range(3), repeat=n - 1), 1):
        rows = (0,) + tail
        columns = []
        for j in range(n):
            forbidden = {(-source[i][j] - rows[i]) % 3 for i in range(n)}
            available = next((x for x in range(3) if x not in forbidden), None)
            if available is None:
                break
            columns.append(available)
        if len(columns) == n:
            result = tuple(
                tuple((source[i][j] + rows[i] + columns[j]) % phase for j in range(n))
                for i in range(n)
            )
            if not is_cgw_exponent_matrix(result, phase, n):
                raise ValueError("the monomial gauge failed to preserve the source Gram identity")
            return result, rows, tuple(columns), checked
    return None


def _order_seven_inventory(block_order):
    # Keep the first fixed-q catalogue exhibiting each effective alphabet.
    # Counts are fixed-alphabet classes, not an arbitrary-torsion completeness claim.
    inventory = {}
    compatible = []
    for ambient in (18, 24, 30, 36, 42, 48):
        relative = f"certificates/finite-butson/bh7-{ambient}.json"
        payload = json.loads((ROOT / relative).read_text())
        if not payload["classification_complete"]:
            raise ValueError("use only a completed fixed-alphabet catalogue for this inventory")
        groups = {}
        for index, cls in enumerate(payload["classes"]):
            matrix = cls["representative"]
            if any(matrix[0]) or any(row[0] for row in matrix):
                raise ValueError("source field degrees require dephased representatives")
            common = reduce(gcd, (x for row in matrix for x in row), ambient)
            effective = ambient // common
            groups.setdefault(effective, []).append((index, matrix, common))
        for effective, points in groups.items():
            if effective in inventory:
                continue
            degree = int(totient(effective))
            admitted = block_order % degree == 0
            inventory[effective] = {
                "effective_phase": effective,
                "cyclotomic_degree": degree,
                "ambient_phase": ambient,
                "fixed_alphabet_class_count": len(points),
                "representation_degree_divides_block_order": admitted,
                "certificate": relative,
            }
            if not admitted:
                continue
            if effective not in (6, 12):
                raise ValueError(
                    "unexpected compatible alphabet; extend the inverse audit explicitly"
                )
            for index, matrix, common in points:
                source = tuple(tuple(x // common for x in row) for row in matrix)
                if not is_cgw_exponent_matrix(source, effective, 7):
                    raise ValueError("invalid stored source torsion point")
                gauge = _avoid_third_subgroup(source, effective)
                if gauge is None:
                    raise ValueError("the compatible point has no required entry-domain gauge")
                adjusted, rows, columns, checked = gauge
                compatible.append(
                    {
                        "name": f"bh7-{ambient}:class-{index}",
                        "certificate": relative,
                        "certificate_sha256": hashlib.sha256(
                            (ROOT / relative).read_bytes()
                        ).hexdigest(),
                        "class_index": index,
                        "phase": effective,
                        "source_order": 7,
                        "row_phase_shifts": rows,
                        "column_phase_shifts": columns,
                        "gauge_row_cases_checked": checked,
                        "exponents": adjusted,
                        "source_gram_exactly_verified": True,
                        "entry_domain": "mu_6 minus mu_2" if effective == 6 else "mu_12 minus mu_4",
                        "source_entries_use_entire_domain": len(
                            {x for row in adjusted for x in row}
                        )
                        == (4 if effective == 6 else 8),
                    }
                )
    return sorted(inventory.values(), key=lambda x: x["effective_phase"]), compatible


def build_report():
    historical, candidates = _snapshot_candidates()
    if 4676 not in candidates:
        raise ValueError("the selected target is no longer a candidate in the saved snapshot")
    source_three = tuple(tuple(4 if i == j else 2 for j in range(3)) for i in range(3))
    if not is_cgw_exponent_matrix(source_three, 6, 3):
        raise ValueError("invalid order-three source")
    inventory, points = _order_seven_inventory(668)
    sixth = [p for p in points if p["phase"] == 6]
    source_nine = triple_unreal_source(source_three)
    source_thirteen = known_order_thirteen_source()
    seed13 = json.loads(BH13_SEED.read_text())
    extra_sources = (
        {
            "name": "unreal-order-nine-by-tripling",
            "source_order": 9,
            "phase": 6,
            "construction": "A'=B tensor J3, B'=A tensor I3+B tensor C3; source order 3",
            "known_family_reference": "https://link.springer.com/article/10.1007/s10623-015-0045-y",
            "exponents": source_nine,
            "source_gram_exactly_verified": True,
            "entry_domain": "mu_6 minus mu_2",
        },
        {
            "name": "catalogue-M13A-unreal-gauge",
            "source_order": 13,
            "phase": 6,
            "source_url": seed13["source_url"],
            "seed": str(BH13_SEED.relative_to(ROOT)),
            "seed_sha256": hashlib.sha256(BH13_SEED.read_bytes()).hexdigest(),
            "row_phase_shifts": seed13["unreal_row_phase_shifts"],
            "column_phase_shifts": seed13["unreal_column_phase_shifts"],
            "exponents": source_thirteen,
            "source_gram_exactly_verified": True,
            "entry_domain": "mu_6 minus mu_2",
        },
    )
    routes = sorted(
        (
            {
                "target_order": n,
                "source_order": s,
                "block_order": n // s,
                "source_available": True,
                "block_pair_constructed": False,
            }
            for n in candidates
            for s in (3, 7, 9, 13)
            if n % s == 0 and (n // s) % 4 == 0
        ),
        key=lambda x: (x["block_order"], x["target_order"]),
    )
    if routes[0] != {
        "target_order": 4676,
        "source_order": 7,
        "block_order": 668,
        "source_available": True,
        "block_pair_constructed": False,
    }:
        raise ValueError("the minimum-block route changed; review the target selection")

    base = known_h668()
    switched = hall_variant(base, 1)
    seed = json.loads(SEED.read_text())
    if binary_digest(switched) != seed["first_hall_switch_binary_sha256"]:
        raise ValueError("the switched matrix differs from its published hash")
    angle = sixth_root_block_pair_conditions(base, switched)
    if not (angle["first_is_hadamard"] and angle["second_is_hadamard"]) or angle["valid"]:
        raise ValueError("unexpected result for the imported known-block check")
    # Hall switches flip an even number of entries in every row. This
    # excludes every pair in this family by the necessary distance d/4=167.
    parities = set()
    variant_hashes = []
    for mask in range(16):
        h = hall_variant(base, mask)
        if not is_sign_hadamard(h):
            raise ValueError("a claimed Hall variant is not Hadamard")
        parities.update(sum(x < 0 for x in row) % 2 for row in h)
        variant_hashes.append(binary_digest(h))
    if len(parities) != 1:
        raise ValueError("Hall-family row-sign parity is not constant")

    c = Matrix(SKEW_CONFERENCE_4)
    p4, q4 = eye(4) + c, -eye(4) + c
    pair4 = sixth_root_block_pair_conditions(p4.tolist(), q4.tolist())
    if not pair4["valid"]:
        raise ValueError("the inverse-block control pair failed")
    controls = []
    for name, source in (
        ("order-three source", source_three),
        *((p["name"], p["exponents"]) for p in sixth),
        *((p["name"], p["exponents"]) for p in extra_sources),
    ):
        output = lift_unreal_sixth_from_pair(source, p4.tolist(), q4.tolist())
        h = Matrix(output)
        if not (set(h) <= {-1, 1} and h * h.T == h.T * h == len(output) * eye(len(output))):
            raise ValueError("the inverse-block control failed its independent integer Gram check")
        controls.append(
            {
                "name": name,
                "real_order": len(output),
                "dense_integer_grams_valid": True,
                "known_order_control": True,
            }
        )

    p12, q12 = known_order_twelve_pair()
    pair12 = sixth_root_block_pair_conditions(p12, q12)
    arithmetic12 = pair_operator_arithmetic(p12, q12)
    bound12 = prime_block_operator_requirements(12)
    if not pair12["valid"] or arithmetic12["zero_entries"]:
        raise ValueError("the dense denominator-six operator control failed")
    if any(
        counts != {"2": 3, "6": 9} for counts in arithmetic12["reduced_denominator_counts_per_row"]
    ):
        raise ValueError("the order-twelve pair does not attain the proved arithmetic bound")
    r12 = Matrix(q12) * Matrix(p12).T / 12
    if not (
        r12 * r12.T == eye(12)
        and r12 + r12.T == eye(12)
        and r12**2 - r12 + eye(12) == Matrix.zeros(12)
    ):
        raise ValueError("the denominator-six control fails its exact field equations")
    for name, source in (
        ("order-three source", source_three),
        (sixth[0]["name"], sixth[0]["exponents"]),
        *((p["name"], p["exponents"]) for p in extra_sources),
    ):
        output = lift_unreal_sixth_from_pair(source, p12, q12)
        h = Matrix(output)
        if not (set(h) <= {-1, 1} and h * h.T == len(output) * eye(len(output))):
            raise ValueError("the denominator-six scalar transfer control failed")
        controls.append(
            {
                "name": name,
                "block_order": 12,
                "real_order": len(output),
                "dense_integer_grams_valid": True,
                "known_order_control": True,
            }
        )

    return {
        "scope": "Target-first inverse lift for exact dense sign outputs",
        "snapshot_sha256": hashlib.sha256(SNAPSHOT.read_bytes()).hexdigest(),
        "historical_order_count": len(historical),
        "candidate_count_after_reported_exclusions": len(candidates),
        "candidate_status_scope": "dated snapshot, not an independent current priority claim",
        "available_source_routes_sorted_by_block_order": routes,
        "verified_unreal_source_orders_used_in_inverse_map": [3, 7, 9, 13],
        "shared_required_block_pairs": [
            {
                "block_order": d,
                "candidate_target_orders": sorted(
                    route["target_order"] for route in routes if route["block_order"] == d
                ),
                "source_orders": sorted(
                    route["source_order"] for route in routes if route["block_order"] == d
                ),
                "pair_constructed": False,
            }
            for d in sorted({route["block_order"] for route in routes})
        ],
        "prime_block_operator_requirements": [
            prime_block_operator_requirements(d)
            for d in sorted({route["block_order"] for route in routes})
            if d % 4 == 0 and d // 4 > 2 and isprime(d // 4)
        ],
        "selected_target": {
            "real_order": 4676,
            "factorization": "7 * 668",
            "source_order": 7,
            "block_order": 668,
            "selected_source_name": sixth[0]["name"],
            "source_matrix_constructed": True,
            "block_matrix_of_required_order_verified": True,
            "required_block_pair_constructed": False,
            "target_matrix_constructed": False,
            "missing_identity": "P Q^T + Q P^T = 668 I, with P,Q real Hadamard blocks",
            "required_operator_arithmetic": prime_block_operator_requirements(668),
        },
        "order_seven_alphabet_inventory": inventory,
        "compatible_torsion_points": points,
        "additional_verified_unreal_sixth_root_sources": extra_sources,
        "twelfth_root_branch_requires_at_least_the_sixth_root_pair": True,
        "inverse_theorem": {
            "necessary_and_sufficient_block_conditions": "PP^T=QQ^T=dI and PQ^T+QP^T=dI, with all entries +/-1",
            "representation": "rho(zeta_6)=Q P^T/d",
            "frame": "Z=P-Q",
            "substitution": {"1": "P", "2": "Q", "4": "-P", "5": "-Q"},
            "output": "real Hadamard of order dn for every admitted BH(n,6)",
        },
        "verified_known_668_blocks": {
            "source_url": seed["source_url"],
            "source_data_sha256": hashlib.sha256(SEED.read_bytes()).hexdigest(),
            "base_binary_sha256": binary_digest(base),
            "switched_binary_sha256": binary_digest(switched),
            "published_hashes_match": True,
            "pair_conditions": angle,
            "hall_variant_binary_sha256": variant_hashes,
            "all_sixteen_variants_exactly_hadamard": True,
            "row_negative_count_parity": next(iter(parities)),
            "all_labelled_hall_pairs_excluded_by_row_parity": True,
            "minimum_one_sided_column_sign_changes_necessary_without_column_permutation": 163,
            "exclusion_scope": "also independent signed row permutations and shared column operations; independent column operations are not ruled out",
        },
        "control_block_pair": {
            "P": [[int(x) for x in row] for row in p4.tolist()],
            "Q": [[int(x) for x in row] for row in q4.tolist()],
            "conditions": pair4,
        },
        "denominator_six_control_block_pair": {
            "P": p12,
            "Q": q12,
            "conditions": pair12,
            "operator_arithmetic": arithmetic12,
            "prime_block_bound": bound12,
            "arithmetic_bound_attained_in_every_row": True,
            "field_equations_exactly_verified": True,
        },
        "known_order_controls": controls,
        "new_real_hadamard_order_claimed": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", type=Path)
    args = parser.parse_args()
    if args.output is not None and args.check is not None:
        parser.error("use either --output or --check")
    report = build_report()
    encoded = json.dumps(report, indent=2) + "\n"
    if args.check is not None:
        if json.loads(encoded) != json.loads(args.check.read_text()):
            raise SystemExit("inverse-lift report replay differs")
    else:
        path = args.output or DEFAULT_REPORT
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(encoded)
    print(
        json.dumps(
            {
                "target": report["selected_target"],
                "compatible_torsion_points": len(report["compatible_torsion_points"]),
                "known_block_pair_valid": report["verified_known_668_blocks"]["pair_conditions"][
                    "valid"
                ],
                "new_order_claimed": False,
            }
        )
    )


if __name__ == "__main__":
    main()

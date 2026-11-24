"""Completed scalar endpoint recipes from the known unreal order-13 family.

The retained finite executor checks sources of orders 39, 117 and 351 and
reuses the certified width-2524 T631 frame. It never allocates a dense endpoint.
Orders 98436 and 295308 are new to the best of the author's knowledge; order
885924 is known by Paley's construction. No family-wide novelty is claimed.
"""
from __future__ import annotations

import argparse
from functools import cache
import hashlib
import json
import math
import operator
from pathlib import Path

import numpy as np

from experiments.scalar_lift.forward_source_registry import load_source
from experiments.scalar_lift.reuse_hadamard_7572_frame import (
    NOVELTY_CONVENTION, ReusedEndpoint, verified_frames,
)

ROOT = Path(__file__).resolve().parents[2]
SEED = ROOT / "experiments/scalar_lift/bh13_6_seed.json"
BUNDLE = ROOT / "artifacts/hadamard-7572"
SOURCE_ORDERS = (13, 39, 117, 351)
WIDTH = 2524
C = np.array(((0, 1, -1), (-1, 0, 1), (1, -1, 0)), dtype=np.int64)
J = np.ones((3, 3), dtype=np.int64)
I3 = np.eye(3, dtype=np.int64)
ROOTS = np.array(((1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)), dtype=np.int64)
UNREAL_FROM_SIGNS = np.array(((4, 2), (5, 1)), dtype=np.int64)


def require(value, message):
    if not value:
        raise ValueError(message)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def exponent_hash(matrix):
    return digest(np.asarray(matrix, dtype="<i8").tobytes())


def source_gram(exponents):
    """Check every row and column Gram entry in the exact basis (1,zeta6)."""
    e = np.asarray(exponents)
    require(e.ndim == 2 and e.shape[0] == e.shape[1]
            and np.issubdtype(e.dtype, np.integer)
            and np.all((0 <= e) & (e < 6)), "invalid sixth-root source")
    n = len(e)
    expected = n * np.eye(n, dtype=np.int64)
    for matrix in (e, e.T):
        a, b = ROOTS[matrix, 0], ROOTS[matrix, 1]
        real = a @ (a + b).T + b @ b.T
        imaginary = b @ a.T - a @ b.T
        require(np.array_equal(real, expected) and not np.any(imaginary),
                "complete exact sixth-root Gram failed")
    return {"order": n, "phase": 6, "source_gram_scale": n,
            "all_row_Gram_entries_checked": n * n,
            "all_column_Gram_entries_checked": n * n,
            "coefficient_basis": ["1", "zeta_6"],
            "exact_row_and_column_Gram": f"{n} I",
            "exponent_sha256": exponent_hash(e)}


def components(unreal):
    require(set(map(int, np.unique(unreal))) <= {1, 2, 4, 5},
            "unreal sixth-root entries required")
    a = np.array((0, 1, -1, 0, -1, 1), dtype=np.int64)[unreal]
    b = np.array((0, 1, 1, 0, -1, -1), dtype=np.int64)[unreal]
    return a, b


def triple(unreal):
    """U=(A+sqrt(-3)B)/2 -> (B tensor J+sqrt(-3)(A tensor I+B tensor C))/2."""
    a, b = components(unreal)
    aa, bb = np.kron(b, J), np.kron(a, I3) + np.kron(b, C)
    require(np.all(np.abs(aa) == 1) and np.all(np.abs(bb) == 1),
            "tripling left the unreal entry alphabet")
    return UNREAL_FROM_SIGNS[(aa + 1) // 2, (bb + 1) // 2]


@cache
def seed():
    record = json.loads(SEED.read_text())
    raw = np.asarray(record["exponents"], dtype=np.int64)
    raw_proof = source_gram(raw)
    rows = np.asarray(record["unreal_row_phase_shifts"], dtype=np.int64)
    columns = np.asarray(record["unreal_column_phase_shifts"], dtype=np.int64)
    unreal = (raw + rows[:, None] + columns[None, :]) % 6
    components(unreal)
    unreal_proof = source_gram(unreal)
    require(np.array_equal((unreal - 1) % 6, load_source("unreal13-q6").exponents),
            "seed differs from the retained complete source")
    require(np.array_equal(C.T, -C)
            and np.array_equal(C @ C.T, 3 * I3 - J)
            and not np.any(C @ J) and not np.any(J @ C)
            and np.array_equal(J @ J, 3 * J), "tripling tensor identities failed")
    proof = {"raw_seed_file": str(SEED.relative_to(ROOT)),
             "raw_seed_file_sha256": digest(SEED.read_bytes()),
             "primary_literal_source_url": record["source_url"],
             "catalogue_url": record["catalogue_url"],
             "catalogue_attribution": record["catalogue_attribution"],
             "raw_exponents": raw.tolist(), "raw_full_Gram": raw_proof,
             "unreal_row_shifts_mod6": rows.tolist(),
             "unreal_column_shifts_mod6": columns.tolist(),
             "unreal_exponent_alphabet": sorted(map(int, np.unique(unreal))),
             "unreal_full_Gram": unreal_proof}
    unreal.setflags(write=False)
    return unreal, proof


@cache
def verified_source(order):
    require(not isinstance(order, (bool, np.bool_)), "Boolean source order")
    order = operator.index(order)
    require(order in SOURCE_ORDERS, "retained executor supports source orders13,39,117,351")
    k = SOURCE_ORDERS.index(order)
    unreal, _ = seed()
    stages = []
    for stage in range(1, k + 1):
        unreal = triple(unreal)
        stages.append({"tripling": stage, **source_gram(unreal)})
    e = (unreal - 1) % 6
    require(set(map(int, np.unique(e))) == {0, 1, 3, 4}, "frame palette failed")
    proof = source_gram(e)
    registry = load_source(f"cubic13-q24-t{k}")
    require(np.array_equal(registry.exponents, (4 * unreal + 2) % 24),
            "literal tripling differs from the complete retained phase24 source")
    require(np.array_equal(4 * e, (registry.exponents.astype(np.int64) - 6) % 24),
            "phase24 gauge and sixth-root extraction disagree")
    proof.update({"construction": "known unreal seed13 followed by exact tensor triplings",
                  "tripling_count": k, "all_intermediate_full_Grams": stages,
                  "raw_unreal_exponent_alphabet": sorted(map(int, np.unique(unreal))),
                  "global_gauge_mod6": 5,
                  "actual_exponent_alphabet": [0, 1, 3, 4],
                  "retained_registry_id": registry.id,
                  "retained_registry_phase": registry.phase,
                  "retained_registry_hash": registry.source_hash,
                  "registry_to_source_gauge": "e6=((e24-6) mod24)/4",
                  "complete_source_constructed_and_verified": True,
                  "source_priority_claimed": False})
    e.setflags(write=False)
    return e, proof


def endpoint(order):
    e, _ = verified_source(order)
    p, q, _ = verified_frames()
    return ReusedEndpoint(f"tripled-unreal13-{order}", e, p, q)


def prime(n):
    return n > 1 and all(n % d for d in range(2, math.isqrt(n) + 1))


def row_controls(matrix):
    selected = [0, 1, WIDTH - 1, WIDTH, matrix.order - 1]
    packed_rows = []
    for i in selected:
        row = matrix.row(i)
        require(row.shape == (matrix.order,) and np.all(np.abs(row) == 1),
                "literal endpoint row entry test failed")
        for j in (0, 1, WIDTH - 1, WIDTH, matrix.order - 1):
            require(matrix.entry(i, j) == int(row[j]), "entry and row recipes disagree")
        packed_rows.append(int.from_bytes(np.packbits(row < 0, bitorder="little").tobytes(), "little"))
    for i, row in enumerate(packed_rows):
        for j in range(i):
            require((row ^ packed_rows[j]).bit_count() == matrix.order // 2,
                    "sampled literal endpoint row Gram failed")
    return {"row_indices": selected, "all_selected_rows_have_only_signs": True,
            "selected_pairwise_integer_Grams_verified": True,
            "scope": "literal interface controls; the complete endpoint Gram is proved algebraically"}


def certificate():
    _, seed_proof = seed()
    _, _, frame_proof = verified_frames()
    original = json.loads((BUNDLE / "certificate.json").read_text())
    records = []
    for n in SOURCE_ORDERS[1:]:
        matrix = endpoint(n)
        _, source_proof = verified_source(n)
        known = n == 351
        if known:
            require(prime(matrix.order - 1) and (matrix.order - 1) % 4 == 3,
                    "Paley known-order control failed")
        records.append({"source": source_proof, "endpoint_order": matrix.order,
                        "source_gram_scale": n, "frame_gram_scale": WIDTH,
                        "endpoint_gram_scale": n * WIDTH,
                        "actual_endpoint_alphabet": [-1, 1],
                        "literal_entry_and_row_recipe_available": True,
                        "endpoint_Gram_proved_exactly": True,
                        "proof": "B(z)=rho(z)P; sum_k B(h_ik)B(h_jk)^T=2524*rho(sum_k h_ik*conjugate(h_jk))=2524*n*delta_ij*I",
                        "dense_endpoint_allocated": False,
                        "dense_endpoint_Gram_multiplication_performed": False,
                        "literal_row_controls": row_controls(matrix),
                        "ordinary_real_Hadamard_Kronecker_product": False,
                        "Kronecker_reason": "order>4 and exactly two factors2; two nontrivial real Hadamard orders cannot have this product",
                        "known_order": known,
                        "prior_construction_identified": known,
                        "novelty_status": "known" if known else "new_to_best_of_authors_knowledge",
                        "known_order_construction": {"method": "Paley I", "prime_power": 885923}
                            if known else None,
                        "new_Hadamard_order_claimed": not known,
                        "new_order_claim": None if known else
                            f"To the best of the author's knowledge, {matrix.order} is a new real Hadamard order.",
                        "new_equivalence_class_claimed": False})
    return {"schema": "completed-tripled-unreal-scalar-endpoint-recipes-v1",
            "seed": seed_proof,
            "known_source_family_primary_reference": {
                "authors": "B. Compton, R. Craigen, W. de Launey",
                "title": "Unreal BH(n,6)'s and Hadamard matrices",
                "journal": "Designs, Codes and Cryptography79(2016),219-229",
                "url": "https://link.springer.com/article/10.1007/s10623-015-0045-y",
                "publisher_stated_family": "n=3^k*h*m, k>=0, h a real Hadamard order, m in{1,7,10,13}",
                "specific_family_used_here": "unreal BH(13*3^k,6), then global sixth-root gauge-1",
                "method_priority_claimed": False},
            "independent_tripling_proof": {
                "input": "U=(A+sqrt(-3)B)/2, A,B sign matrices; AA^T+3BB^T=4nI and AB^T=BA^T",
                "C3": C.tolist(), "J3": J.tolist(),
                "tensor_identities": "C^T=-C, CC^T=3I-J, CJ=JC=0, J^2=3J",
                "recurrence": "A'=B tensor J; B'=A tensor I+B tensor C",
                "output_entries": "A',B' are signs, so(A'+sqrt(-3)B')/2 has only unreal sixth-root entries",
                "Gram_expansion": "A'A'^T+3B'B'^T=3(AA^T+3BB^T) tensor I=12nI",
                "cross_expansion": "A'B'^T-B'A'^T=(BA^T-AB^T) tensor J=0",
                "family_rule": "all source orders13*3^k; endpoint orders32812*3^k, k>=0",
                "actual_complete_source_Gram_executions": [39, 117, 351]},
            "fixed_frame": frame_proof,
            "actual_scalar_action_entry_alphabet": original["scalar_restriction"]["actual_R_entry_alphabet"],
            "actual_intermediate_entry_alphabet": "{0,1,-1} union the recorded rational R-entry alphabet and its negatives",
            "actual_used_scalar_images": {"0": "P", "1": "Q", "3": "-P", "4": "-Q"},
            "immutable_frame_certificate": "artifacts/hadamard-7572/certificate.json",
            "immutable_frame_certificate_sha256": digest((BUNDLE / "certificate.json").read_bytes()),
            "complete_endpoint_recipes": records,
            "novelty_convention": dict(NOVELTY_CONVENTION),
            "new_real_orders_claimed_to_best_of_authors_knowledge": [98436, 295308],
            "novelty_scope": "Completed real Hadamard orders98436 and295308 are new to the best of the author's knowledge. Order885924 is known by Paley I at prime885923. No family-wide order novelty, source-method priority or equivalence-class novelty is claimed.",
            "implementation_sha256": digest(Path(__file__).read_bytes())}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = certificate()
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"completed_endpoint_orders": [x["endpoint_order"] for x in report["complete_endpoint_recipes"]],
                      "source_orders": list(SOURCE_ORDERS[1:]), "new_order_claimed": True,
                      "new_real_orders_claimed_to_best_of_authors_knowledge": report["new_real_orders_claimed_to_best_of_authors_knowledge"],
                      "known_Paley_order": 885924}))


if __name__ == "__main__":
    main()

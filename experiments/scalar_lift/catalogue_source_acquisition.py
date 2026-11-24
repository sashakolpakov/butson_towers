"""Fixed published Butson inputs for subsequent scalar restrictions only.

This imports literal data, replays exact Grams, and audits alphabet gauges.
It neither constructs a larger source nor changes an ongoing restriction tree.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from functools import cache
from pathlib import Path

import numpy as np

from experiments.scalar_lift.forward_butson_descent import exponent_hash, verify_source
from experiments.scalar_lift.forward_source_registry import VerifiedSource
from experiments.scalar_lift.scalar_restriction_obstructions import alphabet_obstruction
from hadamard.finite_equivalence import canonical_finite_butson_key

ROOT = Path(__file__).resolve().parents[2]
DATA = Path(__file__).with_name("catalogue_inputs")
SOURCE_IDS = ("catalogue-p7-q6", "catalogue-p7-q12", "catalogue-d14-q4")


def data(name):
    return json.loads((DATA / f"{name}.json").read_text())


def cubic_gauge(matrix):
    """Exhaust all ternary row gauges, stopping at a witness or a full cover.

    A column gauge can omit one cubic coset iff its current row colours
    occupy at most two cosets. Fixing the first row gauge loses no solution.
    Signs do not affect this criterion; permutations do not affect feasibility.
    """
    matrix = np.asarray(matrix, dtype=np.int64) % 3
    n = len(matrix)
    nodes = 0

    def visit(rows, masks):
        nonlocal nodes
        nodes += 1
        if len(rows) == n:
            columns = [next(c for c in range(3)
                            if all((int(matrix[i, j]) + rows[i] + c) % 3 != 2
                                   for i in range(n))) for j in range(n)]
            return rows, columns
        i = len(rows)
        for shift in ((0,) if i == 0 else range(3)):
            updated = [mask | (1 << ((int(matrix[i, j]) + shift) % 3))
                       for j, mask in enumerate(masks)]
            if 7 in updated:
                continue
            witness = visit(rows + [shift], updated)
            if witness is not None:
                return witness
        return None

    witness = visit([], [0] * n)
    return {"exists": witness is not None, "visited_row_prefixes": nodes,
            "row_shifts": witness[0] if witness else None,
            "column_shifts": witness[1] if witness else None,
            "scope": "All sixth-root row/column gauges and permutations of this fixed matrix"}


def cross_ratio(matrix, phase):
    residues = (matrix - matrix[:, :1] - matrix[:1, :] + matrix[0, 0]) % phase
    for i, j in np.ndindex(matrix.shape):
        exponent = int(residues[i, j])
        if np.gcd(exponent, phase) == 1:
            return {"rows": [0, i], "columns": [0, j], "exponent": exponent,
                    "root_order": phase,
                    "invariant_under_row_column_phase_gauges": True}
    return None


def twelfth_root_provenance():
    """Match the stored literal input to the author's published family."""
    base = np.asarray(data("P7")["exponents"], dtype=np.int64)
    mask = np.zeros((7, 7), dtype=np.int64)
    mask[1:3, 1:3] = 1
    mask[3:5, 3:5] = -1
    published = (2 * base + mask) % 12  # the author's parameter p=1/12
    stored = np.asarray(data("P7_q12")["exponents"], dtype=np.int64)
    verify_source(published, 12)
    verify_source(stored, 12)
    first = canonical_finite_butson_key(published.tolist(), 12)
    second = canonical_finite_butson_key(stored.tolist(), 12)
    if first != second or first[1] != 12:
        raise ValueError("the stored twelfth-root matrix does not match the published family")
    witness = cross_ratio(stored, 12)
    if witness is None:
        raise ValueError("no genuine twelfth-root cross ratio")
    return {"published_parameter": "1/12", "published_specialization_exponents": published.tolist(),
            "exact_monomial_equivalence_verified": True,
            "canonical_key": [first[0], first[1], list(first[2])],
            "primitive_cross_ratio": witness,
            "not_equivalent_to_a_sixth_root_matrix": True}


def block_width_gate(width):
    """Necessary order condition for a normalized uniform sign block."""
    if type(width) is not int or width < 1:
        raise ValueError("a positive integer block width is required")
    return {"width": width, "passes": width in (1, 2) or width % 4 == 0,
            "required_identity": "B_e B_e^T = width*I for every used scalar block",
            "scope": "Uniform orthogonal scalar action and a frame F with FF^T=width*I",
            "general_weighted_or_nonuniform_scalar_maps_excluded": False}


@cache
def load_source(source_id):
    if source_id not in SOURCE_IDS:
        raise KeyError(source_id)
    name = {SOURCE_IDS[0]: "P7", SOURCE_IDS[1]: "P7_q12", SOURCE_IDS[2]: "D14"}[source_id]
    record = data(name)
    matrix = np.asarray(record["exponents"], dtype=np.int64)
    phase = record["phase"]
    acquisition = {"data_file": str((DATA / f"{name}.json").relative_to(ROOT)),
                   "data_sha256": hashlib.sha256((DATA / f"{name}.json").read_bytes()).hexdigest(),
                   "new_source_construction_operation": False,
                   "scope": "This fixed starting input; subsequent operations are scalar restrictions"}
    if name == "P7":
        gauge = cubic_gauge(matrix)
        if not gauge["exists"]:
            raise ValueError("the published P7 lost its admitted gauge")
        matrix = (matrix + np.asarray(gauge["row_shifts"])[:, None]
                  + np.asarray(gauge["column_shifts"])[None, :]) % phase
        acquisition["fixed_starting_gauge"] = gauge
    elif name == "P7_q12":
        acquisition["published_family_equivalence"] = twelfth_root_provenance()
    proof = {**verify_source(matrix, phase), "acquisition": acquisition}
    alphabet = tuple(map(int, np.unique(matrix)))
    proof["alphabet_gate"] = alphabet_obstruction(phase, alphabet)
    if proof["alphabet_gate"]["excluded"]:
        raise ValueError("a retained starting input has an odd scalar relation")
    matrix = matrix.astype(np.uint16)
    matrix.setflags(write=False)
    return VerifiedSource(source_id, len(matrix), phase, exponent_hash(matrix), alphabet,
                          proof, matrix, lambda i, j: int(matrix[i, j]))


def certificate():
    frontier_file = ROOT / "certificates/scalar-lift/forward-novelty-frontier-2026-11-22.json"
    targets = json.loads(frontier_file.read_text())["retained_source_listed_unresolved_orders"]
    retained = []
    for source_id in SOURCE_IDS:
        source = load_source(source_id)
        retained.append({**source.metadata(), "exponents": source.exponents.tolist(),
                         "frontier_factor_routes": [
                             {"target_order": target, "required_total_width": target // source.order,
                              "normalized_uniform_sign_block_gate": block_width_gate(target // source.order),
                              "ambient_field_and_sign_frame_constructed": False}
                             for target in targets if target % source.order == 0]})
    rejected = []
    for name in ("M13B", "M13C"):
        record = data(name)
        proof = verify_source(np.asarray(record["exponents"]), 6)
        gate = alphabet_obstruction(6, range(6))
        gauge = cubic_gauge(record["exponents"])
        if gauge["exists"]:
            raise ValueError("the saved exact gauge exclusion changed")
        rejected.append({"name": name, "source_url": record["source_url"],
                         "full_Gram": proof, "raw_alphabet_gate": gate,
                         "two_cubic_coset_gauge": gauge,
                         "signed_three_class_relation": {
                             "exhaustive_gauge_scope": "Sixth-root gauges only; larger-field gauges not tested",
                             "argument": "Every gauge uses all three exponent classes modulo3. Choosing one root from each class and changing their integer coefficient signs yields a rotated 1+zeta3+zeta3^2=0.",
                             "coefficient_sum_is_odd": True},
                         "general_nonuniform_scalar_maps_excluded": False})
    return {"schema": "catalogue-source-acquisition-v1",
            "catalogue_commit": data("P7")["catalogue_commit"],
            "retained_inputs": retained, "rejected_uniform_alphabet_candidates": rejected,
            "frontier_file": str(frontier_file.relative_to(ROOT)),
            "frontier_sha256": hashlib.sha256(frontier_file.read_bytes()).hexdigest(),
            "largest_matching_frontier_order": max(r["target_order"] for s in retained
                                                   for r in s["frontier_factor_routes"]),
            "new_Hadamard_order_constructed": False,
            "scope": "Exact acquisition and necessary factor/alphabet tests; no target sign frame"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--output", type=Path)
    group.add_argument("--check", type=Path)
    args = parser.parse_args()
    report = certificate()
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    elif json.loads(args.check.read_text()) != report:
        raise ValueError("the catalogue acquisition certificate differs from exact replay")
    print(json.dumps({"verified_new_starting_inputs": len(report["retained_inputs"]),
                      "rejected_gauge_orbits": len(report["rejected_uniform_alphabet_candidates"])}))


if __name__ == "__main__":
    main()

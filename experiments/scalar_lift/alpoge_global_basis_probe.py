"""Global rational two-basis restrictions anchored at published H668/H716.

The output isometry is a product of shuffled Walsh-four layers.  Each layer
acts on every coordinate; no low-rank Householder ansatz or source search is
used.  Candidate partners are evaluated by exact quartet sums only.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from experiments.scalar_lift.alpoge_mixed_basis_probe import (
    binary_rank,
    canonical_fraction,
    entry_statistics,
    frame_inputs,
    twice_root_action,
    weighted_trace_model,
)
from experiments.scalar_lift.ambient_field_candidates import load_ambient_source, period_field

WALSH4 = np.array([[1, 1, 1, 1], [1, -1, 1, -1],
                   [1, 1, -1, -1], [1, -1, -1, 1]], dtype=np.int64)


def apply_layer(values, permutation):
    """Apply a shuffled block-diagonal Walsh4 numerator in O(d * columns)."""
    values = np.asarray(values, dtype=np.int64)
    permutation = np.asarray(permutation, dtype=np.int64)
    order = values.shape[0]
    if order % 4 or permutation.shape != (order,) or not np.array_equal(np.sort(permutation), np.arange(order)):
        raise ValueError("a permutation of a multiple-of-four coordinate space is required")
    if np.max(np.abs(values), initial=0) >= 2**60:
        raise OverflowError("the four-term layer exceeds the exact int64 bound")
    groups = permutation.reshape(-1, 4)
    a, b, c, d = (values[groups[:, j]] for j in range(4))
    result = np.empty_like(values)
    result[groups[:, 0]] = a+b+c+d
    result[groups[:, 1]] = a-b+c-d
    result[groups[:, 2]] = a+b-c-d
    result[groups[:, 3]] = a-b-c+d
    return result


def global_partner(frame, layers):
    """Exact numerator/denominator of Q R0 Q^T F, Q=Q_k ... Q_1."""
    if not layers:
        raise ValueError("at least one global layer is required")
    values = np.asarray(frame, dtype=np.int64)
    if values.ndim != 2 or values.shape[0] != values.shape[1] or not np.all(np.abs(values) == 1):
        raise ValueError("a fixed square sign frame is required")
    for permutation in reversed(layers):
        values = apply_layer(values, permutation)
    values = twice_root_action(values)
    for permutation in layers:
        values = apply_layer(values, permutation)
    return canonical_fraction(values, 2**(2*len(layers)+1))


def global_rank_lower_bound(order, layers):
    """rank_Q(Q-I)>=rank_F2(Q_num), as Q_num-2^k I reduces to Q_num."""
    values = np.eye(order, dtype=np.uint8)
    for permutation in layers:
        groups = np.asarray(permutation).reshape(-1, 4)
        sums = np.bitwise_xor.reduce(values[groups], axis=1)
        result = np.empty_like(values)
        result[groups] = sums[:, None, :]
        values = result
    return binary_rank(values)


def nonproduct_column(order, layers):
    """Exhibit a nonzero tensor minor of a column of the output basis Q^T."""
    r = order//4
    for column in range(min(order, 8)):
        vector = np.zeros(order, dtype=np.int64)
        vector[column] = 1
        for permutation in reversed(layers):
            vector = apply_layer(vector, permutation)
        table = vector.reshape(4, r)
        for a in range(4):
            for b in range(a+1, 4):
                occupied = np.flatnonzero((table[a] != 0) | (table[b] != 0))
                if not len(occupied):
                    continue
                c = int(occupied[0])
                minors = table[a, c]*table[b]-table[b, c]*table[a]
                nonzero = np.flatnonzero(minors)
                if len(nonzero):
                    e = int(nonzero[0])
                    return {"basis_column": column, "local_rows": [a, b],
                            "period_columns": [c, e],
                            "integer_minor_entries": [[int(table[a, c]), int(table[a, e])],
                                                      [int(table[b, c]), int(table[b, e])]],
                            "integer_minor_determinant": int(minors[e]),
                            "basis_column_denominator": 2**len(layers)}
    return None


def recipes(order, per_depth=16, maximum_depth=4):
    for depth in range(1, maximum_depth+1):
        for index in range(per_depth):
            seed = 2026112300+order*10000+depth*100+index
            rng = np.random.default_rng(seed)
            layers = [rng.permutation(order).astype(np.int64) for _ in range(depth)]
            yield {"depth": depth, "index": index, "seed": seed}, layers


def partner_hash(numerator, denominator):
    return hashlib.sha256(numerator.astype("<i8").tobytes()+str(denominator).encode()).hexdigest()


def _census(numerator):
    values, counts = np.unique(numerator, return_counts=True)
    return {str(int(value)): int(count) for value, count in zip(values, counts, strict=True)}


def run_case(order, frame, proof):
    baseline = twice_root_action(frame)
    if np.any(baseline % 2):
        raise ArithmeticError("the fixed frame baseline is not integral")
    required_rank = (binary_rank(baseline//2-1)+1)//2
    sources = [load_ambient_source(identifier) for identifier in
               (("catalogue-p7-q6", "unreal13-q6") if order == 668 else ("unreal13-q6",))]
    records, seen = [], {}
    for recipe, layers in recipes(order):
        numerator, denominator = global_partner(frame, layers)
        digest = partner_hash(numerator, denominator)
        rank = global_rank_lower_bound(order, layers)
        witness = nonproduct_column(order, layers)
        if witness is None:
            raise ArithmeticError("a global basis has no checked nonproduct column")
        record = {**recipe, "layer_permutations": [p.tolist() for p in layers],
                  "partner_sha256": digest, "partner_denominator": denominator,
                  "rank_Q_minus_I_certified_lower_bound": rank,
                  "meets_low_rank_obstruction_threshold": rank >= required_rank,
                  "nonproduct_output_basis_column": witness}
        if digest in seen:
            record["literal_duplicate_of"] = seen[digest]
        else:
            seen[digest] = len(records)
            record["entry_statistics"] = entry_statistics(numerator, denominator)
            record["complete_numerator_value_counts"] = _census(numerator)
            if sum(record["complete_numerator_value_counts"].values()) != order**2:
                raise ArithmeticError("the exact complete matrix census has the wrong size")
        records.append(record)
    distinct = [r for r in records if "entry_statistics" in r]
    best = max(distinct, key=lambda r: r["entry_statistics"]["sign_entries"])
    return {"published_frame_order": order, "published_frame_proof": proof,
            "frame_sha256": hashlib.sha256(frame.astype("<i8").tobytes()).hexdigest(),
            "frame_row_order": "four channels: own border rows, then own core rows",
            "ambient_field": period_field(order//4).proof(),
            "weighted_trace_model": weighted_trace_model(),
            "faithful_auxiliary_generator": "Q (I4 tensor S^-1 M_eta0 S) Q^T",
            "fixed_Butson_sources": [s.metadata() for s in sources],
            "restricted_source_orders": [s.order*order for s in sources],
            "minimum_rank_Q_minus_I_for_any_sign_partner": required_rank,
            "rank_lower_bound_proof": "rank_Q(Q-I)=rank_Q(Qnum-2^depth I)>=rank_F2(Qnum)",
            "candidate_visits": len(records), "distinct_partners": len(seen),
            "literal_duplicates_pruned": len(records)-len(seen), "candidates": records,
            "all_candidates_certifiably_above_low_rank_cutoff": all(r["meets_low_rank_obstruction_threshold"] for r in records),
            "best_recipe": {k: best[k] for k in ("depth", "index", "seed")},
            "best_sign_entries": best["entry_statistics"]["sign_entries"],
            "total_partner_entries": order**2,
            "sign_partner_found": any(r["entry_statistics"]["all_entries_signs"] for r in distinct)}


class SavedGlobalBasisNode:
    """Replay a saved partner and access the full restricted source implicitly."""

    def __init__(self, path, case=0, candidate=0):
        report = json.loads(Path(path).read_text())
        case_record = report["cases"][case]
        record = case_record["candidates"][candidate]
        order, frame, proof = frame_inputs()[case]
        if order != case_record["published_frame_order"] or proof != case_record["published_frame_proof"]:
            raise ValueError("the published fixed frame does not replay")
        if hashlib.sha256(frame.astype("<i8").tobytes()).hexdigest() != case_record["frame_sha256"]:
            raise ValueError("the fixed frame coordinates do not replay")
        if (period_field(order//4).proof() != case_record["ambient_field"]
                or weighted_trace_model() != case_record["weighted_trace_model"]):
            raise ValueError("the faithful ambient field or weighted trace basis does not replay")
        self.frame, self.order = frame, order
        self.layers = [np.asarray(p, dtype=np.int64) for p in record["layer_permutations"]]
        self.numerator, self.denominator = global_partner(frame, self.layers)
        if partner_hash(self.numerator, self.denominator) != record["partner_sha256"]:
            raise ValueError("the exact scalar partner does not replay")
        if "complete_numerator_value_counts" in record and _census(self.numerator) != record["complete_numerator_value_counts"]:
            raise ValueError("the complete partner census does not replay")
        self.sources = {s["id"]: load_ambient_source(s["id"]) for s in case_record["fixed_Butson_sources"]}
        for source in case_record["fixed_Butson_sources"]:
            if self.sources[source["id"]].source_hash != source["source_hash"]:
                raise ValueError("the fixed Butson source does not replay")

    def entry(self, source_id, row, column):
        source = self.sources[source_id]
        total_order = source.order*self.order
        if type(row) is not int or type(column) is not int or not 0 <= row < total_order or not 0 <= column < total_order:
            raise ValueError("restricted matrix coordinate outside its order")
        i, a = divmod(row, self.order)
        j, b = divmod(column, self.order)
        exponent = source.entry(i, j)
        f = int(self.frame[a, b])*self.denominator
        g = int(self.numerator[a, b])
        numerator = (f, g, g-f, -f, -g, f-g)[exponent]
        return {"numerator": numerator, "denominator": self.denominator}


def report():
    if not np.array_equal(WALSH4 @ WALSH4.T, 4*np.eye(4, dtype=np.int64)):
        raise ArithmeticError("the exact Walsh-four isometry identity failed")
    cases = [run_case(*case) for case in frame_inputs()]
    return {"schema": "alpoge-global-two-basis-restriction-v1", "cases": cases,
            "coordinate_map": "B(z)=Q rho0(z) Q^T F; Q=(T_k/2)...(T_1/2)",
            "output_basis": "Q^T times the weighted-trace orthonormal field basis",
            "input_basis": "output basis times the fixed published frame F",
            "weighted_basis_metrics": "G_out=I_d, G_in=F^T F=d I_d",
            "scalar_blocks": ["F", "G", "G-F", "-F", "-G", "F-G"],
            "exact_partner_identities": "G G^T=d I and G F^T+F G^T=d I",
            "ordinary_full_restricted_source_Gram": "n*d I",
            "orthogonality_proof": "Each T_j is a permutation conjugate of blockdiag(Walsh4), T_j T_j^T=4 I",
            "candidate_visits": sum(c["candidate_visits"] for c in cases),
            "distinct_partners": sum(c["distinct_partners"] for c in cases),
            "literal_duplicates_pruned": sum(c["literal_duplicates_pruned"] for c in cases),
            "only_operations": "scalar restriction in declared unequal mixed bases and exact literal partner duplicate pruning",
            "source_construction_or_SAT_used": False, "pending_declared_candidates": 0,
            "new_sign_Hadamard_found": any(c["sign_partner_found"] for c in cases),
            "scope": "Sixteen deterministic global layer recipes at each depth1..4 for each fixed frame; not all unequal bases"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    choice = parser.add_mutually_exclusive_group(required=True)
    choice.add_argument("--output", type=Path)
    choice.add_argument("--check", type=Path)
    args = parser.parse_args()
    result = report()
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2)+"\n")
    elif json.loads(args.check.read_text()) != result:
        raise ValueError("global-basis probe certificate differs from exact replay")
    print(json.dumps({k: result[k] for k in
                      ("candidate_visits", "distinct_partners", "literal_duplicates_pruned", "new_sign_Hadamard_found")}))


if __name__ == "__main__":
    main()

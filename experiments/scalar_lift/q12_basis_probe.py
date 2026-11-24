"""Exact mixed two-basis restrictions of a genuine twelfth-root input.

Every used scalar is checked, including the independent fourth-root action.
The only starting representative change is the recorded global phase -1/12.
The fixed Butson source and published real frames are never searched for.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from functools import cache
from pathlib import Path

import numpy as np
from sympy import Matrix, eye

from experiments.scalar_lift.alpoge_global_basis_probe import (
    apply_layer,
    nonproduct_column,
    partner_hash,
    recipes,
)
from experiments.scalar_lift.alpoge_mixed_basis_probe import (
    canonical_fraction,
    frame_inputs,
    weighted_trace_model,
)
from experiments.scalar_lift.ambient_field_candidates import _quaternion, period_field
from experiments.scalar_lift.catalogue_source_acquisition import load_source
from experiments.scalar_lift.forward_butson_descent import exponent_hash, verify_source
from experiments.scalar_lift.forward_source_registry import VerifiedSource
from experiments.scalar_lift.scalar_restriction_obstructions import alphabet_obstruction

ROOT = Path(__file__).resolve().parents[2]
REPRESENTATIVES = (0, 1, 3, 4)
PALETTE = (0, 1, 3, 4, 6, 7, 9, 10)


@cache
def fixed_source():
    base = load_source("catalogue-p7-q12")
    values = (base.exponents.astype(np.int64)-1) % 12
    proof = verify_source(values, 12)
    used = tuple(map(int, np.unique(values)))
    if used != PALETTE:
        raise ValueError("the fixed twelfth-root representative changed alphabet")
    gate = alphabet_obstruction(12, used)
    if gate["excluded"]:
        raise ArithmeticError("the admitted source acquired an odd scalar relation")
    proof = {**proof, "fixed_starting_representative": {
        "base_source_id": base.id, "base_source_hash": base.source_hash,
        "global_exponent_shift": -1, "phase": 12,
        "operation": "multiply the fixed source by zeta12^-1 before restriction",
        "no_new_source_construction": True}, "alphabet_gate": gate,
        "published_provenance": base.proof["acquisition"]}
    values = values.astype(np.uint16)
    values.setflags(write=False)
    return VerifiedSource("catalogue-p7-q12-global-minus-one", 7, 12,
                          exponent_hash(values), used, proof, values,
                          lambda i, j: int(values[i, j]))


@cache
def twice_local_actions():
    """Return exactly 2*rho(zeta12^e), including the nonsextic actions."""
    identity = np.eye(4, dtype=np.int64)
    d = _quaternion(0, 1, 1, 1)
    j = np.array([[0, -1, 0, 0], [1, 0, 0, 0],
                  [0, 0, 0, 1], [0, 0, -1, 0]], dtype=np.int64)
    root = Matrix((j@(identity-d)).tolist())/2
    if root**3 != Matrix(j.tolist()) or root**2-Matrix(j.tolist())*root-eye(4) != 0*eye(4):
        raise ArithmeticError("the twelfth-root and fourth-root relation failed")
    result = []
    for exponent in range(12):
        twice = 2*root**exponent
        if any(value.q != 1 for value in twice):
            raise ArithmeticError("the local action requires an unexpected denominator")
        values = np.asarray(twice.tolist(), dtype=np.int64)
        if not np.array_equal(values@values.T, 4*identity):
            raise ArithmeticError("a local root action is not orthogonal")
        result.append(values)
    table = np.stack(result)
    if not np.array_equal(table[6:], -table[:6]):
        raise ArithmeticError("opposite roots do not have opposite actions")
    table.setflags(write=False)
    return table


def twice_phase_action(values, exponent):
    """Apply a four-dimensional local action tensor the period identity."""
    values = np.asarray(values, dtype=np.int64)
    width = len(values)
    if width % 4 or np.max(np.abs(values), initial=0) >= 2**59:
        raise ValueError("coordinates outside the certified exact arithmetic domain")
    local = twice_local_actions()[int(exponent) % 12]
    return np.tensordot(local, values.reshape((4, width//4)+values.shape[1:]),
                        axes=(1, 0)).reshape(values.shape)


def scalar_blocks(frame, layers, used_alphabet=PALETTE):
    """All needed opposite-root representatives, for any twelfth-root input."""
    frame = np.asarray(frame, dtype=np.int64)
    if frame.ndim != 2 or frame.shape[0] != frame.shape[1] or not np.all(np.abs(frame) == 1):
        raise ValueError("a fixed square sign frame is required")
    back = frame
    for permutation in reversed(layers):
        back = apply_layer(back, permutation)
    result = {}
    for exponent in sorted({int(e) % 6 for e in used_alphabet}):
        if exponent == 0:
            result[exponent] = (frame.copy(), 1)
            continue
        values = twice_phase_action(back, exponent)
        for permutation in layers:
            values = apply_layer(values, permutation)
        result[exponent] = canonical_fraction(values, 2**(2*len(layers)+1))
    return result


def block_statistics(numerator, denominator):
    signs = np.abs(numerator) == denominator
    zeros = numerator == 0
    integers = numerator % denominator == 0
    count = int(np.count_nonzero(signs))
    failed = None
    if count != numerator.size:
        position = int(np.argmax(~signs))
        i, j = divmod(position, numerator.shape[1])
        failed = {"row": i, "column": j, "numerator": int(numerator[i, j]),
                  "denominator": denominator}
    return {"sign_entries": count, "zero_entries": int(np.count_nonzero(zeros)),
            "other_integer_entries": int(np.count_nonzero(integers & ~signs & ~zeros)),
            "noninteger_entries": int(np.count_nonzero(~integers)),
            "total_entries": int(numerator.size), "all_entries_signs": failed is None,
            "first_failed_entry": failed}


def state_hash(source, frame_hash, blocks):
    actions = {str(e): partner_hash(*blocks[e]) for e in sorted(blocks)}
    identity = {"source_hash": source.source_hash, "frame_hash": frame_hash,
                "representative_scalar_blocks": actions,
                "opposite_roots_are_negative_blocks": True}
    return hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest(), actions


def candidate_record(source, frame_hash, frame, layers, recipe):
    if source.phase != 12:
        raise ValueError("embed the fixed input in twelfth-root coordinates first")
    blocks = scalar_blocks(frame, layers, source.used_alphabet)
    identity, action_hashes = state_hash(source, frame_hash, blocks)
    representative_stats = {e: block_statistics(*blocks[e]) for e in blocks}
    used_stats = {}
    for exponent in source.used_alphabet:
        base, negative = exponent % 6, exponent >= 6
        stats = dict(representative_stats[base])
        if negative and stats["first_failed_entry"] is not None:
            stats["first_failed_entry"] = {**stats["first_failed_entry"],
                                           "numerator": -stats["first_failed_entry"]["numerator"]}
        used_stats[str(exponent)] = stats
    counts = np.bincount(source.exponents.ravel(), minlength=12)
    total_signs = sum(int(counts[e])*used_stats[str(e)]["sign_entries"] for e in source.used_alphabet)
    total = (source.order*len(frame))**2
    passed = all(stats["all_entries_signs"] for stats in used_stats.values())
    first_failed = next(({"exponent": int(e), **used_stats[str(e)]["first_failed_entry"]}
                         for e in source.used_alphabet
                         if used_stats[str(e)]["first_failed_entry"] is not None), None)
    return {**recipe, "layer_permutations": [p.tolist() for p in layers],
            "literal_state_sha256": identity,
            "representative_scalar_sha256": action_hashes,
            "representative_denominators": {str(e): blocks[e][1] for e in sorted(blocks)},
            "nonproduct_output_basis_column": nonproduct_column(len(frame), layers),
            "all_used_scalar_statistics": used_stats,
            "complete_restricted_matrix_sign_entries": total_signs,
            "complete_restricted_matrix_total_entries": total,
            "all_used_scalars_flat": passed, "first_failed_used_scalar": first_failed}


def run_case(order, frame, proof, per_depth=16, maximum_depth=4):
    source = fixed_source()
    frame_hash = hashlib.sha256(frame.astype("<i8").tobytes()).hexdigest()
    records, seen = [], {}
    for recipe, layers in recipes(order, per_depth, maximum_depth):
        record = candidate_record(source, frame_hash, frame, layers, recipe)
        identity = record["literal_state_sha256"]
        if identity in seen:
            record["literal_duplicate_of"] = seen[identity]
        else:
            seen[identity] = len(records)
        if record["nonproduct_output_basis_column"] is None:
            raise ArithmeticError("the mixed basis has no certified nonproduct column")
        records.append(record)
    best = max(records, key=lambda r: r["complete_restricted_matrix_sign_entries"])
    target = source.order*order
    frontier = json.loads((ROOT/"certificates/scalar-lift/forward-novelty-frontier-2026-11-22.json").read_text())
    return {"published_frame_order": order, "published_frame_proof": proof,
            "frame_sha256": frame_hash,
            "frame_row_order": "four channels: own border rows, then own core rows",
            "ambient_field": period_field(order//4).proof(),
            "weighted_trace_model": weighted_trace_model(),
            "faithful_auxiliary_action": "Q (I4 tensor S^-1 M_eta S) Q^T",
            "target_order": target,
            "target_listed_unresolved_in_snapshot": target in frontier["retained_source_listed_unresolved_orders"],
            "candidate_visits": len(records), "literal_distinct_nodes": len(seen),
            "literal_duplicates_pruned": len(records)-len(seen), "candidates": records,
            "best_recipe": {k: best[k] for k in ("depth", "index", "seed")},
            "best_complete_matrix_sign_entries": best["complete_restricted_matrix_sign_entries"],
            "complete_matrix_total_entries": best["complete_restricted_matrix_total_entries"],
            "all_used_roots_checked_in_every_node": True,
            "sign_Hadamard_found": any(r["all_used_scalars_flat"] for r in records)}


def report():
    source = fixed_source()
    twice_local_actions()
    cases = [run_case(*case) for case in frame_inputs()]
    return {"schema": "genuine-q12-mixed-two-basis-probe-v1",
            "fixed_Butson_source": source.metadata(), "literal_source_exponents": source.exponents.tolist(),
            "original_projective_phase": 12, "actual_used_alphabet": list(source.used_alphabet),
            "palette_description": "mu4 union zeta12*mu4; four blocks F,RF,JF,JRF and their negatives",
            "local_actions": {"D": "left(i+j+k)", "J": "right i",
                              "R": "J(I-D)/2", "relations": ["R^3=J", "R^2-JR-I=0", "R^12=I"]},
            "coordinate_map": "B(z)=Q rho0(z) Q^T F; Q=(T_depth/2)...(T_1/2)",
            "input_basis": "output basis times the fixed published F",
            "weighted_basis_metrics": "G_out=I_d, G_in=F^T F=dI_d for the recorded weighted trace",
            "full_restricted_Gram": "n*d I, from exact source Gram and orthogonal scalar action",
            "cases": cases, "candidate_visits": sum(c["candidate_visits"] for c in cases),
            "literal_distinct_nodes": sum(c["literal_distinct_nodes"] for c in cases),
            "literal_duplicates_pruned": sum(c["literal_duplicates_pruned"] for c in cases),
            "pending_declared_candidates": 0,
            "only_operations": "scalar restriction in declared unequal mixed bases; exact literal duplicate pruning",
            "new_source_construction_or_SAT_used": False,
            "new_sign_Hadamard_found": any(c["sign_Hadamard_found"] for c in cases),
            "scope": "Sixteen fixed Walsh-layer recipes at each depth1..4 for both published frames; all eight used roots, not a sextic surrogate; not all bases"}


class SavedQ12BasisNode:
    """Replay all used scalar blocks and access any full restricted entry."""
    def __init__(self, path, case=0, candidate=0):
        saved = json.loads(Path(path).read_text())
        source = fixed_source()
        if saved["fixed_Butson_source"] != source.metadata():
            raise ValueError("the fixed twelfth-root input does not replay")
        case_record = saved["cases"][case]
        record = case_record["candidates"][candidate]
        order, frame, proof = frame_inputs()[case]
        frame_hash = hashlib.sha256(frame.astype("<i8").tobytes()).hexdigest()
        if (order != case_record["published_frame_order"] or proof != case_record["published_frame_proof"]
                or frame_hash != case_record["frame_sha256"]):
            raise ValueError("the published input does not replay")
        if (period_field(order//4).proof() != case_record["ambient_field"]
                or weighted_trace_model() != case_record["weighted_trace_model"]):
            raise ValueError("the faithful field or weighted trace basis does not replay")
        layers = [np.asarray(p, dtype=np.int64) for p in record["layer_permutations"]]
        recipe = {k: record[k] for k in ("depth", "index", "seed")}
        replay = candidate_record(source, frame_hash, frame, layers, recipe)
        expected = {k: v for k, v in record.items() if k != "literal_duplicate_of"}
        if replay != expected:
            raise ValueError("the full used scalar map or exact census does not replay")
        self.source, self.width = source, order
        self.order = source.order*order
        self.blocks = scalar_blocks(frame, layers, source.used_alphabet)

    def entry(self, row, column):
        if type(row) is not int or type(column) is not int or not 0 <= row < self.order or not 0 <= column < self.order:
            raise ValueError("restricted matrix coordinate outside its order")
        i, a = divmod(row, self.width)
        j, b = divmod(column, self.width)
        exponent = self.source.entry(i, j)
        numerator, denominator = self.blocks[exponent % 6]
        sign = -1 if exponent >= 6 else 1
        return {"numerator": sign*int(numerator[a, b]), "denominator": denominator}


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
        raise ValueError("twelfth-root mixed-basis report differs from exact replay")
    print(json.dumps({k: result[k] for k in ("candidate_visits", "literal_distinct_nodes",
                                           "literal_duplicates_pruned", "new_sign_Hadamard_found")}))


if __name__ == "__main__":
    main()

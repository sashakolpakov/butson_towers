"""Share each exact mixed q12 geometry across all compatible fixed sources.

Six antipodal scalar representatives are streamed once per frame/basis.
Every retained source and explicit global anchor is then matched to that
complete action map; no source or frame is constructed or optimized.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from experiments.scalar_lift.alpoge_global_basis_probe import (
    apply_layer,
    nonproduct_column,
    recipes,
)
from experiments.scalar_lift.alpoge_mixed_basis_probe import (
    canonical_fraction,
    weighted_trace_model,
)
from experiments.scalar_lift.ambient_field_candidates import period_field
from experiments.scalar_lift.forward_butson_descent import exponent_hash, verify_source
from experiments.scalar_lift.phase_independent_frontier import (
    all_sources,
    anchor_cases,
    route_records,
)
from experiments.scalar_lift.q12_basis_probe import (
    block_statistics,
    twice_local_actions,
    twice_phase_action,
)

ROOT = Path(__file__).resolve().parents[2]
FRAME_SPECS = (("h668-width4", 668), ("h668-width12", 668), ("h716-width4", 716),
               ("h764-0", 764), ("h764-1", 764), ("h892", 892),
               ("h1132", 1132), ("h1964", 1964))


def json_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def array_hash(values, denominator=None):
    """Hash without the two dense temporary byte copies of astype/tobytes."""
    values = np.ascontiguousarray(values, dtype="<i8")
    digest = hashlib.sha256(memoryview(values).cast("B"))
    if denominator is not None:
        digest.update(b":"+str(denominator).encode())
    return digest.hexdigest()


def _border_order(matrix, border):
    width = len(matrix)
    core, own_border = (width-border)//4, border//4
    permutation = np.concatenate([np.r_[np.arange(i*own_border, (i+1)*own_border),
                                        np.arange(border+i*core, border+(i+1)*core)]
                                  for i in range(4)])
    return np.asarray(matrix, dtype=np.int64)[permutation]


def load_frame(identifier):
    """Load one immutable published frame at a time; do not cache all frames."""
    if identifier == "h668-width4":
        from experiments.scalar_lift.inverse_blocks import SEED, known_h668
        matrix = np.asarray(known_h668(), dtype=np.int64)
        proof = {"source_seed": str(SEED.relative_to(ROOT)),
                 "source_seed_sha256": hashlib.sha256(SEED.read_bytes()).hexdigest(),
                 "source_url": json.loads(SEED.read_text()).get("source_url"),
                 "published_matrix_hash_and_independent_bit_Gram_replayed": True}
        matrix = _border_order(matrix, 4)
    elif identifier == "h668-width12":
        from experiments.scalar_lift.h668_width12 import known_h668_width12
        matrix, proof = known_h668_width12()
        matrix = _border_order(matrix, 12)
    elif identifier == "h716-width4":
        from experiments.scalar_lift.alpoge_published_inputs import known_h716
        matrix, proof = known_h716()
        matrix = _border_order(matrix, 4)
    elif identifier in ("h764-0", "h764-1"):
        from experiments.scalar_lift.published_h764 import known_h764
        matrix, proof = known_h764(int(identifier[-1]))
    elif identifier in {"h892", "h1132", "h1964"}:
        from experiments.scalar_lift.published_alpoge_gs import known_gs_frame
        matrix, proof = known_gs_frame(int(identifier[1:]))
    else:
        raise KeyError(identifier)
    matrix = np.asarray(matrix, dtype=np.int64)
    width = len(matrix)
    if (matrix.shape != (width, width) or not np.all(np.abs(matrix) == 1)
            or not np.array_equal(matrix@matrix.T, width*np.eye(width, dtype=np.int64))):
        raise ArithmeticError("the published frame failed the full exact sign Gram")
    return matrix, proof


def embed_source(source, shift):
    """Exactly embed an explicit source/global-phase recipe in q12 labels."""
    if source.exponents is None:
        raise ValueError("a complete dense fixed source is required")
    original = (source.exponents.astype(np.int64)+shift) % source.phase
    if np.any((12*original) % source.phase):
        raise ValueError("the actual shifted source roots do not embed in q12")
    values = (12*original//source.phase) % 12
    return values


def prepare_bindings(widths, routes=None, source_map=None):
    """Prune only literal embedded matrix duplicates, retaining every origin."""
    routes = route_records() if routes is None else routes
    source_map = {s.id: s for s in all_sources()} if source_map is None else source_map
    groups, selected_routes, uncovered = {}, [], []
    cases = {}
    for route in routes:
        if route["required_total_width"] not in widths or 12 % route["actual_root_generator_order"]:
            uncovered.append(route)
            continue
        selected_routes.append(route)
        source = source_map[route["source_id"]]
        if source.id not in cases:
            cases[source.id] = {case["id"]: case for case in anchor_cases(source)}
        case = cases[source.id][route["alphabet_case"]]
        shifts = {0} if case["contains_raw_representative"] else set()
        shifts.update(recipe["global_exponent_shift"] for recipe in case["anchor_gauge_recipes"])
        for shift in sorted(shifts):
            values = embed_source(source, shift)
            digest = exponent_hash(values)
            identity = json_hash({"order": source.order, "phase": 12, "exponents": digest})
            if identity not in groups:
                used = tuple(map(int, np.unique(values)))
                first = {str(e): list(map(int, np.argwhere(values == e)[0])) for e in used}
                groups[identity] = {"id": identity, "order": source.order, "phase": 12,
                    "source_hash": digest, "literal_exponents": values.tolist(),
                    "full_Gram": verify_source(values, 12), "used_alphabet": list(used),
                    "exponent_frequencies": np.bincount(values.ravel(), minlength=12).tolist(),
                    "first_source_position_for_each_used_root": first,
                    "widths": set(), "origins": {}}
            binding = groups[identity]
            binding["widths"].add(route["required_total_width"])
            key = json_hash({"source_id": source.id, "source_hash": source.source_hash,
                             "phase": source.phase, "global_shift": shift})
            origin = binding["origins"].setdefault(key, {
                "source_id": source.id, "source_hash": source.source_hash,
                "source_proof_sha256": json_hash(source.proof), "declared_phase": source.phase,
                "global_exponent_shift": shift, "alphabet_case": case["id"],
                "embedding": "12*((original_exponent+shift) mod declared_phase)/declared_phase",
                "every_embedded_exponent_checked_integral": True,
                "widths": set(), "preflight_route_ids": set()})
            origin["widths"].add(route["required_total_width"])
            origin["preflight_route_ids"].add(json_hash(route))
    bindings = []
    for identity in sorted(groups):
        binding = groups[identity]
        binding["widths"] = sorted(binding["widths"])
        binding["palette_id"] = json_hash(binding["used_alphabet"])
        origins = []
        for key in sorted(binding["origins"]):
            origin = binding["origins"][key]
            origin["widths"] = sorted(origin["widths"])
            origin["preflight_route_ids"] = sorted(origin["preflight_route_ids"])
            origins.append(origin)
        binding["origins"] = origins
        bindings.append(binding)
    return bindings, selected_routes, uncovered


def replay_binding_origins(binding, source_map=None):
    """Replay literal Gram, every original source, anchor, and root embedding."""
    values = np.asarray(binding["literal_exponents"], dtype=np.int64)
    if (exponent_hash(values) != binding["source_hash"]
            or verify_source(values, 12) != binding["full_Gram"]):
        raise ValueError("a complete fixed source does not replay")
    source_map = {s.id: s for s in all_sources()} if source_map is None else source_map
    for origin in binding["origins"]:
        source = source_map[origin["source_id"]]
        if (source.source_hash != origin["source_hash"]
                or source.phase != origin["declared_phase"]
                or json_hash(source.proof) != origin["source_proof_sha256"]
                or not np.array_equal(embed_source(source, origin["global_exponent_shift"]), values)):
            raise ValueError("an original source, explicit global anchor, or exact q12 embedding does not replay")
    return True


def root_statistics(representative_stats, exponent):
    result = dict(representative_stats[str(exponent % 6)])
    if exponent >= 6 and result["first_failed_entry"] is not None:
        result["first_failed_entry"] = {**result["first_failed_entry"],
                                        "numerator": -result["first_failed_entry"]["numerator"]}
    return result


def geometry_data(frame, layers):
    """Stream all six blocks; retain statistics and hashes, never six arrays."""
    back = frame
    for permutation in reversed(layers):
        back = apply_layer(back, permutation)
    stats, hashes, denominators = {}, {}, {}
    for exponent in range(6):
        if exponent == 0:
            numerator, denominator = frame, 1
        else:
            values = twice_phase_action(back, exponent)
            for permutation in layers:
                values = apply_layer(values, permutation)
            numerator, denominator = canonical_fraction(values, 2**(2*len(layers)+1))
            del values
        stats[str(exponent)] = block_statistics(numerator, denominator)
        hashes[str(exponent)] = array_hash(numerator, denominator)
        denominators[str(exponent)] = denominator
        if exponent:
            del numerator
    return stats, hashes, denominators


def match_palettes(representative_stats, palettes):
    results = {}
    for key, used in palettes.items():
        failures = []
        for exponent in used:
            stat = root_statistics(representative_stats, exponent)
            if stat["first_failed_entry"] is not None:
                failures.append({"exponent": exponent, **stat["first_failed_entry"]})
        results[key] = {"used_alphabet": used, "all_used_roots_flat": not failures,
                        "failure_witness_for_every_failed_used_root": failures}
    return results


def run_geometry(identifier, bindings, per_depth=16, maximum_depth=4):
    frame, proof = load_frame(identifier)
    width = len(frame)
    selected = [b for b in bindings if width in b["widths"]]
    if not selected:
        raise ValueError("the published frame has no compatible retained source route")
    field = period_field(width//4)
    field_proof = field.proof()
    del field
    period_field.cache_clear()
    frame_hash = array_hash(frame)
    palettes = {b["palette_id"]: b["used_alphabet"] for b in selected}
    best = {b["id"]: {"sign_entries": -1, "node": None} for b in selected}
    seen_geometry, seen_source = {}, {b["id"]: {} for b in selected}
    records, source_duplicates = [], 0
    for recipe, layers in recipes(width, per_depth, maximum_depth):
        stats, hashes, denominators = geometry_data(frame, layers)
        identity = json_hash({"frame": frame_hash, "width": width,
                              "antipodal_scalar_hashes": hashes})
        record = {**recipe, "layer_permutations": [p.tolist() for p in layers],
                  "geometry_sha256": identity, "antipodal_scalar_hashes": hashes,
                  "antipodal_scalar_denominators": denominators,
                  "antipodal_scalar_statistics": stats,
                  "nonproduct_output_basis_column": nonproduct_column(width, layers),
                  "palette_matches": match_palettes(stats, palettes)}
        if record["nonproduct_output_basis_column"] is None:
            raise ArithmeticError("a declared mixed basis has no verified nonproduct column")
        if identity in seen_geometry:
            record["literal_geometry_duplicate_of"] = seen_geometry[identity]
        else:
            seen_geometry[identity] = len(records)
        sign_counts, identities, duplicates = [], [], []
        for binding in selected:
            count = sum(binding["exponent_frequencies"][e]*stats[str(e % 6)]["sign_entries"]
                        for e in binding["used_alphabet"])
            sign_counts.append(count)
            source_identity = json_hash({"source": binding["source_hash"],
                "used_scalar_hashes": {str(e): hashes[str(e % 6)] for e in binding["used_alphabet"]}})
            identities.append(source_identity)
            source_seen = seen_source[binding["id"]]
            if source_identity in source_seen:
                duplicates.append({"binding_id": binding["id"], "previous_node": source_seen[source_identity]})
                source_duplicates += 1
            else:
                source_seen[source_identity] = len(records)
            if count > best[binding["id"]]["sign_entries"]:
                best[binding["id"]] = {"sign_entries": count, "node": len(records)}
        record["source_sign_entry_counts_in_binding_order"] = sign_counts
        record["exact_source_used_action_identities_in_binding_order"] = identities
        record["literal_source_output_duplicates"] = duplicates
        records.append(record)
    summaries = [{"binding_id": b["id"], "target_order": b["order"]*width,
                  "total_entries": (b["order"]*width)**2, **best[b["id"]]}
                 for b in selected]
    return {"frame_id": identifier, "width": width, "published_frame_proof": proof,
            "frame_sha256": frame_hash, "ambient_field": field_proof,
            "weighted_trace_model": weighted_trace_model(),
            "faithful_auxiliary_action": "Q (I4 tensor S^-1 M_eta S) Q^T",
            "binding_ids": [b["id"] for b in selected], "palettes": palettes,
            "candidate_geometries": len(records), "literal_distinct_geometries": len(seen_geometry),
            "literal_geometry_duplicates": len(records)-len(seen_geometry),
            "source_matrix_evaluations": len(records)*len(selected),
            "literal_source_output_duplicates": source_duplicates,
            "literal_distinct_source_outputs": sum(len(s) for s in seen_source.values()),
            "palette_evaluations": len(records)*len(palettes),
            "source_binding_best_counts": summaries, "nodes": records,
            "all_sign_output_found": any(s["sign_entries"] == s["total_entries"] for s in summaries)}


def report(per_depth=16, maximum_depth=4, progress=False):
    widths = {width for _identifier, width in FRAME_SPECS}
    all_routes = route_records()
    bindings, selected, uncovered = prepare_bindings(widths, all_routes)
    geometries = []
    for identifier, _width in FRAME_SPECS:
        geometry = run_geometry(identifier, bindings, per_depth, maximum_depth)
        geometries.append(geometry)
        if progress:
            print(json.dumps({"completed_frame": identifier,
                              "geometries": geometry["candidate_geometries"],
                              "source_outputs": geometry["source_matrix_evaluations"],
                              "sign_output_found": geometry["all_sign_output_found"]}), flush=True)
    preflight_file = ROOT/"certificates/scalar-lift/phase-independent-frontier.json"
    return {"schema": "broad-shared-q12-geometry-restriction-v1",
            "preflight_file": str(preflight_file.relative_to(ROOT)),
            "preflight_file_sha256": hashlib.sha256(preflight_file.read_bytes()).hexdigest(),
            "preflight_surviving_routes": len(all_routes), "covered_preflight_routes": len(selected),
            "covered_route_sha256s": sorted(map(json_hash, selected)),
            "uncovered_surviving_routes": uncovered,
            "covered_declared_source_phases": sorted({r["declared_phase"] for r in selected}),
            "covered_source_ids": sorted({r["source_id"] for r in selected}),
            "covered_target_orders": sorted({r["target_order"] for r in selected}),
            "source_bindings": bindings,
            "literal_distinct_embedded_source_bindings": len(bindings),
            "origin_source_global_anchor_recipes": sum(len(b["origins"]) for b in bindings),
            "source_binding_deduplication": "Only identical complete embedded exponent matrices; every source/global-anchor origin retained",
            "geometries": geometries,
            "candidate_geometries": sum(g["candidate_geometries"] for g in geometries),
            "literal_distinct_geometries": sum(g["literal_distinct_geometries"] for g in geometries),
            "source_matrix_evaluations": sum(g["source_matrix_evaluations"] for g in geometries),
            "literal_distinct_source_outputs": sum(g["literal_distinct_source_outputs"] for g in geometries),
            "literal_source_output_duplicates": sum(g["literal_source_output_duplicates"] for g in geometries),
            "source_output_identity_scope": "Distinctness and output duplicate pruning are conservative within each fixed literal embedded binding; cross-binding/Galois output comparisons are not used",
            "palette_evaluations": sum(g["palette_evaluations"] for g in geometries),
            "scalar_representative_evaluations": 6*sum(g["candidate_geometries"] for g in geometries),
            "scalar_block_memory_policy": "one dense root block at a time, with one shared Q^T F; frame providers are not all cached",
            "coordinate_map": "B(z)=Q rho0(z) Q^T F in unequal bases, rho0 the complete q12 action",
            "full_restricted_Gram": "n*d I for every embedded source, by exact source Gram and orthogonal scalar representation",
            "orthogonal_layer_model": "Each shuffled Walsh4 layer is divided by2; all permutations and nonproduct witnesses are saved",
            "all_used_root_failure_witnesses_saved": True,
            "pending_declared_geometries": 0, "new_source_or_frame_search_or_SAT_used": False,
            "only_tree_operations": "scalar restriction in declared unequal mixed bases and exact literal duplicate pruning",
            "new_sign_Hadamard_found": any(g["all_sign_output_found"] for g in geometries),
            "scope": "All necessary-test-surviving q12-embeddable source/global-anchor recipes at the eight supplied fixed frame geometries; fixed finite Walsh bases, not all possible bases or surviving routes"}


def reconstruct_block(frame, layers, exponent):
    if exponent % 12 == 0:
        return frame, 1
    values = frame
    for permutation in reversed(layers):
        values = apply_layer(values, permutation)
    values = twice_phase_action(values, exponent)
    for permutation in layers:
        values = apply_layer(values, permutation)
    return canonical_fraction(values, 2**(2*len(layers)+1))


class SavedBroadBasisNode:
    """Replay a streamed geometry and bind it to any selected actual source."""
    def __init__(self, path, geometry=0, node=0):
        saved = json.loads(Path(path).read_text())
        case, record = saved["geometries"][geometry], saved["geometries"][geometry]["nodes"][node]
        self.frame, proof = load_frame(case["frame_id"])
        if proof != case["published_frame_proof"] or array_hash(self.frame) != case["frame_sha256"]:
            raise ValueError("the published fixed frame does not replay")
        field = period_field(len(self.frame)//4)
        if field.proof() != case["ambient_field"] or weighted_trace_model() != case["weighted_trace_model"]:
            raise ValueError("the faithful field or weighted trace basis does not replay")
        del field
        period_field.cache_clear()
        self.layers = [np.asarray(p, dtype=np.int64) for p in record["layer_permutations"]]
        stats, hashes, denominators = geometry_data(self.frame, self.layers)
        if (stats != record["antipodal_scalar_statistics"] or hashes != record["antipodal_scalar_hashes"]
                or denominators != record["antipodal_scalar_denominators"]):
            raise ValueError("the streamed root actions or complete entry statistics do not replay")
        if nonproduct_column(len(self.frame), self.layers) != record["nonproduct_output_basis_column"]:
            raise ValueError("the mixed-basis witness does not replay")
        self.bindings = {b["id"]: b for b in saved["source_bindings"] if b["id"] in case["binding_ids"]}
        sources = {source.id: source for source in all_sources()}
        for binding in self.bindings.values():
            replay_binding_origins(binding, sources)
        self.width = len(self.frame)
        self.last_exponent, self.last_block = None, None

    def entry(self, binding_id, row, column):
        binding = self.bindings[binding_id]
        size = binding["order"]*self.width
        if type(row) is not int or type(column) is not int or not 0 <= row < size or not 0 <= column < size:
            raise ValueError("restricted matrix coordinate outside its order")
        i, a = divmod(row, self.width)
        j, b = divmod(column, self.width)
        exponent = binding["literal_exponents"][i][j]
        if exponent != self.last_exponent:
            self.last_block = reconstruct_block(self.frame, self.layers, exponent)
            self.last_exponent = exponent
        numerator, denominator = self.last_block
        return {"numerator": int(numerator[a, b]), "denominator": denominator}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    choice = parser.add_mutually_exclusive_group(required=True)
    choice.add_argument("--output", type=Path)
    choice.add_argument("--check", type=Path)
    parser.add_argument("--per-depth", type=int, default=16)
    parser.add_argument("--maximum-depth", type=int, default=4)
    parser.add_argument("--progress", action="store_true")
    args = parser.parse_args()
    if args.per_depth < 1 or not 1 <= args.maximum_depth <= 4:
        raise ValueError("positive recipe count and depth1..4 are required")
    twice_local_actions()
    result = report(args.per_depth, args.maximum_depth, args.progress)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, separators=(",", ":"))+"\n")
    elif json.loads(args.check.read_text()) != result:
        raise ValueError("the broad shared-geometry restriction report differs from exact replay")
    print(json.dumps({k: result[k] for k in ("candidate_geometries", "source_matrix_evaluations",
        "literal_distinct_embedded_source_bindings", "palette_evaluations", "covered_target_orders",
        "new_sign_Hadamard_found")}))


if __name__ == "__main__":
    main()

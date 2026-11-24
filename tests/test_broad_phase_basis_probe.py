import copy
import json
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from experiments.scalar_lift.alpoge_global_basis_probe import WALSH4
from experiments.scalar_lift.broad_phase_basis_probe import (
    SavedBroadBasisNode,
    array_hash,
    embed_source,
    geometry_data,
    match_palettes,
    prepare_bindings,
    reconstruct_block,
    replay_binding_origins,
)
from experiments.scalar_lift.forward_source_registry import load_source
from experiments.scalar_lift.phase_independent_frontier import anchor_cases
from experiments.scalar_lift.q12_basis_probe import (
    PALETTE,
    block_statistics,
    fixed_source,
    scalar_blocks,
)


def route(source, width):
    case = anchor_cases(source)[0]
    return {"source_id": source.id, "source_hash": source.source_hash,
            "declared_phase": source.phase, "required_total_width": width,
            "actual_root_generator_order": case["actual_root_generator_order"],
            "alphabet_case": case["id"], "target_order": source.order*width}


def test_literal_embedded_source_dedup_retains_all_origins_and_replays_them():
    source = fixed_source()
    alias = replace(source, id="same-literal-different-provider")
    sources = {source.id: source, alias.id: alias}
    bindings, selected, uncovered = prepare_bindings({4}, [route(source, 4), route(alias, 4)], sources)
    assert len(selected) == 2 and not uncovered
    assert all(len(binding["origins"]) == 2 for binding in bindings)
    for binding in bindings:
        assert replay_binding_origins(binding, sources)
    bad = copy.deepcopy(bindings[0])
    bad["origins"][0]["source_hash"] = "0"*64
    with pytest.raises(ValueError, match="original source"):
        replay_binding_origins(bad, sources)


def test_declared_phase24_uses_actual_embeddable_roots_not_phase_label():
    source = load_source("cubic13-q24-t0")
    values = embed_source(source, 0)
    assert source.phase == 24 and np.all((2*values-source.exponents.astype(np.int64)) % 24 == 0)
    bad = load_source("fifth-5-q20")
    with pytest.raises(ValueError, match="do not embed"):
        embed_source(bad, 0)


def test_streamed_scalar_maps_equal_dense_api_and_palettes_are_matched_separately():
    frame = WALSH4.copy()
    frame[0] *= -1
    stats, hashes, denominators = geometry_data(frame, [])
    for exponent, (numerator, denominator) in scalar_blocks(frame, [], range(12)).items():
        assert stats[str(exponent)] == block_statistics(numerator, denominator)
        assert hashes[str(exponent)] == array_hash(numerator, denominator)
        assert denominators[str(exponent)] == denominator
    outcomes = match_palettes(stats, {"q12": list(PALETTE), "sextic": [0, 2, 6, 8]})
    assert outcomes["q12"]["all_used_roots_flat"]
    assert not outcomes["sextic"]["all_used_roots_flat"]
    assert {w["exponent"] for w in outcomes["sextic"]["failure_witness_for_every_failed_used_root"]} == {2, 8}


@pytest.mark.parametrize("identifier,twelfth", [("unreal13-q6", False), ("q12", True)])
def test_complete_two_phase_positive_controls_have_correct_expanded_grams(identifier, twelfth):
    source = fixed_source() if twelfth else load_source(identifier)
    exponents = embed_source(source, 0)
    frame = WALSH4.copy()
    if twelfth:
        frame[0] *= -1
    blocks = {int(e): reconstruct_block(frame, [], int(e)) for e in np.unique(exponents)}
    assert all(np.all(np.abs(numerator) == denominator) for numerator, denominator in blocks.values())
    matrix = np.block([[blocks[int(exponents[i, j])][0]//blocks[int(exponents[i, j])][1]
                        for j in range(source.order)] for i in range(source.order)])
    assert np.all(np.abs(matrix) == 1)
    assert np.array_equal(matrix@matrix.T, len(matrix)*np.eye(len(matrix), dtype=np.int64))
    assert len(matrix) == (28 if twelfth else 52)


def test_completed_broad_batch_counts_all_used_roots_and_replays_actual_source_entries():
    path = Path("certificates/scalar-lift/broad-phase-basis-probe.json")
    report = json.loads(path.read_text())
    assert report["candidate_geometries"] == 512
    assert report["source_matrix_evaluations"] == 39104
    assert report["covered_preflight_routes"] == 198
    assert report["literal_distinct_embedded_source_bindings"] == 234
    assert report["covered_declared_source_phases"] == [6, 12, 24]
    assert report["pending_declared_geometries"] == 0
    assert not report["new_sign_Hadamard_found"]
    for case in report["geometries"]:
        for node in case["nodes"]:
            assert set(node["antipodal_scalar_statistics"]) == set(map(str, range(6)))
            assert len(node["source_sign_entry_counts_in_binding_order"]) == len(case["binding_ids"])
    saved = SavedBroadBasisNode(path, geometry=0, node=0)
    binding_id = next(iter(saved.bindings))
    binding = saved.bindings[binding_id]
    for i, j, a, b in ((0, 0, 0, 1), (binding["order"]-1, binding["order"]-1, 667, 667)):
        exponent = binding["literal_exponents"][i][j]
        numerator, denominator = reconstruct_block(saved.frame, saved.layers, exponent)
        assert saved.entry(binding_id, i*668+a, j*668+b) == {
            "numerator": int(numerator[a, b]), "denominator": denominator}
    with pytest.raises(ValueError, match="outside"):
        saved.entry(binding_id, binding["order"]*668, 0)

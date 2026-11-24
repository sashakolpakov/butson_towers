from __future__ import annotations

import json
from pathlib import Path

import pytest

from experiments.scalar_lift import large_butson_factory as factory
from experiments.scalar_lift.large_butson_descendants import SavedButsonMatrix


def test_fixed_catalogue_rejects_source_construction_branches():
    with pytest.raises(ValueError, match="four existing"):
        factory.fixed_source("character-97")


def test_fixed_source_tree_builds_actual_descendants_and_deduplicates(tmp_path):
    tree = factory.RestrictionTree(tmp_path / "tree.sqlite3", seeds=("conference-1000-q8",))
    try:
        first = tree.run(30, max_nodes=1)
        assert first["progress"]["completed_matrices"] == 1
        assert first["pending_recipes"] == 2
        report = tree.run(30)
        assert report["stop_reason"] == "registered_restriction_tree_complete"
        assert report["completed_order_phase_pairs"] == [(1000, 8), (2000, 4), (4000, 2)]
        assert report["progress"]["duplicate_recipes"] == 2
        assert report["real_output_orders"] == [4000]
        assert not report["new_real_Hadamard_order_claimed"]
        assert not report["all_scalar_restriction_methods_exhausted"]
        for node in report["nodes"]:
            assert Path(node["bundle"]).is_file()
            assert json.loads(Path(node["proof_file"]).read_text())["matrix_hash"] == node["matrix_hash"]
            assert node["phase"] == 2 or node["other_field_restrictions_open"]
            saved = SavedButsonMatrix(node["bundle"])
            assert (saved.order, saved.phase) == (node["order"], node["phase"])
            assert 0 <= saved.entry(0, 0) < saved.phase
            assert 0 <= saved.entry(saved.order-1, saved.order-1) < saved.phase
        with pytest.raises(RuntimeError, match="another worker"):
            factory.RestrictionTree(tree.path, seeds=("conference-1000-q8",))
    finally:
        tree.close()
    resumed = factory.RestrictionTree(tmp_path / "tree.sqlite3", seeds=("conference-1000-q8",))
    try:
        assert resumed.run(1)["batch_evaluated_recipes"] == 0
    finally:
        resumed.close()


def test_breadth_first_cursor_places_all_sources_before_descendants(tmp_path):
    tree = factory.RestrictionTree(tmp_path / "tree.sqlite3")
    try:
        report = tree.run(30, max_nodes=1)
        rows = tree.db.execute("SELECT recipe,depth FROM tasks WHERE status='pending' ORDER BY depth,rowid").fetchall()
        assert rows[0][1] == 0
        assert json.loads(rows[0][0])["source"] == "character-241"
        assert any(depth == 1 for _, depth in rows)
        assert report["progress"]["completed_matrices"] == 1
        assert report["nodes"][0]["order"] == 5329
    finally:
        tree.close()


def test_additional_verified_catalogue_input_builds_and_replays_actual_matrices(tmp_path):
    from experiments.scalar_lift.catalogue_source_acquisition import load_source

    def loader(identifier):
        source = load_source(identifier)
        return source, source.proof

    path = tmp_path / "catalogue.sqlite3"
    worker = factory.RestrictionTree(path, seeds=("catalogue-p7-q12",), source_loader=loader)
    try:
        report = worker.run(20)
        assert report["pending_recipes"] == 0
        assert report["completed_order_phase_pairs"] == [(7, 12), (14, 6)]
        for node in report["nodes"]:
            saved = SavedButsonMatrix(node["bundle"])
            assert saved.order == node["order"] and saved.phase == node["phase"]
        assert worker._get("identity")["verified_input_manifest"][0]["phase"] == 12
    finally:
        worker.close()
    resumed = factory.RestrictionTree(path, seeds=("catalogue-p7-q12",), source_loader=loader)
    try:
        assert resumed.run(1)["batch_evaluated_recipes"] == 0
    finally:
        resumed.close()


def test_additional_self_contained_source_replays_and_changed_input_cannot_resume(tmp_path):
    from dataclasses import replace

    import numpy as np

    from experiments.scalar_lift.catalogue_source_acquisition import load_source
    from experiments.scalar_lift.forward_butson_descent import exponent_hash, verify_source

    original = replace(load_source("catalogue-p7-q12"), id="fixed-extra-input")

    def loader(_identifier):
        return original, original.proof

    path = tmp_path / "extra.sqlite3"
    worker = factory.RestrictionTree(path, seeds=(original.id,), source_loader=loader)
    try:
        report = worker.run(20)
        for node in report["nodes"]:
            SavedButsonMatrix(node["bundle"])
    finally:
        worker.close()
    shifted = (original.exponents.astype(np.int64)+1) % original.phase
    changed = replace(original, exponents=shifted, source_hash=exponent_hash(shifted),
                      proof=verify_source(shifted, original.phase))
    with pytest.raises(ValueError, match="another fixed source catalogue"):
        factory.RestrictionTree(path, seeds=(original.id,),
                                source_loader=lambda _identifier: (changed, changed.proof))

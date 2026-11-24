from __future__ import annotations

from functools import cache

import numpy as np
import pytest

from experiments.scalar_lift import weighted_scalar_tree as tree
from experiments.scalar_lift.forward_character_sources import character_source, formal_gram_check
from experiments.scalar_lift.relative_scalar_restriction import edge, restriction_in_basis


def test_declared_field_graph_and_path_counts():
    @cache
    def paths(phase):
        return 1 + sum(paths(q) for q in tree.proper_fields(phase))

    assert tree.proper_fields(24) == [1, 3, 4, 8, 12]
    assert paths(24) == 16 and paths(120) == 88
    assert tree.proper_fields(1) == []


def test_composition_and_metric_congruence_do_not_merge_coordinate_distinct_matrices():
    source = character_source(5, 4)
    parent, relative = edge(12, 4), edge(4, 1)
    basis = tuple(a+3*b for a in parent.basis_exponents for b in relative.basis_exponents)
    child = restriction_in_basis(12, 1, basis)
    composition = tree.verify_composition(parent, relative, child)
    assert composition["every_action_coefficient_composed_exactly"]
    assert composition["metric_trace_transitivity_checked_coefficientwise"]
    assert tree.canonical_basis_certificate(child)["Hermitian_metric_congruence_exact"]
    canonical = edge(12, 1)
    assert not np.array_equal(child.actions.numerators, canonical.actions.numerators)
    assert tree.matrix_identity(source, child) != tree.matrix_identity(source, canonical)


def test_saved_small_endpoint_has_the_complete_weighted_gram_and_literal_zeros(tmp_path):
    source = character_source(5, 4)
    restriction = restriction_in_basis(12, 1, (0, 3, 1, 4))
    path = tmp_path / "weighted.npz"
    proof = {**formal_gram_check(source), "full_Gram_verified": True}
    metadata = tree.save_node(path, source, proof, restriction, [12, 4, 1], {})
    saved = tree.SavedWeightedMatrix(path)
    assert saved.order == 100 and saved.lower_field == 1
    assert all(saved.replay().values())
    values = np.array([[saved.entry(i, j)["numerators"][0] for j in range(saved.order)]
                       for i in range(saved.order)], dtype=np.int64)
    weight = np.kron(np.eye(source.order, dtype=np.int64), saved.inverse_metric[:, :, 0])
    denominator = metadata["action_denominator"]
    assert np.array_equal(values @ weight @ values.T, source.order*denominator**2*weight)
    assert metadata["rational_endpoint"]["complete_matrix_entry_counts"]["0"] == int(np.count_nonzero(values == 0))
    assert not metadata["rational_endpoint"]["literal_entries_all_signs"]
    assert saved.weight_entry(0, saved.width)["numerators"] == [0]
    with pytest.raises(ValueError, match="outside"):
        saved.entry(saved.order, 0)


def test_checkpoint_preserves_weighted_nodes_and_resumes_without_repeating(tmp_path):
    path = tmp_path / "tree.sqlite3"
    worker = tree.WeightedTree(path, seeds=("character-73",))
    try:
        first = worker.run(20, max_nodes=1)
        assert first["progress"]["literal_nodes"] == 1
        assert first["pending_tasks"] == 5
        with pytest.raises(RuntimeError, match="another worker"):
            tree.WeightedTree(path, seeds=("character-73",))
    finally:
        worker.close()
    worker = tree.WeightedTree(path, seeds=("character-73",))
    try:
        report = worker.run(20)
        assert report["stop_reason"] == "registered_divisor_field_tree_complete"
        assert report["pending_tasks"] == 0
        assert report["rational_weighted_orders"] == [42632]
        assert 6 <= report["progress"]["literal_nodes"] <= 16
        for node in report["nodes"]:
            saved = tree.SavedWeightedMatrix(node["bundle"])
            assert saved.order == node["matrix_order"]
            assert saved.entry(0, 0)["lower_field"] == node["lower_field"]
        assert worker.run(1)["batch_nodes"] == 0
    finally:
        worker.close()


def test_existing_dense_source_has_correct_entries_histogram_and_weighted_gram(tmp_path):
    assert len(tree.retained_registry_seeds()) == 18
    source, proof = tree.load_tree_source("unreal3-q6")
    restriction = edge(6, 1)
    path = tmp_path / "dense-source.npz"
    metadata = tree.save_node(path, source, proof, restriction, [6, 1], {})
    saved = tree.SavedWeightedMatrix(path)
    assert saved.order == 6
    assert all(saved.replay().values())
    values = np.array([[saved.entry(i, j)["numerators"][0] for j in range(6)] for i in range(6)], dtype=np.int64)
    weight = np.kron(np.eye(3, dtype=np.int64), saved.inverse_metric[:, :, 0])
    assert np.array_equal(values @ weight @ values.T,
                          3*metadata["action_denominator"]**2*weight)
    assert sum(metadata["rational_endpoint"]["complete_matrix_entry_counts"].values()) == 36
    assert tree.proper_fields(6) == [1]  # Q(zeta3)=Q(zeta6) is not a proper subfield.
    gate = tree.completed_inventory_sign_gate(("unreal3-q6",))
    assert gate["sources"][0]["not_excluded_by_saved_alphabet_gate"]


def test_acquired_twelfth_root_input_runs_and_replays_every_branch(tmp_path):
    worker = tree.WeightedTree(tmp_path / "catalogue.sqlite3", seeds=("catalogue-p7-q12",))
    try:
        report = worker.run(20)
        assert report["pending_tasks"] == 0
        assert report["rational_weighted_orders"] == [28]
        assert report["progress"]["literal_nodes"] > 3
        assert report["uniform_sign_gate_after_inventory"]["sources"][0]["not_excluded_by_saved_alphabet_gate"]
        for node in report["nodes"]:
            saved = tree.SavedWeightedMatrix(node["bundle"])
            assert all(saved.replay().values())
            if saved.lower_field != 1:
                continue
            values = np.array([[saved.entry(i, j)["numerators"][0]
                                for j in range(saved.order)] for i in range(saved.order)])
            weight = np.kron(np.eye(7, dtype=np.int64), saved.inverse_metric[:, :, 0])
            assert np.array_equal(values @ weight @ values.T,
                                  7*saved.metadata["action_denominator"]**2*weight)
            assert not node["rational_endpoint"]["literal_entries_all_signs"]
    finally:
        worker.close()

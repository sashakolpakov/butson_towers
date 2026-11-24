import json
from pathlib import Path

import numpy as np
import pytest

from experiments.scalar_lift.forward_butson_descent import ForwardTower
from experiments.scalar_lift.forward_character_sources import character_source
from experiments.scalar_lift.large_butson_descendants import (
    ImplicitForwardTower,
    SavedButsonMatrix,
)
from experiments.scalar_lift.vanishing_collision_search import integral_roots

ROOT = Path(__file__).resolve().parents[1]


def test_existing_small_character_control_has_complete_descendant_Gram():
    source = character_source(5, 4)
    tower = ImplicitForwardTower(source)
    assert tower.orders == [25, 50]
    assert tower.phases == [12, 6]
    matrix = np.array([[tower.entry(i, j) for j in range(50)] for i in range(50)])
    roots = np.asarray(integral_roots(6), dtype=np.int64)
    gram = roots[(matrix[:, None, :] - matrix[None, :, :]) % 6].sum(axis=2)
    expected = np.zeros_like(gram)
    expected[np.arange(50), np.arange(50), 0] = 50
    assert np.array_equal(gram, expected)
    source_dense = np.array([[source.entry(divmod(i, 5), divmod(j, 5)) for j in range(25)] for i in range(25)])
    assert np.array_equal(matrix, ForwardTower(source_dense, 12).materialize())
    assert tower.used_alphabet == tuple(map(int, np.unique(matrix)))


def test_all_four_actual_large_output_files_reload_with_exact_proofs_and_entries():
    report = json.loads((ROOT / "certificates/scalar-lift/large-butson-descendants.json").read_text())
    orders = []
    for row in report["outputs"]:
        for artifact in row["actual_matrix_files"]:
            path = Path(artifact["path"])
            saved = SavedButsonMatrix(path if path.is_absolute() else ROOT / path)
            orders.append((saved.order, saved.phase))
            assert saved.metadata["matrix_hash"] == artifact["matrix_hash"]
            for i, j, exponent in saved.metadata["stages"][-1]["twenty_literal_entries"]:
                assert saved.entry(i, j) == exponent
            assert saved.metadata["dense_large_matrices_materialized"] is False
    assert orders == [(10658, 12), (21316, 6), (116162, 60), (232324, 30)]


def test_saved_scalar_definition_is_self_contained_and_rejects_corruption(tmp_path):
    source = character_source(73)
    tower = ImplicitForwardTower(source, steps=[(2, "primary")])
    path = tmp_path / "matrix.npz"
    tower.save_bundle(path)
    saved = SavedButsonMatrix(path)
    assert saved.entry(10657, 53) == tower.entry(10657, 53)
    with pytest.raises(ValueError, match="outside"):
        saved.entry(1.5, 2)
    with np.load(path, allow_pickle=False) as stored:
        arrays = {key: stored[key].copy() for key in stored.files}
    arrays["scalar_table_1"][0, 0, 0] ^= 1
    bad = tmp_path / "modified.npz"
    np.savez_compressed(bad, **arrays)
    with pytest.raises(ValueError, match="scalar table was modified"):
        SavedButsonMatrix(bad)


def test_exact_definition_hash_deduplicates_identical_paths_and_keeps_real_alphabet():
    source = character_source(73)
    first = ImplicitForwardTower(source, steps=[(2, "primary"), (2, "primary")])
    second = ImplicitForwardTower(source, steps=[{"prime": 2, "generator_mode": "primary"},
                                                {"prime": 2, "generator_mode": "primary"}])
    assert first.matrix_hash == second.matrix_hash
    assert first.metadata()["stages"][-1]["real_sign_endpoint"] is False
    root = ImplicitForwardTower(source, steps=[])
    assert root.order == 5329 and root.phase == 24
    assert root.entry(17, 5311) == source.entry(divmod(17, 73), divmod(5311, 73))

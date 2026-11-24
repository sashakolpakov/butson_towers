from __future__ import annotations

import json

import pytest

from experiments.scalar_lift import forward_design_frame_search as search
from experiments.scalar_lift.joint_t_sequence import normalize_orbit_scaling
from hadamard.inverse_lift import is_sign_hadamard


def small_control():
    states, _ = normalize_orbit_scaling([0, 1, 1, 2, 6, 4, 6], 7, 1)
    cell = {"prime": 7, "multiplier": 1, "row_sums": [-1, 1, 1, 2],
            "fixed_orbit_states": [[2, states[2]], [3, states[3]]]}
    return cell, states


def test_finite_cells_are_exact_disjoint_partitions_and_sources_are_complete():
    cells = search.cells_for_chart()
    assert len(cells) == 192 and len({search.digest(c) for c in cells}) == 192
    assert len(search.cells_for_chart(491, 7)) == 64
    for sums in {tuple(c["row_sums"]) for c in cells}:
        assert len([c for c in cells if tuple(c["row_sums"]) == sums]) == 64
    for source in search.source_certificate():
        assert source["full_Gram_verified"]
        assert set(source["used_alphabet"]) <= {1, 2, 4, 5}
        assert source["registered_source_exponent_shift"] == 1


def test_small_exact_candidate_constructs_literal_complete_forward_endpoints():
    cell, states = small_control()
    witness = search.verify_candidate(cell, states)
    assert [r["order"] for r in witness["real_outputs"]] == [364, 84]
    for n in (13, 3):
        order = 28 * n
        matrix = [[search.endpoint_entry(witness, i, j, n) for j in range(order)] for i in range(order)]
        assert is_sign_hadamard(matrix)
    bad = dict(cell, fixed_orbit_states=[[2, (states[2]+1) % 8], [3, states[3]]])
    with pytest.raises(ValueError, match="declared cell"):
        search.verify_candidate(bad, states)


def test_native_small_cell_positive_then_checkpoint_replay_does_not_repeat(tmp_path):
    pytest.importorskip("ortools.sat.python.cp_model")
    cell, _ = small_control()
    path = tmp_path / "positive.sqlite3"
    worker = search.DesignFrameSearch(path, cells=[cell], cell_seconds=1)
    try:
        report = worker.run(4, max_attempts=1)
        assert report["cell_statuses"] == {"verified": 1}
        assert report["stop_reason"] == "verified_frame_found"
        assert not report["new_order_priority_established"]
    finally:
        worker.close()
    worker = search.DesignFrameSearch(path, cells=[cell])
    try:
        assert worker.run(1)["batch_attempts"] == 0
        assert worker.retry_open(1) == 0
    finally:
        worker.close()


def test_unknown_attempt_is_durable_open_and_requires_explicit_new_seed(tmp_path):
    pytest.importorskip("ortools.sat.python.cp_model")
    cell, _ = small_control()
    path = tmp_path / "unknown.sqlite3"
    worker = search.DesignFrameSearch(path, cells=[cell], cell_seconds=1e-9)
    try:
        report = worker.run(3, max_attempts=1)
        assert report["cell_statuses"] == {"open": 1}
        assert report["attempt_statuses"] == {"unknown": 1}
        assert worker.run(1)["batch_attempts"] == 0
        assert worker.retry_open(0) == 0
        assert worker.retry_open(9) == 1
        assert worker.retry_open(9) == 0
    finally:
        worker.close()
    worker = search.DesignFrameSearch(path, cells=[cell], cell_seconds=1e-9)
    try:
        assert worker.report()["attempt_statuses"] == {"pending": 1, "unknown": 1}
        with pytest.raises(RuntimeError, match="another worker"):
            search.DesignFrameSearch(path, cells=[cell])
    finally:
        worker.close()


def test_interrupted_attempt_preserves_parameters_and_saved_candidate_is_verified(tmp_path):
    cell, states = small_control()
    path = tmp_path / "interrupted.sqlite3"
    worker = search.DesignFrameSearch(path, cells=[cell])
    key = search.digest(cell)
    with worker.db:
        worker.db.execute("UPDATE attempts SET status='running',result=?", (json.dumps({"parameters": {"random_seed": 0}}),))
        worker.db.execute("UPDATE cells SET status='running'")
    worker.close()
    worker = search.DesignFrameSearch(path, cells=[cell])
    try:
        assert worker.report()["attempt_statuses"] == {"interrupted": 1}
        raw = json.loads(worker.db.execute("SELECT result FROM attempts").fetchone()[0])
        assert raw["parameters"] == {"random_seed": 0}
        assert raw["explicit_new_seed_required"]
        with worker.db:
            worker.db.execute("UPDATE attempts SET status='candidate',result=? WHERE cell_id=?",
                              (json.dumps({"orbit_states": states}), key))
            worker.db.execute("UPDATE cells SET status='candidate'")
    finally:
        worker.close()
    worker = search.DesignFrameSearch(path, cells=[cell])
    try:
        assert worker.run(3)["cell_statuses"] == {"verified": 1}
    finally:
        worker.close()

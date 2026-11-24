import json
from pathlib import Path

import numpy as np
import pytest
from sympy import Matrix, eye, kronecker_product

from experiments.scalar_lift.alpoge_global_basis_probe import WALSH4
from experiments.scalar_lift.q12_basis_probe import (
    PALETTE,
    SavedQ12BasisNode,
    candidate_record,
    fixed_source,
    scalar_blocks,
    twice_local_actions,
)


def dense_layer(permutation):
    values = np.zeros((len(permutation), len(permutation)), dtype=np.int64)
    for group in permutation.reshape(-1, 4):
        values[np.ix_(group, group)] = WALSH4
    return Matrix(values.tolist())/2


def test_all_twelfth_root_actions_equal_dense_conjugation_and_orthogonal_grams():
    frame = np.kron(WALSH4, np.array([[1, 1], [1, -1]], dtype=np.int64))
    layers = [np.array([0, 3, 6, 1, 5, 2, 7, 4]), np.array([6, 2, 0, 7, 5, 3, 1, 4])]
    change = eye(8)
    for permutation in layers:
        change = dense_layer(permutation)*change
    blocks = scalar_blocks(frame, layers, range(12))
    assert set(blocks) == set(range(6))
    actions = twice_local_actions()
    root = Matrix(actions[1].tolist())/2
    J = Matrix(actions[3].tolist())/2
    assert root**3 == J and J**2+eye(4) == 0*eye(4)
    assert root**2-J*root-eye(4) == 0*eye(4)
    assert root**2-root+eye(4) != 0*eye(4)  # this is not a sextic substitute
    for exponent, (numerator, denominator) in blocks.items():
        block = Matrix(numerator.tolist())/denominator
        action = kronecker_product(Matrix(actions[exponent].tolist())/2, eye(2))
        assert block == change*action*change.T*Matrix(frame.tolist())
        assert block*block.T == 8*eye(8)


def test_genuine_q12_positive_control_builds_full_sign_hadamard28():
    source = fixed_source()
    assert source.phase == 12 and source.used_alphabet == PALETTE
    assert source.proof["fixed_starting_representative"]["global_exponent_shift"] == -1
    frame = WALSH4.copy()
    frame[0] *= -1
    blocks = scalar_blocks(frame, [])
    matrix = np.block([[
        (-1 if source.entry(i, j) >= 6 else 1)
        *blocks[source.entry(i, j) % 6][0]//blocks[source.entry(i, j) % 6][1]
        for j in range(source.order)] for i in range(source.order)])
    assert matrix.shape == (28, 28)
    assert np.all(np.abs(matrix) == 1)
    assert np.array_equal(matrix@matrix.T, 28*np.eye(28, dtype=np.int64))
    record = candidate_record(source, "control", frame, [], {"depth": 0, "index": 0, "seed": 0})
    assert record["all_used_scalars_flat"]
    assert record["complete_restricted_matrix_sign_entries"] == 28**2


def test_fourth_root_success_does_not_hide_failed_twelfth_root_blocks():
    record = candidate_record(fixed_source(), "control", WALSH4, [],
                              {"depth": 0, "index": 0, "seed": 0})
    stats = record["all_used_scalar_statistics"]
    assert stats["0"]["all_entries_signs"] and stats["3"]["all_entries_signs"]
    assert not stats["1"]["all_entries_signs"] and not stats["4"]["all_entries_signs"]
    assert not record["all_used_scalars_flat"]
    assert set(map(int, stats)) == set(PALETTE)


def test_completed_q12_batch_and_saved_actual_entries():
    path = Path("certificates/scalar-lift/q12-basis-probe.json")
    report = json.loads(path.read_text())
    assert report["candidate_visits"] == report["literal_distinct_nodes"] == 128
    assert report["pending_declared_candidates"] == 0
    assert not report["new_sign_Hadamard_found"]
    for case in report["cases"]:
        for node in case["candidates"]:
            stats = node["all_used_scalar_statistics"]
            assert set(map(int, stats)) == set(PALETTE)
            for values in stats.values():
                assert sum(values[k] for k in ("sign_entries", "zero_entries", "other_integer_entries", "noninteger_entries")) == values["total_entries"]
    node = SavedQ12BasisNode(path, case=0, candidate=0)
    assert node.order == 4676
    for i, j, a, b in ((0, 0, 0, 1), (1, 2, 17, 19), (6, 6, 667, 667)):
        exponent = node.source.entry(i, j)
        num, den = node.blocks[exponent % 6]
        expected = (-1 if exponent >= 6 else 1)*int(num[a, b])
        assert node.entry(i*668+a, j*668+b) == {"numerator": expected, "denominator": den}
    with pytest.raises(ValueError, match="outside"):
        node.entry(node.order, 0)

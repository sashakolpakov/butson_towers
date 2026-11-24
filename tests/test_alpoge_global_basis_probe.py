import json
from pathlib import Path

import numpy as np
from sympy import Matrix

from experiments.scalar_lift.alpoge_global_basis_probe import (
    WALSH4,
    SavedGlobalBasisNode,
    apply_layer,
    global_partner,
    global_rank_lower_bound,
)
from experiments.scalar_lift.alpoge_mixed_basis_probe import twice_root_action


def dense_layer(permutation):
    result = np.zeros((len(permutation), len(permutation)), dtype=np.int64)
    for group in permutation.reshape(-1, 4):
        result[np.ix_(group, group)] = WALSH4
    return Matrix(result.tolist())


def test_global_layers_equal_independent_dense_rational_conjugation():
    frame = np.kron(WALSH4, np.array([[1, 1], [1, -1]], dtype=np.int64))
    layers = [np.array([0, 3, 6, 1, 5, 2, 7, 4]), np.array([6, 2, 0, 7, 5, 3, 1, 4])]
    change = Matrix.eye(8)
    for layer in layers:
        T = dense_layer(layer)
        assert 4*Matrix.eye(8) == T*T.T
        change = T*change/2
    numerator, denominator = global_partner(frame, layers)
    partner = Matrix(numerator.tolist())/denominator
    root = Matrix(twice_root_action(np.eye(8, dtype=np.int64)).tolist())/2
    F = Matrix(frame.tolist())
    assert partner == change*root*change.T*F
    assert partner*partner.T == 8*Matrix.eye(8)
    assert partner*F.T+F*partner.T == 8*Matrix.eye(8)
    assert global_rank_lower_bound(8, layers) <= (change-Matrix.eye(8)).rank()
    values = np.arange(24, dtype=np.int64).reshape(8, 3)
    assert Matrix(apply_layer(values, layers[0]).tolist()) == dense_layer(layers[0])*Matrix(values.tolist())


def test_width_four_positive_control_accepts_a_sign_partner():
    numerator, denominator = global_partner(WALSH4, [np.arange(4), np.arange(4)])
    assert np.all(np.abs(numerator) == denominator)


def test_completed_global_batch_clears_low_rank_gate_and_replays_source_entries():
    path = Path("certificates/scalar-lift/alpoge-global-basis-probe.json")
    report = json.loads(path.read_text())
    assert report["candidate_visits"] == report["distinct_partners"] == 128
    assert report["pending_declared_candidates"] == 0
    assert not report["new_sign_Hadamard_found"]
    for case in report["cases"]:
        assert case["all_candidates_certifiably_above_low_rank_cutoff"]
        for node in case["candidates"]:
            assert node["rank_Q_minus_I_certified_lower_bound"] >= case["minimum_rank_Q_minus_I_for_any_sign_partner"]
            assert sum(node["complete_numerator_value_counts"].values()) == case["published_frame_order"]**2
            assert node["nonproduct_output_basis_column"]["integer_minor_determinant"]
    node = SavedGlobalBasisNode(path, case=0, candidate=0)
    source = node.sources["unreal13-q6"]
    for i, j, a, b in ((0, 0, 0, 1), (1, 2, 17, 19), (12, 12, 667, 667)):
        scalar = source.entry(i, j)
        f, g = int(node.frame[a, b])*node.denominator, int(node.numerator[a, b])
        expected = (f, g, g-f, -f, -g, f-g)[scalar]
        assert node.entry(source.id, i*668+a, j*668+b) == {"numerator": expected, "denominator": node.denominator}

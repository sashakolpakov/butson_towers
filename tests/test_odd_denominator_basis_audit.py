import json
from itertools import product
from pathlib import Path

import numpy as np

from experiments.scalar_lift.odd_denominator_basis_audit import (
    audit_frame,
    column_factorization,
    exhaustive_small_control,
    kernel_sign_points,
    native,
    support_graph,
)


def test_modular_column_factorization_and_compressed_K_rank_are_exact():
    frame = np.array([[1, 0, 1, 1], [0, 1, 1, -1],
                      [1, 1, -1, 0], [1, -1, 0, -1]], dtype=np.int64)
    columns, V, reduced, pivots = column_factorization(frame, 3)
    T = np.diag([1, 1, -1, -1])
    full = native(frame.T@T@frame, 3)
    compressed = V.transpose()*native(T@columns, 3)
    assert full.rank() == compressed.rank() == 2
    assert native(frame, 3) == V*native(reduced, 3)
    assert pivots == [0, 1]


def test_graph_stabilizers_and_small_kernel_sign_points_match_all16_vectors():
    reduced = np.array([[1, 0, 1, 1], [0, 1, 1, 1]], dtype=np.uint16)
    graph = support_graph(reduced, [0, 1])
    assert graph["component_count"] == 1
    free, points = kernel_sign_points(reduced, [0, 1], 3)
    assert free == [2, 3]
    enumerated = {tuple(values) for _signs, values in points}
    direct = {values for values in product((-1, 1), repeat=4)
              if not np.any((reduced.astype(np.int64)@values) % 3)}
    assert direct == enumerated == {(-1, -1, -1, -1), (1, 1, 1, 1)}
    disconnected = support_graph(np.array([[1, 0, 1, 0], [0, 1, 0, 1]], dtype=np.uint16), [0, 1])
    assert disconnected["component_count"] == 2
    control = exhaustive_small_control()
    assert control["sign_choices_checked"] == 16 and control["all_three_tests_agree"]


def test_all_eight_exact_families_close_without_hiding_rank_defects():
    report = json.loads(Path("certificates/scalar-lift/odd-denominator-basis-audit.json").read_text())
    assert report["frame_count"] == report["families_closed"] == 8
    assert report["pending_declared_audits"] == 0 and report["open_frames"] == []
    for frame in report["frames"]:
        assert frame["modular_F_rank"] == frame["required_half_rank"]
        assert frame["complete_column_factorization_verified"]
        assert frame["kernel_support_graph"]["constant_signs_are_only_stabilizers"]
        assert frame["all_non_two_power_root_families_closed"]
        assert all(case["family_closed"] for case in frame["root_cases"].values())
        assert all(case["family_closed"] for case in frame["related_root_cases"].values())
    for index in (0, 2):
        frame = report["frames"][index]
        assert not frame["root_cases"]["1"]["column_K_equals_kernel_F"]
        for case in (*frame["root_cases"].values(), *frame["related_root_cases"].values()):
            gate = case["small_linear_gate"]
            assert gate["constraint_nullity"] == 3
            assert gate["free_sign_assignments_checked"] == 8
            assert gate["all_kernel_sign_points_exhausted"]
            assert gate["all_sign_points_excluded"]
            assert all(p["baseline_binary_rank"] > 2*p["rank_Q_minus_global_sign_I"]
                       for p in gate["sign_points"])


def test_full_deficient_frame_replay_includes_transposed_primitive_root():
    report = json.loads(Path("certificates/scalar-lift/odd-denominator-basis-audit.json").read_text())
    assert audit_frame("h668-width4") == report["frames"][0]
    transpose = report["frames"][0]["related_root_cases"]["5"]
    assert transpose["small_linear_gate"]["all_sign_points_excluded"]

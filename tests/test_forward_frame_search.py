import itertools

import numpy as np

from experiments.scalar_lift.fifth_source_assembly import exact_column_gram
from experiments.scalar_lift.forward_frame_search import (
    check_frame,
    frame_block,
    frame_options,
    frame_proof,
    global_affine_fit,
    row_column_gauge_search,
)
from experiments.scalar_lift.phase120_common_base import source_fifteen


def test_all_frame_certificates_and_cross_grams():
    for phase in (3, 5, 6, 10, 12, 20, 24, 40, 48, 80, 96, 120):
        for record in frame_options(phase):
            proof = frame_proof(record["id"])
            assert proof["all_cross_Grams_follow"]
            assert proof["admitted"] == record["admitted"]
            frame = frame_block(record["id"], 0)
            for a, b in ((0, 1), (1, phase - 1), (phase // 2, phase // 3)):
                assert np.array_equal(
                    frame_block(record["id"], a) @ frame_block(record["id"], b).T,
                    frame_block(record["id"], a - b) @ frame.T,
                )
            assert check_frame(phase, record["admitted"], record["id"])["fits"]
            assert not check_frame(phase, range(phase), record["id"])["fits"]


def test_regauged_actual_source_reaches_exact_sign_endpoint():
    source = np.asarray(source_fifteen(), dtype=np.int64)
    n = len(source)
    scrambled = (source + np.arange(n)[:, None] + (np.arange(n)**2)[None, :]) % 120
    assert global_affine_fit(120, set(map(int, scrambled.flat)), "common-120-raw")["status"] == "EXHAUSTED"
    result = row_column_gauge_search(scrambled, 120, "common-120-raw", seconds=5)
    assert result["status"] == "FOUND"
    output = np.array(result["exponents"])
    gram = exact_column_gram(output, 120)
    expected = np.zeros_like(gram)
    expected[np.arange(n), np.arange(n), 0] = n
    assert np.array_equal(gram, expected)
    real = np.block([[frame_block("common-120-raw", e) for e in row] for row in output])
    assert np.all(np.abs(real) == 1)
    assert np.array_equal(real @ real.T, len(real) * np.eye(len(real), dtype=np.int64))


def test_finite_gauge_exclusion_and_budget_scope():
    tails = list(itertools.product(range(3), repeat=2))
    colors = np.array([[0] * 9, [a for a, _ in tails], [b for _, b in tails]])
    result = row_column_gauge_search(colors, 24, "cubic-24", seconds=3)
    assert result["status"] == "INFEASIBLE"
    assert row_column_gauge_search(colors, 24, "cubic-24", max_entries=1)["status"] == "NOT_RUN"
    assert row_column_gauge_search(colors, 24, "cubic-24", seconds=0)["status"] == "UNKNOWN"


def test_global_gauge_and_galois_feasibility_equivalence():
    # An invertible color multiplier changes the excluded residue; a global
    # translation restores it. Thus full row/column gauges need only u=1.
    for p in (3, 5):
        for forbidden in range(p):
            for unit in range(1, p):
                shift = forbidden - unit * forbidden
                for value in range(p):
                    assert (value != forbidden) == ((unit * value + shift) % p != forbidden)

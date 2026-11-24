import json
from pathlib import Path

import numpy as np
import pytest

from experiments.scalar_lift.rectangular_petrescu_sources import (
    boundary_problem,
    branch_parameters,
    branch_prefilter,
    certificate,
    metric_certificate,
    root_pair_control,
)

ROOT = Path(__file__).resolve().parents[1]


def test_actual_smaller_sources_and_finite_branch_catalogue_replay():
    report = certificate()
    assert report == json.loads((ROOT / "certificates/scalar-lift/rectangular-petrescu-sources.json").read_text())
    assert (report["finite_branch_count"], report["excluded_finite_branches"],
            report["open_dense_completion_branches"]) == (36, 30, 6)
    assert not report["new_characteristic_zero_source_or_real_order_constructed"]


@pytest.mark.parametrize(("order", "selected", "delta"), [(567, 466, 365), (729, 385, 41)])
def test_independent_admitted_frames_and_actual_length_root_pairs(order, selected, delta):
    descriptor = {"kind": "stride", "step": 2, "shift": 7}
    data = boundary_problem(order, descriptor)
    assert data["T_L"].shape == data["T_R"].shape == (selected, order)
    assert np.array_equal((-data["T_R"].T) % 120, data["T_L"].T)
    assert not np.any(data["T_L"] % 3 == 2)
    assert not np.any(data["T_L"] % 5 == 0)
    row = root_pair_control(order, descriptor)
    assert row["full_two_by_two_Gram_of_length1499_rows_verified"]
    assert row["E_row_squared_norm"] == delta and row["F_row_squared_norm"] == 1499


def test_nonuniform_metric_gate_and_order_parity():
    with pytest.raises(ValueError, match="same parity"):
        branch_parameters(1499, 500)
    result = branch_prefilter(1499, 729, {"step": 1},
                             {"kind": "normalized_small_Butson_blocks", "block_order": 35})
    assert result["status"] == "EXCLUDED_ROOT_SUM_INTEGRAL_NORM"
    assert branch_parameters(1499, 665)["known_input_constructor_available"] is False


def test_new_585_actual_source_row_gate_and_essential_fifth_metric_replay():
    report = certificate((585,))
    saved = json.loads((ROOT / "certificates/scalar-lift/rectangular-petrescu-585.json").read_text())
    assert report == saved
    branch = report["branches"][0]
    assert branch["parameters"]["selected_rows"] == 457
    assert branch["parameters"]["omitted_rows"] == 128
    assert branch["known_input"]["complete_585_Gram_verified_by_tensor_identity"]
    assert branch["known_input"]["factor117_exact_Gram"]["all_row_pairs_checked"] == 6786
    assert branch["actual_root_pair"]["unit_root_sum_positions"] == 329
    assert branch["actual_root_pair"]["zero_root_sum_positions"] == 128
    metric = metric_certificate(585)
    assert metric["selected_scalar_subfield_coefficients"] == [18, 0, 0, 2, 0, -1, 0, 2]
    assert metric["that_automorphism_changes_witness_to_distinct_conjugate"]
    assert metric["phase24_metric_obstruction"]["determinant_norm_prime_valuation"] == 457
    assert (report["finite_branch_count"], report["excluded_finite_branches"],
            report["open_dense_completion_branches"]) == (12, 9, 3)

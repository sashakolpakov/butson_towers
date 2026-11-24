import hashlib
import json
from pathlib import Path

import numpy as np

from experiments.scalar_lift.cubic_core_support_lift import (
    direct_small_control,
    full_bitset_check,
    large_problem,
    phase_coordinates,
    root_candidate_exponents,
    source_layers,
    template,
)
from experiments.scalar_lift.mixed_cubic_ramification import derivative


def test_template_and_full_bitset_derivative_against_direct_field_control():
    record = direct_small_control()
    h0, h1 = np.array(record["H0"]), np.array(record["H1"])
    check = full_bitset_check(h0, h1)
    assert check["all_row_pairs_including_diagonal_checked"] == 23*24//2
    assert check["full_first_Gram_derivative_zero"]
    # An actual changed phase bit must be noticed by both independent checks.
    changed = h1.copy()
    changed[0, 0] ^= int(h0[0, 0])
    assert np.any(derivative(h0, changed))
    try:
        full_bitset_check(h0, changed)
    except ValueError as error:
        assert "derivative" in str(error)
    else:
        raise AssertionError("a broken first Gram derivative was accepted")


def test_large_constructions_match_the_fully_verified_artifacts():
    path = Path(__file__).resolve().parents[1] / "certificates/scalar-lift/cubic-core-support-lift.json"
    report = json.loads(path.read_text())
    assert report["template"] == {key: value.tolist() for key, value in template().items()}
    for route in report["routes"]:
        problem = large_problem(route["source_order"], route["phase"])
        assert problem["template_groups"] == route["template_groups"]
        assert problem["inactive_fillers"] == route["inactive_fillers"]
        assert problem["U_exponents"].tolist() == route["U_exponents"]
        assert problem["D_exponents"].tolist() == route["D_exponents"]
        assert problem["primitive_cross_ratio"] == route["boundary_primitive_cross_ratio"]
        h0, h1 = source_layers(problem)
        for name, matrix in (("H0", h0), ("H1", h1)):
            digest = hashlib.sha256(matrix.astype(np.uint8).tobytes()).hexdigest()
            assert digest == route["full_source_check"][f"{name}_uint8_sha256"]
        assert route["root_phases_unitary_mod_pi_squared"]
        assert not route["characteristic_zero_source_constructed"]


def test_actual_admitted_roots_reconstruct_the_saved_first_phase_layers():
    problem = large_problem(1499, 96)
    candidate = root_candidate_exponents(problem)
    colors, dyadic = phase_coordinates(candidate, 96)
    h0, h1 = source_layers(problem)
    assert np.array_equal(h0, np.ones_like(colors) ^ (3*colors))
    assert np.array_equal(h1, h0*(dyadic % 2))
    assert np.all((colors == 0) | (colors == 1))
    assert candidate.shape == (1499, 1499)

from __future__ import annotations

import numpy as np
import pytest

from experiments.scalar_lift.joint_gs_223 import orbit_data
from experiments.scalar_lift.joint_t_sequence import (
    STATES,
    build_model,
    normalize_orbit_scaling,
    row_sum_types,
    verify_t_states,
)


def test_trivial_multiplier_model_accepts_a_verified_base_sequence_witness():
    cp_model = pytest.importorskip("ortools.sat.python.cp_model")
    # BS(4,3), with signed axis normalization, supplies these T7 states.
    states = [0, 1, 1, 2, 6, 4, 6]
    verify_t_states(states, prime=7, multiplier=1)
    normalized, _scale = normalize_orbit_scaling(states, prime=7, multiplier=1)
    assert (-1, 1, 1, 2) in row_sum_types(7, 1)
    model, variables, _metadata = build_model(7, 1)
    for variable, value in zip(variables, normalized, strict=True):
        model.add(variable == value)
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = 1
    solver.parameters.max_time_in_seconds = 2
    assert solver.solve(model) in (cp_model.FEASIBLE, cp_model.OPTIMAL)


def test_one_hot_control_gives_a_full_exact_sign_pair_and_rejects_a_bad_shift():
    verified = verify_t_states([0, 2, 4], prime=3, multiplier=1)
    pair = verified["full_GS_sign_pair_verification"]["full_sign_block_pair_independently_verified"]
    assert pair["valid"] and pair["block_order"] == 12
    with pytest.raises(ValueError, match="correlation failed"):
        verify_t_states([0, 0, 4], prime=3, multiplier=1)


def test_large_compressed_dot_equations_equal_direct_integer_correlations():
    prime = 223
    orbits, index, shifts, pairs, coefficients = orbit_data(prime, 3)
    states = [(i * i + 3 * i + 1) % 8 for i in range(len(orbits))]
    orbit_vectors = np.array([STATES[s] for s in states], dtype=np.int64)
    dots = np.array([int(orbit_vectors[i] @ orbit_vectors[j]) for i, j in pairs])
    expanded = orbit_vectors[index]
    direct = [int(np.sum(expanded * np.roll(expanded, -shift, axis=0))) for shift in shifts]
    assert np.array_equal(coefficients @ dots, direct)


def test_high_order_trivial_character_types_have_the_exact_norm_and_residues():
    assert row_sum_types(223, 3) == [(-14, 3, 3, 3), (-5, 6, 9, 9), (13, 3, 3, 6)]
    assert row_sum_types(491, 7) == [(1, 0, 7, 21)]
    for prime, multiplier in ((223, 3), (491, 7)):
        for row in row_sum_types(prime, multiplier):
            assert sum(x * x for x in row) == prime
            assert tuple(x % multiplier for x in row) == (1, 0, 0, 0)


@pytest.mark.parametrize("prime,multiplier,states,wanted", [
    (223, 3, [0] + [0] * 33 + [1] * 38 + [2, 4, 6], 2),
    (491, 7, [0] + [0] * 33 + [1] * 33 + [4, 6, 6, 6], 4),
])
def test_orbit_normalization_preserves_sums_and_permutes_every_periodic_correlation(
    prime, multiplier, states, wanted
):
    normalized, scale = normalize_orbit_scaling(states, prime, multiplier)
    _orbits, index, _shifts, _pairs, _coefficients = orbit_data(prime, multiplier)
    before = np.array([STATES[states[index[t]]] for t in range(prime)])
    after = np.array([STATES[normalized[index[t]]] for t in range(prime)])
    assert normalized[0] == 0 and normalized[1] == wanted
    assert np.array_equal(after, before[(scale * np.arange(prime)) % prime])
    assert np.array_equal(after.sum(axis=0), before.sum(axis=0))
    before_paf = [int(np.sum(before * np.roll(before, -s, axis=0))) for s in range(prime)]
    after_paf = [int(np.sum(after * np.roll(after, -s, axis=0))) for s in range(prime)]
    assert after_paf == [before_paf[(scale * s) % prime] for s in range(prime)]

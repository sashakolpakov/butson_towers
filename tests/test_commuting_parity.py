from itertools import product

import numpy as np

from experiments.scalar_lift.commuting_parity import (
    build_commuting_parity_report,
    commuting_factors,
    cube_switch,
    paley_core,
    reconstruct_flat_image,
)


def test_all_local_orientations_and_required_borders():
    report = build_commuting_parity_report()
    assert len(report["local_orientation_checks"]) == 16
    assert len(report["required_directions"]) == 4
    assert len(report["constant_column_controls"]) == 8
    assert not report["required_flat_column_constructed"]
    assert not report["full_frame_constructed"]


def test_cubic_criterion_matches_direct_rational_image_on_controls():
    b, d, e = commuting_factors()
    q = paley_core(179)
    probes = []
    for local in product((-1, 1), repeat=4):
        v = np.repeat(np.array(local)[:, None], 179, axis=1)
        probes.append(v)
        altered = v.copy()
        altered[0, 0] *= -1
        probes.append(altered)
    for v in probes:
        direct_numerator = 179 * v + b @ (v @ q.T) + (d @ v.sum(axis=1))[:, None]
        direct_is_flat = bool(np.all(np.abs(direct_numerator) == 358))
        reconstructed = reconstruct_flat_image(v, e, d, q)
        assert (reconstructed is not None) == direct_is_flat
        if reconstructed is not None:
            assert np.array_equal(358 * reconstructed, direct_numerator)


def test_parity_switch_is_an_involution_but_not_a_linear_cube_map():
    _, d, _ = commuting_factors()
    cube = np.array(list(product((-1, 1), repeat=4)), dtype=np.int64).T
    switched = cube_switch(cube, d)
    assert np.array_equal(cube_switch(switched, d), cube)
    assert np.any(cube.T @ cube != switched.T @ switched)

import numpy as np
import pytest

from experiments.scalar_lift.cubic_core_support_lift import phase_coordinates
from experiments.scalar_lift.fifth_source_assembly import exact_column_gram
from experiments.scalar_lift.variable_cubic_frame_bits import (
    adjacent_gauge_count,
    full_parameter_minor,
    orthogonal_parameter_problem,
)


@pytest.mark.parametrize("order,phase,patterns", [(1499, 24, 2), (2999, 48, 4)])
def test_actual_orthogonal_frames_preserve_physical_pairs_and_root_phases(order, phase, patterns):
    p = orthogonal_parameter_problem(order, phase)
    assert len(p["orthogonal_pattern_classes"]) == patterns
    assert p["boundary_pattern_Gram"] == (patterns*np.eye(patterns, dtype=np.int64)).tolist()
    assert len(p["actual_row_pairs"]) == 93*patterns
    for first, second in p["actual_row_pairs"]:
        assert p["row_parameters"][first] == p["row_parameters"][second]
    colors, dyadic = phase_coordinates(p["U_exponents"], phase)
    assert np.array_equal(colors, p["A"])
    assert np.array_equal(dyadic % (phase//6), p["A"]*p["row_parameters"][:, None])
    assert len(full_parameter_minor(p, phase//6)["normalized_star_minor"]) == patterns+1
    ratio = p["primitive_cross_ratio"]
    a, b = ratio["rows"]
    c, d = ratio["columns"]
    u = p["U_exponents"]
    assert (int(u[a,c])+int(u[b,d])-int(u[a,d])-int(u[b,c])) % phase == ratio["exponent"]
    # A single changed root really breaks the exact frame Gram.
    broken = u.copy()
    broken[0, 0] = (broken[0, 0]+phase//2) % phase
    assert not np.array_equal(exact_column_gram(broken, phase), exact_column_gram(u, phase))


def test_adjacent_channel_obstruction_is_distinct_from_a_binary_phase_control():
    star = np.array([[0, 0, 1], [0, 0, 3], [1, 3, 0]])
    assert adjacent_gauge_count(star) == {"normalized_gauges_checked": 128, "admissible_gauges": 0}
    assert adjacent_gauge_count(star % 2)["admissible_gauges"] > 0

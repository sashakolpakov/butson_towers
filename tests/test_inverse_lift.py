from __future__ import annotations

import json
from pathlib import Path

import pytest
from sympy import Matrix, eye, kronecker_product, symbols

from hadamard.inverse_lift import (
    is_sign_hadamard,
    lift_unreal_sixth_from_pair,
    sixth_root_block_pair_conditions,
)
from hadamard.scalar_lift import SKEW_CONFERENCE_4


def _pair():
    c = Matrix(SKEW_CONFERENCE_4)
    return eye(4) + c, -eye(4) + c


def test_inverse_pair_recovers_the_field_representation_and_dense_frame_images():
    p, q = _pair()
    assert sixth_root_block_pair_conditions(p.tolist(), q.tolist())["valid"]
    r = q * p.T / 4
    z = p - q
    assert r * r.T == eye(4)
    assert r**2 - r + eye(4) == Matrix.zeros(4)
    assert z * z.T == 4 * eye(4)
    assert r * z == p and r**2 * z == q
    w, v = (p - q) / 2, (p + q) / 2
    assert w * w.T == eye(4) and v * v.T == 3 * eye(4)
    assert w * v.T + v * w.T == Matrix.zeros(4)


def test_inverse_pair_composes_with_a_source_and_a_larger_block_degree():
    p, q = _pair()
    h2 = Matrix([[1, 1], [1, -1]])
    p, q = kronecker_product(p, h2), kronecker_product(q, h2)
    source = tuple(tuple(4 if i == j else 2 for j in range(3)) for i in range(3))
    h = Matrix(lift_unreal_sixth_from_pair(source, p.tolist(), q.tolist()))
    assert h.shape == (24, 24) and set(h) == {-1, 1}
    assert h * h.T == h.T * h == 24 * eye(24)


def test_two_hadamard_blocks_do_not_automatically_form_a_scalar_lift():
    p, _ = _pair()
    result = sixth_root_block_pair_conditions(p.tolist(), p.tolist())
    assert result["first_is_hadamard"] and result["second_is_hadamard"]
    assert not result["diagonal_angle_condition"] and not result["valid"]
    with pytest.raises(ValueError, match="angle identity"):
        lift_unreal_sixth_from_pair(((1,),), p.tolist(), p.tolist())


def test_dense_domain_and_source_orthogonality_are_checked():
    p, q = _pair()
    for invalid in (((None,),), ((0,),)):
        with pytest.raises(ValueError, match="dense"):
            lift_unreal_sixth_from_pair(invalid, p.tolist(), q.tolist())
    with pytest.raises(ValueError, match="not an exact"):
        lift_unreal_sixth_from_pair(((1, 1), (1, 1)), p.tolist(), q.tolist())


@pytest.mark.parametrize("bad", ((), ((1, 1),), ((0,),), ((True,),)))
def test_non_sign_blocks_are_rejected(bad):
    with pytest.raises(ValueError):
        is_sign_hadamard(bad)


def test_binary_gram_check_handles_odd_orders_without_rounding():
    assert is_sign_hadamard(((1,),))
    assert not is_sign_hadamard(((1, 1, 1), (1, -1, 1), (1, 1, -1)))


def test_order_twelve_control_is_a_symbolic_full_orthogonal_design():
    from experiments.scalar_lift.inverse_blocks import quaternion_order_twelve_array

    x, y, z, w = symbols("x y z w")
    design = Matrix(quaternion_order_twelve_array(x, y, z, w))
    assert set(design) == {x, -x, y, -y, z, -z, w, -w}
    gram = (design * design.T).applyfunc(lambda v: v.expand())
    assert gram == 3 * (x**2 + y**2 + z**2 + w**2) * eye(12)


def test_dense_sextic_operator_requires_the_odd_prime_denominator():
    from experiments.scalar_lift.inverse_blocks import (
        known_order_twelve_pair,
        pair_operator_arithmetic,
        prime_block_operator_requirements,
    )

    p, q = known_order_twelve_pair()
    assert sixth_root_block_pair_conditions(p, q)["valid"]
    r = Matrix(q) * Matrix(p).T / 12
    assert r**2 - r + eye(12) == Matrix.zeros(12)
    assert r * r.T == eye(12)
    arithmetic = pair_operator_arithmetic(p, q)
    assert arithmetic["zero_entries"] == 0
    assert all(
        counts == {"2": 3, "6": 9} for counts in arithmetic["reduced_denominator_counts_per_row"]
    )
    bound = prime_block_operator_requirements(668)
    assert bound["allowed_reduced_denominators"] == [2, 334]
    assert bound["minimum_denominator_twice_prime_entries_per_row"] == 665
    assert bound["half_integral_operator_excluded"]
    with pytest.raises(ValueError, match="odd prime"):
        prime_block_operator_requirements(36)
    with pytest.raises(ValueError, match="valid block pair"):
        pair_operator_arithmetic(p, p)


def test_tripling_preserves_dense_sixth_root_sources_of_different_orders():
    from experiments.scalar_lift.inverse_blocks import (
        known_order_thirteen_source,
        triple_unreal_source,
    )
    from hadamard.cgw import is_cgw_exponent_matrix

    source_three = tuple(tuple(4 if i == j else 2 for j in range(3)) for i in range(3))
    for source in (((1,),), source_three, known_order_thirteen_source()):
        tripled = triple_unreal_source(source)
        assert len(tripled) == 3 * len(source)
        assert {e for row in tripled for e in row} <= {1, 2, 4, 5}
        assert is_cgw_exponent_matrix(tripled, 6, len(tripled))
    with pytest.raises(ValueError, match="dense"):
        triple_unreal_source(((None,),))
    with pytest.raises(ValueError, match="exact Hadamard"):
        triple_unreal_source(((1, 1), (1, 1)))


def test_target_first_certificate_replays_sources_blocks_and_inverse_equations():
    from scripts.check_inverse_dense_lifts import build_report

    root = Path(__file__).resolve().parents[1]
    expected = json.loads((root / "certificates/scalar-lift/inverse-dense-lifts.json").read_text())
    actual = json.loads(json.dumps(build_report()))
    assert actual == expected
    assert actual["selected_target"]["source_matrix_constructed"]
    assert not actual["selected_target"]["required_block_pair_constructed"]
    assert actual["denominator_six_control_block_pair"]["arithmetic_bound_attained_in_every_row"]
    shared = next(p for p in actual["shared_required_block_pairs"] if p["block_order"] == 668)
    assert shared["candidate_target_orders"] == [4676, 8684]
    assert shared["source_orders"] == [7, 13]
    assert any(
        route["target_order"] == 6876 and route["source_order"] == 9 and route["block_order"] == 764
        for route in actual["available_source_routes_sorted_by_block_order"]
    )

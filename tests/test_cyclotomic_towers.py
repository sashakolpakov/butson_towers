from __future__ import annotations

import json
from math import prod
from pathlib import Path

import pytest
from sympy import I, Matrix, eye, im, kronecker_product, re, simplify

from experiments.scalar_lift.cyclotomic_towers import (
    cyclotomic_orders_with_degree_dividing,
    dyadic_flat_source_reduction,
    gaussian_terminal_degree_test,
    representation_obstructions,
    source_obstructions,
    target_factor_audit,
)
from hadamard.cgw import is_cgw_exponent_matrix
from hadamard.dense_lift import (
    TWELFTH_POWER_QUATERNARY_BLOCKS,
    lift_unreal_tenth_hadamard,
    lift_unreal_twelfth_hadamard,
    unreal_twelfth_to_quaternary,
)


def _clean(matrix):
    return matrix.applyfunc(lambda x: simplify(x.expand()))


def _realify(matrix):
    return Matrix.vstack(
        *(
            Matrix.hstack(*(Matrix([[re(x), -im(x)], [im(x), re(x)]]) for x in row))
            for row in matrix.tolist()
        )
    )


def test_iterated_twelfth_root_transfer_has_the_exact_composite_frame():
    s = Matrix([[I, 1 + I], [-1 + I, -I]])
    t = I * (eye(2) - s) / 2
    r = _realify(t)
    frame = kronecker_product(eye(2), Matrix([[0, 2], [2, 0]]))
    assert _clean(r * r.T) == eye(4)
    assert _clean(r**4 - r**2 + eye(4)) == Matrix.zeros(4)
    assert frame * frame.T == 4 * eye(4)
    for e, stored in TWELFTH_POWER_QUATERNARY_BLOCKS.items():
        block = Matrix([[I**x for x in row] for row in stored])
        assert block == _clean((1 + I) * t**e)
        h = Matrix(lift_unreal_twelfth_hadamard(((e,),)))
        assert _clean(r**e * frame) == h
        assert set(h) <= {-1, 1} and h * h.T == 4 * eye(4)


def test_sixth_root_source_really_passes_through_fourth_roots_before_signs():
    source = tuple(tuple(4 if i == j else 2 for j in range(3)) for i in range(3))
    middle = unreal_twelfth_to_quaternary(source, 6)
    assert len(middle) == 6 and is_cgw_exponent_matrix(middle, 4, 6)
    h = Matrix(lift_unreal_twelfth_hadamard(source, 6))
    assert h.shape == (12, 12) and set(h) == {-1, 1}
    assert h * h.T == 12 * eye(12)


def test_tower_partial_domain_and_source_identity_are_enforced():
    for e in (0, 3, 6, 9):
        with pytest.raises(ValueError, match="forbidden fourth root"):
            unreal_twelfth_to_quaternary(((e,),))
    for source in (((None,),), ((True,),)):
        with pytest.raises(ValueError, match="integer exponents"):
            unreal_twelfth_to_quaternary(source)
    with pytest.raises(ValueError, match="not an exact"):
        unreal_twelfth_to_quaternary(((1, 1), (1, 1)))
    with pytest.raises(ValueError, match="third, sixth or twelfth"):
        unreal_twelfth_to_quaternary(((1,),), 7)


def test_inverse_totient_enumeration_is_finite_and_includes_large_conductors():
    assert cyclotomic_orders_with_degree_dividing(2) == [1, 2, 3, 4, 6]
    assert cyclotomic_orders_with_degree_dividing(4) == [1, 2, 3, 4, 5, 6, 8, 10, 12]
    assert cyclotomic_orders_with_degree_dividing(668) == [1, 2, 3, 4, 5, 6, 8, 10, 12]
    assert 2154 in cyclotomic_orders_with_degree_dividing(716)


def test_standard_metric_filters_include_discriminant_and_dyadic_obstructions():
    assert any(
        o["kind"] == "hermitian_trace_determinant" for o in representation_obstructions(4, 10)
    )
    assert not representation_obstructions(8, 10)
    assert not representation_obstructions(4, 8)
    assert any(o["kind"] == "dyadic_split_conjugation" for o in representation_obstructions(12, 7))
    assert not representation_obstructions(24, 7)
    assert any(
        o["kind"] == "dyadic_split_conjugation" for o in representation_obstructions(716, 2154)
    )
    assert not representation_obstructions(36, 38)
    assert not source_obstructions(191, 38)
    assert gaussian_terminal_degree_test(36, 38)["passes_necessary_degree_test"]
    assert gaussian_terminal_degree_test(36, 76)["passes_necessary_degree_test"]
    assert not gaussian_terminal_degree_test(36, 114)["passes_necessary_degree_test"]
    assert gaussian_terminal_degree_test(36, 114)["compositum_absolute_degree"] == 72
    assert any(o["kind"] == "source_determinant_norm" for o in source_obstructions(191, 18))


def test_prime_factor_audit_excludes_the_higher_degree_branch_by_flat_alphabet():
    a = target_factor_audit(4676)
    assert [
        (b["source_order"], b["block_order"], b["surviving_phases"])
        for b in a["branches"]
        if b["surviving_phases"]
    ] == [
        (2338, 2, [4]),
        (7, 668, [6, 12]),
    ]
    a = target_factor_audit(6876)
    branch = next(b for b in a["branches"] if b["block_order"] == 36)
    assert branch["source_order"] == 191 and branch["surviving_phases"] == []
    for phase in (38, 76, 114):
        field = next(f for f in branch["coefficient_fields"] if f["phase"] == phase)
        assert any(o["kind"] == "dyadic_flat_source_alphabet" for o in field["obstructions"])
    assert a["necessary_only"] and a["source_or_lift_existence_not_asserted"]


def test_odd_denominator_stabilizer_is_exactly_the_quartic_subgroup_in_the_control():
    def integral_at_two(matrix):
        return all(int(simplify(x).q) % 2 for x in matrix)

    def common_row_product(matrix):
        products = [int(prod(row)) for row in matrix.tolist()]
        assert len(set(products)) == 1
        return products[0]

    s = Matrix([[I, 1 + I], [-1 + I, -I]])
    r = _realify(I * (eye(2) - s) / 2)
    assert [e for e in range(12) if integral_at_two(r**e)] == [0, 3, 6, 9]
    blocks = {e: Matrix(lift_unreal_twelfth_hadamard(((e,),))) for e in range(12) if e % 3}
    for e, a in blocks.items():
        for f, b in blocks.items():
            quotient = a * b.T / 4
            assert quotient == _clean(r ** ((e - f) % 12))
            assert integral_at_two(quotient) == (common_row_product(a) == common_row_product(b))
    primitive_third_root = r**4
    assert _clean(primitive_third_root**2 + primitive_third_root + eye(4)) == Matrix.zeros(4)
    assert _clean(primitive_third_root - eye(4)).det() == 9
    assert not integral_at_two(primitive_third_root)


def test_flat_alphabet_reduction_keeps_the_known_sources_and_respects_dimension_boundary():
    for d, n, k in ((4, 3, 6), (668, 7, 6), (668, 13, 12)):
        reduction = dyadic_flat_source_reduction(d, n, k)
        assert reduction["applicable"] and not reduction["excluded"]
    assert dyadic_flat_source_reduction(36, 191, 38)["excluded"]
    assert dyadic_flat_source_reduction(36, 191, 76)["excluded"]
    reduced114 = dyadic_flat_source_reduction(36, 191, 114)
    assert reduced114["available_reduced_alphabets"] == [6]
    assert reduced114["excluded"]
    assert not dyadic_flat_source_reduction(8, 5, 10)["applicable"]
    # Four flat fifth-root classes really do exist in block dimension 8,
    # so extending the dimension-4-mod-8 theorem to that case would be false.
    for e in (1, 2, 3, 4):
        h = Matrix(lift_unreal_tenth_hadamard(((e,),), 5))
        assert set(h) <= {-1, 1} and h * h.T == 8 * eye(8)


def test_saved_tower_report_replays_all_factor_filters_and_actual_dense_controls():
    from scripts.check_cyclotomic_towers import build_report

    root = Path(__file__).resolve().parents[1]
    expected = json.loads((root / "certificates/scalar-lift/cyclotomic-towers.json").read_text())
    report = json.loads(json.dumps(build_report()))
    assert report == expected
    assert report["candidate_orders_audited"] == 85
    assert all(c["real_signs_and_integer_gram_exact"] for c in report["known_two_stage_controls"])
    auxiliary = report["auxiliary_cyclotomic_subfield_degree_witness"]
    assert auxiliary["period_field_degree"] == 167
    assert auxiliary["compositum_with_Q_zeta12_degree"] == 668
    assert auxiliary["basis_trace_gram_exactly_verified"]
    assert not auxiliary["dense_sign_frame_constructed"]
    assert report["higher_degree_branch_audit"]["excluded_for_uniform_faithful_field_scalar_lifts"]
    assert not report["new_real_hadamard_order_claimed"]

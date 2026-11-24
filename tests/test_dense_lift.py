from __future__ import annotations

from itertools import product

import pytest
from sympy import I, Matrix, eye, simplify

from hadamard.cgw import is_cgw_exponent_matrix
from hadamard.dense_lift import (
    FIFTH_POWER_QUATERNARY_BLOCKS,
    QUATERNARY_REAL_BLOCKS,
    lift_quaternary_hadamard,
    lift_unreal_tenth_hadamard,
    unreal_tenth_to_quaternary,
)


def clean(matrix):
    return matrix.applyfunc(lambda x: simplify(x.expand()))


def complex_matrix(exponents):
    return Matrix([[I**e for e in row] for row in exponents])


def chirped_fourier_five():
    # 3(i+j)^2+1 never vanishes mod 5, since 3 is a nonsquare.
    return tuple(tuple((3 * (i + j) ** 2 + 1) % 5 for j in range(5)) for i in range(5))


def test_gaussian_scalar_frame_identity_is_exact_for_the_full_quaternary_alphabet():
    j = Matrix([[0, -1], [1, 0]])
    z = Matrix([[1, 1], [1, -1]])
    for e, stored in enumerate(QUATERNARY_REAL_BLOCKS):
        b = Matrix(stored)
        assert b == j**e * z
        assert set(b) <= {-1, 1}
        assert b * b.T == 2 * eye(2)


def test_all_four_fifth_root_powers_are_flat_and_represent_the_field_with_involution():
    m = complex_matrix(FIFTH_POWER_QUATERNARY_BLOCKS[0]) / 2
    assert clean(m * m.conjugate().T) == eye(4)
    assert clean(sum((m**e for e in range(5)), Matrix.zeros(4))) == Matrix.zeros(4)
    assert clean(m**5) == eye(4)
    for e, stored in enumerate(FIFTH_POWER_QUATERNARY_BLOCKS, 1):
        block = complex_matrix(stored)
        assert block == clean(2 * m**e)
        assert set(block) <= {1, -1, I, -I}
        assert clean(block * block.conjugate().T) == 4 * eye(4)


@pytest.mark.parametrize("phase", (5, 10))
def test_a_nonreal_fourier_input_produces_a_dense_real_hadamard_of_order_forty(phase):
    source = tuple(tuple(phase // 5 * e for e in row) for row in chirped_fourier_five())
    assert is_cgw_exponent_matrix(source, phase, 5)
    q = unreal_tenth_to_quaternary(source, phase)
    assert is_cgw_exponent_matrix(q, 4, 20)
    h = Matrix(lift_unreal_tenth_hadamard(source, phase))
    assert h.shape == (40, 40)
    assert set(h) == {-1, 1}
    assert h * h.T == h.T * h == 40 * eye(40)
    assert all(
        sum(h[i, k] == h[j, k] for k in range(40)) == 20 for i in range(40) for j in range(i)
    )


def test_genuinely_signed_tenth_root_input_also_has_a_dense_lift():
    f = chirped_fourier_five()
    source = tuple(
        tuple((2 * f[i // 2][j // 2] + 5 * (i % 2) * (j % 2)) % 10 for j in range(10))
        for i in range(10)
    )
    assert any(e % 2 for row in source for e in row)
    h = Matrix(lift_unreal_tenth_hadamard(source))
    assert h.shape == (80, 80)
    assert set(h) == {-1, 1}
    assert h * h.T == h.T * h == 80 * eye(80)


def test_every_allowed_tenth_root_has_the_correct_dense_eight_by_eight_block():
    m = complex_matrix(FIFTH_POWER_QUATERNARY_BLOCKS[0]) / 2
    for e in (1, 2, 3, 4, 6, 7, 8, 9):
        q = unreal_tenth_to_quaternary(((e,),))
        assert complex_matrix(q) == clean(2 * (-1) ** e * m ** (3 * e % 5))
        h = Matrix(lift_unreal_tenth_hadamard(((e,),)))
        assert set(h) == {-1, 1}
        assert h * h.T == 8 * eye(8)


@pytest.mark.parametrize("bad", ((), ((0, 1),), ((True,),), ((0.5,),), ((None,),)))
def test_dense_transfers_reject_malformed_inputs_and_zeros(bad):
    with pytest.raises(ValueError):
        lift_unreal_tenth_hadamard(bad)
    with pytest.raises(ValueError):
        lift_quaternary_hadamard(bad)


def test_partial_domain_and_exact_source_gram_are_enforced():
    for e in (0, 5):
        with pytest.raises(ValueError, match="forbidden real"):
            lift_unreal_tenth_hadamard(((e,),))
    with pytest.raises(ValueError, match="not an exact"):
        lift_unreal_tenth_hadamard(((1, 1), (1, 1)))
    with pytest.raises(ValueError, match="not an exact"):
        lift_quaternary_hadamard(((0, 0), (0, 0)))
    with pytest.raises(ValueError, match="fifth or tenth"):
        lift_unreal_tenth_hadamard(((1,),), 6)


def test_multiplicative_projection_has_the_claimed_boundary_counterexamples():
    for p in (3, 5):
        n = 4 * p
        source = tuple(
            tuple(((i // 2) * (j // 2) + p * (i % 2) * (j % 2)) % (2 * p) for j in range(n))
            for i in range(n)
        )
        assert is_cgw_exponent_matrix(source, 2 * p, n)
        s = Matrix([[(-1) ** e for e in row] for row in source])
        assert s.rank() == 4
        assert s * s.T != n * eye(n)
        assert all(int((s * s.T)[i, j]) % p == 0 for i, j in product(range(n), repeat=2))


def test_f9_norm_character_gives_a_dense_hadamard_without_a_completion():
    from scripts.check_dense_lifts import norm_f9_hadamard

    h = norm_f9_hadamard()
    assert h.shape == (20, 20)
    assert set(h) == {-1, 1}
    assert h * h.T == h.T * h == 20 * eye(20)


def test_saved_dense_lift_report_replays_the_field_and_integer_sign_identities():
    import json
    from pathlib import Path

    from scripts.check_dense_lifts import build_report

    root = Path(__file__).resolve().parents[1]
    expected = json.loads((root / "certificates/scalar-lift/dense-lifts.json").read_text())
    assert json.loads(json.dumps(build_report())) == expected

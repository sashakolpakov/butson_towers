import json
from pathlib import Path

import pytest

from experiments.scalar_lift.h668_width12 import decode_signs, known_h668_width12


def test_printed_width12_seed_reconstructs_exact_hadamard668():
    _matrix, report = known_h668_width12()
    root = Path(__file__).resolve().parents[1]
    assert report == json.loads((root / "certificates/scalar-lift/h668-width12.json").read_text())
    assert report["core_row_sums"] == [14, -18, -2, -12]
    assert report["full_integer_and_independent_bit_distance_Gram_checks_passed"]
    assert not report["inequivalence_independently_verified"]
    assert not report["scalar_lift_pair_or_new_order_constructed"]


def test_hex_decoding_preserves_leading_zeros_and_rejects_nonzero_padding():
    assert decode_signs("04", 6).tolist() == [1, 1, 1, 1, 1, -1]
    with pytest.raises(ValueError, match="zero padding"):
        decode_signs("05", 6)

"""Independent dense controls for the additional scalar source convention."""
import copy
import json

import numpy as np
import pytest

from experiments.scalar_lift import near_williamson_scalar_endpoints as source
from experiments.scalar_lift.williamson_scalar_endpoints import (
    ScalarEndpoint, verify_orthogonal_design,
)


@pytest.mark.parametrize("width", [47, 53, 69, 73, 127])
def test_complete_small_endpoints_and_literal_interfaces(width):
    # The order-one T-input makes a small complete endpoint. Multiplying
    # every dense entry independently catches reversal and block mistakes.
    x, _ = verify_orthogonal_design(np.array([[1], [0], [0], [0]]))
    a, _ = source.verified_source(width)
    h = sum(np.kron(xi.astype(np.int64), ai.astype(np.int64))
            for xi, ai in zip(x, a, strict=True))
    assert np.all(abs(h) == 1)
    assert np.array_equal(h @ h.T, len(h)*np.eye(len(h), dtype=np.int64))
    literal = ScalarEndpoint(x, a)
    for i in (0, width-1, width, len(h)-1):
        assert np.array_equal(literal.row(i), h[i])
        for j in (0, width-1, width, len(h)-1):
            assert literal.entry(i, j) == h[i, j]


def test_unreversed_near_williamson_input_is_rejected():
    records = json.loads(source.INPUTS.read_text())["entries"]
    record = copy.deepcopy(next(r for r in records if r["width"] == 53))
    record["first_matrix_reversed"] = False
    with pytest.raises(ValueError, match="symmetric"):
        source.verify_quadruple(record)


def test_changed_literal_is_rejected_before_endpoint_construction():
    records = json.loads(source.INPUTS.read_text())["entries"]
    record = copy.deepcopy(next(r for r in records if r["width"] == 69))
    word = record["first_rows_pm"][0]
    record["first_rows_pm"][0] = ("-" if word[0] == "+" else "+") + word[1:]
    with pytest.raises(ValueError):
        source.verify_quadruple(record)

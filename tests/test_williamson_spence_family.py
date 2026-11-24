import json
from pathlib import Path

import numpy as np
import pytest

from experiments.scalar_lift import williamson_spence_family as family


def test_complete_source_family_certificate_replays_exactly():
    path = Path(__file__).resolve().parents[1]/'certificates/scalar-lift/williamson-spence-family.json'
    report = family.certificate()
    assert report == json.loads(path.read_text())
    assert [x['source_order'] for x in report['first_expanded_members']] == [91,153,1891,5151,7381]
    assert [x['literal_T631_endpoint_order_and_Gram_scale']
            for x in report['first_expanded_members']] == [229684,386172,4772884,13001124,18629644]


def test_deeper_implicit_recurrence_and_exact_convolution_reject_corruption():
    source = family.SpenceSource(13,2)
    rows = source.coefficient_rows()
    family.exact_group_norms(rows, source.group_dimensions, 13, 2)
    assert source.order == 1183
    assert source.entry(2,source.order-1,100) == source.row(2,source.order-1)[100]
    bad = rows.copy()
    bad[0,0] *= -1
    with pytest.raises(ValueError,match='group-ring norm'):
        family.exact_group_norms(bad, source.group_dimensions, 13, 2)


def test_source_domain_and_prime_power_additive_group():
    source = family.SpenceSource(121,1)
    assert source.group_dimensions == (61,11,11)
    assert family.field_subtract(121,0,1) == 10
    assert family.field_subtract(121,0,11) == 110
    assert family.character(121,11) in (-1,1)
    for q,r in ((True,0),(13,False),(13,-1),(73,0)):
        with pytest.raises(ValueError):family.SpenceSource(q,r)
    for channel,i,j in ((True,0,0),(4,0,0),(0,-1,0),(0,0,source.order)):
        with pytest.raises(ValueError):source.entry(channel,i,j)
    assert isinstance(source.entry(np.int64(0),np.int64(0),np.int64(0)),int)


def test_implicit_source_substitution_has_full_exact_H2548_Gram():
    from experiments.scalar_lift.williamson_scalar_endpoints import (
        ImplicitScalarEndpoint, verify_orthogonal_design,
    )
    t = np.zeros((4,7),dtype=np.int64)
    t[0,0] = 1
    t[0,[1,2,4]] = -1
    t[1,3] = t[2,6] = t[3,5] = 1
    x,_ = verify_orthogonal_design(t)
    endpoint = ImplicitScalarEndpoint(x,family.SpenceSource(13,1))
    assert endpoint.order == 2548
    packed = []
    for i in range(endpoint.order):
        row = endpoint.row(i)
        assert np.all(np.abs(row) == 1)
        for j in (0,90,91,endpoint.order-1):
            assert endpoint.entry(i,j) == row[j]
        packed.append(int.from_bytes(np.packbits(row<0,bitorder='little').tobytes(),'little'))
    for i,row in enumerate(packed):
        assert all((row^other).bit_count() == endpoint.order//2 for other in packed[:i])

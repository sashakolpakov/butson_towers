import numpy as np
import pytest

from experiments.scalar_lift import hadamard_frame_search as search
from experiments.scalar_lift import williamson_scalar_endpoints as transfer


@pytest.fixture(scope="module")
def small_design():
    b = search.expand([0,1,2],search.chart(7,3,2))
    return transfer.verify_orthogonal_design(b)


def test_polynomial_design_matches_all_dense_coefficient_grams(small_design):
    x,proof = small_design
    x = x.astype(np.int64)
    for a in range(4):
        assert np.array_equal(x[a]@x[a].T,7*np.eye(28,dtype=np.int64))
        for b in range(a):
            assert not np.any(x[a]@x[b].T+x[b]@x[a].T)
    assert proof["all_polynomial_periodic_Gram_coefficients_checked"] == 70


@pytest.mark.parametrize("w",[5,9,17])
def test_scalar_recipe_has_every_dense_endpoint_row_product(small_design,w):
    x,_ = small_design
    a,_ = transfer.verified_williamson(w)
    endpoint = transfer.ScalarEndpoint(x,a)
    actual = np.array([endpoint.row(i) for i in range(endpoint.order)],dtype=np.int64)
    assert set(np.unique(actual)) == {-1,1}
    assert np.array_equal(actual@actual.T,endpoint.order*np.eye(endpoint.order,dtype=np.int64))
    for i in [0,w-1,w,endpoint.order-1]:
        for j in [0,w,w+1,endpoint.order-1]:
            assert endpoint.entry(i,j) == actual[i,j]


def test_corrupted_T_source_is_rejected():
    b = search.expand([0,0,2],search.chart(7,3,2))
    with pytest.raises(ValueError,match="periodic T Gram"):
        transfer.verify_orthogonal_design(b)


def test_corrupted_Williamson_source_is_rejected():
    a,proof = transfer.verified_williamson(17)
    rows = a[:,0,:].copy()
    rows[0,0] *= -1
    with pytest.raises(ValueError,match="summed Gram"):
        transfer.verify_williamson(rows)


def test_nonintegral_Williamson_signs_are_rejected():
    a,_ = transfer.verified_williamson(17)
    with pytest.raises(ValueError,match="sign words"):
        transfer.verify_williamson(1.9*a[:,0,:])


def test_cyclotomic_algebra_dimension_is_not_a_single_field_degree():
    algebra = transfer.algebra_description(15)
    assert algebra["rational_dimension"] == 15
    assert algebra["number_field_factors"] == [
        {"conductor":1,"field_degree":1,"module_rank":1},
        {"conductor":3,"field_degree":2,"module_rank":1},
        {"conductor":5,"field_degree":4,"module_rank":1},
        {"conductor":15,"field_degree":8,"module_rank":1},
    ]
    assert not algebra["algebra_is_single_number_field"]


def test_noncyclic_group_uses_field_multiplicities():
    algebra = transfer.group_algebra_description([61,11,11])
    assert algebra["rational_dimension"] == 7381
    assert algebra["group_exponent"] == 671
    assert algebra["cyclotomic_field_factors"] == [
        {"conductor":1,"field_degree":1,"multiplicity":1},
        {"conductor":11,"field_degree":10,"multiplicity":12},
        {"conductor":61,"field_degree":60,"multiplicity":1},
        {"conductor":671,"field_degree":600,"multiplicity":12},
    ]

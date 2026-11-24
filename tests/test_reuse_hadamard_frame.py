import numpy as np
import pytest

from experiments.scalar_lift import hadamard_frame_search as search
from experiments.scalar_lift import reuse_hadamard_frame as reuse


def test_completed_generic_reduction_has_full_small_endpoint_Grams(tmp_path):
    record = {'p':7, 'm':3, 'mu':2, 'states':[0,1,2]}
    report = reuse.certificate(record, tmp_path/'verified')
    assert report['output_orders'] == [84,196,252,364]
    assert report['field_degree'] == 2
    assert report['module_rank'] == 14
    assert report['rational_dimension'] == 28
    data = search.chart(7,3,2)
    b = search.expand(record['states'],data)
    a = search.D4 @ b
    P,Q = search.gs(a-b),search.gs(a+b)
    for result in report['outputs']:
        literal = reuse.ReusedEndpoint(result['source_id'],
            np.array(result['source_exponents']),P,Q)
        matrix = np.array([literal.row(i) for i in range(literal.order)],dtype=np.int64)
        assert set(np.unique(matrix)) == {-1,1}
        assert np.array_equal(matrix @ matrix.T,literal.order*np.eye(literal.order,dtype=np.int64))
        assert result['source_gram_scale'] * result['frame_gram_scale'] == literal.order
        assert result['complete_algebraic_endpoint_Gram_verified']
    assert not report['new_orders_claimed']


def test_generic_frame_replay_rejects_a_corrupted_witness(tmp_path):
    with pytest.raises(ValueError,match='coefficient Gram'):
        reuse.replay({'p':7,'m':3,'mu':2,'states':[0,0,2]},tmp_path/'invalid')
    assert not (tmp_path/'invalid').exists()

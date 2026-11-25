"""Check projection optima against finite enumeration and exact T controls."""
from itertools import product
import shutil

import numpy as np
import pytest

from experiments.scalar_lift import hadamard_frame_search as frame
from experiments.scalar_lift.t_sequence_projection import (
    Projections, build_assignment, search, verify,
)


@pytest.fixture(scope='module')
def assignment(tmp_path_factory):
    if not (shutil.which('c++') or shutil.which('clang++')):
        pytest.skip('native numeric assignment requires a C++17 compiler')
    function,library=build_assignment(tmp_path_factory.mktemp('t-projection'))
    yield function
    del library


@pytest.mark.parametrize('heads',[3,None])
def test_discrete_projection_matches_every_admissible_T19_function(assignment,heads):
    data=frame.chart(19,3,7)
    engine=Projections(data,heads,1,1,assignment)
    rng=np.random.default_rng(101)
    # Enumerate all actual signed-axis orbit words independently. This checks
    # the orbit-coordinate convention as well as fixed/free count optimization.
    admissible=[]
    for word in product(range(8),repeat=6):
        if heads is not None and sum(s//2==0 for s in word)!=heads:continue
        if sum((1-2*(s%2)) for s in word if s//2==0)!=1:continue
        if sum((1-2*(s%2)) for s in word if s//2!=0)!=1:continue
        states=[0]+list(word)
        admissible.append(frame.expand(states,data))
    candidates=np.stack(admissible)
    for _ in range(3):
        values=rng.normal(size=(4,19))
        states,discrete=engine.discrete(values)
        brute=np.einsum('nij,ij->n',candidates,values).max()
        assert float(np.sum(discrete*values))==pytest.approx(brute,abs=1e-9)
        assert np.array_equal(discrete.sum(axis=1),engine.target)
        assert states[0]==0


@pytest.mark.parametrize('heads',[5,19,71,None])
def test_fourier_projection_has_correct_norm_DC_and_covariance(assignment,heads):
    data=frame.chart(223,3,183)
    engine=Projections(data,heads,-5,3,assignment)
    states,b=engine.discrete(np.random.default_rng(102).normal(size=(4,223)))
    transformed=np.fft.fft(b,axis=1)
    power=np.sum(np.abs(transformed)**2,axis=0)
    report=verify(states,data,engine.target)
    assert float(np.sum((power-223)**2)/(6*223))==pytest.approx(report['best_energy'],abs=1e-9)
    projected=engine.spectral(transformed)
    assert np.allclose(np.sum(np.abs(projected)**2,axis=0),223,atol=1e-10)
    assert np.array_equal(projected[:,0],engine.target)
    values=np.fft.ifft(projected,axis=1)
    assert np.max(np.abs(values.imag))<1e-12
    assert np.allclose(values[:,183*np.arange(223)%223],values[[0,3,1,2]],atol=1e-12)


@pytest.mark.parametrize('orientation',['discrete-first','spectral-first'])
@pytest.mark.parametrize('heads',[1,None])
def test_DR_exact_integer_control_is_checked_at_every_shift(assignment,tmp_path,orientation,heads):
    data=frame.chart(31,3,5)
    # The seeded solution makes this an exact endpoint and fixed-point test;
    # it does not depend on a stochastic search reaching zero within a timeout.
    planted=[0,1,2,2,7,6,2,6,2,5,3]
    result=search(data,heads,-1,3,assignment,2,20261125,.8,
                  tmp_path/'control.json',initial=planted,orientation=orientation)
    assert result['exact_T_design']
    assert result['best_energy']==0
    assert result['all_expanded_shift_correlations_checked']
    assert result['row_sums']==[-2,3,3,3]
    assert result['maximum_absolute_nonzero_shift']==0
    damaged=planted.copy();damaged[2]=4
    bad=verify(damaged,data,np.array([-2,3,3,3]))
    assert bad['best_energy']>0
    assert not bad['exact_T_design']

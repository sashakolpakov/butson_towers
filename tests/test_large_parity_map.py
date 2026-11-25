"""Exact mathematical controls for the larger determinant-coordinate map."""
import random

from experiments.scalar_lift import hadamard_frame_search as old
from experiments.scalar_lift.large_parity_map import (
    BinaryField, DeterminantSampler, normalize_head, recover_map, validate_map,
    validate_expanded_geometry,
)


def correlation_parities(axes, chart):
    return [sum(count for (i,j,delta),count in row.items()
                if axes[i] == (0 if axes[j] == 0 else 1+(axes[j]-1+delta)%3)) % 2
            for row in chart['rows']]


def test_field_inverse_half_trace_and_irreducibility():
    field=BinaryField(3,0b1011)
    assert field.irreducible()
    assert not BinaryField(3,0b1001).irreducible()
    for x in range(1,8):
        assert field.mul(x,field.inverse(x))==1
    for c in range(8):
        roots=[x for x in range(8) if field.square(x)^x==c]
        y=field.artin_schreier(c)
        assert (y is None)==(not roots)
        if roots:assert y in roots


def test_complete_order19_determinant_census():
    config=recover_map(19,3,7,1,1)
    assert validate_map(config)['coefficients_checked']==79
    assert validate_expanded_geometry(config)['exact_correlations_compared']==79*18
    field=BinaryField(3,config['field_modulus'])
    chart=old.chart(19,3,7)
    admitted=set()
    for p in range(8):
        for q in range(8):
            if not(p or q):continue
            for w in range(8):
                if p:b0,b1=w,field.mul(1^field.mul(q,w),field.inverse(p))
                else:b0,b1=field.inverse(q),w
                packed=p|(q<<3)|(b0<<6)|(b1<<9)
                color=old.binary_linear(config['inverse_columns'],packed)
                axes=[0]+[(color>>(2*j))&3 for j in range(6)]
                assert not any(correlation_parities(axes,chart))
                admitted.add(color)
    actual=set()
    for color in range(1<<12):
        axes=[0]+[(color>>(2*j))&3 for j in range(6)]
        if not any(correlation_parities(axes,chart)):actual.add(color)
    assert admitted==actual
    assert len(actual)==8*(8*8-1)==504


def test_normalization_is_actual_multiplicative_reindexing():
    config=recover_map(19,3,7,1,1)
    chart=old.chart(19,3,7)
    sampler=DeterminantSampler(config,rng=random.Random(901))
    def expand(axes):
        return [0 if axes[int(chart['index'][x])]==0 else
                1+(axes[int(chart['index'][x])]-1+int(chart['phase'][x]))%3
                for x in range(19)]
    checked=0
    for _ in range(100):
        axes=sampler.next_axes()
        if 0 not in axes[1:]:continue
        normalized,scale=normalize_head(axes,chart)
        before,after=expand(axes),expand(normalized)
        assert after==[before[scale*x%19] for x in range(19)]
        assert normalized[1]==0
        assert not any(correlation_parities(normalized,chart))
        checked+=1
    assert checked>10


def test_full_large_map_and_existing_field_controls():
    for p,m,mu,head,tail,expected in (
        (127,9,22,1,1,407),
        (631,21,458,1,1,1831),
        (223,3,183,-5,3,11027),
    ):
        config=recover_map(p,m,mu,head,tail)
        audit=validate_map(config)
        assert audit['all_quadratic_coefficients_verified']
        assert audit['coefficients_checked']==expected
        chart=old.chart(p,m,mu)
        sampler=DeterminantSampler(config,seed=901)
        for _ in range(12):
            axes=sampler.next_axes()
            assert not any(correlation_parities(axes,chart))

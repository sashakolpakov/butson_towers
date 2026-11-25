#!/usr/bin/env python3
"""Check six further scalar endpoints from the H7572 frame.

Place this beside recipe.py, t631.json and frame-reuse.json (as in anc/ of
the source archive). The larger matrices are specified entry by entry;
their complete Grams follow algebraically from the checked source/frame
identities. This script does not allocate the larger endpoint matrices.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np


def require(condition, message):
    if not condition:
        raise ValueError(message)


def source_gram(exponents):
    roots = np.array(((1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)),dtype=np.int64)
    E = np.asarray(exponents)
    require(E.ndim==2 and E.shape[0]==E.shape[1] and np.issubdtype(E.dtype,np.integer),
            'source must be a square integer exponent matrix')
    require(set(map(int,np.unique(E))) <= {0,1,3,4}, 'unadmitted scalar symbol')
    n = len(E)
    for matrix in (E,E.T):
        for i in range(n):
            for j in range(n):
                actual = roots[(matrix[i]-matrix[j])%6].sum(axis=0)
                require(actual.tolist()==([n,0] if i==j else [0,0]), 'complete source Gram failed')


def triple_unreal(gauged_exponents):
    # Undo the global zeta_6^-1 gauge, perform the known unreal tripling,
    # then restore that gauge. All arithmetic remains exact integer data.
    U = (np.asarray(gauged_exponents,dtype=np.int64)+1)%6
    require(set(map(int,np.unique(U))) <= {1,2,4,5},'source is not unreal')
    A = np.where(np.isin(U,(1,5)),1,-1)
    B = np.where(np.isin(U,(1,2)),1,-1)
    C = np.array(((0,1,-1),(-1,0,1),(1,-1,0)),dtype=np.int64)
    Ap = np.kron(B,np.ones((3,3),dtype=np.int64))
    Bp = np.kron(A,np.eye(3,dtype=np.int64))+np.kron(B,C)
    require(np.all(Ap*Ap==1) and np.all(Bp*Bp==1),'tripled coefficients are not signs')
    exponents = np.array(((4,2),(5,1)),dtype=np.int64)[(Ap>0).astype(int),(Bp>0).astype(int)]
    return (exponents-1)%6


def verify(directory=None):
    here = Path(__file__).resolve().parent if directory is None else Path(directory)
    spec = importlib.util.spec_from_file_location('_h7572_core',here/'recipe.py')
    core = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(core)
    a,b = core.expand_sequences()
    P,Q = core.goethals_seidel(a-b),core.goethals_seidel(a+b)
    d = len(P)
    pb,pr = core.pack_signs(P)
    qb,qr = core.pack_signs(Q)
    require(core.sha256(pb.tobytes())==core.P_SHA256,'P differs from certified frame')
    require(core.sha256(qb.tobytes())==core.Q_SHA256,'Q differs from certified frame')
    core.verify_gram(pr,d)
    core.verify_gram(qr,d)
    for i in range(d):
        require((pr[i]^qr[i]).bit_count()==d//4,'cross diagonal failed')
        for j in range(i):
            require((pr[i]^qr[j]).bit_count()+(qr[i]^pr[j]).bit_count()==d,
                    'complete symmetric cross identity failed')
    certificate = json.loads((here/'frame-reuse.json').read_text())
    sources = list(certificate['outputs'])
    E = np.array(next(item['exponents'] for item in sources if item['source_order']==13),dtype=np.int64)
    for _ in range(3):
        E = triple_unreal(E)
        n = len(E)
        sources.append({'exponents':E.tolist(),'source_order':n,'source_gram_scale':n,
                        'frame_width':d,'frame_gram_scale':d,
                        'output_order':d*n,'output_gram_scale':d*n})
    outputs = []
    for item in sources:
        E = np.array(item['exponents'],dtype=np.int64)
        source_gram(E)
        n = len(E)
        require(item['source_order']==n and item['source_gram_scale']==n,'wrong source scale')
        require(item['frame_width']==d and item['frame_gram_scale']==d,'wrong frame scale')
        require(item['output_order']==d*n and item['output_gram_scale']==d*n,'wrong endpoint scale')
        # Literal row/block controls supplement the full algebraic proof.
        selected = sorted({0,d-1,d,d+1,d*n-1})
        rows = []
        for index in selected:
            i,a = divmod(index,d)
            rows.append(np.concatenate([(1 if e<3 else -1)*(P if e%3==0 else Q)[a]
                                        for e in E[i]]))
        rows = np.array(rows,dtype=np.int64)
        require(np.array_equal(rows@rows.T,d*n*np.eye(len(rows),dtype=np.int64)),
                'literal endpoint row controls failed')
        outputs.append({'order':d*n,'source_order':n,'source_scale':n,'frame_scale':d,
                        'all_entry_images_are_signs':True,
                        'complete_algebraic_Gram_verified':True,
                        'full_dense_endpoint_Gram_performed':False,
                        'order_status':('known' if d*n==885924 else 'new'),
                        'order_status_basis':('Paley I with prime 885923'
                                              if d*n==885924 else
                                              "New to the best of the author's knowledge; no earlier construction located in the dated literature search.")})
    return {'schema':'portable-h7572-further-endpoints-audit-v1',
            'frame_width':d,'field':'Q(zeta_6)/Q','field_degree':2,'module_rank':d//2,
            'conjugation':'transpose','canonical_trace_metric':[[2,1],[1,2]],
            'full_P_Q_and_symmetric_cross_integer_identities_verified':True,
            'frame_reuse_certificate_sha256':hashlib.sha256((here/'frame-reuse.json').read_bytes()).hexdigest(),
            'outputs':outputs,'new_orders_claimed':True,
            'new_orders':[item['order'] for item in outputs if item['order_status']=='new'],
            'novelty_convention':"Known where an earlier construction was located; otherwise new to the best of the author's knowledge."}


if __name__=='__main__':
    print(json.dumps(verify(),indent=2))

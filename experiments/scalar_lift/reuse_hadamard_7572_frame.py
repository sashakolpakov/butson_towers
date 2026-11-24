"""Exact completed scalar endpoint recipes reusing the certified H7572 frame.

The resulting real orders 17668, 22716 and 32812 are new to the best of the
author's knowledge. Full algebraic Grams follow from complete exact source Grams,
the independently replayed full frame identities and all used entry tests.
No large endpoint matrix is allocated or claimed checked by dense arithmetic.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from functools import cache
import hashlib
import importlib.util
import json
import operator
from pathlib import Path

import numpy as np

from experiments.scalar_lift.forward_butson_descent import exponent_hash
from experiments.scalar_lift.forward_source_registry import load_source as registry_source
from experiments.scalar_lift.catalogue_source_acquisition import load_source as catalogue_source, data as catalogue_data
from experiments.scalar_lift.broad_source_acquisition import load_source as broad_source, dataset

ROOT = Path(__file__).resolve().parents[2]
BUNDLE = ROOT / 'artifacts/hadamard-7572'
WIDTH = 2524
SOURCES = ('catalogue-p7-q6', 'aalto-9-6-000', 'unreal13-q6')
ROOTS = np.array(((1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)), dtype=np.int64)
NOVELTY_CONVENTION = {
    'known': 'an identified prior construction gives this real Hadamard order',
    'new_to_best_of_authors_knowledge': 'a completed endpoint for which the author knows no prior construction',
    'new_order_claims_are_qualified_by_authors_knowledge': True,
    'absolute_priority_claimed': False,
    'equivalence_class_novelty_claimed': False,
}


def _digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _json_digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def verify_source_gram(exponents):
    values = np.asarray(exponents)
    if values.ndim != 2 or values.shape[0] != values.shape[1] or not np.issubdtype(values.dtype, np.integer):
        raise ValueError('a complete square integer exponent matrix is required')
    if np.any(values < 0) or np.any(values >= 6):
        raise ValueError('source exponents must lie in 0..5')
    expected = np.zeros((*values.shape,2), dtype=np.int64)
    expected[np.arange(len(values)), np.arange(len(values)), 0] = len(values)
    for matrix in (values, values.T):
        a,b = ROOTS[matrix][...,0], ROOTS[matrix][...,1]
        # z^2=z-1 and conjugate(z)=1-z.
        gram = np.stack((a@(a+b).T+b@b.T, b@a.T-a@b.T), axis=2)
        if not np.array_equal(gram, expected):
            raise ValueError('complete exact sixth-root source Gram failed')
    return {'order':len(values), 'phase':6, 'row_and_column_Grams_verified':True,
            'source_gram_scale':len(values), 'exponent_sha256':exponent_hash(values)}


@cache
def verified_frames():
    path = BUNDLE / 'recipe.py'
    spec = importlib.util.spec_from_file_location('_verified_h7572_recipe_for_reuse', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    a,b = module.expand_sequences()
    P,Q = module.goethals_seidel(a-b), module.goethals_seidel(a+b)
    packed_p, rows_p = module.pack_signs(P)
    packed_q, rows_q = module.pack_signs(Q)
    hp,hq = module.sha256(packed_p.tobytes()),module.sha256(packed_q.tobytes())
    if hp != module.P_SHA256 or hq != module.Q_SHA256:
        raise ValueError('frame differs from the immutable completed construction')
    module.verify_gram(rows_p, WIDTH)
    module.verify_gram(rows_q, WIDTH)
    for i in range(WIDTH):
        if (rows_p[i]^rows_q[i]).bit_count() != WIDTH//4:
            raise ValueError('frame cross diagonal failed')
        for j in range(i):
            if (rows_p[i]^rows_q[j]).bit_count()+(rows_q[i]^rows_p[j]).bit_count() != WIDTH:
                raise ValueError('frame symmetric cross identity failed')
    for matrix in (P,Q):
        matrix.setflags(write=False)
    proof = {'rational_dimension':WIDTH, 'field':'Q(zeta_6)/Q', 'field_degree':2,
        'module_rank':WIDTH//2, 'conjugation':'zeta_6 -> 1-zeta_6; scalar-action transpose',
        'canonical_trace_basis':['1','zeta_6'], 'canonical_trace_metric':[[2,1],[1,2]],
        'actual_Euclidean_coordinate_trace_form':'h(u,v)=u^T*v/2-(u^T*S*v)*sqrt(-3)/6, S=2R-I; Tr(h)=u^T*v',
        'scalar_action':'R=Q*P^T/2524', 'scalar_action_identities':['R*R^T=I','R+R^T=I','R^2-R+I=0'],
        'frame':'F=P', 'frame_gram_scale':WIDTH,
        'full_P_integer_Gram_verified':True, 'full_Q_integer_Gram_verified':True,
        'full_symmetric_cross_identity_verified':True,
        'unordered_row_pairs_checked_per_identity':WIDTH*(WIDTH-1)//2,
        'P_packed_sha256':hp, 'Q_packed_sha256':hq,
        'actual_frame_and_used_image_alphabet':[-1,1],
        'flat_entry_images':{'0':'P','1':'Q','3':'-P','4':'-Q'},
        'unused_nonflat_images':{'2':'Q-P','5':'P-Q'},
        'immutable_recipe':'artifacts/hadamard-7572/recipe.py',
        'immutable_recipe_sha256':_digest(path),
        'immutable_sequences_sha256':_digest(BUNDLE/'t631.json')}
    return P,Q,proof


def verified_input(identifier):
    if identifier == 'catalogue-p7-q6':
        source = catalogue_source(identifier)
        raw = np.asarray(catalogue_data('P7')['exponents'], dtype=np.int64)
        gauge = source.proof['acquisition']['fixed_starting_gauge']
        rows,columns,shift = gauge['row_shifts'],gauge['column_shifts'],0
        raw_path = ROOT/'experiments/scalar_lift/catalogue_inputs/P7.json'
    elif identifier == 'unreal13-q6':
        source = registry_source(identifier)
        raw_path = ROOT/'experiments/scalar_lift/bh13_6_seed.json'
        record = json.loads(raw_path.read_text())
        raw = np.asarray(record['exponents'], dtype=np.int64)
        rows,columns = record['unreal_row_phase_shifts'],record['unreal_column_phase_shifts']
        shift = -1
    elif identifier == 'aalto-9-6-000':
        source = broad_source(identifier)
        raw = np.asarray(dataset(9,6)[0], dtype=np.int64)
        gauge = source.proof['monomial_gauge_audit']
        rows,columns = gauge['row_shifts'],gauge['column_shifts']
        shift = -source.proof['acquisition']['global_anchor_subtracted_after_gauge']
        raw_path = ROOT/'experiments/scalar_lift/broad_catalogue_inputs/BH-9-6.txt'
    else:
        raise ValueError('only the three audited complete inputs are admitted')
    final = (raw+np.asarray(rows)[:,None]+np.asarray(columns)[None,:]+shift)%6
    if not np.array_equal(final, source.exponents):
        raise ValueError('recorded literal source gauge does not reconstruct the input')
    if set(map(int,np.unique(final))) != {0,1,3,4}:
        raise ValueError('the actual used palette is not covered by this frame')
    raw_proof, final_proof = verify_source_gram(raw),verify_source_gram(final)
    return source, {'raw_source_data':str(raw_path.relative_to(ROOT)),
        'raw_source_data_sha256':_digest(raw_path),
        'raw_source_full_Gram':raw_proof,
        'row_exponent_shifts_mod6':list(map(int,rows)),
        'column_exponent_shifts_mod6':list(map(int,columns)),
        'global_exponent_shift_mod6':int(shift)%6,
        'literal_reconstruction':'(raw[i,j]+row[i]+column[j]+global_shift) mod6',
        'gauged_source_full_Gram':final_proof,
        'gauged_source_matches_active_verified_registry':True,
        'active_source_proof':source.proof}


@dataclass
class ReusedEndpoint:
    source_id: str
    exponents: np.ndarray = field(repr=False)
    P: np.ndarray = field(repr=False)
    Q: np.ndarray = field(repr=False)

    def __post_init__(self):
        values=np.asarray(self.exponents)
        if (values.ndim!=2 or values.shape[0]!=values.shape[1]
                or not np.issubdtype(values.dtype,np.integer)
                or not set(map(int,np.unique(values)))<={0,1,3,4}):
            raise ValueError('literal endpoint needs a square source in the certified palette')
        if (self.P.ndim!=2 or self.P.shape[0]!=self.P.shape[1]
                or self.P.shape!=self.Q.shape or not self.P.shape[0]):
            raise ValueError('literal endpoint needs two certified square frames of the same width')

    @property
    def width(self):
        return self.P.shape[0]

    @property
    def order(self):
        return self.width*len(self.exponents)

    @staticmethod
    def _index(value, bound):
        if isinstance(value,(bool,np.bool_)):
            raise ValueError('a Boolean is not a matrix index')
        try:
            value=operator.index(value)
        except TypeError as error:
            raise ValueError('matrix indices must be exact integers') from error
        if not 0<=value<bound:
            raise ValueError('matrix index outside the completed endpoint')
        return value

    def entry(self, row, column):
        row,column=self._index(row,self.order),self._index(column,self.order)
        i,a=divmod(row,self.width);j,b=divmod(column,self.width)
        exponent=int(self.exponents[i,j]);matrix=self.P if exponent%3==0 else self.Q
        return (1 if exponent<3 else -1)*int(matrix[a,b])

    def row(self, row):
        row=self._index(row,self.order);i,a=divmod(row,self.width)
        return np.concatenate([(1 if int(e)<3 else -1)*(self.P if int(e)%3==0 else self.Q)[a]
                               for e in self.exponents[i]])


def endpoint(identifier):
    source,_=verified_input(identifier)
    P,Q,_=verified_frames()
    return ReusedEndpoint(identifier,source.exponents,P,Q)


def certificate():
    P,Q,frame_proof=verified_frames()
    outputs=[]
    for identifier in SOURCES:
        source,gauges=verified_input(identifier)
        literal=ReusedEndpoint(identifier,source.exponents,P,Q)
        selected=[0,WIDTH-1,WIDTH,WIDTH+1,literal.order-1]
        rows=np.array([literal.row(i) for i in selected],dtype=np.int64)
        if not np.array_equal(rows@rows.T,literal.order*np.eye(len(selected),dtype=np.int64)):
            raise ValueError('independent literal endpoint row controls failed')
        for i in selected:
            row=literal.row(i)
            for j in selected:
                if literal.entry(i,j)!=int(row[j]):
                    raise ValueError('literal entry and row recipes disagree')
        recipe={'source_id':identifier,'source_order':source.order,'source_phase':6,
            'source_hash':source.source_hash,'exponents':source.exponents.tolist(),
            'source_used_alphabet':list(source.used_alphabet),'source_gram_scale':source.order,
            'source_gauges_and_exact_proofs':gauges,
            'frame_width':WIDTH,'frame_gram_scale':WIDTH,'output_order':literal.order,
            'output_gram_scale':WIDTH*source.order,'actual_endpoint_alphabet':[-1,1],
            'literal_entry_recipe':'divmod(i,2524)=(r,a), divmod(j,2524)=(c,b); exponent=source[r,c]; return P[a,b],Q[a,b],-P[a,b],-Q[a,b] for exponent0,1,3,4',
            'full_endpoint_Gram_proof':'B(z)=rho(z)P; sum_k B(Hik) B(Hjk)^T =2524*rho(sum_k Hik*conjugate(Hjk))=2524*n*delta_ij*I',
            'complete_algebraic_endpoint_Gram_verified':True,
            'all_used_entry_images_verified_sign':True,
            'full_dense_endpoint_Gram_performed':False,
            'literal_endpoint_sampled_rows':selected,
            'all_25_sampled_integer_row_products_verified':True,
            'new_order_status':'new_to_best_of_authors_knowledge',
            'prior_construction_identified':False,
            'new_Hadamard_order_claimed':True,
            'new_order_claim':f'To the best of the author\'s knowledge, {literal.order} is a new real Hadamard order.',
            'new_equivalence_class_claimed':False}
        recipe['literal_recipe_sha256']=_json_digest(recipe)
        outputs.append(recipe)
    return {'schema':'completed-h7572-frame-reuse-v1','frame_proof':frame_proof,
            'outputs':outputs,'output_orders':[r['output_order']for r in outputs],
            'novelty_convention':dict(NOVELTY_CONVENTION),
            'all_three_are_completed_literal_sign_endpoint_recipes':True,
            'large_endpoint_matrices_allocated':False,
            'new_real_orders_claimed':True,
            'new_real_orders_claimed_to_best_of_authors_knowledge':[r['output_order']for r in outputs],
            'new_equivalence_classes_claimed':False}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--report',type=Path)
    group.add_argument('--check',type=Path)
    args=parser.parse_args();report=certificate()
    if args.report is not None:
        args.report.parent.mkdir(parents=True,exist_ok=True)
        args.report.write_text(json.dumps(report,indent=2)+'\n')
    elif json.loads(args.check.read_text())!=report:
        raise ValueError('frame reuse certificate differs from exact replay')
    print(json.dumps({'output_orders':report['output_orders'],
        'completed_literal_sign_recipes':True,'full_algebraic_Grams_verified':True,
        'new_real_orders_claimed':report['new_real_orders_claimed'],
        'new_real_orders_claimed_to_best_of_authors_knowledge':report['new_real_orders_claimed_to_best_of_authors_knowledge']},indent=2))


if __name__=='__main__':main()

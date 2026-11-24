"""Replay a completed twisted T-frame with every admitted sixth-root source.

This performs the full frame and order-three endpoint verification afresh.
Larger endpoints have literal entry/row access and exact algebraic Gram
proofs. Their novelty is a separate, dated literature question.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import tempfile

import numpy as np

from experiments.scalar_lift import hadamard_frame_search as search
from experiments.scalar_lift.reuse_hadamard_7572_frame import (
    ReusedEndpoint, SOURCES, verified_input, verify_source_gram,
)


def witness(record):
    """Read either an endpoint record or a retained positive-prime certificate."""
    if 'known_order_positive_repeat' in record:
        record = record['known_order_positive_repeat']['endpoint']
    for key in ('p', 'm', 'mu', 'states'):
        if key not in record:
            raise ValueError(f'missing completed-witness field: {key}')
    return record


def replay(record, directory):
    record = witness(record)
    data = search.chart(record['p'], record['m'], record['mu'])
    # This verifies all shifts, the actual P/Q integer Grams and symmetric
    # cross identity, the actual action alphabet, and every H(12p) row pair.
    proof = search.verify_endpoint(record['states'], data, directory)
    b = search.expand(record['states'], data)
    a = search.D4 @ b
    P, Q = search.gs(a-b), search.gs(a+b)
    P.setflags(write=False)
    Q.setflags(write=False)
    return P, Q, proof


def complete_sources():
    exponents = np.array(((3,1,1),(1,3,1),(1,1,3)), dtype=np.int64)
    yield 'diagonal-minus-one-BH3-6', exponents, verify_source_gram(exponents)
    for identifier in SOURCES:
        source, proof = verified_input(identifier)
        yield identifier, source.exponents, proof


def certificate(record, directory):
    P, Q, frame = replay(record, directory)
    width = len(P)
    outputs = []
    for identifier, exponents, source_proof in complete_sources():
        literal = ReusedEndpoint(identifier, exponents, P, Q)
        selected = sorted({0, width-1, width, width+1, literal.order-1})
        rows = np.array([literal.row(i) for i in selected], dtype=np.int64)
        if not np.array_equal(rows @ rows.T,
                              literal.order * np.eye(len(selected), dtype=np.int64)):
            raise ValueError('literal endpoint row checks failed')
        for i, row in zip(selected, rows):
            for j in selected:
                if literal.entry(i,j) != int(row[j]):
                    raise ValueError('literal entry and row recipes disagree')
        outputs.append({
            'source_id': identifier,
            'source_order': len(exponents),
            'source_phase': 6,
            'source_exponents': exponents.tolist(),
            'actual_source_exponents': sorted(map(int,np.unique(exponents))),
            'source_complete_exact_Gram_proof': source_proof,
            'source_gram_scale': len(exponents),
            'frame_dimension': width,
            'frame_gram_scale': width,
            'output_order': literal.order,
            'output_gram_scale': literal.order,
            'actual_endpoint_alphabet': [-1,1],
            'all_used_entry_images_verified_sign': True,
            'literal_entry_recipe': 'With d=frame_dimension, divmod(i,d)=(r,a), divmod(j,d)=(c,b); source exponent 0,1,3,4 gives P[a,b],Q[a,b],-P[a,b],-Q[a,b].',
            'complete_algebraic_endpoint_Gram_verified': True,
            'full_Gram_identity': 'B(z)=rho(z)P; sum_k B(Hik)B(Hjk)^T=d*rho(sum_k Hik*conjugate(Hjk))=d*n*delta_ij*I.',
            'complete_endpoint_row_pair_check_performed': len(exponents)==3,
            'sampled_literal_rows': selected,
            'sampled_literal_row_Grams_verified': True,
            'new_order_status': 'unverified; compare this completed endpoint with dated primary records',
            'new_equivalence_class_claimed': False,
        })
    return {
        'schema': 'completed-generic-sixth-root-frame-reuse-v1',
        'frame_proof': frame,
        'field': 'Q(zeta_6)/Q',
        'field_degree': 2,
        'module_rank': width//2,
        'rational_dimension': width,
        'conjugation': 'zeta_6 -> 1-zeta_6; adjoint is transpose',
        'canonical_trace_metric': [[2,1],[1,2]],
        'incoming_trace_form': frame['incoming_form'],
        'action': frame['action'],
        'flat_entry_images': {'0':'P','1':'Q','3':'-P','4':'-Q'},
        'unused_nonflat_images': {'2':'Q-P','5':'P-Q'},
        'outputs': outputs,
        'output_orders': [item['output_order'] for item in outputs],
        'new_orders_claimed': False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--witness', type=Path, required=True)
    parser.add_argument('--report', type=Path, required=True)
    parser.add_argument('--verified-directory', type=Path)
    args = parser.parse_args()
    record = json.loads(args.witness.read_text())
    if args.verified_directory is None:
        with tempfile.TemporaryDirectory(prefix='scalar-frame-replay-') as folder:
            report = certificate(record, Path(folder)/'verified')
    else:
        report = certificate(record, args.verified_directory)
    report['input_record_sha256'] = hashlib.sha256(args.witness.read_bytes()).hexdigest()
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'output_orders':report['output_orders'],
                      'full_frame_checks_passed':True,
                      'all_endpoint_algebraic_Grams_verified':True,
                      'new_orders_claimed':False}))


if __name__ == '__main__':
    main()

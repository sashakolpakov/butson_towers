"""Exact two-basis restrictions anchored at the published H668 and H716.

The output field basis is changed by a dense rational Householder isometry.
The input basis is its product with the fixed published Hadamard frame F.
Thus B(1)=F and B(zeta6)=Q R0 Q^T F, without changing a source matrix.
This is a bounded coordinate probe, not a search for new sign sources.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np

from experiments.scalar_lift.alpoge_published_inputs import known_h716
from experiments.scalar_lift.ambient_field_candidates import (
    _quaternion,
    load_ambient_source,
    period_field,
)
from experiments.scalar_lift.h668_width12 import known_h668_width12


def twice_root_action(values):
    """Apply I+L_(i+j+k), with coordinates (local, period)."""
    values = np.asarray(values, dtype=np.int64)
    d = values.shape[0]
    if d % 4:
        raise ValueError("a multiple of four is required")
    t = np.eye(4, dtype=np.int64) + _quaternion(0, 1, 1, 1)
    return np.tensordot(t, values.reshape((4, d//4) + values.shape[1:]),
                        axes=(1, 0)).reshape(values.shape)


def weighted_trace_model():
    identity = np.eye(4, dtype=np.int64)
    d = _quaternion(0, 1, 1, 1)
    j = np.array([[0, -1, 0, 0], [1, 0, 0, 0],
                  [0, 0, 0, 1], [0, 0, -1, 0]], dtype=np.int64)
    k = d @ j
    theta = 3*identity-k
    alpha = (2*identity, 2*j, identity+d-j+k, -identity+d-j-k)
    if not np.array_equal(k @ k, 3*identity):
        raise ArithmeticError('the real weight generator has the wrong square')
    for a in range(4):
        if not np.array_equal(alpha[a][:, 0], 2*identity[:, a]):
            raise ArithmeticError('field basis identification failed')
        for b in range(4):
            if int(np.trace(theta @ alpha[a].T @ alpha[b])) != 48*(a == b):
                raise ArithmeticError('the weighted trace basis failed')
    return {'local_field': 'Q(zeta12)=Q(D,J), D^2=-3, J^2=-1, DJ=JD',
            'weight': 'theta=(3-DJ)/12; both real embeddings are positive',
            'local_basis': ['1', 'J', '(1+D-J+DJ)/2', '(-1+D-J-DJ)/2'],
            'local_weighted_trace_Gram': 'I4, all sixteen entries checked exactly',
            'auxiliary_basis': 'period basis eta transformed by the certified S, S^T G S=I',
            'full_metric': 'Tr_E/Q(theta * conjugate(x) * y), Gram I_(4r)',
            'unweighted_canonical_trace_Gram_claimed': False}


def reflected_partner(frame, vector=None):
    """Return exact integer numerator and denominator for Q R0 Q^T F.

    Q=I-2vv^T/m, m=v^Tv. The identity v^T(2R0)v=m avoids
    dense matrix products and permits exact int64 arithmetic at these sizes.
    """
    frame = np.asarray(frame, dtype=np.int64)
    if frame.ndim != 2 or frame.shape[0] != frame.shape[1] or not np.all(np.abs(frame) == 1):
        raise ValueError("a square sign frame is required")
    a = twice_root_action(frame)
    if vector is None:
        return a, 2
    v = np.asarray(vector, dtype=np.int64)
    if v.shape != (len(frame),) or np.max(np.abs(v)) > 16 or len(frame) > 1024:
        raise ValueError("vector outside the exact bounded arithmetic domain")
    m = int(v @ v)
    if not m:
        raise ValueError("the reflecting vector is zero")
    b, c, row = twice_root_action(v), v @ frame, v @ a
    numerator = m*a - 2*np.outer(v, row) + 2*np.outer(2*v-b, c)
    return numerator, 2*m


def canonical_fraction(numerator, denominator):
    factor = math.gcd(int(np.gcd.reduce(numerator.ravel())), denominator)
    return numerator//factor, denominator//factor


def tensor_minor(vector):
    """Exhibit a nonproduct column of the new output basis Q^T."""
    v = np.asarray(vector, dtype=np.int64)
    m, r = int(v @ v), len(v)//4
    for p in np.flatnonzero(v):
        column = -2*v*int(v[p])
        column[p] += m
        table = column.reshape(4, r)
        for a in range(4):
            for b in range(a+1, 4):
                nonzero = np.flatnonzero((table[a] != 0) | (table[b] != 0))
                if not len(nonzero):
                    continue
                c = int(nonzero[0])
                determinants = table[a, c]*table[b] - table[b, c]*table[a]
                hits = np.flatnonzero(determinants)
                if len(hits):
                    e = int(hits[0])
                    return {"basis_column": int(p), "local_rows": [a, b],
                            "period_columns": [c, e], "integer_minor": int(determinants[e]),
                            "column_denominator": m}
    return None


def frame_inputs():
    result = []
    for order, loader, border in ((668, known_h668_width12, 12), (716, known_h716, 4)):
        matrix, proof = loader()
        m, b = (order-border)//4, border//4
        permutation = np.concatenate([np.r_[np.arange(i*b, (i+1)*b),
                                     np.arange(border+i*m, border+(i+1)*m)] for i in range(4)])
        result.append((order, matrix[permutation], proof))
    return result


def vectors(frame):
    for j in range(16):
        yield f"published-column-{j}", frame[:, j].copy()
    for j in range(1, 17):
        for coefficient in (1, 2):
            yield f"published-columns-0-{j}-weight-{coefficient}", frame[:, 0] + coefficient*frame[:, j]
    rng = np.random.default_rng(20261123 + len(frame))
    for j in range(16):
        yield f"integer-vector-{j}", rng.integers(-2, 3, size=len(frame), dtype=np.int64)


def entry_statistics(numerator, denominator):
    signs = np.abs(numerator) == denominator
    zeros = numerator == 0
    integers = numerator % denominator == 0
    bad = np.argwhere(~signs)
    witness = None
    if len(bad):
        i, j = map(int, bad[0])
        witness = {"row": i, "column": j, "numerator": int(numerator[i, j]),
                   "denominator": denominator}
    return {"sign_entries": int(np.count_nonzero(signs)), "zero_entries": int(np.count_nonzero(zeros)),
            "other_integer_entries": int(np.count_nonzero(integers & ~signs & ~zeros)),
            "noninteger_entries": int(np.count_nonzero(~integers)),
            "total_entries": numerator.size, "all_entries_signs": not len(bad),
            "first_failed_entry": witness}


def binary_rank(matrix):
    pivots = {}
    for row in np.asarray(matrix, dtype=np.uint8) % 2:
        value = int.from_bytes(np.packbits(row, bitorder='little').tobytes(), 'little')
        while value:
            pivot = value.bit_length()-1
            if pivot not in pivots:
                pivots[pivot] = value
                break
            value ^= pivots[pivot]
    return len(pivots)


def restricted_entry(source, frame, partner_numerator, partner_denominator, row, column):
    """Read any entry of the complete restricted source, without expanding it."""
    width = len(frame)
    if source.phase != 6 or not 0 <= row < source.order*width or not 0 <= column < source.order*width:
        raise ValueError('source or coordinate outside the sextic restriction')
    i, a = divmod(row, width)
    j, b = divmod(column, width)
    f = int(frame[a, b])*partner_denominator
    g = int(partner_numerator[a, b])
    value = (f, g, g-f, -f, -g, f-g)[source.entry(i, j)]
    return Fraction(value, partner_denominator)


def run_case(order, frame, proof):
    records, seen = [], set()
    for name, v in [("identity-output-basis", None), *vectors(frame)]:
        numerator, denominator = canonical_fraction(*reflected_partner(frame, v))
        identity = hashlib.sha256(numerator.astype('<i8').tobytes()
                                  + str(denominator).encode()).hexdigest()
        duplicate = identity in seen
        seen.add(identity)
        records.append({"basis": name, "vector": None if v is None else v.tolist(),
                        "partner_sha256": identity, "partner_denominator": denominator,
                        "literal_duplicate": duplicate,
                        "nonproduct_basis_column": None if v is None else tensor_minor(v),
                        "entry_statistics": entry_statistics(numerator, denominator)})
    best = max(records, key=lambda row: row['entry_statistics']['sign_entries'])
    targets = [4676, 8684] if order == 668 else [9308]
    baseline_twice = twice_root_action(frame)
    if np.any(baseline_twice % 2):
        raise ArithmeticError('the four-term half sums are nonintegral')
    rank = binary_rank(baseline_twice//2-1)
    field = period_field(order//4)
    sources = ['catalogue-p7-q6', 'unreal13-q6'] if order == 668 else ['unreal13-q6']
    return {"published_frame_order": order, "published_frame_proof": proof,
            "ambient_field": field.proof(),
            "normalized_auxiliary_operator": "I4 tensor (S^-1 M_eta S); conjugate by the displayed Q",
            "frame_sha256": hashlib.sha256(frame.astype('<i8').tobytes()).hexdigest(),
            "frame_row_order": "four channels: own border rows, then own core rows",
            "target_orders": targets, "candidates": records, "candidate_visits": len(records),
            "verified_sources": [load_ambient_source(source).metadata() for source in sources],
            "distinct_partners": len(seen), "best_basis": best['basis'],
            "best_sign_entries": best['entry_statistics']['sign_entries'],
            "total_partner_entries": order**2,
            "binary_rank_baseline_minus_all_ones": rank,
            "minimum_number_of_Householder_factors": (rank+1)//2,
            "rank_gate": "rank_F2(R0 F-J) <= rank_Q(G-R0 F) <= 2 rank(Q-I) <= 2k",
            "all_single_Householder_output_bases_excluded": rank > 2,
            "sign_partner_found": any(row['entry_statistics']['all_entries_signs'] for row in records)}


def report():
    cases = [run_case(*case) for case in frame_inputs()]
    return {"schema": "alpoge-two-basis-householder-probe-v1", "cases": cases,
            "coordinate_map": "B(z)=Q rho_0(z) Q^T F; Q=I-2vv^T/(v^Tv)",
            "input_basis": "output_basis times the fixed published F",
            "weighted_trace_model": weighted_trace_model(),
            "basis_metric_identity": "for the stated weighted trace form, G_out=I and G_in=F^T F=d I",
            "ordinary_output_Gram": "n*d I for the complete restricted source",
            "sextic_pair_identity": "G G^T=d I and G F^T+F G^T=d I by exact rational conjugacy",
            "field_faithfulness": "conjugate both known cyclotomic and auxiliary period operators by Q",
            "source_palette": "{1,zeta6,-1,-zeta6}; published BH(7,6) and saved BH(13,6)",
            "no_source_construction_or_SAT": True, "pending_declared_candidates": 0,
            "new_sign_Hadamard_found": any(c['sign_partner_found'] for c in cases),
            "scope": "These fixed F inputs and the 64 displayed dense rational reflections per input; not all unequal bases"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--output', type=Path)
    group.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = report()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    elif json.loads(args.check.read_text()) != result:
        raise ValueError('two-basis probe certificate differs from exact replay')
    print(json.dumps({"candidate_visits": sum(c['candidate_visits'] for c in result['cases']),
                      "distinct_partners": sum(c['distinct_partners'] for c in result['cases']),
                      "new_sign_Hadamard_found": result['new_sign_Hadamard_found']}))


if __name__ == '__main__':
    main()

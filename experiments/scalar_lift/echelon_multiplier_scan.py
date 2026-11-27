"""Replay a bounded, symmetry-first T-input selection audit through order100000.

This is a necessary integer row-norm screen, not a T-input search.  Both odd
coordinate actions are covered.  The orbit bound is explicit: its prime bounds
are proved independently of any list of desired matrix orders.  Surviving
outputs retain the evidence labels of echelon_frontier; most higher outputs
have unassessed literature status.  No job is started by this module.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
import hashlib
import json
from math import gcd, isqrt, lcm
from pathlib import Path

from sympy import divisors, factorint, isprime, primerange
from sympy.ntheory import n_order

from .echelon_frontier import build_report as frontier_report, candidate_status
from .low_gap_target_scan import identity_dc_rows, c3_dc_rows
from .multiplier_target_scan import parameter_candidates, prime_bound
from .williamson_symmetry import select_charts as williamson_charts

ROOT = Path(__file__).resolve().parents[2]
CERTIFICATE = Path('certificates/scalar-lift/echelon-multiplier-selection.json')
SURVEY = 'https://arxiv.org/html/2610.27175v2'
SEBERRY = 'https://documents.uow.edu.au/~jennie/WEBPDF/139_1990.pdf'
MIYAMOTO = 'https://github.com/vatskishan03/Hadamard/blob/main/hadamard/MIYAMOTO_BLOCKS.md'


def knownness_provenance():
    """Attribute each extra positive gate, separately from literal inputs."""
    return {
        'known-T-times-Williamson-type': {
            'product': {'url': 'https://doi.org/10.1017/S0004972700045019',
                        'location': 'Cooper-Wallis Theorem3 and Corollary4, p.274'},
            'literal_Williamson_widths': ['experiments/scalar_lift/williamson_inputs.json',
                                        'experiments/scalar_lift/near_williamson_inputs.json'],
            'additional_Williamson_widths': {
                '193': {'url': SEBERRY, 'location': 'Corollary29(ii), explicit list, printed p.121'},
                '241': {'url': SEBERRY, 'location': 'Corollary29(ii), explicit list, printed p.121'},
                '433': {'url': SEBERRY, 'location': 'Corollary29(ii), explicit list, printed p.121'},
                '397': {'url': SEBERRY+'#page=24', 'location': 'Corollary29(ii), p.120; conference order198'},
                '761': {'url': SEBERRY, 'location': 'Corollary29(ii), explicit list, printed p.121'},
                '883': {'url': SEBERRY, 'location': 'Corollary31 and Remark32, explicit list, printed p.122'}}},
        'known_T_inputs': {
            'odd_at_most100': {'url': 'https://arxiv.org/html/0912.5091#S3',
                               'location': 'Base-sequence-to-T construction (2.3) and section3 existence range'},
            'Golay_and_base_sequence_families': {'url': SURVEY+'#S4',
                'location': 'T-matrix constructions and Table3: Golay lengths, 2g+1, g+h, and TT(n)'},
            'published_special_lengths': {'lengths': [103, 109, 151, 163, 181, 197],
                                          'url': SURVEY+'#S4.T3', 'location': 'Table3'},
            'completed_exact_project_lengths': {'lengths': [127, 271, 547, 631],
                'record': 'certificates/scalar-lift/completed-order-inventory.json'}},
        'known-T-times-Turyn-Williamson': {'url': SURVEY+'#S2.SS1',
            'location': 'Turyn Williamson construction; q=2s-1 a prime power congruent1 modulo4'},
        'known-T-times-commuting-Miyamoto-blocks': {'url': MIYAMOTO,
            'location': 'Explicit commuting 16-block sources of sizes4*269,4*509,4*613,4*809 and T-input substitution',
            'block_sizes': [269, 509, 613, 809]}}


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True,
                                     separators=(',', ':')).encode()).hexdigest()


def identity_parameter_rows(bound):
    """Enumerate (p,m,row) from N=2h+m(h²+a²+b²+c²), independently."""
    result = set()
    radius = isqrt(bound+1)
    for h in range(-1-radius, radius):
        tail_bound = bound-2*h-h*h
        if tail_bound < 0:
            continue
        for a in range(isqrt(tail_bound)+1):
            for b in range(a+1):
                if a*a+b*b > tail_bound:
                    break
                for c in range(b+1):
                    S = h*h+a*a+b*b+c*c
                    if not S or a*a+b*b+c*c > tail_bound:
                        continue
                    for m in range(1, (bound-2*h)//S+1, 2):
                        N = 2*h+m*S
                        if not 2 <= N <= bound or N % 2:
                            continue
                        minimum = abs(h)+a+b+c
                        if minimum > N or (N-minimum) % 2:
                            continue
                        p = 1+m*N
                        if isprime(p):
                            result.add((p, m, 1+m*h, m*a, m*b, m*c))
    return result


def arithmetic_charts(bound=74):
    """Cross-check prime-first against independent parameter-first screens."""
    if type(bound) is not int or bound < 2:
        raise ValueError('orbit bound must be an integer at least2')
    rows, identity, c3 = [], set(), set()
    bounds = {'identity': 1+bound*(bound+2)//2, 'C3': prime_bound(bound)}
    for pp in primerange(3, max(bounds.values())+1):
        p = int(pp)
        for mm in divisors(p-1):
            m = int(mm)
            N = (p-1)//m
            if m % 2 == 0 or N > bound:
                continue
            for action, fn in [('identity', identity_dc_rows), ('C3', c3_dc_rows)]:
                if p > bounds[action] or (action == 'C3' and m % 3):
                    continue
                necessary = fn(p, m)
                if not necessary:
                    continue
                if action == 'identity':
                    identity.update((p, m, *r['row']) for r in necessary)
                else:
                    c3.update((p, m, r['head_sign_sum'], r['tail_sign_sum'])
                              for r in necessary)
                rows.append({'p': p, 'm': m, 'coordinate_action': action,
                             'nonzero_orbits': N, 'necessary_DC_rows': necessary})
    if identity != identity_parameter_rows(bound):
        raise ArithmeticError('independent identity row screens differ')
    if c3 != parameter_candidates(bound):
        raise ArithmeticError('independent C3 row screens differ')
    return sorted(rows, key=lambda r: (r['nonzero_orbits'], r['p'],
                                      r['coordinate_action'], -r['m']))


@lru_cache(None)
def prime_power(n):
    return n > 1 and len(factorint(n)) == 1


def known_t_lengths(limit):
    """Positive lengths from the retained published construction inventory."""
    golay = {2**a*10**b*26**c for a in range(limit.bit_length())
             for b in range(6) for c in range(5)
             if 2**a*10**b*26**c <= limit}
    return ({n for n in range(1, 101, 2)} |
            {3*n-1 for n in range(2, 47, 2)} | golay |
            {2*g+1 for g in golay if 2*g+1 <= limit} |
            {g+h for g in golay for h in golay if g+h <= limit} |
            {103, 109, 151, 163, 181, 197, 127, 271, 547, 631})


def build_report(root=ROOT, *, max_orbits=74, ceiling=100000):
    root = Path(root)
    inputs = {}

    def load(relative):
        raw = (root/relative).read_bytes()
        inputs[str(relative)] = hashlib.sha256(raw).hexdigest()
        return json.loads(raw)

    widths = set()
    for path, key in [('williamson_inputs.json', 'order'),
                      ('near_williamson_inputs.json', 'width')]:
        record = load(Path('experiments/scalar_lift')/path)
        widths.update(row[key] for row in record['entries'])
    # Finite exact sixth-root sources retained in the scalar-endpoint package.
    # The first four widths already occur above;117 and351 are additional.
    butson = {3, 7, 9, 13, 39, 117, 351}
    load(Path('certificates/scalar-lift/hadamard-7572-frame-reuse.json'))
    load(Path('certificates/scalar-lift/tripled-unreal-scalar-endpoints.json'))
    frontier = frontier_report(root, ceiling, include_candidate_pool=True)
    pool = set(frontier['candidate_pool'])
    charts = arithmetic_charts(max_orbits)
    all_charts_digest = digest(charts)
    strongest = {}
    for chart in charts:
        key = chart['p'], chart['coordinate_action']
        if key not in strongest or chart['m'] > strongest[key]['m']:
            strongest[key] = chart
    knownT = known_t_lengths(ceiling//4)
    knownWT = widths | {193, 241, 397, 433, 761, 883}
    commuting_miyamoto = {269, 509, 613, 809}

    @lru_cache(None)
    def positive_product(n):
        quarter = n//4
        for tt in divisors(quarter):
            t = int(tt)
            if t not in knownT:
                continue
            s = quarter//t
            if s in knownWT:
                return {'construction': 'known-T-times-Williamson-type', 't': t, 's': s}
            if s in commuting_miyamoto:
                return {'construction': 'known-T-times-commuting-Miyamoto-blocks', 't': t, 's': s}
            if (2*s-1) % 4 == 1 and prime_power(2*s-1):
                return {'construction': 'known-T-times-Turyn-Williamson', 't': t, 's': s}
        return None

    survivors, excluded = [], []
    for key, original in sorted(strongest.items()):
        row = dict(original)
        p, m, action = row['p'], row['m'], row['coordinate_action']
        N = row['nonzero_orbits']
        outputs = {}
        for width in sorted(widths | butson):
            base = 4*p*width
            order, exponent = base, 0
            base_gate = positive_product(base)
            while order <= ceiling:
                if order in pool:
                    if base_gate:
                        excluded.append([p, m, action, width, exponent, order, base_gate])
                    else:
                        outputs.setdefault(order, {'order': order,
                            'status': candidate_status(order, frontier), 'routes': []})['routes'].append({
                                'width': width, 'source_kind': ('literal-Williamson-or-near-Williamson'
                                    if width in widths else 'finite-exact-sixth-root-Butson'),
                                'tensor_H2_exponent': exponent})
                order *= 2
                exponent += 1
        if not outputs:
            continue
        order2 = int(n_order(2, p))
        H = lcm(order2//gcd(order2, 2), m)
        field = H//m
        row.update({'unfiltered_axis_bits': 2*N,
            'unfiltered_axis_word_count': str(4**N),
            'unfiltered_axis_and_sign_bits': 3*N,
            'frequency_arithmetic': {'ord_p_2': order2, 'H_order': H,
                'extension_index': field, 'half_frequency_orbits': (p-1)//(2*H),
                'scope': 'Factor-size arithmetic only; no parity-algebra validation.'},
            'conditional_generic_fixed_axis_sign_bits_estimate': max(0, N//2-(2 if action == 'C3' else 4)),
            'sign_estimate_scope': 'Conditional generic affine-rank estimate; actual rank depends on each axis word.',
            'retained_outputs': [outputs[n] for n in sorted(outputs)]})
        # Bind existing exact parity maps without treating them as sign witnesses.
        name = (f'identity-parity-p{p}.json' if action == 'identity'
                else f'hadamard-frame-p{p}.json')
        path = Path('certificates/scalar-lift')/name
        if (root/path).exists():
            cert = load(path)
            if cert.get('p') == p and cert.get('m') == m:
                row['retained_parity_map'] = {'path': path.as_posix(),
                    'field_degree': cert.get('field_degree'),
                    'validation': cert.get('parity_map_validation', cert.get('validation')),
                    'scope': 'Retained exact unsigned parity map; no signed T-input claimed.'}
        survivors.append(row)
    survivors.sort(key=lambda r: (r['nonzero_orbits'], r['p'], r['coordinate_action']))
    by_tier = []
    for tier in frontier['tiers']:
        selected = sorted({o['order'] for r in survivors for o in r['retained_outputs']
                           if tier['lower_inclusive'] <= o['order'] <= tier['upper_inclusive']})
        explicit = [n for n in selected if candidate_status(n, frontier) != 'unassessed_candidate']
        by_tier.append({'name': tier['name'], 'selected_order_count': len(selected),
                        'source_listed_selected_orders': explicit,
                        'unassessed_selected_count': len(selected)-len(explicit)})
    # Replay the five comparison rows used in the strategy note. These are
    # dense Williamson sources rather than T-inputs; no coefficient table or
    # all-shift parity audit is inferred from their row-sum screen.
    comparison_lengths = [1, 2, 3, 5, 7, 9, 11, 13, 127, 271, 547, 631]
    comparison = []
    for record in williamson_charts(sorted(pool), comparison_lengths,
                                    max_orbits=80, order_ceiling=ceiling+1):
        if record['p'] not in {463, 919, 1801, 1951, 2269}:
            continue
        record['unfiltered_axis_bits'] = 2*record['nonzero_orbits']
        record['unfiltered_axis_word_count'] = str(4**record['nonzero_orbits'])
        for route in record['target_routes']:
            route['status'] = candidate_status(route['target_order'], frontier)
        comparison.append(record)
    math_path = Path('experiments/scalar_lift/williamson_symmetry.py')
    inputs[str(math_path)] = hashlib.sha256((root/math_path).read_bytes()).hexdigest()
    return {'schema': 'symmetry-first-T-multiplier-echelon-selection-v1',
        'scope': 'Exact necessary integer norm/support screens; no matrix search or new construction.',
        'max_nonzero_orbits': max_orbits, 'max_output_order_inclusive': ceiling,
        'prime_bounds': {'identity': 1+max_orbits*(max_orbits+2)//2, 'C3': prime_bound(max_orbits)},
        'identity_prime_bound_proof': [
            'p=1+mN=(1+mh)^2+m^2(a^2+b^2+c^2), so N=2h+mS.',
            'Odd prime p forces S>=2. N+2-2m=2h+2+m(S-2)>=0.',
            'For h>=-1 this is immediate; for h<=-2 use S>=h^2 and m>=1, giving h(h+2)>=0.',
            'Therefore m<=(N+2)/2 and p<=1+N(N+2)/2.'],
        'independent_parameter_and_prime_screens_agree': True,
        'all_arithmetic_chart_count': len(charts), 'all_arithmetic_charts_sha256': all_charts_digest,
        'strongest_action_chart_count_before_endpoint_gates': len(strongest),
        'source_widths': sorted(widths), 'additional_finite_Butson_widths': sorted(butson-widths),
        'power_two_daughters_included': True, 'input_sha256': inputs,
        'frontier_pool_sha256': frontier['candidate_pool_sha256'],
        'positive_product_exclusion_count': len(excluded),
        'positive_product_exclusions_sha256': digest(excluded),
        'positive_knownness_provenance': knownness_provenance(),
        'positive_knownness_references': [
            'https://arxiv.org/html/2610.27175v2#S4.T3',
            'https://documents.uow.edu.au/~jennie/WEBPDF/139_1990.pdf',
            'https://github.com/vatskishan03/Hadamard/blob/main/hadamard/MIYAMOTO_BLOCKS.md'],
        'successful_control': {'p': 631, 'm': 21, 'coordinate_action': 'C3',
            'nonzero_orbits': 30, 'unfiltered_axis_bits': 60,
            'unsigned_parity_field_degree': 15, 'checked_fixed_axis_free_sign_bits': 13,
            'scope': 'Completed exact source, distinct from the necessary screens above.'},
        'tiers': by_tier, 'strongest_surviving_charts': survivors,
        'dense_Williamson_comparison': {
            'scope': 'Five strategy rows replayed by exact DC/support selection; correlations remain unsolved.',
            'known_T_lengths': comparison_lengths, 'max_nonzero_orbits': 80,
            'rows': comparison,
            'fixed_axis_sign_dimension': 0,
            'sign_dimension_scope': 'Signs are forced by the modulo-four argument for this dense model. Axis choices remain4^N.',
            'full_parity_map_or_sign_witness_claimed_here': False},
        'interpretation': {'survivor_proves_existence': False,
            'absence_from_positive_gates_proves_global_unknownness': False,
            'symmetry_class_exclusion_proves_nonexistence': False,
            'matrix_search_performed': False}}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-orbits', type=int, default=74)
    parser.add_argument('--ceiling', type=int, default=100000)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args(argv)
    record = build_report(max_orbits=args.max_orbits, ceiling=args.ceiling)
    target = args.output or ROOT/CERTIFICATE
    if args.check:
        if json.loads(target.read_text()) != record:
            raise SystemExit('selection certificate differs')
    elif args.output:
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(record, indent=2)+'\n')
    print(json.dumps({'arithmetic_charts': record['all_arithmetic_chart_count'],
                      'surviving_strongest_charts': len(record['strongest_surviving_charts']),
                      'tiers': record['tiers']}, indent=2))


if __name__ == '__main__':
    main()

"""Dated target provenance and positive construction gates through order 100000.

Explorer's complete implementation table stops at 10000. Higher orders left
by this module are unassessed targets, not an assertion that they are unknown.
The positive gates deliberately form an extensible subset of constructions.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from math import isqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
EXPLORER = Path('certificates/scalar-lift/williamson-symmetry-selection.json')
KUMAR = Path('experiments/scalar_lift/frontier_120_snapshot.json')
COMPLETED = Path('certificates/scalar-lift/completed-order-inventory.json')
SNAPSHOT = Path('experiments/scalar_lift/echelon_frontier_snapshot.json')


def digest(value):
    return hashlib.sha256(json.dumps(value, separators=(',', ':')).encode()).hexdigest()


def prime_powers(limit):
    """Map prime powers up to limit to their prime and positive exponent."""
    sieve = bytearray(b'\1') * (limit + 1)
    sieve[:2] = b'\0\0'
    for p in range(2, isqrt(limit) + 1):
        if sieve[p]:
            sieve[p*p:limit+1:p] = b'\0' * ((limit-p*p)//p + 1)
    powers = {}
    for p in range(2, limit + 1):
        if sieve[p]:
            value, exponent = p, 1
            while value <= limit:
                powers[value] = (p, exponent)
                value *= p
                exponent += 1
    return powers


def positive_known_orders(explorer_gaps, completed, ceiling):
    """Known <=10000 seeds, Paley, a recorded Cooper-Wallis case and products.

    Return an exact finite recipe for each positive order. Exclusion from
    this dictionary means only that these gates supplied no recipe.
    """
    if type(ceiling) is not int or not 10000 <= ceiling <= 100000:
        raise ValueError('ceiling must be an integer from 10000 through 100000')
    gaps = set(explorer_gaps)
    known = {1: {'construction': 'trivial'},
             2: {'construction': 'Sylvester-H2'}}
    known.update({n: {'construction': 'Explorer-v61-positive-table'}
                  for n in range(4, 10001, 4) if n not in gaps})
    for n in completed:
        if type(n) is not int or n < 4 or n % 4:
            raise ValueError('completed orders must be positive multiples of four')
        if n <= ceiling:
            known[n] = {'construction': 'Butson-Towers-completed-ledger'}
    for q, (p, exponent) in prime_powers(ceiling - 1).items():
        if q % 4 == 3:
            known.setdefault(q+1, {'construction': 'Paley-I', 'q': q,
                                  'prime': p, 'exponent': exponent})
        elif q % 4 == 1 and 2*(q+1) <= ceiling:
            known.setdefault(2*(q+1), {'construction': 'Paley-II', 'q': q,
                                      'prime': p, 'exponent': exponent})
    # The recorded classical positive case at Williamson-type width397
    # follows from Paley conference198, Seberry-Yamada Corollary29(ii),
    # and Cooper-Wallis Theorem3/Corollary4 for each odd T-length<=19.
    for t in range(1, 20, 2):
        n = 4*397*t
        if n <= ceiling:
            known.setdefault(n, {'construction': 'Cooper-Wallis-W397',
                                 'Williamson_type_width': 397, 'T_length': t,
                                 'conference_order': 198, 'Paley_q': 197})
    # Every nontrivial factor is smaller than n; one increasing pass closes
    # the positive seed set under all products in this finite order range.
    for n in range(4, ceiling + 1, 4):
        if n in known:
            continue
        for factor in range(2, isqrt(n) + 1):
            if n % factor == 0 and factor in known and n//factor in known:
                known[n] = {'construction': 'Kronecker',
                            'factors': [factor, n//factor]}
                break
    return known


def build_report(root=ROOT, ceiling=100000, *, include_candidate_pool=False):
    """Rebuild three tiers without any network request or matrix search."""
    root = Path(root)
    hashes = {}

    def load(path):
        raw = (root/path).read_bytes()
        hashes[path.as_posix()] = hashlib.sha256(raw).hexdigest()
        return json.loads(raw)

    explorer, kumar, ledger = load(EXPLORER), load(KUMAR), load(COMPLETED)
    current = explorer['frontier']['orders']
    if len(current) != 81 or current != sorted(set(current)):
        raise ValueError('the retained Explorer v61 frontier must have 81 orders')
    if digest(current) != explorer['frontier']['ordered_list_sha256']:
        raise ValueError('Explorer order-list hash differs')
    dated = [n for n in kumar['orders'] if n >= 10000]
    if len(dated) != 37 or any(n > 11996 for n in dated):
        raise ValueError('the retained higher Kumar frontier must have 37 orders')
    completed = ledger['completed_new_real_orders'] + ledger['completed_known_real_orders']
    known = positive_known_orders(current, completed, ceiling)
    candidates = [n for n in range(4, ceiling+1, 4) if n not in known]
    dated_survivors = [n for n in dated if n not in known and n <= ceiling]
    unassessed = [n for n in candidates if n not in current and n not in dated]
    tiers = []
    for lower, upper, name in [(4, 9999, 'below10000'),
                               (10000, 50000, '10000through50000'),
                               (50001, 100000, 'above50000through100000')]:
        upper = min(upper, ceiling)
        if lower > upper:
            continue
        def within(rows):
            return [n for n in rows if lower <= n <= upper]
        values = within(candidates)
        tiers.append({'name': name, 'lower_inclusive': lower,
                      'upper_inclusive': upper,
                      'admissible_order_count': len(range((lower+3)//4*4, upper+1, 4)),
                      'current_explorer_gaps': within(current),
                      'source_listed_kumar_gaps': within(dated_survivors),
                      'unassessed_survivor_count': len(within(unassessed)),
                      'candidate_pool_count': len(values),
                      'candidate_pool_sha256': digest(values)})
    positive = [[n, known[n]] for n in sorted(known)]
    return {'schema': 'hadamard-echelon-frontier-v1',
            'scope': ('Three-tier target selection through 100000. Current source-listed gaps '
                      'are distinguished from candidates surviving incomplete positive gates.'),
            'network_recheck_utc_date': '2026-11-27',
            'input_sha256': hashes,
            'primary_records': {
                'explorer': {'url': 'https://www.math74.com/Hadamard/index.html',
                             'unknown_url': 'https://www.math74.com/Hadamard/unknown.html',
                             'version': 61, 'coverage_max': 10000,
                             'rechecked_index_sha256': 'f4a5f857a36f81c0867ba4e032656ce8765ca7d3afb5a82f06bc9d51d465bdfc',
                             'rechecked_unknown_sha256': 'd752d50acaee5e0aecc30df45ed4f9f30f918982f8783b339fbbeb972a0ee331'},
                'kumar': {'url': 'https://github.com/vatskishan03/Hadamard/blob/0a73016a2a6e80bfd4a36b530ad1e5ed02387947/FRONTIER.md',
                          'rechecked_file_sha256': 'e67e681689e3a02a3426d8c70800244bb7e2b6eb0f024987fa7ca0ff63f0e9bd',
                          'scope': 'Author-listed unresolved odd-quarter orders through 11996; same 37 higher rows as retained snapshot.'},
                'paley_and_kronecker': {'url': 'https://arxiv.org/html/2610.27175v2#S2',
                                         'scope': 'Theorem 2 and preceding Sylvester product rule.'},
                'Williamson_type397': {
                    'source_url': 'https://documents.uow.edu.au/~jennie/WEBPDF/139_1990.pdf#page=24',
                    'source_theorem': 'Seberry-Yamada Corollary29(ii), p.120; Paley conference order198.',
                    'endpoint_url': 'https://doi.org/10.1017/S0004972700045019',
                    'endpoint_theorem': 'Cooper-Wallis Theorem3and Corollary4, p.274; all odd T-lengths<=19.',
                    'retained_record': 'docs/HADAMARD_FRAME_TARGET_SCAN.md'}},
            'positive_gates': ['Explorer-v61positive-table-through10000',
                               'Butson-Towers-completed-ledger',
                               'Paley-I-odd-prime-power-q=3mod4',
                               'Paley-II-odd-prime-power-q=1mod4',
                               'recorded-Williamson-type397-times-odd-T-length-at-most19',
                               'full-bounded-Kronecker-closure'],
            'positive_gate_count_including1and2': len(known),
            'positive_recipe_sha256': digest(positive),
            'excluded_source_listed_orders': [n for n in dated if n in known],
            **({'candidate_pool': candidates} if include_candidate_pool else {}),
            'candidate_pool_sha256': digest(candidates),
            'tiers': tiers,
            'interpretation': {
                'current_explorer_gap': 'No construction in Explorer v61 at the retrieval date, within its stated range.',
                'source_listed_kumar_gap': 'Unresolved in rechecked author frontier; no positive recipe from this finite gate set.',
                'unassessed_candidate': 'No recipe from the declared positive gates; literature status remains unassessed.',
                'full_global_unknown_frontier_through100000_claimed': False,
                'matrix_search_performed': False}}


def candidate_status(order, report):
    """Attribute an order in the retained candidate pool to its evidence class."""
    if 'candidate_pool' not in report:
        raise ValueError('request build_report(include_candidate_pool=True) for membership checks')
    if type(order) is not int or order not in report['candidate_pool']:
        raise ValueError('order is outside the candidate pool')
    for tier in report['tiers']:
        if order in tier['current_explorer_gaps']:
            return 'current_explorer_gap'
        if order in tier['source_listed_kumar_gaps']:
            return 'source_listed_kumar_gap'
    return 'unassessed_candidate'


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args(argv)
    report = build_report()
    if args.check:
        if report != json.loads((ROOT/SNAPSHOT).read_text()):
            raise SystemExit('echelon frontier replay differs')
    if args.output:
        args.output.write_text(json.dumps(report, indent=2)+'\n')
    summary = [{key: tier[key] for key in (
        'name', 'candidate_pool_count', 'unassessed_survivor_count',
        'candidate_pool_sha256')} for tier in report['tiers']]
    print(json.dumps({'candidate_pool_count': sum(
        tier['candidate_pool_count'] for tier in report['tiers']),
        'tiers': summary}, indent=2))


if __name__ == '__main__':
    main()

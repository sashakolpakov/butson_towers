"""Independent arithmetic checks for the bounded positive construction gates."""
import json

import pytest

from experiments.scalar_lift.echelon_frontier import (
    ROOT, SNAPSHOT, build_report, candidate_status, positive_known_orders,
    prime_powers,
)


def test_prime_power_recognition():
    expected = {2:(2,1), 3:(3,1), 4:(2,2), 5:(5,1), 7:(7,1), 8:(2,3),
                9:(3,2), 11:(11,1), 13:(13,1), 16:(2,4), 17:(17,1),
                19:(19,1), 23:(23,1), 25:(5,2), 27:(3,3), 29:(29,1),
                31:(31,1), 32:(2,5)}
    assert prime_powers(35) == expected


def test_positive_recipes_and_product_closure():
    record = build_report(include_candidate_pool=True)
    gaps = record['tiers'][0]['current_explorer_gaps']
    ledger = json.loads((ROOT/'certificates/scalar-lift/completed-order-inventory.json').read_text())
    completed = ledger['completed_new_real_orders'] + ledger['completed_known_real_orders']
    known = positive_known_orders(gaps, completed, 100000)
    # Check every generated parameter independently of the sieve.
    def prime(n):
        return n >= 2 and all(n%d for d in range(2, int(n**0.5)+1))
    for n, recipe in known.items():
        kind = recipe['construction']
        if kind in ('Paley-I', 'Paley-II'):
            p, exponent, q = recipe['prime'], recipe['exponent'], recipe['q']
            assert prime(p) and exponent >= 1 and p**exponent == q
            assert (n == q+1 and q%4 == 3) if kind == 'Paley-I' else (n == 2*(q+1) and q%4 == 1)
        elif kind == 'Cooper-Wallis-W397':
            assert n == 4*397*recipe['T_length']
            assert recipe['T_length'] in range(1, 20, 2)
            assert recipe['conference_order'] == 198
            assert recipe['Paley_q'] == 197 and prime(recipe['Paley_q'])
        elif kind == 'Kronecker':
            a, b = recipe['factors']
            assert a < n and b < n and a in known and b in known and a*b == n
    positives = sorted(known)
    for a in positives:
        for b in positives:
            if a*b > 100000:
                break
            assert a*b in known
    assert known[14292]['construction'] == 'Cooper-Wallis-W397'
    assert not set(gaps).intersection(known)
    assert len(record['candidate_pool']) + len(known) - 2 == 25000


def test_evidence_classes_and_snapshot_replay():
    report = build_report(include_candidate_pool=True)
    assert candidate_status(2676, report) == 'current_explorer_gap'
    assert candidate_status(10012, report) == 'source_listed_kumar_gap'
    assert any(candidate_status(n, report) == 'unassessed_candidate'
               for n in report['candidate_pool'] if n > 50000)
    with pytest.raises(ValueError):
        candidate_status(10000, report)
    compact = dict(report)
    del compact['candidate_pool']
    assert compact == json.loads((ROOT/SNAPSHOT).read_text())


def test_invalid_positive_gate_bounds():
    for bound in (True, 9999, 100001, 50000.0):
        with pytest.raises(ValueError):
            positive_known_orders([], [], bound)

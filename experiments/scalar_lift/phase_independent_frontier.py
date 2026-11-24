"""Finite, phase-independent source/target/anchor restriction preflight.

All inputs are already verified literal matrices or exact implicit sources.
This enumerates necessary conditions, not new Hadamard matrices or basis
searches. Actual root fields are computed from used exponents, rather than
assuming that a phase label is a minimal entry-field conductor.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter
from functools import cache
from pathlib import Path

import numpy as np
from sympy import totient

from experiments.scalar_lift.catalogue_source_acquisition import SOURCE_IDS as CATALOGUE_IDS
from experiments.scalar_lift.catalogue_source_acquisition import data as catalogue_data
from experiments.scalar_lift.catalogue_source_acquisition import load_source as catalogue_source
from experiments.scalar_lift.cyclotomic_towers import representation_obstructions
from experiments.scalar_lift.forward_butson_descent import exponent_hash, verify_source
from experiments.scalar_lift.forward_source_registry import VerifiedSource
from experiments.scalar_lift.forward_source_registry import load_source as registry_source
from experiments.scalar_lift.forward_source_registry import sources as registry_descriptors
from experiments.scalar_lift.nine_source_acquisition import load_source as nine_source
from experiments.scalar_lift.scalar_restriction_obstructions import alphabet_obstruction

ROOT = Path(__file__).resolve().parents[2]
FRONTIER = ROOT/'certificates/scalar-lift/forward-novelty-frontier-2026-11-22.json'
CONTROL10 = ROOT/'certificates/scalar-lift/modified-gs-butson-control.json'
CONTROL12_ROWS = (
    (0,1,3,0,3,1,0,0,0,0,2,2), (3,0,1,1,0,3,0,0,0,2,0,2),
    (1,3,0,3,1,0,0,0,0,2,2,0), (2,1,3,0,1,3,0,2,2,2,2,2),
    (3,2,1,3,0,1,2,0,2,2,2,2), (1,3,2,1,3,0,2,2,0,2,2,2),
    (2,2,2,2,0,0,0,1,3,0,3,1), (2,2,2,0,2,0,3,0,1,1,0,3),
    (2,2,2,0,0,2,1,3,0,3,1,0), (2,0,0,0,0,0,2,1,3,0,1,3),
    (0,2,0,0,0,0,3,2,1,3,0,1), (0,0,2,0,0,0,1,3,2,1,3,0),
)


def _hash_json(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def _literal_source(identifier, exponents, phase, provenance):
    array = np.asarray(exponents, dtype=np.int64) % phase
    proof = {**verify_source(array, phase), 'acquisition': provenance,
             'source_construction_search_performed': False}
    digest = exponent_hash(array)
    array = array.astype(np.uint16)
    array.setflags(write=False)
    return VerifiedSource(identifier, len(array), phase, digest,
                          tuple(map(int, np.unique(array))), proof, array,
                          lambda i, j: int(array[i, j]))


@cache
def all_sources(include_expanded=False):
    if include_expanded:
        from experiments.scalar_lift.expanded_gauge_inputs import (
            admissible_sources as expanded_descriptors,
        )
        from experiments.scalar_lift.expanded_gauge_inputs import load_source as expanded_source

        combined = (*all_sources(), *(expanded_source(row['id']) for row in expanded_descriptors()))
        if len({source.id for source in combined}) != len(combined):
            raise ArithmeticError('expanded input identifiers are not distinct')
        return combined
    inputs = [registry_source(row['id']) for row in registry_descriptors()]
    inputs.extend(catalogue_source(identifier) for identifier in CATALOGUE_IDS)
    for name in ('M13B', 'M13C'):
        record = catalogue_data(name)
        inputs.append(_literal_source(f'catalogue-{name.lower()}-q6-raw', record['exponents'], 6,
                                     {'source_url': record['source_url'],
                                      'acquisition': 'previously saved literal alphabet-excluded input'}))
    inputs.append(nine_source())
    control = json.loads(CONTROL10.read_text())
    inputs.append(_literal_source('saved-modified-gs10-q4', control['literal_Butson_exponents'], 4,
                                 {'existing_certificate': str(CONTROL10.relative_to(ROOT)),
                                  'existing_certificate_sha256': hashlib.sha256(CONTROL10.read_bytes()).hexdigest()}))
    inputs.append(_literal_source('saved-williamson12-q4', CONTROL12_ROWS, 4,
                                 {'existing_note': 'docs/COMPLEX_GS_SCALAR_RESTRICTION.md',
                                  'literal_definition': 'Hermitian circulants (1,i,-i),(1,-i,i),(1,1,1),(1,-1,-1); Williamson array'}))
    from experiments.scalar_lift.alpoge_published_inputs import known_h716
    from experiments.scalar_lift.h668_width12 import known_h668_width12
    from experiments.scalar_lift.inverse_blocks import known_h668

    real_inputs = [('published-real668-width4', np.asarray(known_h668(), dtype=np.int64),
                    {'existing_seed': 'experiments/scalar_lift/h668_seed.json'}),
                   ('published-real668-width12', *known_h668_width12()),
                   ('published-real716-width4', *known_h716())]
    for identifier, matrix, provenance in real_inputs:
        inputs.append(_literal_source(identifier, (1-np.asarray(matrix, dtype=np.int64))//2, 2, provenance))
    from experiments.scalar_lift.broad_source_acquisition import load_source as broad_source
    from experiments.scalar_lift.broad_source_acquisition import sources as broad_descriptors

    inputs.extend(broad_source(row['id']) for row in broad_descriptors())
    identifiers = [source.id for source in inputs]
    if len(identifiers) != len(set(identifiers)):
        raise ArithmeticError('a starting input identifier was repeated')
    return tuple(inputs)


def field_from_alphabet(phase, alphabet):
    generator_gcd = math.gcd(phase, *map(int, alphabet))
    root_order = phase//generator_gcd
    return {'actual_root_generator_order': root_order,
            'exponent_gcd_with_phase': generator_gcd,
            'effective_root_exponents': [int(e)//generator_gcd for e in alphabet],
            'entry_field_degree': int(totient(root_order)),
            'field_definition': f'Q(zeta_{root_order})'}


def source_global_phase_hash(source, common_phase=240):
    """Exact literal global-phase equivalence, not row/column equivalence."""
    if common_phase % source.phase:
        raise ValueError('the common root label does not contain this source alphabet')
    if source.exponents is None:
        return _hash_json({'implicit_source_hash': source.source_hash, 'order': source.order})
    dtype = np.uint16 if common_phase < 32768 else np.int64
    values = ((source.exponents.astype(dtype)+source.phase-source.entry(0, 0))
              *(common_phase//source.phase)) % common_phase
    return hashlib.sha256(str(source.order).encode()+b':'+str(common_phase).encode()
                          +values.astype('<i8').tobytes()).hexdigest()


@cache
def _alphabet_gate(phase, alphabet):
    return alphabet_obstruction(phase, alphabet)


def anchor_cases(source):
    """Exact distinct alphabets from every used-root anchor, retaining recipes.

    Equal alphabets share the same necessary tests, not a literal matrix
    identity. No Galois or arbitrary row/column gauge equivalence is assumed.
    """
    raw = tuple(sorted(source.used_alphabet))
    palettes = {raw: {'contains_raw_representative': True, 'used_entry_anchors': []}}
    for exponent in raw:
        palette = tuple(sorted((u-exponent) % source.phase for u in raw))
        record = palettes.setdefault(palette, {'contains_raw_representative': False,
                                              'used_entry_anchors': []})
        record['used_entry_anchors'].append(exponent)
    result = []
    for index, (alphabet, record) in enumerate(palettes.items()):
        recipes = [{'anchor_exponent': anchor, 'global_exponent_shift': -anchor,
                    'source_recipe_sha256': _hash_json({'source_hash': source.source_hash,
                                                       'phase': source.phase, 'global_exponent_shift': -anchor}),
                    'hash_kind': 'exact source plus global phase recipe, not a dense matrix hash'}
                   for anchor in record['used_entry_anchors']]
        result.append({'id': f'{source.id}:alphabet-{index}', **record,
                       'used_alphabet': list(alphabet),
                       **field_from_alphabet(source.phase, alphabet),
                       'anchor_gauge_recipes': recipes,
                       'source_Gram_preserved_by_unit_global_scalar': True,
                       'alphabet_gate': _alphabet_gate(source.phase, alphabet),
                       'literal_matrix_nodes_deduplicated': False})
    return result


def _route(source, case, target):
    width = target//source.order
    degree = case['entry_field_degree']
    errors = []
    if width % degree:
        errors.append({'kind': 'entry_field_degree', 'required_degree_divisor': degree})
    if width not in (1, 2) and width % 4:
        errors.append({'kind': 'sign_block_width', 'required': 'd=1, d=2, or d divisible by4'})
    if case['alphabet_gate']['excluded']:
        errors.append({'kind': 'odd_integer_scalar_relation', 'alphabet_case': case['id']})
    standard = representation_obstructions(width, case['actual_root_generator_order'])
    for failure in standard:
        if failure['kind'] != 'field_degree':
            errors.append({'kind': 'standard_metric_scalar_representation', 'details': failure})
    return {'source_id': source.id, 'source_hash': source.source_hash,
            'source_order': source.order, 'declared_phase': source.phase,
            'alphabet_case': case['id'], 'actual_root_generator_order': case['actual_root_generator_order'],
            'entry_field_degree': degree, 'target_order': target, 'required_total_width': width,
            'effective_root_exponents': case['effective_root_exponents'],
            'necessary_tests_pass': not errors, 'exclusions': errors,
            'actual_basis_or_sign_endpoint_constructed': False}


def route_records(include_excluded=False, source_ids=None, *, include_expanded=False):
    """Return all arithmetic routes, or just necessary-test survivors by default."""
    targets = json.loads(FRONTIER.read_text())['retained_source_listed_unresolved_orders']
    selected = [s for s in all_sources(include_expanded) if source_ids is None or s.id in source_ids]
    rows = [_route(source, case, target) for source in selected for case in anchor_cases(source)
            for target in targets if target % source.order == 0]
    return rows if include_excluded else [row for row in rows if row['necessary_tests_pass']]


def report(include_expanded=False):
    targets = json.loads(FRONTIER.read_text())['retained_source_listed_unresolved_orders']
    sources, records = [], []
    selected_sources = all_sources(include_expanded)
    common_phase = math.lcm(*(source.phase for source in selected_sources))
    for source in selected_sources:
        if source.proof.get('full_Gram_verified') is not True:
            raise ArithmeticError('a starting input has no completed source Gram proof')
        cases = anchor_cases(source)
        sources.append({'id': source.id, 'order': source.order, 'declared_phase': source.phase,
                        'source_hash': source.source_hash,
                        'source_proof_sha256': _hash_json(source.proof),
                        'global_phase_literal_class_sha256': source_global_phase_hash(source, common_phase),
                        'source_Gram_replayed': True, 'dense_exponents_available': source.exponents is not None,
                        'actual_raw_field': field_from_alphabet(source.phase, source.used_alphabet),
                        'anchor_alphabet_cases': cases})
        records.extend(_route(source, case, target) for case in cases for target in targets
                       if target % source.order == 0)
    retained = [r for r in records if r['necessary_tests_pass']]
    retained_ids = {r['source_id'] for r in retained}
    classes = {}
    for source in sources:
        classes.setdefault(source['global_phase_literal_class_sha256'], []).append(source['id'])
    summary = {}
    for target in sorted({r['target_order'] for r in retained}):
        rows = [r for r in retained if r['target_order'] == target]
        summary[str(target)] = {'source_ids': sorted({r['source_id'] for r in rows}),
                                'declared_phases': sorted({r['declared_phase'] for r in rows}),
                                'actual_entry_field_degrees': sorted({r['entry_field_degree'] for r in rows}),
                                'anchor_cases': len(rows),
                                'required_widths': sorted({r['required_total_width'] for r in rows})}
    result = {'schema': 'phase-independent-frontier-preflight-v1',
            'model_scope': 'uniform scalar input/output bases, with every used scalar block having Gram dI; nonuniform source-coordinate maps are outside these exclusions',
            'frontier_file': str(FRONTIER.relative_to(ROOT)),
            'frontier_sha256': hashlib.sha256(FRONTIER.read_bytes()).hexdigest(),
            'frontier_is_attributed_bounded_snapshot': True,
            'sources': sources, 'source_count': len(sources), 'target_count': len(targets),
            'common_root_label_for_literal_global_phase_comparison': common_phase,
            'literal_global_phase_input_class_count': len(classes),
            'literal_global_phase_duplicate_input_classes': [ids for ids in classes.values() if len(ids)>1],
            'source_target_divisibility_checks': len(sources)*len(targets),
            'declared_phases_checked': sorted({s['declared_phase'] for s in sources}),
            'actual_raw_root_orders_checked': sorted({s['actual_raw_field']['actual_root_generator_order'] for s in sources}),
            'distinct_anchor_alphabet_cases': sum(len(s['anchor_alphabet_cases']) for s in sources),
            'used_entry_anchor_recipes_checked': sum(len(c['used_entry_anchors']) for s in sources for c in s['anchor_alphabet_cases']),
            'divisible_source_target_pairs': sum(N % s['order'] == 0 for s in sources for N in targets),
            'anchor_target_routes_checked': len(records), 'routes': records,
            'retained_anchor_target_routes': len(retained),
            'retained_source_binding_count': len(retained_ids),
            'retained_global_phase_literal_source_class_count': len({
                s['global_phase_literal_class_sha256'] for s in sources if s['id'] in retained_ids}),
            'retained_target_order_count': len(summary),
            'retained_routes': retained, 'retained_target_summary': summary,
            'exclusion_reason_counts': dict(Counter(e['kind'] for r in records for e in r['exclusions'])),
            'retained_target_orders_above20000': sorted({r['target_order'] for r in retained if r['target_order'] > 20000}),
            'raw_and_all_used_entry_global_anchor_gauges_covered': True,
            'arbitrary_row_column_gauges_or_nonuniform_maps_covered': False,
            'Galois_equivalence_used_for_pruning': False,
            'only_deduplication': 'identical scalar alphabets share necessary tests; literal matrices are not merged',
            'remaining_work': 'construct actual mixed input/output bases for retained routes and verify every sign entry and expanded Gram',
            'new_Hadamard_found': False, 'source_construction_or_SAT_search_run': False,
            'all_declared_preflight_cases_finished': True}
    if include_expanded:
        expansion = ROOT/'certificates/scalar-lift/expanded-gauge-inputs.json'
        result['optional_expanded_gauge_inputs'] = {
            'base_source_count': len(all_sources()),
            'additional_fixed_representatives': len(selected_sources)-len(all_sources()),
            'provider': 'experiments/scalar_lift/expanded_gauge_inputs.py',
            'acquisition_certificate': str(expansion.relative_to(ROOT)),
            'acquisition_certificate_sha256': hashlib.sha256(expansion.read_bytes()).hexdigest(),
            'original_default_report_changed': False}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--output', type=Path)
    group.add_argument('--check', type=Path)
    parser.add_argument('--expanded-gauges', action='store_true',
                        help='add the separately acquired higher-root fixed representatives; default stays the original555 inputs')
    args = parser.parse_args()
    result = report(args.expanded_gauges)
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    elif json.loads(args.check.read_text()) != result:
        raise ValueError('phase-independent frontier certificate differs from exact replay')
    print(json.dumps({k: result[k] for k in ('source_count', 'target_count', 'declared_phases_checked',
                                           'distinct_anchor_alphabet_cases', 'anchor_target_routes_checked',
                                           'retained_target_orders_above20000', 'new_Hadamard_found')}))


if __name__ == '__main__':
    main()

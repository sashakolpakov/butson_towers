"""Bind the completed broad input, matrix-tree, and basis audits together.

This summarizes existing exact run reports; it does not launch a search or
merge matrix-equivalence classes across independently saved trees.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DIRECTORY = ROOT/'certificates/scalar-lift'


def report():
    references = {}

    def read(name):
        data = (DIRECTORY/f'{name}.json').read_bytes()
        references[name] = {'path': f'certificates/scalar-lift/{name}.json',
                            'sha256': hashlib.sha256(data).hexdigest()}
        return json.loads(data)

    acquired = read('broad-source-acquisition')
    inputs = {'published_matrices_exactly_verified': acquired['source_count'],
              'own_alphabet_admissible_representatives': acquired['retained_fixed_representative_count']}
    expanded = read('expanded-gauge-inputs')
    inputs.update({'q6_inputs_checked_in_q12': expanded['original_published_input_count'],
                   'newly_rescued_representatives': expanded['rescued_input_count'],
                   'q12_gauged_representatives': expanded['total_admissible_twelfth_root_representatives'],
                   'new_projective_Butson_classes_claimed': False})
    frontier = read('phase-independent-expanded-frontier')
    preflight = {key: frontier[key] for key in
                 ('source_count', 'target_count', 'declared_phases_checked',
                  'source_target_divisibility_checks', 'distinct_anchor_alphabet_cases',
                  'anchor_target_routes_checked', 'retained_target_orders_above20000')}
    preflight['necessary_test_survivors'] = len(frontier['retained_routes'])
    preflight['retained_target_orders'] = sorted(map(int, frontier['retained_target_summary']))
    preflight['survivors_are_not_sign_matrix_constructions'] = True

    trees = []
    for name in ('broad-root-flat-all-inputs-tree', 'expanded-gauge-root-flat-tree'):
        tree = read(name)
        assert tree['pending_recipes'] == 0
        assert len(tree['nodes']) == tree['progress']['completed_matrices']
        trees.append({'report': name+'.json', 'fixed_input_bindings': len(tree['fixed_existing_seeds']),
                      'completed_literal_nodes': len(tree['nodes']),
                      'duplicates_pruned': tree['progress']['duplicate_recipes'],
                      'pending': tree['pending_recipes'], 'storage_bytes': tree['storage_bytes'],
                      'real_orders': tree['real_output_orders']})
    partial = read('partial-fifth-descent')
    assert partial['pending_nodes'] == 0
    trees.append({'report': 'partial-fifth-descent.json', 'fixed_input_bindings': partial['source_count'],
                  'completed_literal_nodes': partial['literal_node_count'],
                  'duplicates_pruned': partial['exact_literal_duplicates_pruned'],
                  'pending': partial['pending_nodes'], 'storage_bytes': partial['literal_storage_bytes'],
                  'real_orders': partial['real_orders'], 'real_nodes': partial['real_literal_nodes'],
                  'partial_third_edges': partial['partial_third_edges'],
                  'partial_fifth_edges': partial['partial_fifth_edges']})
    assert partial['all_terminal_sign_Grams_independently_checked']

    basis = read('broad-phase-basis-probe')
    assert basis['pending_declared_geometries'] == 0
    assert basis['preflight_file_sha256'] == hashlib.sha256(
        (ROOT/basis['preflight_file']).read_bytes()).hexdigest()
    basis_summary = {key: basis[key] for key in
                     ('candidate_geometries', 'source_matrix_evaluations', 'palette_evaluations',
                      'covered_preflight_routes', 'covered_target_orders', 'new_sign_Hadamard_found')}
    gate = read('normalization-denominator-gate')
    odd = read('odd-denominator-basis-audit')
    assert odd['pending_declared_audits'] == 0
    return {'schema': 'broad-scalar-restriction-status-v1', 'inputs': inputs,
            'preflight': preflight, 'completed_matrix_trees': trees,
            'tree_counts_are_per_run_not_global_equivalence_classes': True,
            'basis_experiments': basis_summary,
            'basis_family_theorems': {
                'dyadic_gate_report_schema': gate['schema'],
                'dyadic_families_excluded_for_odd_ratio_at_width4_mod8_above4': True,
                'odd_denominator_family': odd['basis_family'],
                'odd_denominator_frame_families_audited': odd['frame_count'],
                'odd_denominator_families_closed': odd['families_closed'],
                'odd_denominator_open_frames': odd['open_frames'],
                'all_rational_bases_excluded': False},
            'all_declared_runs_finished': True, 'new_unknown_Hadamard_order_found': False,
            'only_matrix_tree_child_operations': 'scalar restriction in declared input/output bases and exact literal duplicate pruning',
            'references': references}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    choice = parser.add_mutually_exclusive_group(required=True)
    choice.add_argument('--output', type=Path)
    choice.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = report()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    elif json.loads(args.check.read_text()) != result:
        raise ValueError('the broad status snapshot differs from its referenced reports')
    print(json.dumps({'verified_input_bindings': result['preflight']['source_count'],
                      'newly_rescued_representatives': result['inputs']['newly_rescued_representatives'],
                      'completed_matrix_nodes_per_tree': [tree['completed_literal_nodes']
                                                          for tree in result['completed_matrix_trees']],
                      'basis_experiments': result['basis_experiments']['source_matrix_evaluations'],
                      'new_unknown_Hadamard_order_found': False}))


if __name__ == '__main__':
    main()

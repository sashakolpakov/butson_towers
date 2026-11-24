import json
from pathlib import Path

import numpy as np
import pytest
from sympy import Matrix, diag

from experiments.scalar_lift.ambient_field_candidates import (
    SavedPeriodBundle,
    family,
    load_ambient_source,
    period_field,
)
from experiments.scalar_lift.pure_period_restriction import period_towers, tower_report


def test_literal_two_towers_and_complete_weighted_source_Gram_without_frame():
    field = period_field(3, 31, (2, 3, 1))
    source = load_ambient_source('unreal3-q6')
    first, second = period_towers(field, source)
    assert [s.order for s in first] == [9, 18, 36]
    assert [s.order for s in second] == [6, 12, 36]
    for path in (first, second):
        stage = path[-1]
        def block(e, current=stage):
            return Matrix(current.scalar_table(e)[:, :, 0].tolist())
        source_matrix = Matrix.vstack(*(Matrix.hstack(*(block(source.entry(i, j)) for j in range(3)))
                                        for i in range(3)))
        G = Matrix(stage.metric_table()[:, :, 0].tolist())
        weight = diag(G.inv(), G.inv(), G.inv())
        assert source_matrix*weight*source_matrix.T == 3*weight
        assert source_matrix*source_matrix.T != 36*Matrix.eye(36)
        assert 0 in source_matrix
        for i, j in ((0, 0), (3, 9), (17, 35), (35, 1)):
            value = stage.entry(i, j)
            assert value['denominator'] == 1
            assert value['numerators'] == [source_matrix[i, j]]
    assert first[-1].literal_hash() != second[-1].literal_hash()
    report = tower_report(field, source)
    assert report['source_Gram_scale_unchanged'] == 3
    assert report['product_basis_sign_endpoint_excluded']
    assert not report['nonproduct_joint_bases_or_unequal_input_output_frames_excluded']


def test_acquired_twelfth_source_and_preexisting_nine_have_the_actual_large_paths():
    for r, source_id, target in ((167, 'catalogue-p7-q12', 4676),
                                 (191, 'saved-unreal9-q6', 6876)):
        field, source = period_field(r), load_ambient_source(source_id)
        report = tower_report(field, source)
        assert report['paths'][0][-1]['order'] == target
        assert report['paths'][1][-1]['order'] == target
        assert report['literal_state_count'] == report['literal_distinct_state_count'] == 7
        for path in report['paths']:
            assert all(not stage['right_frame_inserted'] for stage in path)
            assert all(stage['Gram_scale'] == source.order for stage in path)
        assert report['initial_source']['coefficient_field_degree'] == 4*r


def test_bundle_replays_pure_stages_and_rejects_array_corruption(tmp_path):
    path = tmp_path/'source7.npz'
    report = family(167, 'catalogue-p7-q6', path)
    saved = SavedPeriodBundle(path)
    assert saved.order == 4676
    assert saved.report == report
    assert saved.entry(0, 4)['numerators'] == [0]
    with pytest.raises(ValueError, match='indices'):
        saved.entry(1.5, 0)
    with np.load(path, allow_pickle=False) as data:
        arrays = {key: np.array(data[key]) for key in data.files}
    arrays['period_multiplication'][0, 0] += 1
    np.savez_compressed(path, **arrays)
    with pytest.raises(ValueError, match='arrays'):
        SavedPeriodBundle(path)


def test_executed_finite_report_contains_eight_targets_and_no_pending_product_search():
    path = 'certificates/scalar-lift/ambient-field-candidates.json'
    report = json.loads(Path(path).read_text())
    summary = report['scalar_tree_execution']
    assert summary['recipe_count'] == 12
    assert summary['retained_target_orders'] == [4676, 6876, 7924, 8684, 9308, 9716, 10052, 11596]
    assert summary['pending_declared_recipes'] == 0
    assert summary['product_basis_family_excluded_by_structural_zeros']

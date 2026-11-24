import json
import hashlib
import itertools
import shutil
import subprocess
import sys

import numpy as np

import pytest

from experiments.scalar_lift import hadamard_frame_search as search


def test_complete_parity_polynomial_identity_and_a_corrupted_map(monkeypatch):
    report = search.validate_parity_map()
    assert report['color_map_rank'] == 60
    assert report['determinant_map_rank'] == 15
    assert report['full_quadratic_coefficients_checked'] == 1831
    columns = list(search.PARITY_INVERSE_COLUMNS_631)
    columns[0], columns[1] = columns[1], columns[0]
    monkeypatch.setattr(search, 'PARITY_INVERSE_COLUMNS_631', tuple(columns))
    with pytest.raises(ValueError):
        search.validate_parity_map()


def test_complete_native_small_control_and_full_target_filter_set(tmp_path):
    if not any(shutil.which(name) for name in ('c++', 'clang++', 'g++')):
        pytest.skip('C++17 compiler unavailable')
    binary = search.build_native(tmp_path)
    report = search.self_test(tmp_path, binary)
    assert report['complete_T19_axis_patterns'] == 90
    assert report['complete_T19_solution_set_size'] == 150
    assert report['small_endpoint_integer_Gram_verified']
    assert report['target_affine_masks_checked'] == 8192
    assert report['target_complete_mod8_survivor_set_size'] == 1
    assert report['planted_T631_retained']
    assert report['independent_prime_repeat_endpoint_order'] == 3252
    assert report['independent_prime_repeat_full_integer_Gram_verified']


def test_larger_unresolved_target_does_not_silently_use_small_native_words(tmp_path):
    data = search.chart(223, 3)
    assert data['n'] == 75
    with pytest.raises(ValueError, match='at most63'):
        search.write_tables(tmp_path / 'oversized.bin', data, -5, 3)


def test_planted_T547_replays_complete_H6564_without_stochastic_search(tmp_path):
    record = json.loads((search.ROOT / 'certificates/scalar-lift/hadamard-frame-p547.json').read_text())
    audit = search.validate_parity_map(config=record)
    assert audit['full_quadratic_coefficients_checked'] == 1379
    endpoint = record['known_order_positive_repeat']['endpoint']
    data = search.chart(record['p'], record['m'], record['mu'])
    assert search.expand(endpoint['states'], data).sum(axis=1).tolist() == [-20, 7, 7, 7]
    replay = search.verify_endpoint(endpoint['states'], data, tmp_path / 'verified')
    assert replay['endpoint_order'] == 6564
    assert replay['source_gram_scale'] == 3
    assert replay['frame_gram_scale'] == 2188
    assert replay['module_rank'] == 1094
    assert replay['all_coefficient_pair_cross_and_endpoint_integer_Grams_verified']
    assert replay['packed_sha256'] == endpoint['packed_sha256']
    assert replay['actual_R_entry_numerators'] == endpoint['actual_R_entry_numerators']
    assert not record['known_order_positive_repeat']['new_order_claimed']


@pytest.mark.parametrize('prime', [307, 919])
def test_degree17_sampler_and_complete_high_word_filter(prime, tmp_path):
    if not any(shutil.which(name) for name in ('c++', 'clang++', 'g++')):
        pytest.skip('C++17 compiler unavailable')
    config = json.loads((search.ROOT / f'certificates/scalar-lift/hadamard-frame-p{prime}.json').read_text())
    assert search.validate_parity_map(config=config)['full_quadratic_coefficients_checked'] == 2347
    data = search.chart(prime, config['m'], config['mu'])
    binary = search.build_native(tmp_path)
    table, columns, output = tmp_path / 'table.bin', tmp_path / 'columns.bin', tmp_path / 'samples.json'
    search.write_tables(table, data, config['head_net'], config['tail_net'])
    search.write_parity_map(columns, config)
    subprocess.run([str(binary), 'search', str(table), str(columns), str(output),
                    '0', '2026112601', '256'], check=True, stdout=subprocess.DEVNULL)
    samples = json.loads(output.read_text())['parameter_examples']
    chosen = None
    assert any(row['color_bits'] >> 64 for row in samples)
    for row in samples:
        p, q, u, v = [row[key] for key in ('p_field', 'q_field', 'b0', 'b1')]
        assert search.field_multiply(p, v, 17, config['field_modulus']) ^ search.field_multiply(
            q, u, 17, config['field_modulus']) == 1
        spectral = p | q << 17 | u << 34 | v << 51
        assert spectral == row['spectral']
        colors = search.binary_linear(config['inverse_columns'], spectral)
        assert colors == row['color_bits']
        axes = [0] + [(colors >> (2*j)) & 3 for j in range(34)]
        if chosen is None and axes[1] == 0:
            try:
                masks, correlations = search.affine_signs(
                    axes, data, config['head_net'], config['tail_net'], 1)
            except ValueError:
                continue
            chosen = axes, masks, correlations
    assert chosen is not None
    axes, masks, correlations = chosen
    assert len(masks) == 32768
    expected = {value for value, row in zip(masks, correlations) if np.all(row % 8 == 0)}
    pool, result = tmp_path / 'pool.bin', tmp_path / 'lift.json'
    search.write_pool(pool, data, [axes], config['head_net'], config['tail_net'], 1)
    subprocess.run([str(binary), 'lift', str(table), str(pool), str(result),
                    '0', '1', '1', '1'], check=True, stdout=subprocess.DEVNULL)
    actual = json.loads(result.read_text())['patterns'][0]['mod8_negative_masks']
    assert set(actual) == expected


def test_completed_T127_and_full_scalar_endpoint(tmp_path):
    record = json.loads((search.ROOT / 'certificates/scalar-lift/hadamard-frame-p127.json').read_text())
    assert search.validate_parity_map(config=record)['full_quadratic_coefficients_checked'] == 407
    completed = record['completed_T_matrices']
    data = search.chart(127, 9, 22)
    b = search.expand(completed['states'], data)
    assert np.array_equal(b, completed['b_sequences'])
    indices = np.arange(127)
    matrices = [row[(indices[None, :] - indices[:, None]) % 127] for row in b]
    assert np.all(sum(np.abs(t) for t in matrices) == 1)
    assert np.array_equal(sum(t @ t.T for t in matrices), 127 * np.eye(127, dtype=np.int64))
    replay = search.verify_endpoint(completed['states'], data, tmp_path / 'verified')
    assert replay['endpoint_order'] == 1524
    assert replay['packed_sha256'] == record['known_order_positive_repeat']['endpoint']['packed_sha256']


def test_degree22_sampler_inverse_and_million_mask_filter(tmp_path):
    config = json.loads((search.ROOT / 'certificates/scalar-lift/hadamard-frame-p397.json').read_text())
    assert search.validate_parity_map(config=config)['full_quadratic_coefficients_checked'] == 3917
    binary = search.build_native(tmp_path)
    data = search.chart(397, 9, 79)
    table, columns, output = tmp_path / 'table.bin', tmp_path / 'columns.bin', tmp_path / 'samples.json'
    search.write_tables(table, data, -2, 2)
    search.write_parity_map(columns, config)
    subprocess.run([str(binary), 'search', str(table), str(columns), str(output),
                    '0', '20261126397', '256'], check=True, stdout=subprocess.DEVNULL)
    samples = json.loads(output.read_text())['parameter_examples']
    for row in samples:
        p, q, u, v = [row[k] for k in ('p_field', 'q_field', 'b0', 'b1')]
        assert search.field_multiply(p, v, 22, config['field_modulus']) ^ search.field_multiply(
            q, u, 22, config['field_modulus']) == 1
        assert row['color_bits'] == search.binary_linear(config['inverse_columns'], row['spectral'])
    colors = samples[36]['color_bits']
    axes = [0] + [(colors >> (2*j)) & 3 for j in range(44)]
    system = search.affine_sign_system(axes, data, -2, 2, 1)
    particular, basis = system['particular'], system['basis']
    dimension = len(basis)
    assert dimension == 20

    def parity_word(z):
        negative = particular
        for i, vector in enumerate(basis):
            if z >> i & 1:
                negative ^= vector
        b = search.expand([2*a + ((negative >> i) & 1) for i, a in enumerate(axes)], data)
        correlations = [int(np.sum(b * np.roll(b, -s, axis=1))) for s in data['shifts']]
        assert all(c % 4 == 0 for c in correlations)
        return sum(((c // 4) & 1) << i for i, c in enumerate(correlations))

    constant = parity_word(0)
    linear = [parity_word(1 << i) ^ constant for i in range(dimension)]
    z = np.arange(1 << dimension, dtype=np.uint32)
    values = np.full(len(z), constant, dtype=np.uint32)
    for i in range(dimension):
        values ^= ((z >> i) & 1) * np.uint32(linear[i])
        for j in range(i):
            coefficient = parity_word((1 << i) | (1 << j)) ^ linear[i] ^ linear[j] ^ constant
            values ^= ((z >> i) & (z >> j) & 1) * np.uint32(coefficient)
    expected = set()
    for z0 in np.flatnonzero(values == 0).tolist():
        negative = particular
        for i, vector in enumerate(basis):
            if z0 >> i & 1:
                negative ^= vector
        expected.add(negative)
    pool, result = tmp_path / 'pool.bin', tmp_path / 'lift.json'
    search.write_pool(pool, data, [axes], -2, 2, 1)
    subprocess.run([str(binary), 'lift', str(table), str(pool), str(result),
                    '0', '1', '1', '1'], check=True, stdout=subprocess.DEVNULL)
    assert set(json.loads(result.read_text())['patterns'][0]['mod8_negative_masks']) == expected
    assert len(expected) == 2


def test_headless_small_conditional_class_is_complete_and_positive(tmp_path):
    config = json.loads((search.ROOT / 'certificates/scalar-lift/hadamard-frame-p13-headless.json').read_text())
    binary = search.build_native(tmp_path)
    data = search.chart(13, 3, config['mu'])
    table, columns, result = tmp_path / 'table.bin', tmp_path / 'columns.bin', tmp_path / 'result.json'
    search.write_tables(table, data, 0, 2)
    search.write_parity_map(columns, config)
    subprocess.run([str(binary), 'headless', str(table), str(columns), str(result),
                    '0', '0', 'exhaustive'], check=True, stdout=subprocess.DEVNULL)
    native = json.loads(result.read_text())
    expected_axes, expected_solutions = set(), set()
    for rest in itertools.product((1, 2, 3), repeat=3):
        axes = [0, 1, *rest]
        unsigned = search.expand([2*a for a in axes], data)
        if all(int(np.sum(unsigned * np.roll(unsigned, -s, axis=1))) % 2 == 0 for s in range(1, 13)):
            expected_axes.add(tuple(axes))
        for position in (2, 3, 4):
            states = tuple(2*a + int(i == position) for i, a in enumerate(axes))
            b = search.expand(list(states), data)
            if all(np.sum(b * np.roll(b, -s, axis=1)) == 0 for s in range(1, 13)):
                expected_solutions.add(states)
    axes_found, solutions = set(), set()
    for row in native['records']:
        axes = [0] + [(row['colors'] >> (2*j)) & 3 for j in range(4)]
        axes_found.add(tuple(axes))
        solutions.update(tuple(2*a + ((negative >> i) & 1) for i, a in enumerate(axes))
                         for negative in row['negative_masks'])
    assert axes_found == expected_axes and len(axes_found) == len(native['records']) == 12
    assert solutions == expected_solutions and len(solutions) == 18
    assert search.verify_endpoint(list(sorted(solutions)[0]), data, tmp_path / 'verified')['endpoint_order'] == 156


def test_headless_large_conditional_filter_checks_every_w(tmp_path):
    config = json.loads((search.ROOT / 'certificates/scalar-lift/hadamard-frame-p1453.json').read_text())
    binary = search.build_native(tmp_path)
    data = search.chart(1453, 33, config['mu'])
    table, columns, result = tmp_path / 'table.bin', tmp_path / 'columns.bin', tmp_path / 'result.json'
    search.write_tables(table, data, 0, 2)
    search.write_parity_map(columns, config)
    subprocess.run([str(binary), 'headless', str(table), str(columns), str(result),
                    '0', '0', 'checkpair', '4112976', '4149651'], check=True, stdout=subprocess.DEVNULL)
    record = json.loads(result.read_text())
    assert record['all_w_checked'] == 4194304
    assert record['complete_w_sets_equal'] and record['headless_w_count'] == 9


def test_run_binds_native_source_before_build(tmp_path, monkeypatch):
    source = tmp_path / 'kernel.cpp'
    source.write_text('before compilation')
    expected = hashlib.sha256(source.read_bytes()).hexdigest()
    monkeypatch.setattr(search, 'NATIVE_SOURCE', source)

    def build(workdir):
        source.write_text('changed after source capture')
        return tmp_path / 'binary'

    def run(command, **kwargs):
        (tmp_path / 'seed-1.json').write_text('{}')

    monkeypatch.setattr(search, 'build_native', build)
    monkeypatch.setattr(search.subprocess, 'run', run)
    monkeypatch.setattr(search, 'validate_parity_map', lambda *a, **k: {})
    monkeypatch.setattr(search, 'write_tables', lambda *a, **k: None)
    monkeypatch.setattr(search, 'write_parity_map', lambda *a, **k: None)
    monkeypatch.setattr(sys, 'argv', ['frame-search', '--workdir', str(tmp_path)])
    search.main()
    report = json.loads((tmp_path / 'seed-1-run.json').read_text())
    assert report['native_source_sha256'] == expected
    assert report['source_capture_phase'] == 'before_native_build_and_search'

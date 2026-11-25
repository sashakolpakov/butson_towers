"""Exact controls for head-count-changing local moves."""
import json
from pathlib import Path
import shutil
import subprocess

import pytest

from experiments.scalar_lift import hadamard_frame_search as frame
from experiments.scalar_lift import t_sequence_local_search as local


def test_balanced_head_tail_conversions_and_planted_checkpoint(tmp_path):
    if shutil.which('c++') is None:
        pytest.skip('C++17 compiler unavailable')
    data = frame.chart(31, 3, 5)
    table = tmp_path/'coefficients.bin'
    local.write_tables(table, data, -1, 3)
    binary = local.build(tmp_path)
    output = tmp_path/'checkpoint.json'
    checked = subprocess.run(
        [str(binary), str(table), str(output), '1', '223', '1', 'check',
         '3', '5', '-', '1'], check=True, capture_output=True, text=True)
    records = [json.loads(line) for line in checked.stdout.splitlines()]
    assert records[0] == {'incremental_random_move_checks': 500, 'passed': True}
    assert records[1] == {'balanced_conversion_checks': 300, 'increases': 150,
                          'decreases': 150, 'exact_DC_preserved': True}
    states = [0, 1, 2, 2, 7, 6, 2, 6, 2, 5, 3]
    initial = tmp_path/'initial.txt'
    initial.write_text(' '.join(map(str, states))+'\n')
    subprocess.run([str(binary), str(table), str(output), '0.01', '223', '1',
                    'anneal', '3', '5', str(initial), '1'], check=True,
                   capture_output=True, text=True)
    record = json.loads(output.read_text())
    assert record['free_heads']
    assert record['head_orbits'] == 1
    assert record['best_energy'] == record['search_cost'] == 0
    proof = local.verify_checkpoint(record, data, -1, 3)
    assert proof['T_witness_found'] and proof['all_expanded_correlations_checked']

"""The campaign retains a witness only after exact expanded verification."""
import json
import shutil

import pytest

from experiments.scalar_lift.low_order_campaign import run_projection


def test_seeded_small_projection_campaign_writes_exact_witness(tmp_path):
    if not (shutil.which('c++') or shutil.which('clang++')):
        pytest.skip('native assignment requires C++17')
    states = [0,1,2,2,7,6,2,6,2,5,3]
    seed = tmp_path/'initial.json'
    seed.write_text(json.dumps({'states':states}))
    output = tmp_path/'campaign'
    report = run_projection({'p':31,'m':3,'mu':5,'head_net':-1,'tail_net':3},
                            output, seconds=2, trial_seconds=1, initial=seed)
    assert report['finished'] and report['exact_T_design']
    witness = json.loads((output/'verified-witness.json').read_text())
    assert witness['best_energy'] == 0
    assert witness['all_expanded_shift_correlations_checked']
    assert witness['maximum_absolute_nonzero_shift'] == 0
    assert witness['row_sums'] == [-2,3,3,3]
    assert len(report['trials']) == 1

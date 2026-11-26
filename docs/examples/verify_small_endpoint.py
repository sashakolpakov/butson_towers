"""Reconstruct a complete known-order H84 control from a literal T7 input."""
from pathlib import Path
import json
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from experiments.scalar_lift.hadamard_frame_search import expand
from experiments.scalar_lift.unrestricted_parity_map import periodic_chart
from experiments.scalar_lift.unrestricted_t_sign_model import verify_states
from experiments.scalar_lift.williamson_scalar_endpoints import (
    ScalarEndpoint, verify_orthogonal_design, verified_williamson,
)

STATES = [0, 6, 0, 4, 4, 5, 2]
TARGET_SUMS = (2, 1, 1, 1)


def verify_control():
    """Prove the input/source identities and independently check the dense H84."""
    data = periodic_chart(7)
    t_proof = verify_states(STATES, data, TARGET_SUMS)
    rows = expand(STATES, data)
    coefficients, od_proof = verify_orthogonal_design(rows)
    source, source_proof = verified_williamson(3)
    endpoint = ScalarEndpoint(coefficients, source)
    dense = np.stack([endpoint.row(i) for i in range(endpoint.order)]).astype(np.int64)
    signs = bool(np.all(np.abs(dense) == 1))
    gram = bool(np.array_equal(dense @ dense.T,
                              endpoint.order*np.eye(endpoint.order, dtype=np.int64)))
    if not signs or not gram:
        raise ArithmeticError('independent dense endpoint check failed')
    return {
        'scope': 'known-order control, not a new-order claim',
        'T_length': 7,
        'T_row_sums': t_proof['sequence_row_sums'],
        'T_integer_shifts_checked': t_proof['shift_count'],
        'orthogonal_design_order': od_proof['orthogonal_design_order'],
        'orthogonal_design_proved': od_proof['full_polynomial_Gram_verified'],
        'source_width': source_proof['order'],
        'order': endpoint.order,
        'all_entries_are_signs': signs,
        'full_integer_Gram_verified': gram,
    }


if __name__ == '__main__':
    print(json.dumps(verify_control(), indent=2, sort_keys=True))

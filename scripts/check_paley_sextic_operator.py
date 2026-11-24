"""Replay the dense Paley sextic operator and its exact no-frame proof."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from experiments.scalar_lift.paley_sextic_operator import build_paley_sextic_operator_report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = build_paley_sextic_operator_report()
    if args.check is not None and json.loads(args.check.read_text()) != report:
        raise SystemExit("saved operator report differs from exact replay")
    if args.output is not None:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(
        json.dumps(
            {
                "order": report["order"],
                "exact_integer_operator_squares_verified": True,
                "excluded_orientation_flat_vector_count": report["flat_vector_count"],
                "excluded_orientation_flat_span_dimension": report["dimension_of_flat_vector_span"],
                "opposite_chirality_operator": report["compatible_opposite_chirality_candidate"],
                "required_regular_block_sums_M": report["compatible_candidate_regular_block_sum_audit"][
                    "required_block_sum_matrix_M"
                ],
                "new_real_hadamard_order_claimed": report["new_real_hadamard_order_claimed"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

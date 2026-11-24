"""Replay the complete order-668 tensor-operator exclusion and its certificates."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from experiments.scalar_lift.paley_frame_family import build_paley_frame_family_report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--check", type=Path)
    group.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = build_paley_frame_family_report()
    encoded = json.dumps(report, indent=2) + "\n"
    if args.check is not None and json.loads(args.check.read_text()) != json.loads(encoded):
        raise SystemExit("saved family report differs from exact replay")
    if args.output is not None:
        args.output.write_text(encoded)
    print(
        json.dumps(
            {
                "normalized_B_count": report["normalized_B_count"],
                "D_count": report["D_count"],
                "operator_count": report["operator_count"],
                "local_direction_count": report["local_antipodal_direction_count"],
                "all_operators_excluded": report["all_operators_excluded_by_exact_modular_certificates"],
                "new_real_order_claimed": False,
            }
        )
    )


if __name__ == "__main__":
    main()

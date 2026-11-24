"""Replay the exact integer obstruction to the fixed Paley sextic frame."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from experiments.scalar_lift.paley_frame_moments import build_paley_frame_moment_report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.check is not None and args.output is not None:
        parser.error("use only one of --check and --output")
    report = build_paley_frame_moment_report()
    encoded = json.dumps(report, indent=2) + "\n"
    if args.check is not None and json.loads(args.check.read_text()) != json.loads(encoded):
        raise SystemExit("saved moment certificate differs from exact replay")
    if args.output is not None:
        args.output.write_text(encoded)
    print(
        json.dumps(
            {
                "order": report["order"],
                "local_state_count": len(report["moment_system"]["local_states"]),
                "integer_contradiction": report["left_null_product_with_rhs"],
                "full_sign_frame_excluded_for_this_fixed_operator": True,
                "new_real_hadamard_order_claimed": False,
            }
        )
    )


if __name__ == "__main__":
    main()

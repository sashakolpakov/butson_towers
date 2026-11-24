"""Replay all coefficient cases in the order-716 constant tensor audit."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from experiments.scalar_lift.order_716_family import build_order_716_family_report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = build_order_716_family_report()
    encoded = json.dumps(report, indent=2) + "\n"
    if args.check is not None and json.loads(args.check.read_text()) != json.loads(encoded):
        raise SystemExit("saved order-716 family report differs from exact replay")
    if args.output is not None:
        args.output.write_text(encoded)
    print(json.dumps({key: report[key] for key in (
        "order", "case_count", "excluded_case_count", "surviving_case_count",
        "survivors_are_exactly_the_commuting_B_D_cases", "flat_frame_constructed",
    )}))


if __name__ == "__main__":
    main()

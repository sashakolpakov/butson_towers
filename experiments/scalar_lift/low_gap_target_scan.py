"""Finite row-sum screen for multiplier charts aimed at unresolved low orders.

A prospective periodic T-input has value +e0 at zero. A nonzero orbit has
odd length m. Up to signed coordinate permutations, an odd coordinate
action is either identity or a three-cycle fixing coordinate zero.
The underlying permutation has order1 or3. Odd order forces a positive
sign product around every cycle; diagonal sign conjugation removes those
cycle signs. This accounts for both coordinate actions used here.

This module checks the necessary row-sum norm and signed orbit counts.
It constructs no parity map, signed T-input, frame or Hadamard endpoint.
The frontier is an attributed, pinned target list, not a novelty proof.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from math import isqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CERTIFICATE = Path("certificates/scalar-lift/low-gap-target-scan.json")
FRONTIER = Path("experiments/scalar_lift/frontier_120_snapshot.json")
LEDGER = Path("certificates/scalar-lift/completed-order-inventory.json")
SOURCES = (
    (Path("experiments/scalar_lift/williamson_inputs.json"), "order"),
    (Path("experiments/scalar_lift/near_williamson_inputs.json"), "width"),
)
# Recorded independent comparison with the author's current public Markdown.
# Replay checks the retained ordered list; it does not perform a network fetch.
PRIMARY_URL = (
    "https://raw.githubusercontent.com/vatskishan03/Hadamard/main/FRONTIER.md"
)
PRIMARY_SHA256 = (
    "e67e681689e3a02a3426d8c70800244bb7e2b6eb0f024987fa7ca0ff63f0e9bd"
)
ORDER_LIST_SHA256 = (
    "e287da90e67171dc9274ffb8c8c8a810ba5d03775f6aeb8375a8c015dc02e3ec"
)


def _positive_integer(value, name):
    if type(value) is not int or value < 1:
        raise ValueError(f"{name} must be a positive integer")
    return value


def is_prime(value):
    """Complete trial division suffices for the bounded target factors."""
    return (type(value) is int and value >= 2 and
            all(value % divisor for divisor in range(2, isqrt(value) + 1)))


def _chart_size(p, m):
    if not is_prime(p) or p == 2:
        raise ValueError("an odd prime length is required")
    _positive_integer(m, "multiplier order")
    if m % 2 == 0 or (p - 1) % m:
        raise ValueError("multiplier order must be an odd divisor of p-1")
    return (p - 1) // m


def identity_dc_rows(p, m):
    """All identity-action row sums modulo tail signs and permutations.

    The sums are (1+m*h, m*b, m*c, m*d), with b>=c>=d>=0. A signed
    count vector is realizable on N nonzero orbits exactly when its
    absolute sum is at most N and the remaining number of orbits is even.
    Allocate the demanded signed orbits first, then opposite-sign pairs.
    This is complete for these row sums, not for periodic correlations.
    """
    count = _chart_size(p, m)
    result = []
    for a in range(-isqrt(p), isqrt(p) + 1):
        if (a - 1) % m or (p - a * a) % (m * m):
            continue
        head = (a - 1) // m
        square_sum = (p - a * a) // (m * m)
        for small in range(isqrt(square_sum) + 1):
            rest = square_sum - small * small
            for middle in range(small, isqrt(rest) + 1):
                large_squared = rest - middle * middle
                large = isqrt(large_squared)
                if large < middle or large * large != large_squared:
                    continue
                minimum = abs(head) + small + middle + large
                if minimum <= count and (count - minimum) % 2 == 0:
                    result.append({
                        "row": [a, m * large, m * middle, m * small],
                        "orbit_sign_sums": [head, large, middle, small],
                    })
    return sorted(result, key=lambda row: row["row"])


def c3_dc_rows(p, m):
    """All three-cycle row sums, with the common tail sum positive.

    Each tail orbit contributes m/3 to each of coordinates1,2,3. Thus the
    sums are (1+m*h, (m/3)*k, (m/3)*k, (m/3)*k). Negative tail sums may
    be reversed simultaneously. A zero tail sum would make p a square.
    """
    count = _chart_size(p, m)
    if m % 3:
        raise ValueError(
            "a three-cycle coordinate action requires 3 to divide m")
    unit = m // 3
    result = []
    for a in range(-isqrt(p), isqrt(p) + 1):
        if (a - 1) % m or (p - a * a) % 3:
            continue
        tail_squared = (p - a * a) // 3
        tail = isqrt(tail_squared)
        if tail == 0 or tail * tail != tail_squared or tail % unit:
            continue
        head, signed_tail = (a - 1) // m, tail // unit
        minimum = abs(head) + signed_tail
        if minimum <= count and (count - minimum) % 2 == 0:
            result.append({
                "row": [a, tail, tail, tail],
                "head_sign_sum": head,
                "tail_sign_sum": signed_tail,
            })
    return sorted(result, key=lambda row: row["row"])


def scan_targets(orders, widths, completed_orders=(), *, max_orbits=100,
                 max_order=12000):
    """Exhaust the declared prime/odd-multiplier/action possibilities."""
    _positive_integer(max_orbits, "max_orbits")
    _positive_integer(max_order, "max_order")
    orders, widths = tuple(orders), tuple(widths)
    if any(type(n) is not int or n < 4 or n % 4 for n in orders):
        raise ValueError("target orders must be positive multiples of four")
    for width in widths:
        _positive_integer(width, "source width")
    orders, widths = sorted(set(orders)), sorted(set(widths))
    completed = set(completed_orders)
    charts = []
    for order in orders:
        if order > max_order or order in completed:
            continue
        for width in widths:
            if order % (4 * width):
                continue
            p = order // (4 * width)
            if not is_prime(p) or p == 2:
                continue
            for m in range(1, p, 2):
                if (p - 1) % m or (p - 1) // m > max_orbits:
                    continue
                actions = [("identity", identity_dc_rows(p, m))]
                if m % 3 == 0:
                    actions.append(("C3", c3_dc_rows(p, m)))
                for action, rows in actions:
                    if rows:
                        charts.append({
                            "target": order, "p": p, "width": width, "m": m,
                            "nonorigin_orbits": (p - 1) // m,
                            "action": action, "necessary_DC_rows": rows,
                        })
    return sorted(charts, key=lambda row: (
        row["target"], row["nonorigin_orbits"], row["action"],
        -row["m"], row["width"]))


def build_report(root=ROOT, *, max_orbits=100, max_order=12000):
    """Bind the scan to the frontier, completed ledger and source tables."""
    root = Path(root)
    hashes = {}

    def load(relative):
        raw = (root / relative).read_bytes()
        hashes[relative.as_posix()] = hashlib.sha256(raw).hexdigest()
        return json.loads(raw)

    frontier, ledger = load(FRONTIER), load(LEDGER)
    orders = frontier["orders"]
    order_bytes = json.dumps(orders, separators=(",", ":")).encode()
    order_digest = hashlib.sha256(order_bytes).hexdigest()
    if (orders != sorted(set(orders)) or len(orders) != 120 or
            order_digest != ORDER_LIST_SHA256):
        raise ValueError(
            "retained frontier differs from the cross-checked 120-order list")
    completed = sorted(set(ledger["completed_new_real_orders"]) |
                       set(ledger["completed_known_real_orders"]))
    if completed != sorted(row["order"] for row in ledger["orders"]
                           if row["construction_complete"]):
        raise ValueError(
            "completed-order lists disagree with the positive ledger")
    widths, provenance = set(), {}
    for path, field in SOURCES:
        entries = load(path)["entries"]
        for index, entry in enumerate(entries):
            width = _positive_integer(entry[field], "source width")
            widths.add(width)
            provenance.setdefault(str(width), []).append({
                "path": path.as_posix(), "selector": f"entries[{index}]",
            })
    charts = scan_targets(orders, widths, completed, max_orbits=max_orbits,
                          max_order=max_order)
    screened = [order for order in orders
                if order <= max_order and order not in completed]
    by_cost = sorted(range(len(charts)), key=lambda i: (
        charts[i]["nonorigin_orbits"], charts[i]["target"],
        charts[i]["action"]))
    return {
        "schema": "low-gap-prime-multiplier-necessary-arithmetic-v1",
        "producer": "experiments/scalar_lift/low_gap_target_scan.py",
        "scope": {
            "object": (
                "Necessary row-sum arithmetic for prime periodic T-inputs "
                "with an odd multiplier."),
            "normalization": (
                "Value at zero is +e0; tail signs and coordinate "
                "permutations are normalized."),
            "coordinate_actions": [
                "identity",
                "C3 fixing coordinate0 and cycling coordinates1,2,3",
            ],
            "action_coverage": (
                "Every odd-order signed permutation of four coordinates "
                "is conjugate to one of these actions; even multiplier "
                "orders are outside this scan."),
            "complete_only_for_declared_DC_and_orbit_count_conditions": True,
            "parity_map_constructed": False,
            "signed_correlations_verified": False,
            "endpoint_constructed": False,
            "existence_nonexistence_or_novelty_claim": False,
        },
        "limits": {
            "max_target_order_inclusive": max_order,
            "max_nonorigin_orbits": max_orbits,
        },
        "frontier": {
            "path": FRONTIER.as_posix(),
            "snapshot_source_url": frontier["source_url"],
            "snapshot_retrieved_utc_date": frontier["retrieved_utc_date"],
            "scope": frontier["scope"], "rows": len(orders),
            "ordered_list_sha256": order_digest,
            "recorded_primary_cross_check": {
                "url": PRIMARY_URL, "source_sha256": PRIMARY_SHA256,
                "same_ordered_120_rows": True,
                "replay_scope": (
                    "The recorded comparison is retained; offline replay "
                    "does not fetch the live URL."),
            },
        },
        "completed_order_ledger": LEDGER.as_posix(),
        "excluded_completed_orders": sorted(set(orders) & set(completed)),
        "source_widths": sorted(widths), "source_width_provenance": provenance,
        "summary": {
            "completed_orders_in_ledger": len(completed),
            "frontier_orders_screened": len(screened),
            "surviving_target_orders": sorted({
                row["target"] for row in charts
            }),
            "target_count": len({row["target"] for row in charts}),
            "chart_rows": len(charts),
            "distinct_prime_multiplier_action_charts": len({
                (row["p"], row["m"], row["action"]) for row in charts
            }),
            "action_counts": dict(Counter(row["action"] for row in charts)),
        },
        "charts": charts,
        "chart_indices_by_nonorigin_orbit_count": by_cost,
        "source_sha256": dict(sorted(hashes.items())),
    }


def certificate_bytes(report):
    return (json.dumps(report, indent=2, ensure_ascii=False) + "\n").encode()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / CERTIFICATE)
    parser.add_argument(
        "--check", action="store_true",
        help="require an exact match to the existing certificate")
    parser.add_argument("--max-orbits", type=int, default=100)
    parser.add_argument("--max-order", type=int, default=12000)
    args = parser.parse_args(argv)
    try:
        report = build_report(max_orbits=args.max_orbits,
                              max_order=args.max_order)
        expected = certificate_bytes(report)
        if args.check:
            if args.output.read_bytes() != expected:
                raise ValueError(
                    "certificate differs from current inputs; "
                    "regenerate without --check")
        else:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_bytes(expected)
    except (OSError, ValueError, KeyError, TypeError) as exc:
        parser.exit(1, f"Low-gap target audit failed: {exc}\n")
    print(json.dumps({"ok": True, "mode": "check" if args.check else "write",
                      **report["summary"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

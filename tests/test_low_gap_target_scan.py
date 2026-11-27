"""Independent finite signed-orbit controls for the low-gap arithmetic scan."""
from pathlib import Path
import json

import pytest

from experiments.scalar_lift import low_gap_target_scan as scan

ROOT = Path(__file__).resolve().parents[1]


def brute_rows(p, m, action):
    """Add actual orbit row-sum contributions, without the support formula."""
    count = (p - 1) // m
    if action == "identity":
        positive = [tuple(m if j == i else 0 for j in range(4))
                    for i in range(4)]
    else:
        positive = [(m, 0, 0, 0), (0, m // 3, m // 3, m // 3)]
    contributions = positive + [tuple(-x for x in row) for row in positive]
    reachable = {(1, 0, 0, 0)}
    for _ in range(count):
        reachable = {tuple(x + y for x, y in zip(row, add))
                     for row in reachable for add in contributions}
    result = set()
    for row in reachable:
        if sum(x * x for x in row) != p:
            continue
        if action == "identity":
            normalized = (row[0], *sorted(map(abs, row[1:]), reverse=True))
        else:
            normalized = (row[0], *(abs(x) for x in row[1:]))
        result.add(normalized)
    return result


def test_dc_rows_equal_independent_signed_orbit_enumeration():
    for p in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31):
        for m in range(1, p, 2):
            if (p - 1) % m or (p - 1) // m > 8:
                continue
            actual = {tuple(row["row"]) for row in scan.identity_dc_rows(p, m)}
            assert actual == brute_rows(p, m, "identity"), (p, m)
            if m % 3 == 0:
                actual = {tuple(row["row"]) for row in scan.c3_dc_rows(p, m)}
                assert actual == brute_rows(p, m, "C3"), (p, m)


def test_small_target_scan_exhausts_both_actions_and_respects_limits():
    widths = (1, 3)
    orders = {4 * p * w for p in (3, 5, 7, 11, 13, 19, 23) for w in widths}
    completed = {28}
    observed = scan.scan_targets(orders, widths, completed,
                                 max_orbits=8, max_order=160)
    expected = set()
    for order in orders - completed:
        if order > 160:
            continue
        for width in widths:
            if order % (4 * width):
                continue
            p = order // (4 * width)
            if p < 3 or any(p % k == 0 for k in range(2, p)):
                continue
            for m in range(1, p, 2):
                if (p - 1) % m or (p - 1) // m > 8:
                    continue
                for action in ("identity", "C3"):
                    if action == "C3" and m % 3:
                        continue
                    for row in brute_rows(p, m, action):
                        expected.add((order, p, width, m, action, row))
    actual = {(r["target"], r["p"], r["width"], r["m"], r["action"],
               tuple(row["row"]))
              for r in observed for row in r["necessary_DC_rows"]}
    assert actual == expected
    assert all(r["target"] != 28 for r in observed)


def test_public_certificate_and_promising_low_order_charts_replay():
    report = scan.build_report()
    saved = (ROOT / scan.CERTIFICATE).read_bytes()
    assert saved == scan.certificate_bytes(report)
    assert report["excluded_completed_orders"] == [7572]
    assert report["summary"]["completed_orders_in_ledger"] == 65
    assert report["summary"]["frontier_orders_screened"] == 119
    assert report["summary"]["chart_rows"] == 15
    assert report["summary"]["target_count"] == 9
    assert report["summary"]["action_counts"] == {"C3": 5, "identity": 10}
    by_key = {(r["target"], r["p"], r["m"], r["action"]): r
              for r in report["charts"]}
    row = by_key[5892, 491, 7, "identity"]
    assert row["nonorigin_orbits"] == 70
    assert row["necessary_DC_rows"] == [
        {"row": [1, 21, 7, 0], "orbit_sign_sums": [0, 3, 1, 0]}]
    assert by_key[9260, 463, 11, "identity"]["nonorigin_orbits"] == 42
    assert by_key[2676, 223, 3, "identity"]["nonorigin_orbits"] == 74
    assert by_key[2956, 739, 9, "C3"]["nonorigin_orbits"] == 82
    assert not report["scope"]["endpoint_constructed"]
    assert not report["scope"]["existence_nonexistence_or_novelty_claim"]


def test_every_surviving_row_has_the_exact_norm_and_count_parity():
    for chart in scan.build_report()["charts"]:
        p, m, count = chart["p"], chart["m"], chart["nonorigin_orbits"]
        assert scan.is_prime(p) and m % 2 and (p - 1) // m == count <= 100
        for record in chart["necessary_DC_rows"]:
            assert sum(x * x for x in record["row"]) == p
            if chart["action"] == "identity":
                nets = record["orbit_sign_sums"]
                expected = [1 + m * nets[0]] + [m * x for x in nets[1:]]
            else:
                nets = [record["head_sign_sum"], record["tail_sign_sum"]]
                assert m % 3 == 0
                expected = [1 + m * nets[0]] + [m // 3 * nets[1]] * 3
            assert record["row"] == expected
            minimum = sum(map(abs, nets))
            assert minimum <= count and (count - minimum) % 2 == 0


def test_limit_changes_and_completed_order_exclusion_are_effective():
    assert not scan.scan_targets([2676], [3], max_orbits=73)
    assert len(scan.scan_targets([2676], [3], max_orbits=74)) == 2
    assert not scan.scan_targets([2676], [3], [2676])
    assert not scan.scan_targets([2676], [3], max_order=2675)


@pytest.mark.parametrize("p,m", [(2, 1), (9, 1), (7, 2), (7, 5), (7, 0)])
def test_invalid_chart_parameters_are_rejected(p, m):
    with pytest.raises(ValueError):
        scan.identity_dc_rows(p, m)


def test_invalid_action_and_boolean_parameters_are_rejected():
    with pytest.raises(ValueError, match="3 to divide"):
        scan.c3_dc_rows(11, 5)
    with pytest.raises(ValueError):
        scan.scan_targets([2676], [True])
    with pytest.raises(ValueError):
        scan.scan_targets([2676], [3], max_orbits=False)


def test_check_mode_rejects_changed_certificate(tmp_path):
    path = tmp_path / "scan.json"
    assert scan.main(["--output", str(path)]) == 0
    assert scan.main(["--output", str(path), "--check"]) == 0
    report = json.loads(path.read_text())
    report["scope"]["endpoint_constructed"] = True
    path.write_text(json.dumps(report))
    with pytest.raises(SystemExit) as error:
        scan.main(["--output", str(path), "--check"])
    assert error.value.code == 1

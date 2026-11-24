"""Bounded operator-free OD(4p;3p,p) search in odd-multiplier charts.

Four sequences have entries +/-x or +/-y. Coefficientwise complementary
autocorrelations produce a Goethals--Seidel orthogonal design, so both
sextic sign blocks remain unknown. This is an exact restricted chart,
not a repair objective and not a prescribed cyclotomic operator.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import signal
import subprocess
import sys
import time
from contextlib import suppress
from itertools import product
from math import isqrt
from pathlib import Path

import numpy as np
from sympy import primitive_root

from hadamard.inverse_lift import sixth_root_block_pair_conditions

PRIME = 223
MULTIPLIER = 3
SYMBOLS = ((1, 0), (-1, 0), (0, 1), (0, -1))


def goethals_seidel_array(sequences):
    """The signed four-circulant array, valid over the coefficient ring."""
    sequences = np.asarray(sequences, dtype=np.int64)
    if sequences.ndim != 2 or sequences.shape[0] != 4 or sequences.shape[1] == 0:
        raise ValueError("four nonempty equally long sequences are required")
    length = sequences.shape[1]
    indices = np.arange(length)
    a, b, c, d = [row[(indices[None, :] - indices[:, None]) % length] for row in sequences]
    reversal = np.eye(length, dtype=np.int64)[::-1]
    return np.block([
        [a, b @ reversal, c @ reversal, d @ reversal],
        [-b @ reversal, a, -d.T @ reversal, c.T @ reversal],
        [-c @ reversal, d.T @ reversal, a, -b.T @ reversal],
        [-d @ reversal, -c.T @ reversal, b.T @ reversal, a],
    ])


def orbit_data(prime=PRIME, multiplier=MULTIPLIER):
    if (prime - 1) % multiplier:
        raise ValueError("the multiplier order must divide p-1")
    generator = pow(int(primitive_root(prime)), (prime - 1) // multiplier, prime)
    subgroup = sorted({pow(generator, k, prime) for k in range(multiplier)})
    remaining = set(range(1, prime))
    orbits = [(0,)]
    while remaining:
        orbit = tuple(sorted(a * min(remaining) % prime for a in subgroup))
        orbits.append(orbit)
        remaining.difference_update(orbit)
    orbit_index = [0] * prime
    for index, orbit in enumerate(orbits):
        for value in orbit:
            orbit_index[value] = index
    shifts, covered = [], set()
    for shift in range(1, prime):
        if shift not in covered:
            shifts.append(shift)
            covered.update(sign * shift * a % prime for sign in (-1, 1) for a in subgroup)
    pairs = [(i, j) for i in range(len(orbits)) for j in range(i, len(orbits))]
    pair_index = {pair: k for k, pair in enumerate(pairs)}
    coefficients = np.zeros((len(shifts), len(pairs)), dtype=np.int64)
    for row, shift in enumerate(shifts):
        for x in range(prime):
            pair = tuple(sorted((orbit_index[x], orbit_index[(x + shift) % prime])))
            coefficients[row, pair_index[pair]] += 1
    if any(sum(row) != prime for row in coefficients):
        raise ValueError("an orbit difference table has incorrect mass")
    return orbits, orbit_index, shifts, pairs, coefficients


def row_sum_types(prime=PRIME, multiplier=MULTIPLIER):
    """All trivial-character types after exact origin and permutation gauges."""
    p = prime
    a_bound, b_bound = isqrt(3 * p), isqrt(p)
    a_domains = [[x for x in range(-a_bound, a_bound + 1) if x % multiplier == (0 if i == 0 else 1)] for i in range(4)]
    b_domains = [[x for x in range(-b_bound, b_bound + 1) if x % multiplier == (1 if i == 0 else 0)] for i in range(4)]
    first = [row for row in product(*a_domains) if sum(x * x for x in row) == 3 * p]
    second = [row for row in product(*b_domains) if sum(x * x for x in row) == p and row[1] <= row[2] <= row[3]]
    return [(*a, *b) for a in first for b in second
            if all((x + y) % 2 == 1 for x, y in zip(a, b, strict=True))
            and sum(x * y for x, y in zip(a, b, strict=True)) == 0]


def build_model(prime=PRIME, multiplier=MULTIPLIER):
    from ortools.sat.python import cp_model

    if multiplier < 3 or multiplier % 2 == 0:
        raise ValueError("this origin gauge requires an odd multiplier at least three")
    orbits, _index, shifts, pairs, coefficients = orbit_data(prime, multiplier)
    count = len(orbits)
    model = cp_model.CpModel()
    states = [[model.new_int_var(0, 3, f"s{row}_{i}") for i in range(count)] for row in range(4)]
    a, b, channel = [], [], []
    for row in range(4):
        aa, bb, yy = [], [], []
        for i in range(count):
            av = model.new_int_var(-1, 1, f"a{row}_{i}")
            bv = model.new_int_var(-1, 1, f"b{row}_{i}")
            yv = model.new_bool_var(f"y{row}_{i}")
            model.add_allowed_assignments([states[row][i], av, bv, yv],
                                          [[s, x, y, int(y != 0)] for s, (x, y) in enumerate(SYMBOLS)])
            aa.append(av)
            bb.append(bv)
            yy.append(yv)
        a.append(aa)
        b.append(bb)
        channel.append(yy)
        model.add(states[row][0] == (2 if row == 0 else 0))
    sizes = list(map(len, orbits))
    model.add(sum(sizes[i] * channel[row][i] for row in range(4) for i in range(count)) == prime)
    a_sums = [model.new_int_var(-isqrt(3 * prime), isqrt(3 * prime), f"A{row}") for row in range(4)]
    b_sums = [model.new_int_var(-isqrt(prime), isqrt(prime), f"B{row}") for row in range(4)]
    for row in range(4):
        model.add(a_sums[row] == cp_model.LinearExpr.weighted_sum(a[row], sizes))
        model.add(b_sums[row] == cp_model.LinearExpr.weighted_sum(b[row], sizes))
    types = row_sum_types(prime, multiplier)
    model.add_allowed_assignments([*a_sums, *b_sums], types)
    product_table = [[s, t, x * xx, y * yy, x * yy + y * xx]
                     for s, (x, y) in enumerate(SYMBOLS)
                     for t, (xx, yy) in enumerate(SYMBOLS)]
    products = [[], [], []]
    for row in range(4):
        local = [[], [], []]
        for i, j in pairs:
            if i == j:
                values = [1 - channel[row][i], channel[row][i], 0]
            else:
                values = [model.new_int_var(-1, 1, f"c{k}_{row}_{i}_{j}") for k in range(3)]
                model.add_allowed_assignments([states[row][i], states[row][j], *values], product_table)
            for k in range(3):
                local[k].append(values[k])
        for k in range(3):
            products[k].extend(local[k])
    for coefficient_row in coefficients:
        weights = coefficient_row.tolist() * 4
        for k in range(3):
            model.add(cp_model.LinearExpr.weighted_sum(products[k], weights) == 0)
    error = model.validate()
    if error:
        raise ValueError(error)
    return model, states, {
        "prime": prime, "multiplier_order": multiplier,
        "block_order": 4 * prime, "design_type": [3 * prime, prime],
        "algebraic_source_3_and_13_output_orders": [12 * prime, 52 * prime],
        "orbit_count": count, "four_state_orbit_variables": 4 * count,
        "independent_nonzero_shift_classes": len(shifts),
        "coefficient_correlation_equations": 3 * len(shifts),
        "trivial_character_types": len(types),
        "total_model_variables": len(model.proto.variables),
        "total_model_constraints": len(model.proto.constraints),
        "scalar_operator_prescribed": False,
        "all_four_design_sequences_unknown": True,
        "objective": "exact feasibility; every residual constrained to zero",
    }


def verify_sequences(states, prime=PRIME, multiplier=MULTIPLIER):
    """Independently expand and check every coefficient at every cyclic shift."""
    orbits, index, _shifts, _pairs, _coefficients = orbit_data(prime, multiplier)
    if len(states) != 4 or any(len(row) != len(orbits) for row in states):
        raise ValueError("four complete orbit-symbol rows are required")
    if any(type(value) is not int or not 0 <= value < 4 for row in states for value in row):
        raise ValueError("invalid design symbols")
    expanded = [[SYMBOLS[row[index[x]]] for x in range(prime)] for row in states]
    a = [[pair[0] for pair in row] for row in expanded]
    b = [[pair[1] for pair in row] for row in expanded]
    if sum(x * x for row in a for x in row) != 3 * prime or sum(x * x for row in b for x in row) != prime:
        raise ValueError("incorrect orthogonal-design weights")
    for shift in range(prime):
        aa = sum(row[x] * row[(x + shift) % prime] for row in a for x in range(prime))
        bb = sum(row[x] * row[(x + shift) % prime] for row in b for x in range(prime))
        mixed = sum(a[r][x] * b[r][(x + shift) % prime] + b[r][x] * a[r][(x + shift) % prime]
                    for r in range(4) for x in range(prime))
        if (aa, bb, mixed) != ((3 * prime, prime, 0) if shift == 0 else (0, 0, 0)):
            raise ValueError(f"coefficient autocorrelation failed at shift {shift}")
    first = goethals_seidel_array(np.array(a) - np.array(b))
    second = goethals_seidel_array(np.array(a) + np.array(b))
    pair_check = sixth_root_block_pair_conditions(first, second)
    if not pair_check["valid"]:
        raise ValueError("the independent full sign-block pair verification failed")
    return {"orbit_symbols": states, "a_sequences": a, "b_sequences": b,
            "every_expanded_integer_autocorrelation_verified": True,
            "design_order": 4 * prime, "design_type": [3 * prime, prime],
            "full_sign_block_pair_independently_verified": pair_check,
            "first_block_binary_sha256": hashlib.sha256(((first + 1) // 2).astype(np.uint8).tobytes()).hexdigest(),
            "second_block_binary_sha256": hashlib.sha256(((second + 1) // 2).astype(np.uint8).tobytes()).hexdigest()}


def solve_worker(seconds, prime=PRIME, multiplier=MULTIPLIER):
    from ortools.sat.python import cp_model

    start = time.monotonic()

    def emit(stage, **values):
        print(json.dumps({"stage": stage, "elapsed_seconds": time.monotonic() - start, **values}), flush=True)

    emit("worker_started", pid=os.getpid())
    model, states, metadata = build_model(prime, multiplier)
    emit("model_ready", **metadata)
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.max_presolve_iterations = 1
    solver.parameters.num_search_workers = 1
    solver.parameters.random_seed = 0
    status = solver.solve(model)
    witness = None
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        witness = verify_sequences([[solver.value(value) for value in row] for row in states], prime, multiplier)
    emit("solver_finished", solver_status=solver.status_name(status),
         solver_wall_seconds=solver.wall_time, conflicts=solver.num_conflicts,
         branches=solver.num_branches, exact_design_constructed=witness is not None, witness=witness)


def bounded_probe(wall_seconds=60, solver_seconds=45, prime=PRIME, multiplier=MULTIPLIER):
    if wall_seconds <= 0 or solver_seconds <= 0:
        raise ValueError("positive native and wall limits are required")
    root = Path(__file__).resolve().parents[2]
    environment = os.environ.copy()
    environment["PYTHONPATH"] = os.pathsep.join((str(root), str(root / "src")))
    command = [sys.executable, "-m", "experiments.scalar_lift.joint_gs_223", "--worker", "--seconds", str(solver_seconds),
               "--prime", str(prime), "--multiplier", str(multiplier)]
    start = time.monotonic()
    worker = subprocess.Popen(command, cwd=root, env=environment, stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, text=True, start_new_session=True)

    def stop():
        if worker.poll() is None:
            with suppress(ProcessLookupError):
                os.killpg(worker.pid, signal.SIGKILL)

    timed_out = False
    try:
        output, _ = worker.communicate(timeout=wall_seconds)
    except subprocess.TimeoutExpired:
        timed_out = True
        stop()
        output, _ = worker.communicate(timeout=5)
    except BaseException:
        stop()
        worker.communicate(timeout=5)
        raise
    events, diagnostics = [], []
    for line in output.splitlines():
        try:
            events.append(json.loads(line))
        except json.JSONDecodeError:
            diagnostics.append(line)
    final = next((event for event in reversed(events) if event.get("stage") == "solver_finished"), None)
    witness = None if final is None else final.get("witness")
    if witness is not None and verify_sequences(witness["orbit_symbols"], prime, multiplier) != witness:
        raise ValueError("the independent supervisor replay disagrees")
    return {
        "scope": "One bounded exact feasibility search for OD(4p;3p,p) in a four-sequence odd-multiplier Goethals-Seidel chart",
        "prime": prime, "multiplier_order": multiplier,
        "algebraic_source_3_and_13_output_orders": [12 * prime, 52 * prime], "block_order": 4 * prime,
        "native_limit_seconds": solver_seconds, "wall_limit_seconds": wall_seconds,
        "elapsed_seconds": time.monotonic() - start,
        "worker_pid": worker.pid, "worker_exit_code": worker.returncode,
        "worker_reaped": worker.returncode is not None,
        "external_watchdog_triggered": timed_out,
        "events": events, "diagnostics": diagnostics,
        "solver_status": None if final is None else final["solver_status"],
        "exact_design_constructed": witness is not None,
        "inconclusive": timed_out or worker.returncode != 0 or final is None or final["solver_status"] not in ("OPTIMAL", "FEASIBLE", "INFEASIBLE"),
        "mathematical_nonexistence_certificate_produced": False,
        "new_real_hadamard_order_claimed": False,
        "witness": witness,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--prime", type=int, default=PRIME)
    parser.add_argument("--multiplier", type=int, default=MULTIPLIER)
    parser.add_argument("--seconds", type=float, default=45)
    parser.add_argument("--wall-seconds", type=float, default=60)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.worker:
        solve_worker(args.seconds, args.prime, args.multiplier)
        return
    report = bounded_probe(args.wall_seconds, args.seconds, args.prime, args.multiplier)
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({key: report[key] for key in (
        "prime", "multiplier_order", "algebraic_source_3_and_13_output_orders", "solver_status", "elapsed_seconds", "worker_reaped",
        "external_watchdog_triggered", "exact_design_constructed", "inconclusive",
    )}))


if __name__ == "__main__":
    main()

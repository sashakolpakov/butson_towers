"""Exact high-order joint-design search using periodic one-hot quadruples.

At every cyclic coordinate b is one of +/-e_0,...,+/-e_3. Setting a=Db
for the skew conference matrix of order four makes the mixed equations
identically zero and the a autocorrelation three times that of b. Only
the perfect summed b autocorrelation remains unknown. The order-4p
scalar action is recovered from a completed sign pair, not prescribed.
"""

from __future__ import annotations

import argparse
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

from experiments.scalar_lift.joint_gs_223 import SYMBOLS, orbit_data, verify_sequences
from hadamard.scalar_lift import SKEW_CONFERENCE_4

STATES = tuple(tuple(sign if i == axis else 0 for i in range(4))
               for axis in range(4) for sign in (1, -1))


def row_sum_types(prime, multiplier):
    """All norm types after b(0)=e0 and signed permutations of axes1,2,3."""
    bound = isqrt(prime)
    tails = [tail for tail in product(range(0, bound + 1, multiplier), repeat=3)
             if tail[0] <= tail[1] <= tail[2]]
    return [(head, *tail) for head in range(-bound, bound + 1)
            if head % multiplier == 1 % multiplier for tail in tails
            if head * head + sum(value * value for value in tail) == prime]


def normal_form_state(row_sums):
    """Positive nonorigin-axis sum guarantees an orbit with this positive state."""
    return 2 * next(axis for axis in range(1, 4) if row_sums[axis] > 0)


def normalize_orbit_scaling(states, prime=223, multiplier=3):
    """Use t -> c*t to place a guaranteed positive axis on the orbit of one."""
    orbits, index, _shifts, _pairs, _coefficients = orbit_data(prime, multiplier)
    if len(states) != len(orbits) or any(type(s) is not int or not 0 <= s < 8 for s in states):
        raise ValueError("a complete eight-state orbit assignment is required")
    if states[0] != 0:
        raise ValueError("the origin must already be e0")
    expanded = np.array([STATES[states[index[t]]] for t in range(prime)])
    sums = tuple(int(x) for x in expanded.sum(axis=0))
    if sums not in row_sum_types(prime, multiplier):
        raise ValueError("the row sums are not in the signed-axis normal form")
    chosen = normal_form_state(sums)
    orbit = next(i for i in range(1, len(orbits)) if states[i] == chosen)
    scale = orbits[orbit][0]
    normalized = [states[index[(scale * representatives[0]) % prime]] for representatives in orbits]
    if normalized[0] != 0 or normalized[1] != chosen:
        raise ValueError("the multiplicative orbit normalization failed")
    return normalized, scale


def build_model(prime=223, multiplier=3, orbit_normal_form=True):
    from ortools.sat.python import cp_model

    orbits, _index, shifts, pairs, coefficients = orbit_data(prime, multiplier)
    model = cp_model.CpModel()
    states = [model.new_int_var(0, 7, f"s{i}") for i in range(len(orbits))]
    model.add(states[0] == 0)
    coordinates = [[model.new_int_var(-1, 1, f"b{axis}_{i}") for i in range(len(orbits))]
                   for axis in range(4)]
    for i in range(len(orbits)):
        model.add_allowed_assignments([states[i], *(coordinates[axis][i] for axis in range(4))],
                                      [[s, *vector] for s, vector in enumerate(STATES)])
    sizes = list(map(len, orbits))
    sums = [model.new_int_var(-isqrt(prime), isqrt(prime), f"sum{axis}") for axis in range(4)]
    for axis in range(4):
        model.add(sums[axis] == cp_model.LinearExpr.weighted_sum(coordinates[axis], sizes))
    types = row_sum_types(prime, multiplier)
    model.add_allowed_assignments(sums, types)
    if orbit_normal_form:
        # Rescaling t -> c*t preserves the orbit chart, all sums and all PAFs.
        model.add_allowed_assignments([*sums, states[1]],
                                      [[*row, normal_form_state(row)] for row in types])
    table = [[s, t, sum(a * b for a, b in zip(v, w, strict=True))]
             for s, v in enumerate(STATES) for t, w in enumerate(STATES)]
    dot_products = []
    for i, j in pairs:
        if i == j:
            dot_products.append(1)
        else:
            dot = model.new_int_var(-1, 1, f"dot{i}_{j}")
            model.add_allowed_assignments([states[i], states[j], dot], table)
            dot_products.append(dot)
    for coefficient_row in coefficients:
        model.add(cp_model.LinearExpr.weighted_sum(dot_products, coefficient_row.tolist()) == 0)
    error = model.validate()
    if error:
        raise ValueError(error)
    return model, states, {
        "prime": prime, "multiplier_order": multiplier,
        "block_order": 4 * prime, "design_type": [3 * prime, prime],
        "algebraic_source_3_and_13_output_orders": [12 * prime, 52 * prime],
        "one_hot_eight_state_orbit_variables": len(orbits),
        "independent_nonzero_shift_classes": len(shifts),
        "row_sum_types": [list(row) for row in types],
        "multiplicative_orbit_normal_form": orbit_normal_form,
        "normal_form_first_nonzero_orbit_states": sorted({normal_form_state(row) for row in types}) if orbit_normal_form else None,
        "total_model_variables": len(model.proto.variables),
        "total_model_constraints": len(model.proto.constraints),
        "all_one_hot_sequences_unknown": True,
        "generator_sequence_relation": "a=D*b with fixed skew conference D4",
        "order_4p_scalar_operator_prescribed": False,
        "objective": "exact feasibility; every autocorrelation residual is zero",
    }


def verify_t_states(states, prime=223, multiplier=3):
    """Expand all cyclic correlations before invoking the full GS verifier."""
    orbits, index, _shifts, _pairs, _coefficients = orbit_data(prime, multiplier)
    if len(states) != len(orbits) or any(type(s) is not int or not 0 <= s < 8 for s in states):
        raise ValueError("a complete eight-state orbit assignment is required")
    b = np.array([STATES[states[index[x]]] for x in range(prime)], dtype=np.int64).T
    d = np.array(SKEW_CONFERENCE_4, dtype=np.int64)
    if np.any(d + d.T) or not np.array_equal(d.T @ d, 3 * np.eye(4, dtype=np.int64)):
        raise ValueError("the conference relation failed")
    a = d @ b
    if np.any(a * b) or not np.all(a * a + b * b == 1):
        raise ValueError("the two colored sequence supports do not complement")
    for shift in range(prime):
        bb = int(np.sum(b * np.roll(b, -shift, axis=1)))
        aa = int(np.sum(a * np.roll(a, -shift, axis=1)))
        mixed = int(np.sum(a * np.roll(b, -shift, axis=1) + b * np.roll(a, -shift, axis=1)))
        if bb != (prime if shift == 0 else 0) or aa != 3 * bb or mixed:
            raise ValueError(f"expanded coefficient correlation failed at shift {shift}")
    symbol_index = {symbol: i for i, symbol in enumerate(SYMBOLS)}
    orbit_b = np.array([STATES[s] for s in states], dtype=np.int64).T
    orbit_a = d @ orbit_b
    design_states = [[symbol_index[int(orbit_a[row, i]), int(orbit_b[row, i])]
                      for i in range(len(orbits))] for row in range(4)]
    gs = verify_sequences(design_states, prime, multiplier)
    return {"one_hot_orbit_states": states,
            "every_expanded_periodic_correlation_verified": True,
            "sequence_row_sums": b.sum(axis=1).tolist(),
            "full_GS_sign_pair_verification": gs}


def solve_worker(seconds, prime, multiplier, orbit_normal_form=True):
    from ortools.sat.python import cp_model

    started = time.monotonic()

    def emit(stage, **values):
        print(json.dumps({"stage": stage, "elapsed_seconds": time.monotonic() - started, **values}), flush=True)

    emit("worker_started", pid=os.getpid())
    model, states, metadata = build_model(prime, multiplier, orbit_normal_form)
    emit("model_ready", **metadata)
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.max_presolve_iterations = 1
    solver.parameters.num_search_workers = 1
    solver.parameters.random_seed = 0
    status = solver.solve(model)
    witness = None
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        witness = verify_t_states([solver.value(state) for state in states], prime, multiplier)
    emit("solver_finished", solver_status=solver.status_name(status),
         solver_wall_seconds=solver.wall_time, branches=solver.num_branches,
         conflicts=solver.num_conflicts, witness=witness)


def bounded_probe(wall_seconds=60, solver_seconds=45, prime=223, multiplier=3, orbit_normal_form=True):
    if wall_seconds <= 0 or solver_seconds <= 0:
        raise ValueError("positive native and wall limits are required")
    root = Path(__file__).resolve().parents[2]
    environment = os.environ.copy()
    environment["PYTHONPATH"] = os.pathsep.join((str(root), str(root / "src")))
    command = [sys.executable, "-m", "experiments.scalar_lift.joint_t_sequence", "--worker",
               "--prime", str(prime), "--multiplier", str(multiplier), "--seconds", str(solver_seconds)]
    if not orbit_normal_form:
        command.append("--no-orbit-normal-form")
    started = time.monotonic()
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
    if witness is not None and verify_t_states(witness["one_hot_orbit_states"], prime, multiplier) != witness:
        raise ValueError("the independent supervisor replay disagrees")
    return {
        "scope": "One bounded exact periodic one-hot quadruple search; a restricted joint GS design chart",
        "prime": prime, "multiplier_order": multiplier, "block_order": 4 * prime,
        "multiplicative_orbit_normal_form": orbit_normal_form,
        "algebraic_source_3_and_13_output_orders": [12 * prime, 52 * prime],
        "native_limit_seconds": solver_seconds, "wall_limit_seconds": wall_seconds,
        "elapsed_seconds": time.monotonic() - started,
        "worker_pid": worker.pid, "worker_exit_code": worker.returncode,
        "worker_reaped": worker.returncode is not None, "external_watchdog_triggered": timed_out,
        "events": events, "diagnostics": diagnostics,
        "solver_status": None if final is None else final["solver_status"],
        "exact_design_constructed": witness is not None,
        "inconclusive": timed_out or worker.returncode != 0 or final is None or final["solver_status"] not in ("OPTIMAL", "FEASIBLE", "INFEASIBLE"),
        "mathematical_nonexistence_certificate_produced": False,
        "new_real_Hadamard_order_claimed": False, "witness": witness,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--prime", type=int, default=223)
    parser.add_argument("--multiplier", type=int, default=3)
    parser.add_argument("--seconds", type=float, default=45)
    parser.add_argument("--wall-seconds", type=float, default=60)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--no-orbit-normal-form", action="store_true")
    args = parser.parse_args()
    if args.worker:
        solve_worker(args.seconds, args.prime, args.multiplier, not args.no_orbit_normal_form)
        return
    report = bounded_probe(args.wall_seconds, args.seconds, args.prime, args.multiplier, not args.no_orbit_normal_form)
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({key: report[key] for key in (
        "prime", "multiplier_order", "block_order", "algebraic_source_3_and_13_output_orders",
        "solver_status", "elapsed_seconds", "worker_reaped", "external_watchdog_triggered",
        "exact_design_constructed", "inconclusive",
    )}))


if __name__ == "__main__":
    main()

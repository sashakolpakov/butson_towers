"""Exact fixed-frame endpoints and finite root-gauge searches.

Failure concerns the stated frames/gauges only. Weighted scalar restrictions
and other frames remain open. A source Gram proof is a separate prerequisite.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from functools import cache
from pathlib import Path

import numpy as np

from experiments.scalar_lift.fifth_source_assembly import root_vectors


def frame_options(phase):
    """JSON metadata, without constructing any potentially large endpoint."""
    phase = int(phase)
    records = []
    if phase in (3, 6, 12):
        records.append((f"twelfth-{phase}", 4, [(3, 0)],
                        "src/hadamard/dense_lift.py"))
    if phase in (5, 10):
        records.append((f"fifth-{phase}-raw", 8, [(5, 0)],
                        "src/hadamard/dense_lift.py"))
    if phase in (24, 48, 96):
        records.append((f"cubic-{phase}", phase // 3, [(3, 2)],
                        "certificates/scalar-lift/dyadic-partial-lift.json"))
    if phase in (20, 40, 80):
        c = phase // 5
        records.extend((f"fifth-{phase}-{kind}", 2 * phase // 5, [(5, forbidden)],
                        "certificates/scalar-lift/fifth-dyadic-tower.json")
                       for kind, forbidden in (("raw", 0), ("shifted", -c % 5)))
    if phase == 120:
        records.extend((f"common-120-{kind}", 32, [(3, 2), (5, forbidden)],
                        "certificates/scalar-lift/phase120-common-base.json")
                       for kind, forbidden in (("raw", 0), ("shifted", 1)))
    return [{"id": name, "phase": phase, "width": width,
             "admitted": [e for e in range(phase)
                          if all(e % p != f for p, f in domains)],
             "coefficient_domains": [{"modulus": p, "forbidden": f}
                                     for p, f in domains],
             "certificate_path": certificate}
            for name, width, domains, certificate in records]


def _record(frame_id):
    parts = frame_id.split("-")
    try:
        phase = int(parts[1])
    except (IndexError, ValueError) as error:
        raise ValueError("unknown fixed frame") from error
    for record in frame_options(phase):
        if record["id"] == frame_id:
            return record
    raise ValueError("unknown fixed frame")


@cache
def _blocks(frame_id):
    record = _record(frame_id)
    phase = record["phase"]
    if frame_id.startswith("twelfth-"):
        blocks = _twelfth_blocks()[::12 // phase]
    elif frame_id.startswith("cubic-"):
        from experiments.scalar_lift.dyadic_partial_lift import phase_blocks
        blocks = phase_blocks(phase // 6)
    elif frame_id.startswith("fifth-"):
        from experiments.scalar_lift.fifth_dyadic_tower import phase_blocks
        blocks = (phase_blocks(20)[::20 // phase] if phase in (5, 10)
                  else phase_blocks(phase, frame_id.endswith("shifted")))
    else:
        from experiments.scalar_lift.phase120_common_base import phase_blocks
        blocks = phase_blocks(frame_id.endswith("shifted"))
    result = np.asarray(blocks, dtype=np.int64)
    result.setflags(write=False)
    return result


@cache
def _twelfth_blocks():
    from hadamard.dense_lift import QUATERNARY_REAL_BLOCKS, TWELFTH_POWER_QUATERNARY_BLOCKS

    known = {e: np.block([[np.array(QUATERNARY_REAL_BLOCKS[v], dtype=np.int64)
                          for v in row] for row in matrix])
             for e, matrix in TWELFTH_POWER_QUATERNARY_BLOCKS.items()}
    generator = known[2] @ known[1].T  # four times the orthogonal scalar action
    initial = generator.T @ known[1]
    if np.any(initial % 4):
        raise ArithmeticError("twelfth frame not integral")
    current = initial // 4
    blocks = []
    for _ in range(12):
        blocks.append(current)
        product = generator @ current
        if np.any(product % 4):
            raise ArithmeticError("twelfth orbit not integral")
        current = product // 4
    if not np.array_equal(current, blocks[0]):
        raise ArithmeticError("twelfth orbit did not close")
    if any(not np.array_equal(blocks[e], block) for e, block in known.items()):
        raise ArithmeticError("twelfth reconstruction differs from literal transfer")
    return blocks


def frame_block(frame_id, exponent):
    """Read-only integer block; excluded phases need not have sign entries."""
    return _blocks(frame_id)[int(exponent) % _record(frame_id)["phase"]]


@cache
def frame_proof(frame_id):
    """Certify the cyclotomic representation and frame using integer identities.

    If B_e = rho(zeta^e) F and FF^T=dI, C=B_1 F^T=d rho(zeta).
    C B_e=d B_(e+1), CC^T=d²I, and all cyclotomic linear relations
    imply every source zero sum maps to zero and conjugation to transpose.
    This proves all cross-Grams without enumerating q² matrix products.
    """
    record = _record(frame_id)
    blocks = _blocks(frame_id)
    phase, width = record["phase"], record["width"]
    identity = np.eye(width, dtype=np.int64)
    frame = blocks[0]
    if not np.array_equal(frame @ frame.T, width * identity):
        raise ArithmeticError("incorrect frame Gram")
    numerator = blocks[1] @ frame.T
    if not np.array_equal(numerator @ numerator.T, width**2 * identity):
        raise ArithmeticError("scalar generator is not orthogonal")
    for e in range(phase):
        if not np.array_equal(numerator @ blocks[e], width * blocks[(e + 1) % phase]):
            raise ArithmeticError("scalar orbit recurrence failed")
    vectors = root_vectors(phase)
    reconstructed = np.tensordot(vectors, blocks[:vectors.shape[1]], axes=1)
    if not np.array_equal(reconstructed, blocks):
        raise ArithmeticError("a cyclotomic relation failed")
    admitted = [e for e, block in enumerate(blocks) if np.all(np.abs(block) == 1)]
    if admitted != record["admitted"]:
        raise ArithmeticError("declared flat alphabet differs from the actual frame")
    return {"frame_id": frame_id, "phase": phase, "width": width,
            "exact_frame_Gram": True, "exact_scalar_recurrence": True,
            "exact_cyclotomic_relations": True, "exact_conjugation": True,
            "all_cross_Grams_follow": True, "admitted": admitted,
            "block_table_sha256": hashlib.sha256(blocks.astype("<i8").tobytes()).hexdigest(),
            "source_Gram_is_separate_prerequisite": True}


def check_frame(phase, used_alphabet, frame_id):
    record = _record(frame_id)
    if int(phase) != record["phase"]:
        raise ValueError("source phase differs from frame phase")
    used = sorted({int(e) % phase for e in used_alphabet})
    excluded = sorted(set(used) - set(record["admitted"]))
    return {"status": "ADMITTED" if not excluded else "EXCLUDED_ALPHABET",
            "fits": not excluded, "frame_id": frame_id, "width": record["width"],
            "used_alphabet": used, "excluded_used_exponents": excluded,
            "proof": frame_proof(frame_id),
            "scope": "this fixed frame and actual source alphabet only"}


def global_affine_fit(phase, used_alphabet, frame_id):
    """Complete finite search of global Galois automorphisms and root phases."""
    record = _record(frame_id)
    if int(phase) != record["phase"]:
        raise ValueError("source phase differs from frame phase")
    used = sorted({int(e) % phase for e in used_alphabet})
    allowed = set(record["admitted"])
    trials = 0
    for multiplier in range(1, phase):
        if math.gcd(multiplier, phase) != 1:
            continue
        for shift in range(phase):
            trials += 1
            transformed = sorted({(multiplier * e + shift) % phase for e in used})
            if set(transformed) <= allowed:
                return {"status": "FOUND", "frame_id": frame_id,
                        "multiplier": multiplier, "shift": shift,
                        "used_alphabet": transformed, "trials": trials}
    return {"status": "EXHAUSTED", "frame_id": frame_id, "trials": trials,
            "scope": "all global Galois/root-phase transformations for this frame"}


def _prime_gauges(colors, modulus, forbidden, seconds):
    """Bounded exact CP-SAT; an UNKNOWN result is never an exclusion."""
    from ortools.sat.python import cp_model

    rows, columns = colors.shape
    model = cp_model.CpModel()
    left = [model.new_int_var(0, modulus - 1, f"r{i}") for i in range(rows)]
    right = [model.new_int_var(0, modulus - 1, f"c{j}") for j in range(columns)]
    model.add(left[0] == 0)
    tables = {e: [(a, b) for a in range(modulus) for b in range(modulus)
                  if (e + a + b) % modulus != forbidden] for e in range(modulus)}
    for i in range(rows):
        for j in range(columns):
            model.add_allowed_assignments([left[i], right[j]], tables[int(colors[i, j])])
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = 1
    solver.parameters.max_time_in_seconds = max(0.001, seconds)
    status = solver.solve(model)
    report = {"modulus": modulus, "forbidden": forbidden,
              "status": solver.status_name(status), "wall_seconds": solver.wall_time,
              "workers": 1, "variables": rows + columns, "entry_constraints": rows * columns}
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        report["row_gauges"] = [solver.value(v) for v in left]
        report["column_gauges"] = [solver.value(v) for v in right]
        values = (colors + np.array(report["row_gauges"])[:, None]
                  + np.array(report["column_gauges"])[None, :]) % modulus
        if np.any(values == forbidden):
            raise ArithmeticError("solver gauge failed independent entry verification")
    return report


def row_column_gauge_search(exponents, phase, frame_id, *, seconds=5.0, max_entries=40000):
    """Search the complete root-row/column-gauge family of a dense source.

    Domains are complements of one color at each coprime odd prime. CRT
    therefore separates the constraints. Galois multipliers add no new
    feasibility: multiplication rescales gauges and a global translation
    restores the forbidden color. No source Gram assertion is made here.
    """
    record = _record(frame_id)
    source = np.asarray(exponents, dtype=np.int64)
    if int(phase) != record["phase"] or source.ndim != 2 or not all(source.shape):
        raise ValueError("a nonempty exponent matrix at the frame phase is required")
    if source.size > max_entries:
        return {"status": "NOT_RUN", "reason": "dense entry budget exceeded",
                "entries": source.size, "max_entries": max_entries}
    started = time.monotonic()
    reports = []
    for domain in record["coefficient_domains"]:
        remaining = seconds - (time.monotonic() - started)
        if remaining <= 0:
            return {"status": "UNKNOWN", "prime_reports": reports, "reason": "time budget"}
        report = _prime_gauges(source % domain["modulus"], domain["modulus"],
                               domain["forbidden"], remaining)
        reports.append(report)
        if report["status"] == "INFEASIBLE":
            return {"status": "INFEASIBLE", "prime_reports": reports,
                    "scope": "all row/column root gauges and Galois automorphisms into this frame"}
        if "row_gauges" not in report:
            return {"status": "UNKNOWN", "prime_reports": reports}
    odd = math.prod(report["modulus"] for report in reports)
    dyadic = phase // odd
    left = np.zeros(source.shape[0], dtype=np.int64)
    right = np.zeros(source.shape[1], dtype=np.int64)
    for report in reports:
        p = report["modulus"]
        crt = phase // p * pow(phase // p, -1, p)
        left += crt * np.array(report["row_gauges"], dtype=np.int64)
        right += crt * np.array(report["column_gauges"], dtype=np.int64)
    left %= phase
    right %= phase
    transformed = (source + left[:, None] + right[None, :]) % phase
    if not set(map(int, transformed.flat)) <= set(record["admitted"]):
        raise ArithmeticError("CRT gauge reconstruction failed")
    if np.any(left % dyadic) or np.any(right % dyadic):
        raise ArithmeticError("CRT gauge unexpectedly changes the dyadic coordinate")
    return {"status": "FOUND", "multiplier": 1, "row_gauges": left.tolist(),
            "column_gauges": right.tolist(), "exponents": transformed.tolist(),
            "used_alphabet": sorted(set(map(int, transformed.flat))),
            "prime_reports": reports, "source_Gram_preserved_algebraically": True,
            "Galois_choices_redundant_for_full_gauge_feasibility": True}


def character_cubic_obstruction(source):
    """Reuse a complete eight-point obstruction for the implicit high sources."""
    from experiments.scalar_lift.forward_character_sources import (
        cubic_gauge_census,
        literal_gauge_template_check,
    )
    checked = literal_gauge_template_check(source)
    census = cubic_gauge_census()
    gamma = checked["normalized_interior_cubic_color"]
    witnesses = [r for r in census if r["normalized_interior_cubic_color"] == gamma
                 and r["feasible_row_gauges"] == 0]
    if not witnesses:
        raise ArithmeticError("no complete template exclusion for this actual source")
    return {"status": "INFEASIBLE", "literal_source_check": checked,
            "finite_obstructions": witnesses,
            "scope": "all two-cubic-coset fixed frames under root row/column gauges and Galois"}


def verified_source_gauge_search(source, frame_id, *, seconds=5.0, max_entries=40000):
    """Adapter for forward_source_registry.VerifiedSource, including implicit inputs."""
    record = _record(frame_id)
    if source.phase != record["phase"] or source.proof.get("full_Gram_verified") is not True:
        raise ValueError("a verified source at the frame phase is required")
    binding = {"source_id": source.id, "source_hash": source.source_hash,
               "frame_id": frame_id, "order": source.order, "phase": source.phase}
    affine = global_affine_fit(source.phase, source.used_alphabet, frame_id)
    if affine["status"] == "FOUND":
        return {**binding, **affine, "gauge_kind": "global_affine",
                "frame_proof": frame_proof(frame_id)}
    if source.exponents is not None:
        return {**binding, **row_column_gauge_search(
            source.exponents, source.phase, frame_id, seconds=seconds, max_entries=max_entries
        ), "gauge_kind": "row_column"}
    cubic = any(d["modulus"] == 3 for d in record["coefficient_domains"])
    if cubic and source.id in ("character-73", "character-241"):
        from experiments.scalar_lift.forward_character_sources import character_source

        prime = int(source.id.split("-")[1])
        original = character_source(prime)
        for shape in ((4, 2), (2, 4)):
            points = [(i, j) for i in range(shape[0]) for j in range(shape[1])]
            for left in points:
                for right in points:
                    actual = source.entry(left[0] * prime + left[1], right[0] * prime + right[1])
                    if actual != original.entry(left, right):
                        raise ValueError("implicit source no longer matches its literal obstruction")
        return {**binding, **character_cubic_obstruction(original),
                "gauge_kind": "row_column", "dense_matrix_materialized": False}
    return {**binding, "status": "NOT_RUN", "gauge_kind": "row_column",
            "reason": "no certified implicit gauge reduction registered",
            "global_affine_result": affine}


def report_digest(report):
    return hashlib.sha256(json.dumps(report, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def certificate():
    from experiments.scalar_lift.fifth_source_assembly import exact_column_gram
    from experiments.scalar_lift.phase120_common_base import source_fifteen

    frames = [frame_proof(record["id"])
              for phase in (3, 5, 6, 10, 12, 20, 24, 40, 48, 80, 96, 120)
              for record in frame_options(phase)]
    source = np.asarray(source_fifteen(), dtype=np.int64)
    n = len(source)
    left, right = np.arange(n), np.arange(n)**2
    scrambled = (source + left[:, None] + right[None, :]) % 120
    global_result = global_affine_fit(120, scrambled.flat, "common-120-raw")
    if global_result["status"] != "EXHAUSTED":
        raise ArithmeticError("the row/column control no longer needs nonconstant gauges")
    gram = exact_column_gram(scrambled, 120)
    expected = np.zeros_like(gram)
    expected[np.arange(n), np.arange(n), 0] = n
    if not np.array_equal(gram, expected):
        raise ArithmeticError("scrambled source Gram failed")
    real = np.block([[frame_block("common-120-raw", e) for e in row] for row in source])
    if not np.all(np.abs(real) == 1) or not np.array_equal(
        real @ real.T, len(real) * np.eye(len(real), dtype=np.int64)
    ):
        raise ArithmeticError("actual complete sign endpoint failed")
    return {"schema": "forward-fixed-frames-v1", "frames": frames,
            "regauged_source_control": {
                "source_order": n, "phase": 120, "actual_source_Gram_verified": True,
                "global_affine_search": global_result,
                "explicit_inverse_row_gauges": ((-left) % 120).tolist(),
                "explicit_inverse_column_gauges": ((-right) % 120).tolist(),
                "source_sha256": hashlib.sha256(scrambled.astype("<i8").tobytes()).hexdigest(),
                "real_order": len(real), "full_integer_sign_Gram_verified": True},
            "scope": {
                "source_existence_is_not_proved_by_frame_admission": True,
                "unknown_solver_status_is_not_an_exclusion": True,
                "other_frames_and_weighted_nonflat_paths_remain_open": True,
                "new_unknown_Hadamard_order_claimed": False}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    result = certificate()
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"frame_count": len(result["frames"]),
                      "actual_sign_control_order": result["regauged_source_control"]["real_order"],
                      "report_sha256": report_digest(result)}, sort_keys=True))


if __name__ == "__main__":
    main()

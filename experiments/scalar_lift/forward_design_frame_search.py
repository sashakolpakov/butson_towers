"""Persistent final-frame search from an actual unreal BH(13,6).

The default finite chart partitions multiplier-invariant periodic T223
quadruples into 192 cells. A witness supplies two order-892 sign blocks,
and hence an actual forward scalar lift to order 11596. The source is
already complete; only the compatible scalar action/frame is unknown.
"""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import math
import sqlite3
import time
import zlib
from functools import cache
from itertools import product
from pathlib import Path

from experiments.scalar_lift.inverse_blocks import known_order_thirteen_source
from experiments.scalar_lift.joint_t_sequence import build_model, row_sum_types, verify_t_states
from hadamard.cgw import is_cgw_exponent_matrix

ROOT = Path(__file__).resolve().parents[2]
SCHEMA = "hadamard.forward-design-frame-search.v1"
REVISION = "fixed-unreal13-T223-row-sums-orbits2-3-v1"


def encoded(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def digest(value):
    return hashlib.sha256(encoded(value).encode()).hexdigest()


@cache
def source_certificate():
    """Replay both complete sources, including the optional order-3 output."""
    from experiments.scalar_lift.forward_source_registry import load_source

    registered = [load_source("unreal13-q6"), load_source("unreal3-q6")]
    matrices = [tuple(tuple((source.entry(i, j) + 1) % 6 for j in range(source.order))
                      for i in range(source.order)) for source in registered]
    if matrices[0] != known_order_thirteen_source():
        raise ValueError("the registry's inverse gauge differs from the independently verified source")
    records = []
    for matrix, registered_source in zip(matrices, registered, strict=True):
        if any(e not in (1, 2, 4, 5) for row in matrix for e in row):
            raise ValueError("the source is outside the actual pair's unreal alphabet")
        if not is_cgw_exponent_matrix(matrix, 6, len(matrix)):
            raise ValueError("the complete source failed its exact cyclotomic Gram")
        records.append({"order": len(matrix), "phase": 6, "exponents": matrix,
                        "source_sha256": digest(matrix), "full_Gram_verified": True,
                        "registered_source_id": registered_source.id,
                        "registered_source_sha256": registered_source.source_hash,
                        "registered_source_exponent_shift": 1,
                        "used_alphabet": sorted({e for row in matrix for e in row}),
                        "source_complete_before_frame_search": True})
    return records


def cells_for_chart(prime=223, multiplier=3):
    """Disjoint, exhaustive cells of the existing signed-axis normal form."""
    if 1 + (prime - 1) // multiplier < 4:
        raise ValueError("the orbit partition requires at least four orbits")
    return [{"prime": prime, "multiplier": multiplier, "row_sums": list(sums),
             "fixed_orbit_states": [[2, a], [3, b]]}
            for sums in row_sum_types(prime, multiplier) for a, b in product(range(8), repeat=2)]


@cache
def _base_model(prime, multiplier):
    model, states, metadata = build_model(prime, multiplier)
    sums = {v.name: i for i, v in enumerate(model.proto.variables) if v.name.startswith("sum")}
    if set(sums) != {f"sum{i}" for i in range(4)}:
        raise ValueError("the existing T model changed its row-sum interface")
    return model, tuple(v.index for v in states), sums, metadata


def build_cell(cell):
    base, indices, sums, metadata = _base_model(cell["prime"], cell["multiplier"])
    if tuple(cell["row_sums"]) not in row_sum_types(cell["prime"], cell["multiplier"]):
        raise ValueError("the cell has an invalid normalized row-sum type")
    model = base.clone()
    states = [model.get_int_var_from_proto_index(i) for i in indices]
    for i, value in enumerate(cell["row_sums"]):
        model.add(model.get_int_var_from_proto_index(sums[f"sum{i}"]) == value)
    if len(cell["fixed_orbit_states"]) != 2 or [i for i, _ in cell["fixed_orbit_states"]] != [2, 3]:
        raise ValueError("the partition must specify orbits two and three")
    for i, value in cell["fixed_orbit_states"]:
        if type(value) is not int or not 0 <= value < 8:
            raise ValueError("invalid orbit state")
        model.add(states[i] == value)
    error = model.validate()
    if error:
        raise ValueError(error)
    return model, states, metadata


def verify_candidate(cell, states):
    """Independent expanded correlations and complete sign-pair verification."""
    proof = verify_t_states(states, cell["prime"], cell["multiplier"])
    if proof["sequence_row_sums"] != cell["row_sums"] or any(states[i] != v for i, v in cell["fixed_orbit_states"]):
        raise ValueError("the verified witness does not belong to its declared cell")
    design = proof["full_GS_sign_pair_verification"]
    width = 4 * cell["prime"]
    frontier = json.loads((ROOT / "experiments/scalar_lift/frontier_120_snapshot.json").read_text())
    outputs = [{"order": source["order"] * width, "source_order": source["order"],
                "source_sha256": source["source_sha256"], "frame_width": width,
                "literal_sign_entries": True, "full_Gram_inherited_exactly": True,
                "listed_in_dated_primary_frontier": source["order"] * width in frontier["orders"],
                "new_order_priority_established": False} for source in source_certificate()]
    return {"cell": cell, "one_hot_orbit_states": states, "verification": proof,
            "a_sequences": design["a_sequences"], "b_sequences": design["b_sequences"],
            "sources": source_certificate(), "real_outputs": outputs,
            "representation": "rho(zeta6)=QP^T/d; frame=P-Q; blocks at exponents 1,2,4,5 are P,Q,-P,-Q",
            "P_definition": "GS(a-b)", "Q_definition": "GS(a+b)",
            "scalar_operator_and_frame_recovered_from_verified_pair": True}


def endpoint_entry(witness, row, column, source_order=13):
    """Literal sign access without allocating the order-11596 endpoint."""
    source = next(s for s in witness["sources"] if s["order"] == source_order)
    p = witness["cell"]["prime"]
    width, order = 4 * p, 4 * p * source_order
    if not 0 <= row < order or not 0 <= column < order:
        raise ValueError("endpoint coordinate outside its matrix")
    i, r = divmod(row, width)
    j, c = divmod(column, width)
    exponent = source["exponents"][i][j]
    outer_r, inner_r = divmod(r, p)
    outer_c, inner_c = divmod(c, p)
    # mode 0 is a circulant, -1 is B R, and +1 is B^T R.
    layout = (((0, 0, 1), (1, -1, 1), (2, -1, 1), (3, -1, 1)),
              ((1, -1, -1), (0, 0, 1), (3, 1, -1), (2, 1, 1)),
              ((2, -1, -1), (3, 1, 1), (0, 0, 1), (1, 1, -1)),
              ((3, -1, -1), (2, 1, -1), (1, 1, 1), (0, 0, 1)))
    sequence, mode, sign = layout[outer_r][outer_c]
    index = ((inner_c - inner_r) if mode == 0 else mode * (inner_r + inner_c + 1)) % p
    a = witness["a_sequences"][sequence][index]
    b = witness["b_sequences"][sequence][index]
    return sign * ({1: a-b, 2: a+b, 4: -a+b, 5: -a-b}[exponent])


class DesignFrameSearch:
    def __init__(self, database, max_db_mib=256, *, include_491=False, cell_seconds=15,
                 charts=None, cells=None):
        import ortools

        if not math.isfinite(max_db_mib) or max_db_mib < 2:
            raise ValueError("at least 2 MiB of database budget is required")
        if not math.isfinite(cell_seconds) or not 0 < cell_seconds <= 60:
            raise ValueError("each solver attempt must have a positive budget at most 60 seconds")
        self.path = Path(database)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.lock = Path(str(self.path) + ".lock").open("a+")  # noqa: SIM115 -- held until close()
        try:
            fcntl.flock(self.lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            self.lock.close()
            raise RuntimeError("another worker owns this design-frame checkpoint") from None
        self.db = sqlite3.connect(self.path)
        self.db.execute("PRAGMA journal_mode=DELETE")
        self.db.execute("PRAGMA cache_size=-2048")
        self.db.execute("PRAGMA temp_store=FILE")
        self.db.executescript("""
            CREATE TABLE IF NOT EXISTS meta(key TEXT PRIMARY KEY,value TEXT NOT NULL);
            CREATE TABLE IF NOT EXISTS cells(id TEXT PRIMARY KEY,recipe TEXT NOT NULL,status TEXT NOT NULL);
            CREATE TABLE IF NOT EXISTS attempts(cell_id TEXT,seed INTEGER,status TEXT NOT NULL,result TEXT,
                PRIMARY KEY(cell_id,seed));
            CREATE TABLE IF NOT EXISTS models(prime INTEGER,multiplier INTEGER,sha256 TEXT,proto_zlib BLOB,
                PRIMARY KEY(prime,multiplier));
            CREATE TABLE IF NOT EXISTS witnesses(cell_id TEXT PRIMARY KEY,record TEXT NOT NULL);
        """)
        page = self.db.execute("PRAGMA page_size").fetchone()[0]
        self.db.execute(f"PRAGMA max_page_count={int(max_db_mib * 1024**2) // (2 * page)}")
        charts = charts if charts is not None else [(223, 3)] + ([(491, 7)] if include_491 else [])
        catalogue = cells if cells is not None else [cell for p, h in charts for cell in cells_for_chart(p, h)]
        files = (Path(__file__), Path(__file__).with_name("joint_t_sequence.py"),
                 Path(__file__).with_name("joint_gs_223.py"))
        identity = {"schema": SCHEMA, "revision": REVISION, "cells": catalogue,
                    "sources": source_certificate(), "ortools_version": ortools.__version__,
                    "implementation_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}
        old = self._get("identity")
        if old is not None and encoded(old) != encoded(identity):
            self.close()
            raise ValueError("checkpoint belongs to another exact chart, source or implementation")
        with self.db:
            if old is None:
                self._put("identity", identity)
                self._put("progress", {"evaluations": 0, "last_attempt": None})
                for cell in catalogue:
                    key = digest(cell)
                    self.db.execute("INSERT INTO cells VALUES (?,?,'pending')", (key, encoded(cell)))
                    self.db.execute("INSERT INTO attempts VALUES (?,0,'pending',NULL)", (key,))
            # A killed CP process loses its native search tree. Retain the fact;
            # do not silently repeat the identical seed on resume.
            for key, seed, raw in self.db.execute("SELECT cell_id,seed,result FROM attempts WHERE status='running'").fetchall():
                result = {**(json.loads(raw) if raw else {}), "native_search_state_retained": False,
                          "explicit_new_seed_required": True}
                self.db.execute("UPDATE cells SET status='open' WHERE id=?", (key,))
                self.db.execute("UPDATE attempts SET status='interrupted',result=? WHERE cell_id=? AND seed=?",
                                (encoded(result), key, seed))
        self.cell_seconds = cell_seconds
        self.stop_reason = "not_started"

    def close(self):
        self.db.close()
        fcntl.flock(self.lock, fcntl.LOCK_UN)
        self.lock.close()

    def _get(self, key):
        row = self.db.execute("SELECT value FROM meta WHERE key=?", (key,)).fetchone()
        return json.loads(row[0]) if row else None

    def _put(self, key, value):
        self.db.execute("INSERT OR REPLACE INTO meta VALUES (?,?)", (key, encoded(value)))

    def retry_open(self, seed):
        """An explicit, recorded new seed; never retry a decided empty cell."""
        if type(seed) is not int or not 0 <= seed < 2**31:
            raise ValueError("CP-SAT seed must be an integer in 0..2^31-1")
        before = self.db.total_changes
        with self.db:
            self.db.execute("INSERT OR IGNORE INTO attempts SELECT id,?,'pending',NULL FROM cells WHERE status='open'", (seed,))
        return self.db.total_changes - before

    def _verify_saved_candidate(self, key, recipe, seed, result):
        witness = verify_candidate(recipe, result["orbit_states"])
        with self.db:
            self.db.execute("INSERT OR REPLACE INTO witnesses VALUES (?,?)", (key, encoded(witness)))
            self.db.execute("UPDATE cells SET status='verified' WHERE id=?", (key,))
            self.db.execute("UPDATE attempts SET status='verified' WHERE cell_id=? AND seed=?", (key, seed))
            self._put("progress", {**self._get("progress"),
                                   "cell_statuses": dict(self.db.execute("SELECT status,COUNT(*) FROM cells GROUP BY status"))})

    def _attempt(self, key, recipe, seed, seconds):
        from ortools.sat.python import cp_model

        model, states, metadata = build_cell(recipe)
        model_hash = hashlib.sha256(str(model.proto).encode()).hexdigest()
        parameters = {"max_time_in_seconds": seconds, "num_search_workers": 1,
                      "random_seed": seed, "max_presolve_iterations": 1}
        with self.db:
            if self.db.execute("SELECT 1 FROM models WHERE prime=? AND multiplier=?",
                               (recipe["prime"], recipe["multiplier"])).fetchone() is None:
                base = _base_model(recipe["prime"], recipe["multiplier"])[0]
                text = str(base.proto).encode()
                self.db.execute("INSERT INTO models VALUES (?,?,?,?)",
                                (recipe["prime"], recipe["multiplier"], hashlib.sha256(text).hexdigest(), zlib.compress(text)))
            self.db.execute("UPDATE attempts SET status='running',result=? WHERE cell_id=? AND seed=?",
                            (encoded({"parameters": parameters, "cell_model_sha256": model_hash}), key, seed))
            self.db.execute("UPDATE cells SET status='running' WHERE id=?", (key,))
            self._put("progress", {**self._get("progress"), "last_attempt": {"cell_id": key, "seed": seed, "cell": recipe}})
        solver = cp_model.CpSolver()
        for name, value in parameters.items():
            setattr(solver.parameters, name, value)
        status = solver.solve(model)
        name = solver.status_name(status)
        if status == cp_model.MODEL_INVALID:
            raise ValueError("CP-SAT rejected the validated cell model")
        result = {"solver_status": name, "solver_wall_seconds": solver.wall_time,
                  "conflicts": solver.num_conflicts, "branches": solver.num_branches,
                  "parameters": parameters, "cell_model_sha256": model_hash, "model_metadata": metadata,
                  "standalone_UNSAT_proof_exported": False}
        found = status in (cp_model.FEASIBLE, cp_model.OPTIMAL)
        if found:
            result["orbit_states"] = [int(solver.value(v)) for v in states]
        attempt_status = "candidate" if found else "infeasible" if status == cp_model.INFEASIBLE else "unknown"
        cell_status = "candidate" if found else "infeasible" if status == cp_model.INFEASIBLE else "open"
        with self.db:
            self.db.execute("UPDATE attempts SET status=?,result=? WHERE cell_id=? AND seed=?",
                            (attempt_status, encoded(result), key, seed))
            self.db.execute("UPDATE cells SET status=? WHERE id=?", (cell_status, key))
            progress = self._get("progress")
            progress["evaluations"] += 1
            progress["cell_statuses"] = dict(self.db.execute("SELECT status,COUNT(*) FROM cells GROUP BY status"))
            progress["last_solver_status"] = name
            self._put("progress", progress)
        if found:
            # Raw candidate is durable before the independent full verifier.
            self._verify_saved_candidate(key, recipe, seed, result)

    def run(self, seconds=30, max_attempts=None):
        if not math.isfinite(seconds) or seconds < 0 or (max_attempts is not None and max_attempts < 0):
            raise ValueError("finite nonnegative budgets are required")
        started, attempted = time.monotonic(), 0
        self.stop_reason = "time_budget"
        while time.monotonic() - started < seconds:
            if max_attempts is not None and attempted >= max_attempts:
                self.stop_reason = "attempt_budget"
                break
            row = self.db.execute("""SELECT c.id,c.recipe,a.seed,a.status,a.result FROM attempts a
                JOIN cells c ON c.id=a.cell_id WHERE a.status IN ('pending','candidate')
                AND c.status NOT IN ('verified','infeasible')
                ORDER BY CASE a.status WHEN 'candidate' THEN 0 ELSE 1 END,a.rowid LIMIT 1""").fetchone()
            if row is None:
                self.stop_reason = "scheduled_attempts_finished"
                break
            key, raw, seed, status, result = row
            try:
                if status == "candidate":
                    self._verify_saved_candidate(key, json.loads(raw), seed, json.loads(result))
                else:
                    remaining = seconds - (time.monotonic() - started)
                    self._attempt(key, json.loads(raw), seed, min(self.cell_seconds, remaining))
            except sqlite3.OperationalError as error:
                if getattr(error, "sqlite_errorcode", None) != sqlite3.SQLITE_FULL:
                    raise
                self.stop_reason = "database_quota"
                break
            attempted += 1
            if self.db.execute("SELECT COUNT(*) FROM witnesses").fetchone()[0]:
                self.stop_reason = "verified_frame_found"
                break
        return {**self.report(), "batch_attempts": attempted, "batch_seconds": time.monotonic()-started}

    def report(self):
        counts = dict(self.db.execute("SELECT status,COUNT(*) FROM cells GROUP BY status"))
        outputs = [output for (raw,) in self.db.execute("SELECT record FROM witnesses")
                   for output in json.loads(raw)["real_outputs"]]
        return {"schema": SCHEMA, "revision": REVISION, "stop_reason": self.stop_reason,
                "progress": self._get("progress"), "cell_statuses": counts,
                "attempt_statuses": dict(self.db.execute("SELECT status,COUNT(*) FROM attempts GROUP BY status")),
                "sources": source_certificate(), "verified_real_outputs": outputs,
                "scheduled_cells": sum(counts.values()), "unknown_cells_remain_open": True,
                "declared_chart_empty_by_solver": bool(counts) and set(counts) == {"infeasible"},
                "standalone_UNSAT_proof_exported": False,
                "all_possible_frames_excluded": False, "new_order_priority_established": False,
                "scope": "fixed complete BH13 source; multiplier-invariant one-hot T quadruples with exact row-sum and two-orbit partition; no general frame nonexistence claim",
                "sqlite_bytes": self.db.execute("PRAGMA page_count").fetchone()[0]*self.db.execute("PRAGMA page_size").fetchone()[0]}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--seconds", type=float, default=3600)
    parser.add_argument("--max-db-mib", type=int, default=256)
    parser.add_argument("--cell-seconds", type=float, default=15)
    parser.add_argument("--include-491", action="store_true")
    parser.add_argument("--retry-seed", type=int)
    parser.add_argument("--attempts", type=int)
    args = parser.parse_args(argv)
    search = DesignFrameSearch(args.database, args.max_db_mib, include_491=args.include_491,
                               cell_seconds=args.cell_seconds)
    try:
        if args.retry_seed is not None:
            search.retry_open(args.retry_seed)
        report = search.run(args.seconds, args.attempts)
    finally:
        search.close()
    args.report.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.report.with_suffix(args.report.suffix + ".tmp")
    temporary.write_text(json.dumps(report, indent=2) + "\n")
    temporary.replace(args.report)
    print(json.dumps({k: report[k] for k in ("stop_reason", "cell_statuses", "attempt_statuses", "verified_real_outputs")}))


if __name__ == "__main__":
    main()

"""Persistent source-first search of certified scalar descents and sign frames.

Only verified source matrices enter the queue. Small scalar blocks and lazy
recipes are retained; descendant matrices are never stored as a frontier.
General weighted-frame problems remain open records after registered maps
finish, and are never treated as nonexistence results.
"""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import math
import sqlite3
import time
from pathlib import Path

import numpy as np
from sympy import mobius, totient

from experiments.scalar_lift.forward_butson_descent import (
    available_steps,
    scalar_blocks,
    scalar_certificate,
)
from experiments.scalar_lift.vanishing_collision_search import integral_roots

SCHEMA = "forward-descent-search-v1"
REVISION = "verified-sources-certified-edges-frames-gauges-v2"
ROOT = Path(__file__).resolve().parents[2]


def encoded(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def digest(value):
    return hashlib.sha256(encoded(value).encode()).hexdigest()


def alphabet_after(alphabet, phase, prime, mode):
    blocks = scalar_blocks(phase, prime, mode)
    return sorted(map(int, np.unique(blocks[list(alphabet)])))


def affine_fit(alphabet, phase, admitted):
    """Find a global phase/Galois gauge; a failure is scoped to this family."""
    admitted = set(admitted)
    examined = 0
    for unit in range(1, phase):
        if math.gcd(unit, phase) != 1:
            continue
        for shift in range(phase):
            examined += 1
            transformed = sorted({(unit * e + shift) % phase for e in alphabet})
            if set(transformed) <= admitted:
                return {"found": True, "unit": unit, "shift": shift,
                        "transformed_alphabet": transformed, "gauges_examined": examined}
    return {"found": False, "gauges_examined": examined,
            "scope": "global root rotation and Galois automorphism for this fixed frame"}


def weighted_metric_record(phase):
    """Exact trace metric for an unnormalised regular representation."""
    degree = int(totient(phase))

    def ramanujan(value):
        quotient = phase // math.gcd(phase, value)
        return int(mobius(quotient)) * degree // int(totient(quotient))

    metric = [[ramanujan(i - j) for j in range(degree)] for i in range(degree)]
    roots = np.asarray(integral_roots(phase), dtype=np.int64)
    companion = roots[np.arange(1, degree + 1) % phase].T
    gram = np.asarray(metric, dtype=np.int64)
    if not np.array_equal(companion.T @ gram @ companion, gram):
        raise ArithmeticError("the actual regular action failed its trace metric identity")
    return {"field_phase": phase, "rational_restriction_degree": degree,
            "basis": "1,zeta,...,zeta^(phi(q)-1)", "trace_metric": metric,
            "metric_sha256": digest(metric),
            "adjoint_law": "rho(z)^T G = G rho(conjugate(z))",
            "weighted_Gram_law": "rho(H)(I tensor G^-1)rho(H)^T = n(I tensor G^-1)",
            "flat_frame_constructed_by_this_record": False}


def signed_triangle(alphabet, phase):
    """A signed odd relation excludes every uniform linear sign frame on it."""
    roots = integral_roots(phase)
    lookup = {roots[e]: e for e in alphabet}
    for a in alphabet:
        for b in alphabet:
            value = tuple(x + y for x, y in zip(roots[a], roots[b], strict=True))
            if value in lookup:
                return {"exponents": [a, b, lookup[value]], "coefficients": [1, 1, -1],
                        "exact_relation": "zeta^a+zeta^b-zeta^c=0",
                        "uniform_linear_raw_alphabet_sign_frame_excluded": True,
                        "proof": "an odd signed sum of literal signs cannot vanish",
                        "row_column_gauges_or_factor_dependent_maps_excluded": False}
            negative = tuple(-v for v in value)
            if negative in lookup:
                return {"exponents": [a, b, lookup[negative]], "coefficients": [1, 1, 1],
                        "exact_relation": "zeta^a+zeta^b+zeta^c=0",
                        "uniform_linear_raw_alphabet_sign_frame_excluded": True,
                        "proof": "an odd signed sum of literal signs cannot vanish",
                        "row_column_gauges_or_factor_dependent_maps_excluded": False}
    return None


def read_frontier():
    comparison = ROOT / "certificates/scalar-lift/forward-novelty-frontier-2026-11-22.json"
    if comparison.exists():
        record = json.loads(comparison.read_text())
        return {"orders": record["retained_source_listed_unresolved_orders"],
                "scope": record["scope"], "sources": record["sources"],
                "retrieved_utc_date": record["retrieved_utc_date"]}
    path = ROOT / "experiments/scalar_lift/frontier_120_snapshot.json"
    if not path.exists():
        return {"orders": [], "scope": "no comparison snapshot available"}
    return json.loads(path.read_text())


class ForwardSearch:
    def __init__(self, database, max_db_mib=256, catalogue=None):
        from experiments.scalar_lift.forward_frame_search import frame_options
        from experiments.scalar_lift.forward_source_registry import sources

        if not math.isfinite(max_db_mib) or max_db_mib < 2:
            raise ValueError("at least 2 MiB of database budget is required")
        self.path = Path(database)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.lock = Path(str(self.path) + ".lock").open("a+")  # noqa: SIM115
        try:
            fcntl.flock(self.lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            self.lock.close()
            raise RuntimeError("another forward worker owns this checkpoint") from None
        self.db = sqlite3.connect(self.path)
        self.db.execute("PRAGMA journal_mode=DELETE")
        self.db.execute("PRAGMA cache_size=-2048")
        self.db.execute("PRAGMA temp_store=FILE")
        self.db.executescript("""
            CREATE TABLE IF NOT EXISTS meta(key TEXT PRIMARY KEY,value TEXT NOT NULL);
            CREATE TABLE IF NOT EXISTS sources(id TEXT PRIMARY KEY,descriptor TEXT NOT NULL,proof TEXT);
            CREATE TABLE IF NOT EXISTS states(id TEXT PRIMARY KEY,recipe TEXT NOT NULL,status TEXT NOT NULL,result TEXT);
            CREATE TABLE IF NOT EXISTS endpoints(id TEXT PRIMARY KEY,state_id TEXT NOT NULL,record TEXT NOT NULL);
            CREATE TABLE IF NOT EXISTS edges(parent TEXT NOT NULL,child TEXT NOT NULL,PRIMARY KEY(parent,child));
        """)
        catalogue = sources() if catalogue is None else catalogue
        identity = {"schema": SCHEMA, "revision": REVISION, "catalogue": catalogue,
                    "frame_catalogue": [option for phase in range(2, 121)
                                        for option in frame_options(phase)]}
        if self._get("identity") not in (None, identity):
            self.close()
            raise ValueError("forward search checkpoint has another recipe/catalogue identity")
        if self._get("identity") is None:
            with self.db:
                self._put("identity", identity)
                self._put("progress", {"evaluations": 0, "verified_sources": 0,
                                       "registered_branches_finished": False})
                for source in catalogue:
                    self.db.execute("INSERT INTO sources VALUES (?,?,NULL)", (source["id"], encoded(source)))
                    self._enqueue({"kind": "source", "source_id": source["id"], "steps": []})
        page = self.db.execute("PRAGMA page_size").fetchone()[0]
        self.db.execute(f"PRAGMA max_page_count={int(max_db_mib * 1024**2) // (2 * page)}")
        self.frontier = read_frontier()
        self.stop_reason = "not_started"

    def close(self):
        self.db.close()
        fcntl.flock(self.lock, fcntl.LOCK_UN)
        self.lock.close()

    def _put(self, key, value):
        self.db.execute("INSERT OR REPLACE INTO meta VALUES (?,?)", (key, encoded(value)))

    def _get(self, key):
        row = self.db.execute("SELECT value FROM meta WHERE key=?", (key,)).fetchone()
        return json.loads(row[0]) if row else None

    def _enqueue(self, recipe, parent=None):
        identifier = digest(recipe)
        self.db.execute("INSERT OR IGNORE INTO states VALUES (?,?,'pending',NULL)", (identifier, encoded(recipe)))
        if parent:
            self.db.execute("INSERT OR IGNORE INTO edges VALUES (?,?)", (parent, identifier))
        return identifier

    def source(self, identifier):
        from experiments.scalar_lift.forward_source_registry import load_source

        descriptor = json.loads(self.db.execute("SELECT descriptor FROM sources WHERE id=?", (identifier,)).fetchone()[0])
        return load_source(descriptor)

    def _source_state(self, recipe, identifier):
        source = self.source(recipe["source_id"])
        if source.proof.get("full_Gram_verified") is not True:
            raise ValueError("a source without a Gram proof cannot enter descent")
        proof = {"source_id": recipe["source_id"], "order": source.order, "phase": source.phase,
                 "source_hash": source.source_hash, "used_alphabet": list(source.used_alphabet),
                 "proof": source.proof}
        self.db.execute("UPDATE sources SET proof=? WHERE id=?", (encoded(proof), recipe["source_id"]))
        self._enqueue({"kind": "descent", "source_id": recipe["source_id"], "steps": [],
                       "phase": source.phase, "order": source.order,
                       "alphabet": list(source.used_alphabet), "source_hash": source.source_hash}, identifier)
        return "expanded", proof

    def _endpoint(self, state_id, recipe, frame_id, gauge, width, frame_proof):
        order = recipe["order"] * width
        row = {"source_id": recipe["source_id"], "source_hash": recipe["source_hash"],
               "steps": recipe["steps"], "frame_id": frame_id, "affine_gauge": gauge,
               "order": order, "frame_width": width, "literal_sign_entries": True,
               "exact_Gram_inherited_from_verified_source_and_scalar_edges": True,
               "frame_proof": frame_proof,
               "matches_comparison_frontier": order in self.frontier.get("orders", []),
               "novelty_status": "requires literature audit; verified construction alone is not priority evidence",
               "matrix_representation": "lazy source/edge/frame entry recipe"}
        self.db.execute("INSERT OR IGNORE INTO endpoints VALUES (?,?,?)", (digest(row), state_id, encoded(row)))
        return {"frame_id": frame_id, "real_order": order, "gauge": gauge}

    def _descent_state(self, recipe, identifier):
        from experiments.scalar_lift.forward_frame_search import (
            frame_options,
            frame_proof,
            verified_source_gauge_search,
        )

        phase, alphabet = recipe["phase"], recipe["alphabet"]
        accepted, attempted = [], []
        if phase == 2:
            accepted.append(self._endpoint(identifier, recipe, "identity_sign", {}, 1,
                                          {"proof": "phase2 consists of the literal signs"}))
        for option in frame_options(phase):
            gauge = affine_fit(alphabet, phase, option["admitted"])
            attempted.append({"frame_id": option["id"], **gauge})
            if gauge["found"]:
                proof = frame_proof(option["id"])
                accepted.append(self._endpoint(identifier, recipe, option["id"], gauge, option["width"], proof))
            elif not recipe["steps"]:
                regauging = verified_source_gauge_search(self.source(recipe["source_id"]), option["id"], seconds=3)
                # Only gauge certificates belong in the graph, not dense transformed matrices.
                regauging.pop("exponents", None)
                attempted[-1]["row_column_gauge_search"] = regauging
                if regauging["status"] == "FOUND":
                    gauge = {"kind": "row_column", "unit": 1, "shift": 0,
                             "row_gauges": regauging["row_gauges"],
                             "column_gauges": regauging["column_gauges"]}
                    accepted.append(self._endpoint(identifier, recipe, option["id"], gauge,
                                                    option["width"], frame_proof(option["id"])))
        generated = []
        for prime in available_steps(phase):
            seen_blocks = set()
            for mode in ("primary", "primitive"):
                blocks = scalar_blocks(phase, prime, mode)
                block_hash = hashlib.sha256(blocks.tobytes()).hexdigest()
                if block_hash in seen_blocks:
                    continue
                seen_blocks.add(block_hash)
                certificate = scalar_certificate(phase, prime, mode)
                edge = {"upper_phase": phase, "prime": prime, "mode": mode,
                        "scalar_certificate_sha256": digest(certificate)}
                child = {**recipe, "steps": [*recipe["steps"], edge], "phase": phase // prime,
                         "order": recipe["order"] * prime,
                         "alphabet": alphabet_after(alphabet, phase, prime, mode)}
                self._enqueue(child, identifier)
                generated.append(edge)
        if phase > 2:
            self._enqueue({"kind": "weighted_open", "source_id": recipe["source_id"],
                           "steps": recipe["steps"], "phase": phase, "order": recipe["order"],
                           "alphabet": alphabet, "source_hash": recipe["source_hash"]}, identifier)
        return "expanded", {"partial_frames_examined": attempted, "real_endpoints": accepted,
                            "certified_edges_generated": generated,
                            "other_scalar_bases_or_frames_excluded": False}

    def _weighted_state(self, recipe):
        return "open", {"metric": weighted_metric_record(recipe["phase"]),
                        "signed_odd_relation_gate": signed_triangle(recipe["alphabet"], recipe["phase"]),
                        "status": "requires_general_final_frame_or_another_source_gauge",
                        "weighted_intermediates_retained": True,
                        "general_nonexistence_claimed": False}

    def evaluate(self, recipe, identifier):
        if recipe["kind"] == "source":
            return self._source_state(recipe, identifier)
        if recipe["kind"] == "descent":
            return self._descent_state(recipe, identifier)
        if recipe["kind"] == "weighted_open":
            return self._weighted_state(recipe)
        raise ValueError("unregistered forward recipe")

    def run(self, seconds=30, max_steps=None):
        if not math.isfinite(seconds) or seconds < 0 or (max_steps is not None and max_steps < 0):
            raise ValueError("finite nonnegative budgets are required")
        started, steps = time.monotonic(), 0
        self.stop_reason = "time_budget"
        while time.monotonic() - started < seconds:
            if max_steps is not None and steps >= max_steps:
                self.stop_reason = "step_budget"
                break
            row = self.db.execute("SELECT id,recipe FROM states WHERE status='pending' ORDER BY rowid LIMIT 1").fetchone()
            if row is None:
                with self.db:
                    progress = self._get("progress")
                    progress["registered_branches_finished"] = True
                    self._put("progress", progress)
                self.stop_reason = "registered_branches_finished"
                break
            identifier, raw = row
            recipe = json.loads(raw)
            try:
                with self.db:
                    status, result = self.evaluate(recipe, identifier)
                    self.db.execute("UPDATE states SET status=?,result=? WHERE id=?", (status, encoded(result), identifier))
                    progress = self._get("progress")
                    progress["evaluations"] += 1
                    progress["verified_sources"] += int(recipe["kind"] == "source")
                    progress["last_evaluated"] = {"kind": recipe["kind"], "source_id": recipe["source_id"]}
                    self._put("progress", progress)
            except sqlite3.OperationalError as error:
                if getattr(error, "sqlite_errorcode", None) != sqlite3.SQLITE_FULL:
                    raise
                self.stop_reason = "database_quota"
                break
            steps += 1
        return {**self.report(), "batch_steps": steps, "batch_seconds": round(time.monotonic() - started, 6)}

    def report(self):
        endpoints = [json.loads(row[0]) for row in self.db.execute("SELECT record FROM endpoints")]
        unique = sorted({row["order"] for row in endpoints})
        return {"schema": SCHEMA, "revision": REVISION, "stop_reason": self.stop_reason,
                "progress": self._get("progress"),
                "states": dict(self.db.execute("SELECT status,COUNT(*) FROM states GROUP BY status")),
                "registered_sources": self.db.execute("SELECT COUNT(*) FROM sources").fetchone()[0],
                "verified_endpoint_recipes": len(endpoints), "distinct_real_orders": unique,
                "comparison_frontier_intersection": sorted(set(unique) & set(self.frontier.get("orders", []))),
                "comparison_frontier_source": {k: v for k, v in self.frontier.items() if k != "orders"},
                "endpoint_examples": endpoints[:16],
                "all_endpoint_recipes_stored_in": "SQLite endpoints table",
                "sqlite_bytes": self.db.execute("PRAGMA page_count").fetchone()[0] * self.db.execute("PRAGMA page_size").fetchone()[0],
                "new_order_priority_established": False,
                "general_weighted_frame_search_exhausted": False,
                "branch_scope": "registered verified sources, repeated-prime edges in both distinct generator conventions, registered partial frames under global Galois and phase gauges"}


def endpoint_entry(database, endpoint_id, row, column):
    """Read a literal sign from a persisted endpoint without dense expansion."""
    from experiments.scalar_lift.forward_butson_descent import block_entry
    from experiments.scalar_lift.forward_frame_search import frame_block, frame_proof
    from experiments.scalar_lift.forward_source_registry import load_source

    with sqlite3.connect(Path(database).resolve().as_uri() + "?mode=ro", uri=True) as db:
        record = json.loads(db.execute("SELECT record FROM endpoints WHERE id=?", (endpoint_id,)).fetchone()[0])
        descriptor = json.loads(db.execute("SELECT descriptor FROM sources WHERE id=?", (record["source_id"],)).fetchone()[0])
    if not 0 <= row < record["order"] or not 0 <= column < record["order"]:
        raise ValueError("endpoint coordinate outside matrix")
    width = record["frame_width"]
    row, r = divmod(row, width)
    column, s = divmod(column, width)
    gauge_row, gauge_column = row, column
    digits = []
    for edge in reversed(record["steps"]):
        row, a = divmod(row, edge["prime"])
        column, b = divmod(column, edge["prime"])
        digits.append((edge, a, b))
    source = load_source(descriptor)
    if source.source_hash != record["source_hash"]:
        raise ValueError("source identity changed since endpoint certification")
    exponent, phase = source.entry(row, column), source.phase
    for edge, a, b in reversed(digits):
        certificate = scalar_certificate(phase, edge["prime"], edge["mode"])
        if digest(certificate) != edge["scalar_certificate_sha256"]:
            raise ValueError("scalar edge changed since endpoint certification")
        exponent = block_entry(exponent, a, b, phase, edge["prime"], edge["mode"])
        phase //= edge["prime"]
    if record["frame_id"] == "identity_sign":
        return 1 - 2 * exponent
    if frame_proof(record["frame_id"]) != record["frame_proof"]:
        raise ValueError("frame changed since endpoint certification")
    gauge = record["affine_gauge"]
    exponent = (gauge["unit"] * exponent + gauge["shift"]) % phase
    if gauge.get("kind") == "row_column":
        exponent = (exponent + gauge["row_gauges"][gauge_row] + gauge["column_gauges"][gauge_column]) % phase
    return int(frame_block(record["frame_id"], exponent)[r, s])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--seconds", type=float, default=30)
    parser.add_argument("--max-db-mib", type=int, default=256)
    parser.add_argument("--steps", type=int)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args(argv)
    search = ForwardSearch(args.database, args.max_db_mib)
    try:
        report = search.run(args.seconds, args.steps)
    finally:
        search.close()
    args.report.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.report.with_suffix(args.report.suffix + ".tmp")
    temporary.write_text(json.dumps(report, indent=2) + "\n")
    temporary.replace(args.report)
    print(json.dumps({k: report[k] for k in ("progress", "states", "distinct_real_orders", "comparison_frontier_intersection")}))


if __name__ == "__main__":
    main()

"""Finite breadth-first scalar restrictions of fixed verified sources.

No source-construction operation or solver appears in this tree. Every node
is a completed matrix, represented by a self-contained literal-entry bundle.
Repeated-prime scalar restrictions are the only child operations.
The default catalogue has four seeds; an explicit loader admits further
fixed inputs and binds every input hash into the checkpoint identity.
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

from experiments.scalar_lift.forward_butson_descent import available_steps

SEEDS = ("character-73", "character-241", "conference-1000-q8", "conference-1000-q16")
SCHEMA = "hadamard.completed-restriction-tree.v1"


def encoded(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def digest(value):
    return hashlib.sha256(encoded(value).encode()).hexdigest()


def fixed_source(identifier):
    """Load an already verified seed; this is not a construction search."""
    if identifier not in SEEDS:
        raise ValueError("this tree is restricted to its four existing sources")
    if identifier.startswith("character-"):
        from experiments.scalar_lift.forward_character_sources import (
            character_source,
            formal_gram_check,
        )

        source = character_source(int(identifier.split("-")[1]))
        proof = {**formal_gram_check(source), "full_Gram_verified": True,
                 "order": source.order, "phase": source.phase,
                 "coefficient_array_sha256": source.metadata()["coefficient_array_sha256"]}
    else:
        from experiments.scalar_lift.forward_source_registry import load_source

        source = load_source(identifier)
        proof = source.proof
    if proof.get("full_Gram_verified") is not True:
        raise ValueError("an incomplete source cannot enter the restriction tree")
    return source, proof


def write_json(path, record):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(record, indent=2) + "\n")
    temporary.replace(path)


class RestrictionTree:
    def __init__(self, database, max_db_mib=256, *, seeds=SEEDS, source_loader=None):
        if not math.isfinite(max_db_mib) or max_db_mib < 4:
            raise ValueError("at least four MiB of storage budget is required")
        if (not seeds or len(set(seeds)) != len(seeds)
                or (source_loader is None and any(seed not in SEEDS for seed in seeds))):
            raise ValueError("use a nonempty subset of the fixed existing seeds")
        self.source_loader = fixed_source if source_loader is None else source_loader
        self.path = Path(database).resolve()
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.artifacts = Path(str(self.path) + ".matrices")
        self.artifacts.mkdir(exist_ok=True)
        self.lock = Path(str(self.path) + ".lock").open("a+")  # noqa: SIM115 -- checkpoint lifetime
        try:
            fcntl.flock(self.lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            self.lock.close()
            raise RuntimeError("another worker owns this restriction tree") from None
        self.db = sqlite3.connect(self.path)
        self.db.execute("PRAGMA journal_mode=DELETE")
        self.db.execute("PRAGMA cache_size=-2048")
        self.db.execute("PRAGMA temp_store=FILE")
        self.db.executescript("""
            CREATE TABLE IF NOT EXISTS meta(key TEXT PRIMARY KEY,value TEXT NOT NULL);
            CREATE TABLE IF NOT EXISTS tasks(id TEXT PRIMARY KEY,recipe TEXT NOT NULL,depth INTEGER,status TEXT NOT NULL);
            CREATE TABLE IF NOT EXISTS nodes(id TEXT PRIMARY KEY,record TEXT NOT NULL);
            CREATE TABLE IF NOT EXISTS edges(task_id TEXT PRIMARY KEY,parent TEXT,child TEXT,label TEXT,duplicate INTEGER);
        """)
        page = self.db.execute("PRAGMA page_size").fetchone()[0]
        self.db.execute(f"PRAGMA max_page_count={int(max_db_mib * 1024**2) // (4 * page)}")
        self.max_bytes = int(max_db_mib * 1024**2)
        self.source_cache = {}
        identity = {"schema": SCHEMA, "seeds": list(seeds),
                    "children": "all repeated-prime restrictions, primary and primitive generator conventions",
                    "source_construction_branches": False}
        if source_loader is not None:
            from experiments.scalar_lift.large_butson_descendants import _prepare_source

            manifest = []
            for identifier in seeds:
                source, proof = source_loader(identifier)
                _, descriptor, _, source_hash = _prepare_source(source, proof)
                manifest.append({"id": identifier, "order": descriptor["order"],
                                 "phase": descriptor["phase"], "source_hash": source_hash})
                self.source_cache[identifier] = source, proof
            identity["verified_input_manifest"] = manifest
        previous = self._get("identity")
        if previous is not None and previous != identity:
            self.close()
            raise ValueError("this checkpoint has another fixed source catalogue")
        with self.db:
            if previous is None:
                self._put("identity", identity)
                self._put("progress", {"evaluated_recipes": 0, "completed_matrices": 0, "duplicate_recipes": 0})
                for identifier in seeds:
                    self._enqueue({"source": identifier, "steps": [], "parent": None})
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

    def _enqueue(self, recipe):
        self.db.execute("INSERT OR IGNORE INTO tasks VALUES (?,?,?,'pending')",
                        (digest(recipe), encoded(recipe), len(recipe["steps"])))

    def storage_bytes(self):
        total = sum(p.stat().st_size for p in self.artifacts.rglob("*") if p.is_file())
        total += sum(p.stat().st_size for p in self.path.parent.glob(self.path.name + "*")
                     if p.is_file() and not p.name.endswith(".lock"))
        return total

    def _evaluate(self, task_id, recipe):
        from experiments.scalar_lift.large_butson_descendants import ImplicitForwardTower

        identifier = recipe["source"]
        if identifier not in self.source_cache:
            self.source_cache[identifier] = self.source_loader(identifier)
        source, proof = self.source_cache[identifier]
        tower = ImplicitForwardTower(source, steps=recipe["steps"], source_proof=proof)
        key = tower.matrix_hash
        duplicate = self.db.execute("SELECT 1 FROM nodes WHERE id=?", (key,)).fetchone() is not None
        if not duplicate:
            bundle = self.artifacts / (key + ".npz")
            proof_path = self.artifacts / (key + ".json")
            tower.save_bundle(bundle)
            record = {"matrix_hash": key, "source": identifier, "source_hash": tower.source_hash,
                      "order": tower.order, "phase": tower.phase, "depth": len(recipe["steps"]),
                      "parent": recipe["parent"], "steps": recipe["steps"],
                      "used_alphabet": list(tower.used_alphabet),
                      "literal_entry_definition": tower.metadata(), "exact_Gram_proof": tower.proof(),
                      "bundle": str(bundle), "bundle_sha256": hashlib.sha256(bundle.read_bytes()).hexdigest(),
                      "bundle_bytes": bundle.stat().st_size, "proof_file": str(proof_path),
                      "dense_matrix_allocated": False, "real_sign_matrix": tower.phase == 2,
                      "other_field_restrictions_open": tower.phase > 2,
                      "new_Butson_order_claimed": False, "new_real_order_claimed": False,
                      "uniform_raw_final_sign_obstruction_from_source": identifier.startswith("character-")}
            write_json(proof_path, record)
            if self.storage_bytes() >= self.max_bytes:
                bundle.unlink(missing_ok=True)
                proof_path.unlink(missing_ok=True)
                raise OverflowError("completed-node files would exceed storage budget")
        with self.db:
            if not duplicate:
                self.db.execute("INSERT INTO nodes VALUES (?,?)", (key, encoded(record)))
                for prime in available_steps(tower.phase):
                    for mode in ("primary", "primitive"):
                        child = {"source": identifier, "parent": key,
                                 "steps": [*recipe["steps"], {"prime": prime, "generator_mode": mode}]}
                        self._enqueue(child)
            label = "existing verified source" if not recipe["steps"] else encoded(recipe["steps"][-1])
            self.db.execute("INSERT INTO edges VALUES (?,?,?,?,?)", (task_id, recipe["parent"], key, label, int(duplicate)))
            self.db.execute("UPDATE tasks SET status='complete' WHERE id=?", (task_id,))
            progress = self._get("progress")
            progress["evaluated_recipes"] += 1
            progress["completed_matrices"] += int(not duplicate)
            progress["duplicate_recipes"] += int(duplicate)
            progress["last_matrix"] = {"order": tower.order, "phase": tower.phase, "depth": len(recipe["steps"])}
            self._put("progress", progress)

    def run(self, seconds=60, max_nodes=None):
        if not math.isfinite(seconds) or seconds < 0 or (max_nodes is not None and max_nodes < 0):
            raise ValueError("finite nonnegative budgets are required")
        started, evaluated = time.monotonic(), 0
        self.stop_reason = "time_budget"
        while time.monotonic() - started < seconds:
            if max_nodes is not None and evaluated >= max_nodes:
                self.stop_reason = "node_budget"
                break
            row = self.db.execute("SELECT id,recipe FROM tasks WHERE status='pending' ORDER BY depth,rowid LIMIT 1").fetchone()
            if row is None:
                self.stop_reason = "registered_restriction_tree_complete"
                break
            if self.storage_bytes() >= self.max_bytes:
                self.stop_reason = "storage_budget"
                break
            try:
                self._evaluate(row[0], json.loads(row[1]))
            except OverflowError:
                self.stop_reason = "storage_budget"
                break
            except sqlite3.OperationalError as error:
                if getattr(error, "sqlite_errorcode", None) != sqlite3.SQLITE_FULL:
                    raise
                self.stop_reason = "storage_budget"
                break
            evaluated += 1
        return {**self.report(), "batch_evaluated_recipes": evaluated, "batch_seconds": time.monotonic()-started}

    def report(self):
        nodes = [json.loads(raw) for (raw,) in self.db.execute("SELECT record FROM nodes ORDER BY rowid")]
        return {"schema": SCHEMA, "stop_reason": self.stop_reason, "progress": self._get("progress"),
                "fixed_existing_seeds": list(self._get("identity")["seeds"]),
                "only_child_operations": "verified repeated-prime scalar restrictions",
                "coordinate_convention": "normalized scalar maps in distinct input/output bases, with the norm-degree Fourier frame recorded in each edge proof",
                "same_basis_unital_restriction_claimed": False,
                "nodes": nodes,
                "edges": [{"parent": a, "child": b, "operation": c, "duplicate_pruned": bool(d)}
                          for a, b, c, d in self.db.execute("SELECT parent,child,label,duplicate FROM edges ORDER BY rowid")],
                "pending_recipes": self.db.execute("SELECT COUNT(*) FROM tasks WHERE status='pending'").fetchone()[0],
                "completed_order_phase_pairs": sorted({(r["order"], r["phase"]) for r in nodes}),
                "real_output_orders": sorted({r["order"] for r in nodes if r["real_sign_matrix"]}),
                "storage_bytes": self.storage_bytes(),
                "duplicate_scope": "identical complete source identity and composed scalar block table; not full Hadamard equivalence",
                "new_source_construction_search": False, "SAT_search_performed": False,
                "all_scalar_restriction_methods_exhausted": False,
                "unimplemented_open_operations": "other cyclotomic subfields and retained-metric, non-Butson intermediate matrices; absence of a repeated-prime child does not exclude these paths",
                "new_real_Hadamard_order_claimed": False}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--seconds", type=float, default=60)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--max-db-mib", type=int, default=256)
    parser.add_argument("--nodes", type=int)
    inputs = parser.add_mutually_exclusive_group()
    inputs.add_argument("--all-verified", action="store_true",
                        help="use the complete fixed input catalogue from the phase-independent audit")
    inputs.add_argument("--expanded-only", action="store_true",
                        help="use the 94 separately acquired higher-alphabet gauge representatives")
    inputs.add_argument("--source", action="append", dest="sources",
                        help="select a fixed verified input; repeat for several inputs")
    args = parser.parse_args(argv)
    options = {}
    if args.all_verified or args.expanded_only or args.sources:
        from experiments.scalar_lift.phase_independent_frontier import all_sources

        available = {source.id: source for source in
                     all_sources(include_expanded=args.expanded_only or bool(args.sources))}
        chosen = (tuple(available) if args.all_verified else
                  tuple(identifier for identifier in available if identifier.startswith("expanded-"))
                  if args.expanded_only else tuple(args.sources))
        if any(identifier not in available for identifier in chosen):
            parser.error("an input is absent from the fixed verified catalogue")

        def load_selected(identifier):
            source = available[identifier]
            return source, source.proof

        options = {"seeds": chosen, "source_loader": load_selected}
    tree = RestrictionTree(args.database, args.max_db_mib, **options)
    try:
        report = tree.run(args.seconds, args.nodes)
    finally:
        tree.close()
    write_json(args.report, report)
    print(json.dumps({k: report[k] for k in ("stop_reason", "progress", "completed_order_phase_pairs", "real_output_orders", "storage_bytes")}))


if __name__ == "__main__":
    main()

"""Pure scalar restrictions through cyclotomic divisor fields, with metrics.

Starting matrices are fixed, verified inputs. Different
tower-coordinate bases survive unless their literal action tables and metrics
are identical. A rational weighted endpoint is not called a sign Hadamard.
"""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import math
import sqlite3
import time
from collections import Counter
from fractions import Fraction
from functools import cache
from pathlib import Path

import numpy as np

from experiments.scalar_lift.large_butson_factory import digest, encoded, fixed_source, write_json
from experiments.scalar_lift.relative_scalar_restriction import (
    RationalArray,
    conjugate_transpose,
    edge,
    regular_matrix,
    restriction_in_basis,
)
from experiments.scalar_lift.vanishing_collision_search import integral_roots

SEEDS = ("character-73", "character-241")
SCHEMA = "hadamard.weighted-scalar-tree.v1"
ROOT = Path(__file__).resolve().parents[2]


def proper_fields(phase):
    """Canonical full cyclotomic divisor fields; phase 1 denotes Q."""
    return [q for q in range(1, phase) if phase % q == 0 and q % 4 != 2
            and len(integral_roots(q)[0]) < len(integral_roots(phase)[0])]


def _array_hash(array):
    array = np.ascontiguousarray(array, dtype="<i8")
    return hashlib.sha256(array.tobytes()).hexdigest()


def retained_registry_seeds():
    path = ROOT / "certificates/scalar-lift/scalar-restriction-all-sources.json"
    return tuple(row["source_id"] for row in json.loads(path.read_text())["source_audits"]
                 if not row["alphabet_gate"]["excluded"])


def additional_input_seeds():
    from experiments.scalar_lift.catalogue_source_acquisition import SOURCE_IDS

    return (*SOURCE_IDS, "saved-unreal9-q6")


@cache
def load_tree_source(identifier):
    if identifier in SEEDS:
        return fixed_source(identifier)
    if identifier.startswith("catalogue-"):
        from experiments.scalar_lift.catalogue_source_acquisition import load_source
    elif identifier == "saved-unreal9-q6":
        from experiments.scalar_lift.nine_source_acquisition import load_source
    else:
        from experiments.scalar_lift.forward_source_registry import load_source

    source = load_source(identifier)
    if source.exponents is None or source.proof.get("full_Gram_verified") is not True:
        raise ValueError("the retained source must be a complete verified registry matrix")
    if source.proof.get("exponent_sha256") != _array_hash(source.exponents):
        raise ValueError("the dense source proof is not bound to its literal exponents")
    return source, source.proof


def source_array(source):
    return source.coefficients if hasattr(source, "coefficients") else source.exponents


def source_metadata(source):
    if hasattr(source, "coefficients"):
        return source.metadata()
    return {"kind": "dense_exponents", "order": source.order, "phase": source.phase,
            "coefficient_array_sha256": _array_hash(source.exponents),
            "registry_source_id": source.id, "registry_source_hash": source.source_hash}


def _equal_rationals(first, first_denominator, second, second_denominator):
    bound = max(int(np.max(np.abs(first), initial=0))*second_denominator,
                int(np.max(np.abs(second), initial=0))*first_denominator)
    if bound >= 2**62:
        first, second = first.astype(object), second.astype(object)
    return np.array_equal(first*second_denominator, second*first_denominator)


def _exact_dot(first, second, axes):
    size = first.shape[axes[0]]
    bound = int(np.max(np.abs(first), initial=0))*int(np.max(np.abs(second), initial=0))*size
    if bound >= 2**62:
        raise OverflowError("exact tensor contraction exceeds the declared int64 bound")
    return np.tensordot(first, second, axes=axes)


@cache
def coordinate_restriction(top, lower, basis):
    return restriction_in_basis(top, lower, tuple(basis))


def verify_composition(parent, relative, child):
    """Check every scalar coefficient and the actual transitive trace metric."""
    a, b = parent.actions, relative.actions
    product = _exact_dot(a.numerators, b.numerators[:parent.lower_degree], (3, 0))
    # e,parent-row,parent-column,new-row,new-column,lower-coefficient.
    composed = product.transpose(0, 1, 3, 2, 4, 5).reshape(child.actions.numerators.shape)
    if not _equal_rationals(composed, a.denominator*b.denominator,
                            child.actions.numerators, child.actions.denominator):
        raise ValueError("literal scalar restriction differs from the child action table")
    traces = [relative.trace_exponent(k) for k in range(relative.upper_phase)]
    if len({t.denominator for t in traces}) != 1:
        raise ValueError("trace table denominators are inconsistent")
    trace_array = np.stack([t.numerators for t in traces])
    d, r, h = parent.width, relative.width, relative.lower_degree
    total = np.zeros((d, r, d, r, h), dtype=np.int64)
    for i in range(r):
        for j in range(r):
            shifts = (np.arange(parent.lower_degree)+j-i) % relative.upper_phase
            total[:, i, :, j, :] = _exact_dot(parent.metric.numerators, trace_array[shifts], (2, 0))
    if not _equal_rationals(total.reshape(child.metric.numerators.shape),
                            parent.metric.denominator*traces[0].denominator,
                            child.metric.numerators, child.metric.denominator):
        raise ValueError("the child metric is not the transitive relative trace metric")
    return {"every_action_coefficient_composed_exactly": True,
            "all_original_root_powers_checked": child.upper_phase,
            "metric_trace_transitivity_checked_coefficientwise": True,
            "width_multiplied_by_relative_degree": child.width == parent.width*relative.width,
            "operation": "entrywise scalar restriction, with no Fourier frame or source modification"}


@cache
def canonical_basis_certificate(restriction):
    """Record equivalence without treating it as literal-matrix equality."""
    canonical = edge(restriction.upper_phase, restriction.lower_phase)
    change = canonical.inverse_basis_matrix * restriction.basis_matrix
    d, h = restriction.width, restriction.lower_degree
    coefficients = [[[change[i*h+b, j*h] for b in range(h)] for j in range(d)] for i in range(d)]
    c = RationalArray.exact(coefficients)
    if regular_matrix(c, restriction.lower_phase) != change:
        raise ValueError("the basis change is not linear over the lower field")
    if canonical.basis_matrix*change != restriction.basis_matrix:
        raise ValueError("the literal basis-change identity failed")
    generator = regular_matrix(restriction.action(1), restriction.lower_phase)
    canonical_generator = regular_matrix(canonical.action(1), restriction.lower_phase)
    if canonical_generator*change != change*generator:
        raise ValueError("the basis change fails scalar conjugacy")
    c_star = regular_matrix(conjugate_transpose(c, restriction.lower_phase), restriction.lower_phase)
    if c_star*regular_matrix(canonical.metric, restriction.lower_phase)*change != regular_matrix(restriction.metric, restriction.lower_phase):
        raise ValueError("the basis change fails the Hermitian metric congruence")
    return {"canonical_basis_exponents": list(canonical.basis_exponents),
            "path_basis_exponents": list(restriction.basis_exponents),
            "change_over_lower_field": c.to_json(),
            "basis_identity_exact": True, "scalar_conjugacy_exact": True,
            "Hermitian_metric_congruence_exact": True,
            "general_basis_equivalence_used_to_prune": False}


def source_identity(source):
    metadata = source_metadata(source)
    if metadata.get("kind") == "dense_exponents":
        return digest({key: metadata[key] for key in ("kind", "order", "phase", "coefficient_array_sha256")})
    return digest({"formula": "H[a,b]=coefficients[(b-a) mod prime] on F_prime^2",
                   "metadata": metadata})


def matrix_identity(source, restriction):
    used = list(map(int, np.unique(source_array(source))))
    return digest({"source": source_identity(source), "lower_field": restriction.lower_phase,
                   "order": source.order*restriction.width,
                   "used_phase_actions_sha256": _array_hash(restriction.actions.numerators[used]),
                   "action_denominator": restriction.actions.denominator,
                   "metric_sha256": _array_hash(restriction.metric.numerators),
                   "metric_denominator": restriction.metric.denominator})


def rational_endpoint_statistics(source, restriction):
    if restriction.lower_phase != 1:
        return None
    phase_counts = np.bincount(source_array(source).ravel(), minlength=source.phase)
    development_repetitions = source.order if hasattr(source, "coefficients") else 1
    counts = Counter()
    for exponent, frequency in enumerate(phase_counts):
        if not frequency:
            continue
        values, multiplicities = np.unique(restriction.actions.numerators[exponent], return_counts=True)
        for value, multiplicity in zip(values, multiplicities, strict=True):
            counts[str(Fraction(int(value), restriction.actions.denominator))] += int(frequency)*development_repetitions*int(multiplicity)
    order = source.order*restriction.width
    if sum(counts.values()) != order**2:
        raise ValueError("the exact endpoint entry census has the wrong total")
    return {"matrix_order": order, "all_entries_integral": all(Fraction(v).denominator == 1 for v in counts),
            "literal_entries_all_signs": set(counts) <= {"1", "-1"},
            "complete_matrix_entry_counts": dict(sorted(counts.items())),
            "ordinary_Euclidean_Hadamard_claimed": False,
            "weighted_Gram_scale": source.order,
            "weighted_Gram_scale_is_matrix_order": source.order == order}


def completed_inventory_sign_gate(source_ids):
    """Replay the existing exact alphabet gate; never prune this inventory."""
    path = ROOT / "certificates/scalar-lift/scalar-restriction-all-sources.json"
    if not path.exists():
        return {"available": False, "weighted_inventory_pruned": False}
    certificate = json.loads(path.read_text())
    audits = {row["source_id"]: row for row in certificate["source_audits"]}
    records = []
    for identifier in source_ids:
        source, _ = load_tree_source(identifier)
        audit = audits.get(identifier)
        if audit is None:
            from experiments.scalar_lift.scalar_restriction_obstructions import alphabet_obstruction

            audit = {"phase": source.phase, "order": source.order,
                     "alphabet_gate": alphabet_obstruction(source.phase, np.unique(source_array(source)))}
        gate = audit["alphabet_gate"]
        used = set(map(int, np.unique(source_array(source))))
        if audit["phase"] != source.phase or audit["order"] != source.order or used != set(gate["used_alphabet"]):
            raise ValueError("the saved alphabet gate differs from the actual matrix alphabet")
        if not gate["excluded"]:
            records.append({"source": identifier, "phase": source.phase, "source_order": source.order,
                            "coefficient_array_sha256": source_metadata(source)["coefficient_array_sha256"],
                            "not_excluded_by_saved_alphabet_gate": True, "sign_flattenability_proved": False})
            continue
        terms = gate["relation"]["terms"]
        roots = np.asarray(integral_roots(source.phase), dtype=np.int64)
        if (audit["phase"] != source.phase or audit["order"] != source.order
                or any(term["exponent"] not in used for term in terms)
                or sum(term["coefficient"] for term in terms) % 2 != 1
                or np.any(sum((term["coefficient"]*roots[term["exponent"]] for term in terms), np.zeros(roots.shape[1], dtype=np.int64)))):
            raise ValueError("the saved uniform sign obstruction does not bind this actual source")
        records.append({"source": identifier, "phase": source.phase, "source_order": source.order,
                        "coefficient_array_sha256": source_metadata(source)["coefficient_array_sha256"],
                        "relation": terms, "exact_integer_relation_replayed": True,
                        "all_uniform_Q_linear_sign_endpoints_excluded": True,
                        "scope": gate["scope"],
                        "source_modifications_or_nonuniform_maps_excluded": False})
    return {"available": True, "certificate": str(path.relative_to(ROOT)),
            "certificate_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "sources": records, "weighted_inventory_pruned": False}


def save_node(path, source, source_proof, restriction, field_path, composition):
    basis_proof = canonical_basis_certificate(restriction)
    key = matrix_identity(source, restriction)
    descriptor = {"matrix_hash": key, "source": source_metadata(source), "source_identity": source_identity(source),
                  "source_Gram_proof": source_proof, "original_phase": source.phase,
                  "lower_field": restriction.lower_phase, "field_path": list(field_path),
                  "basis_exponents": list(restriction.basis_exponents),
                  "width": restriction.width, "matrix_order": source.order*restriction.width,
                  "lower_degree": restriction.lower_degree,
                  "action_denominator": restriction.actions.denominator,
                  "metric_denominator": restriction.metric.denominator,
                  "inverse_metric_denominator": restriction.inverse_metric.denominator,
                  "action_sha256": _array_hash(restriction.actions.numerators),
                  "metric_sha256": _array_hash(restriction.metric.numerators),
                  "inverse_metric_sha256": _array_hash(restriction.inverse_metric.numerators),
                  "restriction_proof": restriction.proof(), "composition_proof": composition,
                  "canonical_basis_equivalence": basis_proof,
                  "global_weight_matrix": "I_source_order tensor inverse_metric",
                  "exact_weighted_Gram": "A W A*=source_order W",
                  "no_frame_inserted": True, "rational_endpoint": rational_endpoint_statistics(source, restriction),
                  "new_real_Hadamard_claimed": False,
                  "field_scope": "full cyclotomic fields with canonical conductors dividing the original phase; not all abelian subfields"}
    arrays = {"source_coefficients": source_array(source).astype(np.uint8),
              "actions": restriction.actions.numerators,
              "metric": restriction.metric.numerators,
              "inverse_metric": restriction.inverse_metric.numerators,
              "metadata_json": np.frombuffer(encoded(descriptor).encode(), dtype=np.uint8)}
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("wb") as handle:
        np.savez_compressed(handle, **arrays)
    temporary.replace(path)
    return descriptor


class SavedWeightedMatrix:
    def __init__(self, path):
        with np.load(path, allow_pickle=False) as arrays:
            self.metadata = json.loads(arrays["metadata_json"].tobytes())
            self.coefficients = arrays["source_coefficients"].copy()
            self.actions = arrays["actions"].copy()
            self.metric = arrays["metric"].copy()
            self.inverse_metric = arrays["inverse_metric"].copy()
        for array, label in ((self.actions, "action"), (self.metric, "metric"), (self.inverse_metric, "inverse_metric")):
            if _array_hash(array) != self.metadata[label + "_sha256"]:
                raise ValueError("a saved weighted action or metric was modified")
        if _array_hash(self.coefficients) != self.metadata["source"]["coefficient_array_sha256"]:
            raise ValueError("the saved source coefficient array was modified")
        self.order = self.metadata["matrix_order"]
        self.width = self.metadata["width"]
        self.lower_field = self.metadata["lower_field"]
        top = self.metadata["original_phase"]
        h = len(integral_roots(self.lower_field)[0])
        if (top != self.metadata["source"]["phase"]
                or self.width*h != len(integral_roots(top)[0])
                or self.order != self.metadata["source"]["order"]*self.width
                or self.actions.shape != (top, self.width, self.width, h)
                or self.metric.shape != (self.width, self.width, h)
                or self.inverse_metric.shape != self.metric.shape
                or self.metadata["lower_degree"] != h):
            raise ValueError("saved matrix dimensions do not match the exact field degree")

    def replay(self):
        """Reconstruct the source and every scalar edge, not just saved hashes."""
        metadata = self.metadata["source"]
        if metadata.get("kind") == "dense_exponents":
            source, _ = load_tree_source(metadata["registry_source_id"])
        else:
            from experiments.scalar_lift.forward_character_sources import (
                character_source,
                formal_gram_check,
            )

            source = character_source(metadata["prime"], metadata["character_orders"][1],
                                      first_order=metadata["character_orders"][0], gamma=metadata["gamma_exponent"])
            formal_gram_check(source)
        if not np.array_equal(source_array(source), self.coefficients):
            raise ValueError("the saved matrix does not reconstruct its verified source")
        field_path = self.metadata["field_path"]
        if field_path[0] != source.phase or field_path[-1] != self.lower_field:
            raise ValueError("the saved scalar field path has inconsistent endpoints")
        basis = (0,)
        current = coordinate_restriction(source.phase, source.phase, basis)
        for lower in field_path[1:]:
            relative = edge(current.lower_phase, lower)
            basis = tuple(a+(source.phase//current.lower_phase)*b for a in basis for b in relative.basis_exponents)
            child = coordinate_restriction(source.phase, lower, basis)
            verify_composition(current, relative, child)
            current = child
        if basis != tuple(self.metadata["basis_exponents"]):
            raise ValueError("the saved coordinate basis differs from its scalar path")
        for saved, exact, denominator in ((self.actions, current.actions, self.metadata["action_denominator"]),
                                          (self.metric, current.metric, self.metadata["metric_denominator"]),
                                          (self.inverse_metric, current.inverse_metric, self.metadata["inverse_metric_denominator"])):
            if not _equal_rationals(saved, denominator, exact.numerators, exact.denominator):
                raise ValueError("the saved weighted matrix or metric failed exact scalar replay")
        canonical_basis_certificate(current)
        if matrix_identity(source, current) != self.metadata["matrix_hash"]:
            raise ValueError("the saved matrix identity failed exact source and scalar replay")
        return {"source_Gram_replayed": True, "every_scalar_edge_replayed": True,
                "actual_actions_and_metrics_replayed": True,
                "matrix_hash_matches": True}

    def entry(self, row, column):
        if type(row) is not int or type(column) is not int or not 0 <= row < self.order or not 0 <= column < self.order:
            raise ValueError("matrix coordinate outside its order")
        a, r = divmod(row, self.width)
        b, s = divmod(column, self.width)
        if self.metadata["source"].get("kind") == "dense_exponents":
            exponent = int(self.coefficients[a, b])
        else:
            p = self.metadata["source"]["prime"]
            a0, a1 = divmod(a, p)
            b0, b1 = divmod(b, p)
            exponent = int(self.coefficients[(b0-a0) % p, (b1-a1) % p])
        return {"lower_field": self.lower_field,
                "numerators": self.actions[exponent, r, s].tolist(),
                "denominator": self.metadata["action_denominator"]}

    def weight_entry(self, row, column):
        if not 0 <= row < self.order or not 0 <= column < self.order:
            raise ValueError("metric coordinate outside its order")
        a, r = divmod(row, self.width)
        b, s = divmod(column, self.width)
        values = self.inverse_metric[r, s] if a == b else np.zeros(self.metadata["lower_degree"], dtype=np.int64)
        return {"lower_field": self.lower_field, "numerators": values.tolist(),
                "denominator": self.metadata["inverse_metric_denominator"]}


class WeightedTree:
    def __init__(self, database, max_db_mib=256, *, seeds=SEEDS):
        if not math.isfinite(max_db_mib) or max_db_mib < 4:
            raise ValueError("at least four MiB of storage are required")
        allowed = set(SEEDS) | set(retained_registry_seeds()) | set(additional_input_seeds())
        if not seeds or any(seed not in allowed for seed in seeds) or len(set(seeds)) != len(seeds):
            raise ValueError("only the fixed verified source inputs are allowed")
        self.path = Path(database).resolve()
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.artifacts = Path(str(self.path) + ".matrices")
        self.artifacts.mkdir(exist_ok=True)
        self.lock = Path(str(self.path) + ".lock").open("a+")  # noqa: SIM115 -- checkpoint lifetime
        try:
            fcntl.flock(self.lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            self.lock.close()
            raise RuntimeError("another worker owns this weighted tree") from None
        self.db = sqlite3.connect(self.path)
        self.db.execute("PRAGMA journal_mode=DELETE")
        self.db.execute("PRAGMA temp_store=FILE")
        self.db.executescript("""
            CREATE TABLE IF NOT EXISTS meta(key TEXT PRIMARY KEY,value TEXT);
            CREATE TABLE IF NOT EXISTS tasks(id TEXT PRIMARY KEY,recipe TEXT,depth INTEGER,status TEXT);
            CREATE TABLE IF NOT EXISTS nodes(id TEXT PRIMARY KEY,record TEXT);
            CREATE TABLE IF NOT EXISTS edges(id TEXT PRIMARY KEY,parent TEXT,child TEXT,record TEXT);
        """)
        page = self.db.execute("PRAGMA page_size").fetchone()[0]
        self.db.execute(f"PRAGMA max_page_count={int(max_db_mib*1024**2)//(4*page)}")
        self.max_bytes = int(max_db_mib*1024**2)
        self.sources = {}
        identity = {"schema": SCHEMA, "seeds": list(seeds), "deduplication": "exact used actions and trace metric only",
                    "source_or_frame_modifications": False}
        previous = self._get("identity")
        if previous is not None and previous != identity:
            self.close()
            raise ValueError("another weighted source catalogue owns this checkpoint")
        with self.db:
            if previous is None:
                self._put("identity", identity)
                self._put("progress", {"evaluated_edges": 0, "literal_nodes": 0, "duplicates_pruned": 0})
                for identifier in seeds:
                    if identifier in SEEDS:
                        top = 24 if identifier == "character-73" else 120
                    else:
                        top = load_tree_source(identifier)[0].phase
                    self._enqueue({"source": identifier, "path": [top], "basis": [0], "parent": None})
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
                        (digest(recipe), encoded(recipe), len(recipe["path"])-1))

    def storage_bytes(self):
        return sum(p.stat().st_size for p in self.artifacts.rglob("*") if p.is_file()) + sum(
            p.stat().st_size for p in self.path.parent.glob(self.path.name + "*") if p.is_file() and not p.name.endswith(".lock"))

    def _evaluate(self, task, recipe):
        if recipe["source"] not in self.sources:
            self.sources[recipe["source"]] = load_tree_source(recipe["source"])
        source, source_proof = self.sources[recipe["source"]]
        top, lower = source.phase, recipe["path"][-1]
        restriction = coordinate_restriction(top, lower, tuple(recipe["basis"]))
        composition = {"existing_verified_source_node": True}
        if recipe["parent"] is not None:
            parent_record = json.loads(self.db.execute("SELECT record FROM nodes WHERE id=?", (recipe["parent"],)).fetchone()[0])
            parent = coordinate_restriction(top, parent_record["lower_field"], tuple(parent_record["basis_exponents"]))
            composition = verify_composition(parent, edge(parent.lower_phase, lower), restriction)
        key = matrix_identity(source, restriction)
        duplicate = self.db.execute("SELECT 1 FROM nodes WHERE id=?", (key,)).fetchone() is not None
        basis_proof = canonical_basis_certificate(restriction)
        if not duplicate:
            bundle = self.artifacts / (key + ".npz")
            descriptor = save_node(bundle, source, source_proof, restriction, recipe["path"], composition)
            record = {"id": key, "source": recipe["source"], "lower_field": lower,
                      "basis_exponents": recipe["basis"], "field_path": recipe["path"],
                      "matrix_order": descriptor["matrix_order"], "width": restriction.width,
                      "bundle": str(bundle), "bundle_sha256": hashlib.sha256(bundle.read_bytes()).hexdigest(),
                      "rational_endpoint": descriptor["rational_endpoint"],
                      "composition_proof": composition, "canonical_basis_equivalence": basis_proof,
                      "global_weighted_Gram_scale": source.order, "new_real_Hadamard_claimed": False}
            if self.storage_bytes() >= self.max_bytes:
                bundle.unlink(missing_ok=True)
                raise OverflowError("completed weighted bundle would exceed storage budget")
        edge_record = {"path": recipe["path"], "literal_duplicate_pruned": duplicate,
                       "basis_equivalence": basis_proof, "composition_proof": composition,
                       "basis_equivalence_alone_used_for_pruning": False}
        with self.db:
            if not duplicate:
                self.db.execute("INSERT INTO nodes VALUES (?,?)", (key, encoded(record)))
                for child_field in proper_fields(lower):
                    relative = edge(lower, child_field)
                    basis = [a + (top//lower)*b for a in recipe["basis"] for b in relative.basis_exponents]
                    self._enqueue({"source": recipe["source"], "path": [*recipe["path"], child_field],
                                   "basis": basis, "parent": key})
            self.db.execute("INSERT INTO edges VALUES (?,?,?,?)", (task, recipe["parent"], key, encoded(edge_record)))
            self.db.execute("UPDATE tasks SET status='complete' WHERE id=?", (task,))
            progress = self._get("progress")
            progress["evaluated_edges"] += 1
            progress["literal_nodes"] += int(not duplicate)
            progress["duplicates_pruned"] += int(duplicate)
            progress["last_matrix"] = {"order": source.order*restriction.width, "lower_field": lower,
                                       "depth": len(recipe["path"])-1}
            self._put("progress", progress)

    def run(self, seconds=60, max_nodes=None):
        if not math.isfinite(seconds) or seconds < 0 or (max_nodes is not None and max_nodes < 0):
            raise ValueError("finite nonnegative budgets are required")
        started, completed = time.monotonic(), 0
        self.stop_reason = "time_budget"
        while time.monotonic()-started < seconds:
            if max_nodes is not None and completed >= max_nodes:
                self.stop_reason = "node_budget"
                break
            task = self.db.execute("SELECT id,recipe FROM tasks WHERE status='pending' ORDER BY depth,rowid LIMIT 1").fetchone()
            if task is None:
                self.stop_reason = "registered_divisor_field_tree_complete"
                break
            if self.storage_bytes() >= self.max_bytes:
                self.stop_reason = "storage_budget"
                break
            try:
                self._evaluate(task[0], json.loads(task[1]))
            except OverflowError:
                self.stop_reason = "storage_budget"
                break
            except sqlite3.OperationalError as error:
                if getattr(error, "sqlite_errorcode", None) != sqlite3.SQLITE_FULL:
                    raise
                self.stop_reason = "storage_budget"
                break
            completed += 1
        return {**self.report(), "batch_nodes": completed, "batch_seconds": time.monotonic()-started}

    def report(self):
        nodes = [json.loads(raw) for (raw,) in self.db.execute("SELECT record FROM nodes ORDER BY rowid")]
        depth_counts = Counter(len(node["field_path"])-1 for node in nodes)
        field_counts = Counter((node["source"], node["lower_field"], node["matrix_order"]) for node in nodes)
        return {"schema": SCHEMA, "stop_reason": self.stop_reason, "progress": self._get("progress"),
                "nodes": nodes,
                "edges": [json.loads(raw) for (raw,) in self.db.execute("SELECT record FROM edges ORDER BY rowid")],
                "pending_tasks": self.db.execute("SELECT COUNT(*) FROM tasks WHERE status='pending'").fetchone()[0],
                "rational_weighted_orders": sorted({n["matrix_order"] for n in nodes if n["lower_field"] == 1}),
                "node_counts_by_depth": dict(sorted(depth_counts.items())),
                "source_field_order_counts": [{"source": s, "field": q, "order": n, "nodes": count}
                                              for (s, q, n), count in sorted(field_counts.items())],
                "uniform_sign_gate_after_inventory": completed_inventory_sign_gate(self._get("identity")["seeds"]),
                "only_operations": "scalar restriction and exact literal duplicate pruning",
                "basis_distinct_candidates_retained": True,
                "field_scope": "full cyclotomic divisor fields; other abelian subfields are not implemented",
                "flatness_or_Euclidean_normalization_assumed": False,
                "all_scalar_restriction_methods_exhausted": False,
                "new_real_Hadamard_claimed": False, "storage_bytes": self.storage_bytes()}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--seconds", type=float, default=60)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--max-db-mib", type=int, default=256)
    parser.add_argument("--nodes", type=int)
    inputs = parser.add_mutually_exclusive_group()
    inputs.add_argument("--retained-registry", action="store_true",
                        help="use the 18 existing sources retained by the saved uniform alphabet gate")
    inputs.add_argument("--additional-inputs", action="store_true",
                        help="use the three acquired catalogue inputs and saved order-nine input")
    inputs.add_argument("--source", action="append", dest="sources",
                        help="select a fixed verified source; repeat for multiple sources")
    args = parser.parse_args(argv)
    seeds = (retained_registry_seeds() if args.retained_registry else additional_input_seeds()
             if args.additional_inputs else tuple(args.sources) if args.sources else SEEDS)
    tree = WeightedTree(args.database, args.max_db_mib, seeds=seeds)
    try:
        report = tree.run(args.seconds, args.nodes)
    finally:
        tree.close()
    write_json(args.report, report)
    print(json.dumps({key: report[key] for key in ("stop_reason", "progress", "pending_tasks", "rational_weighted_orders", "storage_bytes")}))


if __name__ == "__main__":
    main()

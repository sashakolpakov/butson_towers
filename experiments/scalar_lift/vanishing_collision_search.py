"""Disk-backed meet-in-the-middle search for minimal positive root sums.

Modular values are a candidate index only. Every reported relation is checked
over Z[x]/Phi_q and for absence of any proper vanishing submultiset. This is a
finite, resumable substrate search, not a Butson matrix construction.
"""

from __future__ import annotations

import argparse
import fcntl
import json
import math
import sqlite3
import struct
import time
from collections import Counter
from functools import cache
from pathlib import Path

from sympy import Poly, Symbol, cyclotomic_poly, factorint, isprime, primitive_root

SCHEMA = "vanishing-collision-v1"


@cache
def _primes(conductor):
    return tuple(map(int, factorint(conductor)))


@cache
def integral_roots(conductor):
    """All root powers in an exact integral power basis, without floats."""
    x = Symbol("x")
    polynomial = tuple(reversed([int(c) for c in Poly(cyclotomic_poly(conductor, x), x).all_coeffs()]))
    degree = len(polynomial) - 1
    values, current = [], (1,) + (0,) * (degree - 1)
    for _ in range(conductor):
        values.append(current)
        high = current[-1]
        shifted = (0,) + current[:-1]
        current = tuple(a - high * b for a, b in zip(shifted, polynomial[:-1], strict=True))
    if current != values[0]:
        raise ArithmeticError("root powers failed the conductor identity")
    return tuple(values)


def exact_sum(exponents, conductor):
    roots = integral_roots(conductor)
    answer = [0] * len(roots[0])
    for exponent in exponents:
        for i, value in enumerate(roots[exponent % conductor]):
            answer[i] += value
    return tuple(answer)


def canonical_relation(exponents, conductor):
    """Remove common phase and Galois action, retaining multiplicities."""
    terms = tuple(int(e) % conductor for e in exponents)
    return min(tuple(sorted((unit * (e - anchor)) % conductor for e in terms))
               for unit in range(1, conductor) if math.gcd(unit, conductor) == 1
               for anchor in set(terms))


def prime_subsum(exponents, conductor, *, proper=False):
    """Cheap exact rejection; absence of a polygon does not prove minimality."""
    counts = Counter(exponents)
    for prime in _primes(conductor):
        if prime > len(exponents) or (proper and prime == len(exponents)):
            continue
        step = conductor // prime
        for offset in {e % step for e in counts}:
            polygon = tuple(offset + j * step for j in range(prime))
            if all(counts[e] for e in polygon):
                return polygon
    return None


def _subset_sums(terms, roots):
    """Gray-code iteration holds a single degree-dimensional running sum."""
    value = [0] * len(roots[0])
    previous = 0
    yield 0, tuple(value)
    for number in range(1, 1 << len(terms)):
        mask = number ^ (number >> 1)
        changed = mask ^ previous
        index = changed.bit_length() - 1
        sign = 1 if mask & changed else -1
        for i, coefficient in enumerate(roots[terms[index]]):
            value[i] += sign * coefficient
        yield mask, tuple(value)
        previous = mask


def minimality_certificate(exponents, conductor):
    """Exact subset MITM: O(2^(w/2)) records, bounded here to w <= 24."""
    terms = tuple(int(e) % conductor for e in exponents)
    if not 2 <= len(terms) <= 24:
        raise ValueError("exact minimality check supports weights 2 through 24")
    if any(exact_sum(terms, conductor)):
        return {"exact_zero": False, "minimal": False}
    polygon = prime_subsum(terms, conductor, proper=True)
    if polygon:
        return {"exact_zero": True, "minimal": False, "proper_subsum": list(polygon)}
    split = len(terms) // 2
    roots, left = integral_roots(conductor), {}
    for mask, value in _subset_sums(terms[:split], roots):
        if mask and not any(value):
            return {"exact_zero": True, "minimal": False,
                    "proper_subsum": [terms[i] for i in range(split) if mask >> i & 1]}
        left.setdefault(value, mask)
    full = (1 << len(terms)) - 1
    for mask, value in _subset_sums(terms[split:], roots):
        other = left.get(tuple(-v for v in value))
        if other is None:
            continue
        combined = other | (mask << split)
        if combined not in (0, full):
            return {"exact_zero": True, "minimal": False,
                    "proper_subsum": [e for i, e in enumerate(terms) if combined >> i & 1]}
    return {"exact_zero": True, "minimal": True,
            "method": "exact integral subset meet-in-the-middle", "subset_half_bound": 1 << max(split, len(terms) - split)}


class Fingerprint:
    """Two ring homomorphisms into split prime fields; no tolerance window."""

    def __init__(self, conductor, primes=None):
        if primes is None:
            found = []
            candidate = ((1_000_000 // conductor) + 1) * conductor + 1
            while len(found) < 2:
                if isprime(candidate):
                    found.append(candidate)
                candidate += conductor
            primes = found
        self.primes = tuple(int(p) for p in primes)
        self.powers = []
        for prime in self.primes:
            if not isprime(prime) or (prime - 1) % conductor:
                raise ValueError("fingerprint primes must split the conductor")
            root = pow(int(primitive_root(prime)), (prime - 1) // conductor, prime)
            self.powers.append(tuple(pow(root, e, prime) for e in range(conductor)))

    def value(self, terms):
        return tuple(sum(powers[e] for e in terms) % p
                     for powers, p in zip(self.powers, self.primes, strict=True))

    def key(self, terms, negate=False):
        sign = -1 if negate else 1
        return struct.pack("<" + "I" * len(self.primes),
                           *(sign * v % p for v, p in zip(self.value(terms), self.primes, strict=True)))


def _next_multiset(terms, conductor, maximum):
    values = list(terms)
    for i in range(len(values) - 1, -1, -1):
        if values[i] < conductor - 1:
            values[i:] = [values[i] + 1] * (len(values) - i)
            return values
    return [0] * (len(values) + 1) if len(values) < maximum else None


def _encode_terms(terms):
    return struct.pack("<" + "H" * len(terms), *terms)


def _decode_terms(blob):
    return struct.unpack("<" + "H" * (len(blob) // 2), blob)


class CollisionSearch:
    """Streaming half-table and persistent collision cursor, one writer per DB.

    All multisets of sizes 1..half_weight are considered. Halves containing a
    prime polygon are safely skipped for minimal-sum discovery. Matching halves
    may have unequal lengths, so odd weights are included. Nothing is sampled or
    dropped at a crowded hash bucket; the cursor can pause inside any bucket.
    """

    def __init__(self, database, conductor=30, half_weight=4, max_db_mib=256):
        if not 2 <= conductor <= 2310 or not 1 <= half_weight <= 12:
            raise ValueError("supported conductor is 2..2310; half-weight is 1..12")
        if max_db_mib < 1:
            raise ValueError("the database cap must be at least 1 MiB")
        self.path = Path(database)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.lock = self.path.with_suffix(self.path.suffix + ".lock").open("a+")
        try:
            fcntl.flock(self.lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            self.lock.close()
            raise RuntimeError("another worker owns this checkpoint") from None
        self.db = sqlite3.connect(self.path)
        self.db.execute("PRAGMA cache_size=-2048")
        self.db.execute("PRAGMA temp_store=FILE")
        self.db.execute("PRAGMA journal_mode=DELETE")
        self.db.executescript("""
            CREATE TABLE IF NOT EXISTS meta (key TEXT PRIMARY KEY, value TEXT NOT NULL);
            CREATE TABLE IF NOT EXISTS halves (id INTEGER PRIMARY KEY, fingerprint BLOB NOT NULL, terms BLOB NOT NULL);
            CREATE INDEX IF NOT EXISTS half_lookup ON halves(fingerprint, id);
            CREATE TABLE IF NOT EXISTS relations (canonical BLOB PRIMARY KEY, record TEXT NOT NULL);
            CREATE TABLE IF NOT EXISTS seen (canonical BLOB PRIMARY KEY, verdict TEXT NOT NULL);
        """)
        self.conductor, self.half_weight = conductor, half_weight
        self.fingerprint = Fingerprint(conductor)
        identity = {"schema": SCHEMA, "conductor": conductor, "half_weight": half_weight,
                    "fingerprint_primes": self.fingerprint.primes}
        identity = json.dumps(identity, sort_keys=True)
        old = self.db.execute("SELECT value FROM meta WHERE key='identity'").fetchone()
        if old and old[0] != identity:
            self.close()
            raise ValueError("checkpoint belongs to a different search recipe")
        if not old:
            with self.db:
                self.db.execute("INSERT INTO meta VALUES ('identity',?)", (identity,))
                self._save({"next_terms": [0], "pending": None, "finished": False,
                            "halves_visited": 0, "halves_stored": 0, "polygon_halves_skipped": 0,
                            "fingerprint_matches": 0, "exact_zero_matches": 0,
                            "modular_false_positives": 0, "decomposable_matches": 0,
                            "duplicate_relations": 0, "minimal_relations": 0})
        self.state = self._load()
        page_size = self.db.execute("PRAGMA page_size").fetchone()[0]
        self.db.execute(f"PRAGMA max_page_count={int(max_db_mib * 1024 * 1024) // page_size}")
        self.stop_reason = "not_started"

    def close(self):
        self.db.close()
        fcntl.flock(self.lock, fcntl.LOCK_UN)
        self.lock.close()

    def _save(self, state):
        self.db.execute("INSERT OR REPLACE INTO meta VALUES ('state',?)", (json.dumps(state),))

    def _load(self):
        return json.loads(self.db.execute("SELECT value FROM meta WHERE key='state'").fetchone()[0])

    def _candidate(self, terms):
        self.state["fingerprint_matches"] += 1
        if any(exact_sum(terms, self.conductor)):
            self.state["modular_false_positives"] += 1
            return
        self.state["exact_zero_matches"] += 1
        if prime_subsum(terms, self.conductor, proper=True):
            self.state["decomposable_matches"] += 1
            return
        canonical = canonical_relation(terms, self.conductor)
        key = _encode_terms(canonical)
        if self.db.execute("SELECT 1 FROM seen WHERE canonical=?", (key,)).fetchone():
            self.state["duplicate_relations"] += 1
            return
        certificate = minimality_certificate(canonical, self.conductor)
        self.db.execute("INSERT INTO seen VALUES (?,?)", (key, json.dumps(certificate)))
        if not certificate["minimal"]:
            self.state["decomposable_matches"] += 1
            return
        record = {"conductor": self.conductor, "exponents": canonical, "weight": len(canonical),
                  "equivalence": "common rotation and Galois automorphism",
                  "relative_conductor": self.conductor // math.gcd(self.conductor, *canonical),
                  **certificate, "complete_matrix": False}
        self.db.execute("INSERT INTO relations VALUES (?,?)", (key, json.dumps(record)))
        self.state["minimal_relations"] += 1

    def step(self):
        pending = self.state["pending"]
        if pending:
            terms = pending["terms"]
            match = self.db.execute(
                "SELECT id,terms FROM halves WHERE fingerprint=? AND id>? AND id<=? ORDER BY id LIMIT 1",
                (self.fingerprint.key(terms, negate=True), pending["after"], pending["id"])).fetchone()
            if match:
                self._candidate(tuple(terms) + _decode_terms(match[1]))
                pending["after"] = match[0]
                return
            self.state["pending"] = None
        terms = self.state["next_terms"]
        if terms is None:
            self.state["finished"] = True
            return
        self.state["next_terms"] = _next_multiset(terms, self.conductor, self.half_weight)
        self.state["halves_visited"] += 1
        if prime_subsum(terms, self.conductor):
            self.state["polygon_halves_skipped"] += 1
            return
        row = self.db.execute("INSERT INTO halves(fingerprint,terms) VALUES (?,?)",
                              (self.fingerprint.key(terms), _encode_terms(terms)))
        self.state["halves_stored"] += 1
        self.state["pending"] = {"terms": terms, "id": row.lastrowid, "after": 0}

    def run(self, seconds=30, steps=None):
        if seconds < 0 or (steps is not None and steps < 0):
            raise ValueError("work budgets must be nonnegative")
        started, count = time.monotonic(), 0
        self.stop_reason = "time_budget"
        try:
            while not self.state["finished"] and time.monotonic() - started < seconds:
                if steps is not None and count >= steps:
                    self.stop_reason = "step_budget"
                    break
                self.step()
                count += 1
                if count % 128 == 0:
                    self._save(self.state)
                    self.db.commit()
            self._save(self.state)
            self.db.commit()
        except sqlite3.OperationalError as error:
            self.db.rollback()
            self.state = self._load()
            if "full" not in str(error).lower():
                raise
            self.stop_reason = "database_cap"
        except BaseException:
            self.db.rollback()
            self.state = self._load()
            raise
        if self.state["finished"]:
            self.stop_reason = "finite_enumeration_finished"
        return {**self.report(), "batch_steps": count, "batch_seconds": round(time.monotonic() - started, 6)}

    def report(self):
        return {"schema": SCHEMA, "conductor": self.conductor, "half_weight": self.half_weight,
                "weight_bound": 2 * self.half_weight, "stop_reason": self.stop_reason,
                "state": self.state, "raw_half_count": sum(math.comb(self.conductor + k - 1, k)
                    for k in range(1, self.half_weight + 1)),
                "sqlite_page_bytes": self.db.execute("PRAGMA page_count").fetchone()[0]
                    * self.db.execute("PRAGMA page_size").fetchone()[0],
                "sqlite_cache_kib": 2048, "complete_butson_found": False,
                "search_scope": "all root multisets of weights 2..2k at the specified conductor, up to rotation and Galois",
                "examples": [json.loads(row[0]) for row in self.db.execute("SELECT record FROM relations ORDER BY canonical LIMIT 20")]}


def read_status(database):
    path = Path(database).resolve()
    with sqlite3.connect(path.as_uri() + "?mode=ro", uri=True) as db:
        identity = json.loads(db.execute("SELECT value FROM meta WHERE key='identity'").fetchone()[0])
        state = json.loads(db.execute("SELECT value FROM meta WHERE key='state'").fetchone()[0])
        return {**identity, "state": state, "sqlite_file_bytes": path.stat().st_size,
                "complete_butson_found": False}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--conductor", type=int, default=30)
    parser.add_argument("--half-weight", type=int, default=4)
    parser.add_argument("--seconds", type=float, default=30)
    parser.add_argument("--steps", type=int)
    parser.add_argument("--max-db-mib", type=int, default=256)
    parser.add_argument("--status", action="store_true")
    parser.add_argument("--report", type=Path)
    args = parser.parse_args(argv)
    if args.status:
        result = read_status(args.database)
    else:
        search = CollisionSearch(args.database, args.conductor, args.half_weight, args.max_db_mib)
        try:
            result = search.run(args.seconds, args.steps)
        finally:
            search.close()
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.report.with_suffix(args.report.suffix + ".tmp")
        temporary.write_text(json.dumps(result, indent=2) + "\n")
        temporary.replace(args.report)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()

"""Only complete, verified inputs for the forward scalar-restriction engine."""

from __future__ import annotations

import hashlib
import json
from collections.abc import Callable
from dataclasses import dataclass, field
from functools import cache
from pathlib import Path

import numpy as np

from experiments.scalar_lift.forward_butson_descent import (
    ForwardTower,
    exponent_hash,
    verify_source,
)

ROOT = Path(__file__).resolve().parents[2]


@dataclass(frozen=True)
class VerifiedSource:
    id: str
    order: int
    phase: int
    source_hash: str
    used_alphabet: tuple[int, ...]
    proof: dict
    exponents: np.ndarray | None = field(repr=False, compare=False)
    _entry: Callable[[int, int], int] = field(repr=False, compare=False)

    @property
    def hash(self):
        return self.source_hash

    @property
    def used_exponents(self):
        return self.used_alphabet

    def entry(self, row, column):
        if int(row) != row or int(column) != column or not 0 <= row < self.order or not 0 <= column < self.order:
            raise ValueError("source entry lies outside the verified order")
        return int(self._entry(int(row), int(column))) % self.phase

    def metadata(self):
        return {"id": self.id, "order": self.order, "phase": self.phase,
                "source_hash": self.source_hash, "used_alphabet": list(self.used_alphabet),
                "dense_exponents_available": self.exponents is not None, "proof": self.proof}


def sources():
    """Compact descriptors; listing does not allocate or verify large arrays."""
    result = [{"id": f"conference-{n}-q{q}", "family": "conference",
               "order": n, "phase": q} for n in (1000, 2000) for q in (8, 16)]
    result += [{"id": f"unreal{n}-q6", "family": "sextic-shifted",
                "order": n, "phase": 6, "global_exponent_shift": -1} for n in (3, 13)]
    result += [{"id": f"common120-{n}", "family": "common120", "order": n, "phase": 120}
               for n in (15, 35)]
    result += [{"id": f"cubic13-q24-t{k}", "family": "cubic13-tripling",
                "order": 13 * 3**k, "phase": 24, "triplings": k} for k in range(5)]
    result += [{"id": "fifth-5-q20", "family": "fifth-chirp", "order": 5, "phase": 20}]
    result += [{"id": f"fifth-40-q{q}", "family": "fifth-tensor", "order": 40, "phase": q}
               for q in (40, 80)]
    result += [{"id": f"fifth-{n}-q{q}", "family": "fifth-descendant", "order": n, "phase": q,
                "starting_source": "fifth-40-q80", "steps": k}
               for n, q, k in ((80, 40, 1), (160, 20, 2))]
    result += [{"id": f"character-{p}", "family": "character", "order": p*p,
                "phase": q, "prime": p} for p, q in ((73, 24), (241, 120))]
    return result


def _dense(descriptor, array, proof=None):
    phase, n = descriptor["phase"], descriptor["order"]
    matrix = np.asarray(array, dtype=np.int64) % phase
    if n > 2000 or matrix.shape != (n, n):
        raise ValueError("dense registry sources must match their descriptor and have order at most 2000")
    digest = exponent_hash(matrix)
    if proof is None:
        proof = verify_source(matrix, phase)
    else:
        proof = {**proof, "order": n, "phase": phase, "exponent_sha256": digest}
    if proof.get("full_Gram_verified") is not True:
        raise ValueError("an unfinished source cannot enter the registry")
    matrix = matrix.astype(np.uint16)
    matrix.setflags(write=False)
    return VerifiedSource(descriptor["id"], n, phase, digest, tuple(map(int, np.unique(matrix))),
                          proof, matrix, lambda i, j: int(matrix[i, j]))


@cache
def _conference_proof(phase):
    from experiments.scalar_lift.even_petrescu_conference_family import formal_family_certificate
    from experiments.scalar_lift.even_petrescu_lattice import lattice_data

    _, _, s, _, _ = lattice_data()  # Replays the actual full order-500 identities.
    c = s - np.eye(500, dtype=np.int64)
    if not np.array_equal(c.T, -c):
        raise ValueError("the actual conference skew identity failed")
    return {"full_Gram_verified": True, "construction": "actual skew conference plus exact four-block polynomial Gram",
            "formal_Gram": formal_family_certificate(phase),
            "conference_order": 500, "conference_weight": 499}


@cache
def _tripling_source(repetitions):
    from experiments.scalar_lift.inverse_blocks import known_order_thirteen_source
    from experiments.scalar_lift.rectangular_petrescu_sources import CYCLE, _components, _triple

    source = np.asarray(known_order_thirteen_source(), dtype=np.int64)
    seed_proof = verify_source(source, 6)
    a, b = _components(source)
    if (not np.array_equal(a @ a.T + 3 * b @ b.T, 52 * np.eye(13, dtype=np.int64))
            or not np.array_equal(a @ b.T, b @ a.T)):
        raise ValueError("the exact sixth-root coefficient seed equations failed")
    j, identity = np.ones((3, 3), dtype=np.int64), np.eye(3, dtype=np.int64)
    c = np.asarray(CYCLE, dtype=np.int64)
    if (not np.array_equal(c.T, -c) or not np.array_equal(c @ c.T, 3 * identity - j)
            or np.any(c @ j) or np.any(j @ c) or not np.array_equal(j @ j, 3 * j)):
        raise ValueError("the tripling tensor identities failed")
    stages = [{"order": 13, "sixth_exponent_sha256": exponent_hash(source)}]
    for _ in range(repetitions):
        aa, bb = _components(source)
        source = _triple(source)
        expected_a, expected_b = np.kron(bb, j), np.kron(aa, identity) + np.kron(bb, c)
        actual_a, actual_b = _components(source)
        if not np.array_equal(actual_a, expected_a) or not np.array_equal(actual_b, expected_b):
            raise ValueError("the actual expanded source disagrees with the proved tensor construction")
        stages.append({"order": len(source), "sixth_exponent_sha256": exponent_hash(source)})
    return (4 * source + 2) % 24, {
        "full_Gram_verified": True, "construction": "exact unreal-sixth-root tensor tripling, then phase-24 embedding and global root",
        "seed_proof": seed_proof, "stages": stages,
        "coefficient_recurrence": "A'=B tensor J3; B'=A tensor I3+B tensor C3",
        "inductive_Gram_identity": "AA^T+3BB^T=4nI and AB^T=BA^T imply the corresponding identities at 3n",
        "all_literal_tensor_expansions_and_root_coefficients_verified": True,
        "full_final_cubic_matrix_multiplication_needed": False}


def _character(descriptor):
    from experiments.scalar_lift.forward_character_sources import (
        character_source,
        formal_gram_check,
    )

    source = character_source(descriptor["prime"])
    metadata = source.metadata()
    proof = formal_gram_check(source)
    path = ROOT / "certificates/scalar-lift/forward-character-sources.json"
    old = next(record for record in json.loads(path.read_text())["known_source_constructions"]
               if record["prime"] == descriptor["prime"])
    if metadata["coefficient_array_sha256"] != old["coefficient_array_sha256"]:
        raise ValueError("the implicit source differs from its existing exact certificate")
    identity = {"formula": "H[i,j]=coefficients[(j0-i0) mod p,(j1-i1) mod p]",
                "linear_index_convention": "divmod(index,p)", "order": source.order,
                "phase": source.phase, "coefficient_array_sha256": metadata["coefficient_array_sha256"]}
    digest = hashlib.sha256(json.dumps(identity, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    proof = {"order": source.order, "phase": source.phase, "full_Gram_verified": True,
             "hash_kind": "exact implicit formula and complete coefficient table, not a hash of an allocated dense matrix",
             "exponent_sha256": digest, "implicit_identity": identity, "formal_Gram": proof,
             "prior_certificate": str(path.relative_to(ROOT)),
             "prior_certificate_sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
    return VerifiedSource(descriptor["id"], source.order, source.phase, digest,
                          tuple(metadata["used_exponents"]), proof, None,
                          lambda i, j: source.entry(divmod(i, source.prime), divmod(j, source.prime)))


@cache
def _load(identifier):
    descriptor = next((row for row in sources() if row["id"] == identifier), None)
    if descriptor is None:
        raise ValueError("unknown complete source; unfinished targets are not registry entries")
    family = descriptor["family"]
    if family == "conference":
        from experiments.scalar_lift.even_petrescu_conference_family import source_exponents

        n, q = descriptor["order"], descriptor["phase"]
        matrix = source_exponents(n, q)
        path = ROOT / "certificates/scalar-lift/even-petrescu-conference-family.json"
        old = next(row for row in json.loads(path.read_text())["complete_known_sources"]
                   if row["order"] == n and row["phase"] == q)
        if hashlib.sha256(matrix.astype(np.uint8).tobytes()).hexdigest() != old["normalized_exponent_sha256"]:
            raise ValueError("the complete conference source differs from its existing certificate")
        proof = {**_conference_proof(q), "outer_H2_tensor": n == 2000,
                 "primitive_cross_ratio_order": q, "prior_certificate": str(path.relative_to(ROOT)),
                 "prior_certificate_sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
        return _dense(descriptor, matrix, proof)
    if family == "common120":
        from experiments.scalar_lift.phase120_common_base import source_fifteen, source_thirty_five

        matrix = source_fifteen() if descriptor["order"] == 15 else source_thirty_five()
        return _dense(descriptor, matrix)
    if family == "sextic-shifted":
        from experiments.scalar_lift.inverse_blocks import known_order_thirteen_source

        original = (2 + 2 * np.eye(3, dtype=np.int64) if descriptor["order"] == 3
                    else np.asarray(known_order_thirteen_source(), dtype=np.int64))
        matrix = (original - 1) % 6
        if not set(map(int, np.unique(matrix))).issubset({0, 1, 3, 4}):
            raise ValueError("the shifted complete sixth-root source left its declared alphabet")
        proof = verify_source(matrix, 6)
        proof.update({"global_exponent_shift_from_unreal_source": -1,
                      "inverse_shift_for_original_unreal_substitution": 1,
                      "declared_allowed_exponents": [0, 1, 3, 4]})
        return _dense(descriptor, matrix, proof)
    if family == "cubic13-tripling":
        matrix, proof = _tripling_source(descriptor["triplings"])
        return _dense(descriptor, matrix, proof)
    if family in ("fifth-chirp", "fifth-tensor"):
        from experiments.scalar_lift.fifth_source_assembly import complete_tensor_source, fifth_seed

        matrix = fifth_seed(40) // 2 if family == "fifth-chirp" else complete_tensor_source(descriptor["phase"])
        return _dense(descriptor, matrix)
    if family == "fifth-descendant":
        source = _load(descriptor["starting_source"])
        tower = ForwardTower(source.exponents, source.phase)
        matrix = tower.materialize(descriptor["steps"])
        # Independent complete Gram check also binds this particular descendant.
        return _dense(descriptor, matrix)
    if family == "character":
        return _character(descriptor)
    raise ValueError("the source family has no exact constructor")


def load_source(descriptor_or_id):
    """Load a canonical descriptor/id and replay its exact source proof once."""
    if isinstance(descriptor_or_id, str):
        identifier = descriptor_or_id
    elif isinstance(descriptor_or_id, dict):
        identifier = descriptor_or_id.get("id")
        canonical = next((row for row in sources() if row["id"] == identifier), None)
        if canonical is None or descriptor_or_id != canonical:
            raise ValueError("source descriptor differs from its verified registry entry")
    else:
        raise ValueError("use a canonical descriptor or source id")
    return _load(identifier)

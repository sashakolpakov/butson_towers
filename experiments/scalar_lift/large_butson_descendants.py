"""Literal large Butson descendants of fixed verified sources; scalar steps only."""

from __future__ import annotations

import argparse
import hashlib
import json
from functools import cache
from pathlib import Path

import numpy as np
from sympy import Poly, Symbol, cyclotomic_poly

from experiments.scalar_lift.forward_butson_descent import (
    _monomial_power,
    available_steps,
    edge_parameters,
    exponent_hash,
    scalar_blocks,
)
from experiments.scalar_lift.vanishing_collision_search import integral_roots


def _json_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


@cache
def fast_scalar_certificate(phase, prime, mode="primary"):
    """Exact algebra proves all cross-Grams without quadratic pair enumeration."""
    edge = edge_parameters(phase, prime, mode)
    q, p = edge["lower_phase"], edge["prime"]
    roots = np.asarray(integral_roots(q), dtype=np.int64)
    blocks = scalar_blocks(phase, p, mode)
    fourier = np.array([[(q // p) * r * s % q for s in range(p)] for r in range(p)])
    for r in range(p):
        for s in range(p):
            value = roots[(fourier[r] - fourier[s]) % q].sum(axis=0)
            expected = np.zeros(roots.shape[1], dtype=np.int64)
            expected[0] = p if r == s else 0
            if not np.array_equal(value, expected):
                raise ValueError("the complete lower Fourier Gram failed")
    one = _monomial_power(1, phase, p, mode)
    for e in range(phase):
        action = _monomial_power(e, phase, p, mode)
        following = _monomial_power(e + 1, phase, p, mode)
        inverse = _monomial_power(-e, phase, p, mode)
        for r, c, a in action:
            _, next_column, next_phase = one[c]
            if following[r] != (r, next_column, (a + next_phase) % q):
                raise ValueError("the exact scalar multiplication recurrence failed")
            if inverse[c] != (c, r, -a % q):
                raise ValueError("the scalar action does not preserve adjoints")
            if not np.array_equal(blocks[e, r], (a + fourier[c]) % q):
                raise ValueError("the literal dense block differs from the action/frame product")
    x = Symbol("x")
    coefficients = list(reversed([int(c) for c in Poly(cyclotomic_poly(phase, x), x).all_coeffs()]))
    value = np.zeros((p, p, roots.shape[1]), dtype=np.int64)
    for e, coefficient in enumerate(coefficients):
        for r, c, a in _monomial_power(e, phase, p, mode):
            value[r, c] += coefficient * roots[a]
    if np.any(value):
        raise ValueError("the represented primitive root fails its cyclotomic polynomial")
    return {**edge, "Phi_upper_of_represented_primitive_root_exact_zero": True,
            "all_scalar_multiplication_and_adjoint_identities_checked": phase,
            "complete_Fourier_Gram_checked": True, "all_literal_root_blocks_checked": phase,
            "cross_Gram_identity": "B_e B_f* = p rho(zeta_upper^(e-f))",
            "cross_Grams_proved": phase**2, "quadratic_pair_enumeration_needed": False,
            "block_table_sha256": exponent_hash(blocks)}


@cache
def _character_proof(prime, first_order, second_order, gamma):
    from experiments.scalar_lift.forward_character_sources import (
        character_source,
        formal_gram_check,
    )

    source = character_source(prime, second_order, first_order=first_order, gamma=gamma)
    return source, {"order": source.order, "phase": source.phase, "full_Gram_verified": True,
                    "source_metadata": source.metadata(), "formal_Gram": formal_gram_check(source)}


def _prepare_source(source, supplied_proof=None):
    if hasattr(source, "coefficients") and hasattr(source, "prime"):
        verified, proof = _character_proof(source.prime, source.first_order, source.second_order, source.gamma)
        if not np.array_equal(source.coefficients, verified.coefficients):
            raise ValueError("the supplied convolution coefficients differ from the verified source")
        data = np.asarray(source.coefficients, dtype=np.int64)
        descriptor = {"kind": "character_difference_table", "prime": source.prime,
                      "order": source.order, "phase": source.phase,
                      "first_order": source.first_order, "second_order": source.second_order,
                      "gamma": source.gamma, "array_sha256": exponent_hash(data),
                      "entry_formula": "coefficient[(column0-row0) mod prime,(column1-row1) mod prime]"}
        # Recipe parameters do not enter the exact matrix key when the literal
        # coefficient table and its indexing are identical.
        identity = {key: descriptor[key] for key in ("kind", "prime", "order", "phase", "array_sha256", "entry_formula")}
        digest = _json_hash(identity)
    elif getattr(source, "exponents", None) is not None:
        data = np.asarray(source.exponents, dtype=np.int64)
        proof = source.proof
        if (proof.get("full_Gram_verified") is not True or proof.get("order") != source.order
                or proof.get("phase") != source.phase or proof.get("exponent_sha256") != exponent_hash(data)):
            raise ValueError("the complete dense-source proof is not bound to these entries")
        descriptor = {"kind": "dense_exponents", "order": source.order, "phase": source.phase,
                      "array_sha256": exponent_hash(data), "registry_id": getattr(source, "id", None)}
        digest = _json_hash({key: descriptor[key] for key in ("kind", "order", "phase", "array_sha256")})
    elif getattr(source, "exponents", 0) is None and str(getattr(source, "id", "")).startswith("character-"):
        from experiments.scalar_lift.forward_character_sources import character_source

        reconstructed = character_source(int(source.id.split("-")[-1]))
        if reconstructed.metadata()["coefficient_array_sha256"] != source.proof["implicit_identity"]["coefficient_array_sha256"]:
            raise ValueError("the registry source differs from its implicit coefficient definition")
        return _prepare_source(reconstructed, supplied_proof)
    else:
        raise ValueError("use an existing complete CharacterSource or verified dense registry source")
    # It is supplementary: the source is independently replayed above.
    if supplied_proof is not None and (supplied_proof.get("order") != descriptor["order"]
                                       or supplied_proof.get("phase") != descriptor["phase"]):
        raise ValueError("the supplementary source proof has different dimensions or phase")
    return data, descriptor, proof, digest


def _base_entry(data, descriptor, row, column):
    if descriptor["kind"] == "dense_exponents":
        return int(data[row, column])
    p = descriptor["prime"]
    r0, r1 = divmod(row, p)
    c0, c1 = divmod(column, p)
    return int(data[(c0-r0) % p, (c1-r1) % p])


class ImplicitForwardTower:
    def __init__(self, source, *, steps=None, source_proof=None, generator_mode="primary"):
        self.source_data, self.source_descriptor, self.source_proof, self.source_hash = _prepare_source(source, source_proof)
        self.original_order = self.source_descriptor["order"]
        self.original_phase = self.source_descriptor["phase"]
        self.original_alphabet = tuple(map(int, np.unique(self.source_data)))
        if steps is None:
            steps, phase = [], self.original_phase
            while available_steps(phase):
                p = min(available_steps(phase))
                steps.append((p, generator_mode))
                phase //= p
        self.steps, self.edges = [], []
        self.tables = [np.arange(self.original_phase, dtype=np.int64)[:, None, None]]
        self.orders, self.phases, self.alphabets = [self.original_order], [self.original_phase], [self.original_alphabet]
        for step in steps:
            p, mode = ((step["prime"], step.get("generator_mode", generator_mode)) if isinstance(step, dict)
                       else step if isinstance(step, (tuple, list)) else (step, generator_mode))
            edge = fast_scalar_certificate(self.phases[-1], p, mode)
            table, width = self.tables[-1], self.tables[-1].shape[1]
            if width * p > 64:
                raise ValueError("the exact composed table exceeds the declared width bound of 64")
            expanded = scalar_blocks(self.phases[-1], p, mode)[table]
            table = expanded.transpose(0, 1, 3, 2, 4).reshape(self.original_phase, width*p, width*p)
            self.steps.append((int(p), mode))
            self.edges.append(edge)
            self.tables.append(table)
            self.orders.append(self.orders[-1] * p)
            self.phases.append(edge["lower_phase"])
            self.alphabets.append(tuple(map(int, np.unique(table[list(self.original_alphabet)]))))

    @property
    def order(self):
        return self.orders[-1]

    @property
    def phase(self):
        return self.phases[-1]

    @property
    def used_alphabet(self):
        return self.alphabets[-1]

    def entry(self, row, column, stage=None):
        stage = len(self.steps) if stage is None else int(stage)
        if (not 0 <= stage < len(self.tables) or int(row) != row or int(column) != column
                or not 0 <= row < self.orders[stage] or not 0 <= column < self.orders[stage]):
            raise ValueError("entry lies outside the selected exact matrix stage")
        width = self.tables[stage].shape[1]
        a, r = divmod(int(row), width)
        b, s = divmod(int(column), width)
        e = _base_entry(self.source_data, self.source_descriptor, a, b)
        return int(self.tables[stage][e, r, s])

    def stage_hash(self, stage=None):
        stage = len(self.steps) if stage is None else int(stage)
        return _json_hash({"source_hash": self.source_hash, "order": self.orders[stage],
                           "phase": self.phases[stage], "original_used_alphabet": self.original_alphabet,
                           "used_block_table_sha256": exponent_hash(self.tables[stage][list(self.original_alphabet)])})

    @property
    def matrix_hash(self):
        return self.stage_hash()

    def stage_metadata(self, stage=None):
        stage = len(self.steps) if stage is None else int(stage)
        n = self.orders[stage]
        samples = [[i * (n-1) // 19, (i*i + 7*i + 11) % n] for i in range(20)]
        samples = [[i, j, self.entry(i, j, stage)] for i, j in samples]
        return {"order": n, "phase": self.phases[stage], "width": self.tables[stage].shape[1],
                "used_alphabet": list(self.alphabets[stage]), "matrix_hash": self.stage_hash(stage),
                "complete_block_table_sha256": exponent_hash(self.tables[stage]),
                "twenty_literal_entries": samples, "sample_sha256": _json_hash(samples),
                "complete_Gram": f"H H*=H* H={n} I", "Gram_entries_proved": n*n,
                "real_sign_endpoint": self.phases[stage] == 2}

    def metadata(self):
        return {"schema": "large-butson-descendant-v1", "source_descriptor": self.source_descriptor,
                "source_hash": self.source_hash, "steps": [list(step) for step in self.steps],
                "stages": [self.stage_metadata(i) for i in range(len(self.tables))],
                "matrix_hash": self.matrix_hash, "order": self.order, "phase": self.phase,
                "dense_large_matrices_materialized": False, "only_scalar_restrictions_used": True,
                "new_Hadamard_order_claimed": False}

    def proof(self):
        return {"source": self.source_proof, "scalar_edges": self.edges,
                "stage_Grams_inherited_exactly": True,
                "every_phase_actually_used_audited_through_literal_tables": True,
                "matrix_hash": self.matrix_hash, "order": self.order, "phase": self.phase}

    def save_bundle(self, path, stage=None):
        stage = len(self.steps) if stage is None else int(stage)
        if not 0 <= stage < len(self.tables):
            raise ValueError("invalid bundle stage")
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        metadata = self.metadata()
        metadata.update({"steps": metadata["steps"][:stage], "stages": metadata["stages"][:stage+1],
                         "matrix_hash": self.stage_hash(stage), "order": self.orders[stage], "phase": self.phases[stage]})
        proof = {"source": self.source_proof, "scalar_edges": self.edges[:stage],
                 "stage_Grams_inherited_exactly": True}
        arrays = {"source_coefficients": self.source_data.astype(np.uint16),
                  "metadata_json": np.frombuffer(json.dumps(metadata, sort_keys=True).encode(), dtype=np.uint8),
                  "proof_json": np.frombuffer(json.dumps(proof, sort_keys=True).encode(), dtype=np.uint8)}
        arrays.update({f"scalar_table_{i}": table.astype(np.uint16) for i, table in enumerate(self.tables[:stage+1])})
        with path.open("wb") as handle:
            np.savez_compressed(handle, **arrays)
        return {"path": str(path), "bytes": path.stat().st_size,
                "file_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                "matrix_hash": self.stage_hash(stage), "order": self.orders[stage], "phase": self.phases[stage]}


class SavedButsonMatrix:
    def __init__(self, path):
        with np.load(path, allow_pickle=False) as stored:
            self.metadata = json.loads(stored["metadata_json"].tobytes())
            self.proof = json.loads(stored["proof_json"].tobytes())
            self.data = stored["source_coefficients"].copy()
            self.tables = [stored[f"scalar_table_{i}"].copy() for i in range(len(self.metadata["steps"])+1)]
        self.order, self.phase = self.metadata["order"], self.metadata["phase"]
        descriptor = self.metadata["source_descriptor"]
        if exponent_hash(self.data) != descriptor["array_sha256"]:
            raise ValueError("the saved complete source coefficient table was modified")
        if descriptor["kind"] == "character_difference_table":
            verified, source_proof = _character_proof(descriptor["prime"], descriptor["first_order"],
                                                      descriptor["second_order"], descriptor["gamma"])
            if not np.array_equal(self.data, verified.coefficients) or self.proof["source"] != source_proof:
                raise ValueError("the stored complete character Gram proof failed replay")
        elif descriptor.get("registry_id"):
            identifier = descriptor["registry_id"]
            from experiments.scalar_lift.catalogue_source_acquisition import SOURCE_IDS
            from experiments.scalar_lift.forward_source_registry import sources

            known = {row["id"] for row in sources()} | set(SOURCE_IDS) | {"saved-unreal9-q6"}
            if identifier in known:
                from experiments.scalar_lift.weighted_scalar_tree import load_tree_source

                verified, _ = load_tree_source(identifier)
                if not np.array_equal(self.data, verified.exponents) or self.proof["source"] != verified.proof:
                    raise ValueError("the stored complete dense-source proof failed replay")
            else:
                # An additional fixed catalogue input is self-contained. Replay
                # its entire cyclotomic Gram rather than trusting a loader name.
                from experiments.scalar_lift.forward_butson_descent import verify_source

                exact = verify_source(self.data, descriptor["phase"])
                if any(self.proof["source"].get(key) != value for key, value in exact.items()):
                    raise ValueError("the additional fixed source failed full Gram replay")
        if not np.array_equal(self.tables[0], np.arange(descriptor["phase"])[:, None, None]):
            raise ValueError("the source scalar table is not the identity")
        for i, table in enumerate(self.tables):
            if exponent_hash(table) != self.metadata["stages"][i]["complete_block_table_sha256"]:
                raise ValueError("a saved scalar table was modified")
            if i:
                p, mode = self.metadata["steps"][i-1]
                previous, width = self.tables[i-1], self.tables[i-1].shape[1]
                q = self.metadata["stages"][i-1]["phase"]
                expanded = scalar_blocks(q, p, mode)[previous].transpose(0, 1, 3, 2, 4).reshape(len(previous), width*p, width*p)
                if not np.array_equal(table, expanded) or self.proof["scalar_edges"][i-1] != fast_scalar_certificate(q, p, mode):
                    raise ValueError("the saved forward restriction cannot be replayed")
        for i, j, value in self.metadata["stages"][-1]["twenty_literal_entries"]:
            if self.entry(i, j) != value:
                raise ValueError("a saved literal matrix entry failed replay")

    def entry(self, row, column):
        if int(row) != row or int(column) != column or not 0 <= row < self.order or not 0 <= column < self.order:
            raise ValueError("entry outside the saved matrix")
        width = self.tables[-1].shape[1]
        a, r = divmod(int(row), width)
        b, s = divmod(int(column), width)
        e = _base_entry(self.data, self.metadata["source_descriptor"], a, b)
        return int(self.tables[-1][e, r, s])


def build_fixed_outputs(directory):
    from experiments.scalar_lift.forward_character_sources import character_source

    directory = Path(directory)
    rows = []
    for prime in (73, 241):
        tower = ImplicitForwardTower(character_source(prime))
        artifacts = [tower.save_bundle(directory / f"bh-{tower.orders[i]}-q{tower.phases[i]}.npz", i)
                     for i in range(1, len(tower.tables))]
        for artifact in artifacts:
            SavedButsonMatrix(artifact["path"])
        rows.append({"tree": tower.metadata(), "proof": tower.proof(), "actual_matrix_files": artifacts})
    return {"schema": "fixed-large-butson-descendants-v1", "outputs": rows,
            "only_existing_sources_and_scalar_restrictions": True, "new_real_Hadamard_claimed": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    result = build_fixed_outputs(args.output_directory)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"actual_outputs": [[a["order"], a["phase"]] for row in result["outputs"] for a in row["actual_matrix_files"]]}))


if __name__ == "__main__":
    main()

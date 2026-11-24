"""Executed scalar restrictions in nonproduct rational coordinate bases.

The same mixed basis is used in the domain and codomain.  Only exact
coordinate conjugation and literal action/metric duplicate pruning occur.
This finite experiment does not cover unequal input/output basis pairs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from itertools import combinations
from pathlib import Path

import numpy as np
from sympy import Matrix, diag

from experiments.scalar_lift.ambient_field_candidates import load_ambient_source, period_field
from experiments.scalar_lift.pure_period_restriction import PeriodRestrictionStage

ROOT = Path(__file__).resolve().parents[2]


def _hash(array):
    return hashlib.sha256(np.asarray(array, dtype="<i8").tobytes()).hexdigest()


def shear_conjugate(array, shear):
    """In place: (I-tE_pq) A (I+tE_pq), without dense multiplication."""
    p, q, t = map(int, shear)
    if p == q or min(p, q) < 0 or max(p, q) >= len(array):
        raise ValueError("invalid off-diagonal shear")
    if int(np.max(np.abs(array), initial=0)) * (1 + abs(t))**2 >= 2**62:
        raise OverflowError("conjugation exceeds its exact int64 bound")
    array[:, q] += t*array[:, p]
    array[p, :] -= t*array[q, :]


def shear_metric(array, shear):
    """In place: (I+tE_pq)^T G (I+tE_pq)."""
    p, q, t = map(int, shear)
    if p == q or min(p, q) < 0 or max(p, q) >= len(array):
        raise ValueError("invalid off-diagonal shear")
    if int(np.max(np.abs(array), initial=0)) * (1 + abs(t))**2 >= 2**62:
        raise OverflowError("metric congruence exceeds its exact int64 bound")
    array[:, q] += t*array[:, p]
    array[q, :] += t*array[p, :]


def recipes():
    """A fixed breadth-first batch with deliberate inverse-shear duplicates."""
    first = [[(p, q, t)] for p, q, values in
             ((0, 5, (-2, -1, 1, 2)), (1, 10, (-2, -1, 1, 2)), (4, 11, (-1, 1)))
             for t in values]
    second = [path + [(5, 10, t)] for path in first[:4] for t in (-1, 1)]
    returns = [path + [(path[-1][0], path[-1][1], -path[-1][2])] for path in first]
    third = [path + [(10, 0, t)] for path in second for t in (-1, 1)]
    return [[], *first, *second, *returns, *third]


class MixedRestriction:
    """An exact implicit matrix reconstructed from a verified source and basis."""

    def __init__(self, field, source, recipe=()):
        self.field, self.source = field, source
        self.recipe = tuple(tuple(map(int, shear)) for shear in recipe)
        self.stage = PeriodRestrictionStage(field, source, 1, True, "period-local")
        self.width, self.order = self.stage.width, self.stage.order
        self.actions = {e: self.stage.scalar_table(e)[:, :, 0].copy()
                        for e in source.used_alphabet}
        self.metric = self.stage.metric_table()[:, :, 0].copy()
        for shear in self.recipe:
            p, q, _ = shear
            if p//4 == q//4 or p % 4 == q % 4:
                raise ValueError("every declared shear must mix both tensor factors")
            for action in self.actions.values():
                shear_conjugate(action, shear)
            shear_metric(self.metric, shear)

    def entry(self, row, column):
        if type(row) is not int or type(column) is not int or not 0 <= row < self.order or not 0 <= column < self.order:
            raise ValueError("matrix coordinate outside its order")
        i, a = divmod(row, self.width)
        j, b = divmod(column, self.width)
        return int(self.actions[self.source.entry(i, j)][a, b])

    def literal_hash(self):
        digest = hashlib.sha256(self.source.source_hash.encode())
        for exponent in self.source.used_alphabet:
            digest.update(np.asarray(self.actions[exponent], dtype="<i8").tobytes())
        digest.update(np.asarray(self.metric, dtype="<i8").tobytes())
        return digest.hexdigest()

    def nonproduct_basis_witness(self):
        """A nonzero 2x2 tensor minor of one actual mixed basis vector."""
        columns = {j: {j: 1} for p, q, _ in self.recipe for j in (p, q)}
        for p, q, t in self.recipe:
            updated = Counter(columns[q])
            for row, value in columns[p].items():
                updated[row] += t*value
            columns[q] = {row: value for row, value in updated.items() if value}
        for column, vector in columns.items():
            periods = sorted({row//4 for row in vector})
            locals_ = sorted({row % 4 for row in vector})
            for rows in combinations(periods, 2):
                for cols in combinations(locals_, 2):
                    minor = [[vector.get(4*r+c, 0) for c in cols] for r in rows]
                    determinant = minor[0][0]*minor[1][1]-minor[0][1]*minor[1][0]
                    if determinant:
                        return {"basis_column": column, "period_rows": list(rows),
                                "local_columns": list(cols), "minor": minor,
                                "nonzero_minor_determinant": determinant}
        return None

    def census(self):
        frequencies = Counter(map(int, self.source.exponents.flat))
        counts = Counter()
        for exponent, frequency in frequencies.items():
            values, multiplicities = np.unique(self.actions[exponent], return_counts=True)
            for value, multiplicity in zip(values, multiplicities, strict=True):
                counts[str(int(value))] += frequency*int(multiplicity)
        if sum(counts.values()) != self.order**2:
            raise ArithmeticError("complete entry census has the wrong matrix size")
        return dict(sorted(counts.items(), key=lambda item: int(item[0])))

    def replay_inverse_shears(self):
        """Check all action and metric arrays, not only stored hashes."""
        restored = {e: a.copy() for e, a in self.actions.items()}
        metric = self.metric.copy()
        for p, q, t in reversed(self.recipe):
            for action in restored.values():
                shear_conjugate(action, (p, q, -t))
            shear_metric(metric, (p, q, -t))
        if (not np.array_equal(metric, self.stage.metric_table()[:, :, 0])
                or any(not np.array_equal(a, self.stage.scalar_table(e)[:, :, 0])
                       for e, a in restored.items())):
            raise ArithmeticError("mixed-basis inverse replay failed")
        return True


def independent_small_verification(mixed):
    """Dense rational arithmetic is used only for the order-84 control."""
    if mixed.order > 100:
        raise ValueError("independent dense verification is limited to the small control")
    change = Matrix.eye(mixed.width)
    for p, q, t in mixed.recipe:
        shear = Matrix.eye(mixed.width)
        shear[p, q] = t
        change *= shear
    inverse = change.inv()
    if change.det() != 1:
        raise ArithmeticError("the declared basis does not have determinant one")
    for exponent, action in mixed.actions.items():
        expected = inverse*Matrix(mixed.stage.scalar_table(exponent)[:, :, 0].tolist())*change
        if expected != Matrix(action.tolist()):
            raise ArithmeticError("independent full scalar conjugation failed")
    G = Matrix(mixed.metric.tolist())
    if change.T*Matrix(mixed.stage.metric_table()[:, :, 0].tolist())*change != G:
        raise ArithmeticError("independent full metric congruence failed")
    matrix = Matrix([[mixed.entry(i, j) for j in range(mixed.order)] for i in range(mixed.order)])
    weight = diag(*([G.inv()]*mixed.source.order))
    if matrix*weight*matrix.T != mixed.source.order*weight:
        raise ArithmeticError("independent complete weighted Gram failed")
    if abs(matrix.det()) != mixed.source.order**(mixed.order//2):
        raise ArithmeticError("independent determinant check failed")
    return {"full_basis_conjugation_checked": True, "complete_weighted_Gram_checked": True,
            "complete_matrix_determinant_checked": True, "basis_determinant": 1,
            "matrix_order": mixed.order}


def run_case(field, source):
    seen = {}
    records = []
    for path in recipes():
        mixed = MixedRestriction(field, source, path)
        path_json = [list(shear) for shear in path]
        mixed.replay_inverse_shears()
        key = mixed.literal_hash()
        duplicate = key in seen
        if duplicate:
            records.append({"basis_shears": path_json, "literal_duplicate_of": seen[key],
                            "action_metric_sha256": key})
            continue
        seen[key] = len(records)
        counts = mixed.census()
        if set(counts) <= {"-1", "1"}:
            raise ArithmeticError("an all-sign matrix contradicts the checked odd determinant invariant")
        witness = mixed.nonproduct_basis_witness()
        if path and witness is None:
            raise ArithmeticError("a retained declared mixed basis has no nonproduct vector")
        records.append({"basis_shears": path_json, "action_metric_sha256": key,
                        "nonproduct_basis_vector_witness": witness,
                        "metric_sha256": _hash(mixed.metric), "complete_entry_counts": counts,
                        "all_entries_signs": False, "inverse_shears_replayed_exactly": True})
    return {"source_id": source.id, "source_hash": source.source_hash,
            "source_Gram_replayed": source.proof["full_Gram_verified"],
            "period_degree": field.degree, "period_conductor": field.conductor,
            "period_norm_witness": list(field.norm_witness),
            "faithful_period_multiplication_sha256": _hash(field.multiplication),
            "auxiliary_generator_in_mixed_basis": "S^-1 (M_eta0 tensor I_4) S",
            "matrix_order": source.order*4*field.degree, "basis_width": 4*field.degree,
            "node_visits": len(records), "distinct_action_metric_nodes": len(seen),
            "literal_duplicates_pruned": len(records)-len(seen), "nodes": records,
            "determinant_absolute_value": f"{source.order}^{source.order*4*field.degree//2}",
            "determinant_is_odd": bool(source.order % 2),
            "all_same_basis_restrictions_excluded_from_sign_endpoints": bool(source.order % 2),
            "all_same_basis_restrictions_excluded_from_Hadamard_Gram": True,
            "unequal_domain_codomain_bases_excluded": False}


def certificate():
    small_field = period_field(3, 31, (2, 3, 1))
    small_source = load_ambient_source("catalogue-p7-q12")
    cases = [run_case(small_field, small_source),
             run_case(period_field(179), load_ambient_source("unreal13-q6"))]
    verification = independent_small_verification(MixedRestriction(small_field, small_source, recipes()[-1]))
    return {"schema": "mixed-same-basis-scalar-restriction-v1", "cases": cases,
            "finite_declared_batch_finished": True, "pending_recipes": 0,
            "same_domain_codomain_basis": True, "right_frame_or_normalization_inserted": False,
            "only_operations": "scalar restriction in a declared mixed basis and exact action/metric duplicate pruning",
            "basis_definition": "S=product(I+t E_pq); every shear changes both period and cyclotomic coordinates",
            "entry_accessor": "MixedRestriction(field,source,basis_shears).entry(i,j)",
            "small_independent_verification": verification,
            "node_visits": sum(c["node_visits"] for c in cases),
            "distinct_nodes": sum(c["distinct_action_metric_nodes"] for c in cases),
            "literal_duplicates_pruned": sum(c["literal_duplicates_pruned"] for c in cases),
            "new_sign_Hadamard_found": False,
            "scope": "Finite nonproduct same-basis tests; determinant proof closes all same-basis variants, not unequal basis pairs"}


def load_saved_node(path, case=0, node=-1):
    """Reconstruct and verify a saved literal matrix, including duplicates."""
    record = json.loads(Path(path).read_text())["cases"][case]
    source = load_ambient_source(record["source_id"])
    field = period_field(record["period_degree"], record["period_conductor"],
                         tuple(record["period_norm_witness"]))
    if (source.source_hash != record["source_hash"]
            or _hash(field.multiplication) != record["faithful_period_multiplication_sha256"]):
        raise ValueError("saved source or field does not replay")
    definition = record["nodes"][node]
    matrix = MixedRestriction(field, source, definition["basis_shears"])
    if matrix.literal_hash() != definition["action_metric_sha256"]:
        raise ValueError("saved mixed scalar actions or metric do not replay")
    if "complete_entry_counts" in definition and matrix.census() != definition["complete_entry_counts"]:
        raise ValueError("saved complete matrix entry census does not replay")
    return matrix


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    choice = parser.add_mutually_exclusive_group(required=True)
    choice.add_argument("--output", type=Path)
    choice.add_argument("--check", type=Path)
    args = parser.parse_args()
    report = certificate()
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    elif json.loads(args.check.read_text()) != report:
        raise ValueError("the mixed-basis certificate differs from exact replay")
    print(json.dumps({key: report[key] for key in
                      ("node_visits", "distinct_nodes", "literal_duplicates_pruned", "pending_recipes", "new_sign_Hadamard_found")}))


if __name__ == "__main__":
    main()
